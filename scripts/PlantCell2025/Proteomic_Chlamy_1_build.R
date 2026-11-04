# =============================================================================
# Proteomic_Chlamy_1_build.R
# TMTpro proteome of Dupuis, Ojeda et al. 2025, Plant Cell 37: koaf086
# (MassIVE MSV000095622): from MS-GF+ identifications (mzTab) and MASIC
# reporter ions to a filtered PSM table, plus reconstruction of the
# TMT channel -> sample design by correlation with Supplementary Data Set 2.
#
# Inputs (BASE folder):
#   Result_Files/  72 *.mzTab        (MS-GF+, one per plex x fraction)
#   MASIC_Results/ 72 *_ReporterIons.txt (+ *_SICstats.txt, optional)
#   supplementary/ SupplementaryDataSet1_OmicsSampleMetadata.xlsx
#                  SupplementaryDataSet2_mRNA_protein_abundances.xlsx
# Outputs (BASE/processed/):
#   psm_tmt.tsv.gz            filtered PSMs with the 16 reporter intensities
#   qc_files.tsv              PSMs and decoys per file
#   channel_design.tsv        reconstructed channel -> sample table
#   channel_correlation.png   correlation of every channel with every sample
# =============================================================================

# ---- 0. Settings -------------------------------------------------------------
BASE     <- Sys.getenv("PROT_BASE", path.expand("~/BMSi/PlantCell2025/TMT"))
SUPP_DIR <- file.path(BASE, "supplementary")
OUT_DIR  <- file.path(BASE, "processed")
FDR      <- 0.01   # PSM-level q-value, recomputed over all 72 files
MIN_INTERFERENCE <- 0.9  # MASIC interference score: 1 = clean precursor; keep >= this (only if SICstats present)

for (p in c("data.table", "readxl", "clue", "pheatmap"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(data.table); library(readxl)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# TMTpro 16-plex: MASIC column masses, in order
CHANNELS <- c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N",
              "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N")

parse_run <- function(f) {
  m <- regmatches(basename(f), regexec("plex_(\\d+)_f(\\d+)", basename(f)))[[1]]
  list(plex = as.integer(m[2]), fraction = as.integer(m[3]))
}

# ---- 1. Identifications: mzTab PSM section -----------------------------------
mz_files <- list.files(file.path(BASE, "Result_Files"), "\\.mzTab$", full.names = TRUE, recursive = TRUE)
mz_files <- mz_files[!duplicated(basename(mz_files))]
stopifnot("expected 72 mzTab files (3 plexes x 24 fractions)" = length(mz_files) == 72 || Sys.getenv("PROT_BASE") != "")
cat("mzTab files:", length(mz_files), "\n")
stopifnot(length(mz_files) > 0)

# The meaning of search_engine_score[n] is declared in each file's MTD header and
# its order is NOT the same in every file: find the column by its CV accession.
score_column <- function(f, cv) {
  mtd <- grep("^MTD\tpsm_search_engine_score\\[", readLines(f, n = 300), value = TRUE)
  hit <- grep(cv, mtd, fixed = TRUE, value = TRUE)
  if (length(hit) != 1) stop("No ", cv, " declared in ", basename(f))
  sub("^MTD\tpsm_(search_engine_score\\[\\d+\\]).*", "\\1", hit)
}

read_mztab_psm <- function(f) {
  d <- fread(f, skip = "PSH", sep = "\t", quote = "", fill = TRUE, showProgress = FALSE)
  d <- d[PSH == "PSM"]
  r <- parse_run(f)
  ev_col <- score_column(f, "MS:1002052")   # MS-GF:SpecEValue (lower = better)
  ev <- as.numeric(d[[ev_col]])
  if (any(ev < 0 | ev > 1, na.rm = TRUE)) stop("SpecEValue out of [0, 1] in ", basename(f))
  data.table(
    plex      = r$plex,
    fraction  = r$fraction,
    scan      = as.integer(sub(".*scan=(\\d+).*", "\\1", d$spectra_ref)),
    sequence  = d$sequence,
    mods      = d$modifications,
    charge    = as.integer(d$charge),
    accession = d$accession,
    spec_evalue = ev,
    decoy     = as.integer(d$`opt_global_cv_MS:1002217_decoy_peptide`) == 1L |
                grepl("^XXX_", d$accession)
  )
}
psm <- rbindlist(lapply(mz_files, read_mztab_psm))

# One row per spectrum: a peptide shared by several proteins appears once per protein
psm[, gene := sub("_4532\\..*$", "", sub("^XXX_", "", accession))]
spec <- psm[, .(sequence = sequence[1], mods = mods[1], charge = charge[1],
                spec_evalue = min(spec_evalue),
                decoy = all(decoy),
                genes = paste(sort(unique(gene[!decoy])), collapse = ";")),
            by = .(plex, fraction, scan)]

# Global target-decoy FDR over the whole experiment (MS-GF+ q-values are per file)
setorder(spec, spec_evalue, na.last = TRUE)
spec[, fdr := cumsum(decoy) / pmax(cumsum(!decoy), 1)]
spec[, qvalue := rev(cummin(rev(fdr)))]
qc <- spec[, .(spectra = .N, decoys = sum(decoy),
               targets_q01 = sum(!decoy & qvalue <= FDR)), by = .(plex, fraction)]
setorder(qc, plex, fraction)
fwrite(qc, file.path(OUT_DIR, "qc_files.tsv"), sep = "\t")
cat(sprintf("Identified spectra: %d targets + %d decoys; at q <= %.2f: %d PSMs\n",
            spec[decoy == FALSE, .N], spec[decoy == TRUE, .N], FDR, spec[decoy == FALSE & qvalue <= FDR, .N]))

spec <- spec[decoy == FALSE & qvalue <= FDR]
spec[, n_genes := lengths(strsplit(genes, ";"))]
cat(sprintf("Peptides unique to one gene: %.1f%% of PSMs (shared ones are kept, flagged)\n",
            100 * spec[, mean(n_genes == 1)]))

# ---- 2. Quantification: MASIC reporter ions ----------------------------------
ri_files <- list.files(file.path(BASE, "MASIC_Results"), "_ReporterIons\\.txt$", full.names = TRUE, recursive = TRUE)
cat("MASIC ReporterIons files:", length(ri_files), "\n")
stopifnot(length(ri_files) > 0)

read_masic <- function(f) {
  d <- fread(f, showProgress = FALSE)
  ion_cols <- grep("^Ion_1[0-9]{2}\\.[0-9]+$", names(d), value = TRUE)  # corrected intensities
  stopifnot(length(ion_cols) == 16)
  r <- parse_run(f)
  out <- d[, c("ScanNumber", ion_cols), with = FALSE]
  setnames(out, c("scan", CHANNELS))
  out[, `:=`(plex = r$plex, fraction = r$fraction)]
  # optional interference score from the matching SICstats file
  sic <- sub("_ReporterIons\\.txt$", "_SICstats.txt", f)
  if (file.exists(sic)) {
    s <- fread(sic, showProgress = FALSE)
    icol <- grep("interference", names(s), ignore.case = TRUE, value = TRUE)[1]
    scol <- intersect(c("FragScanNumber", "Frag_Scan_Number", "FragScan"), names(s))[1]
    if (!is.na(icol) && !is.na(scol)) {
      s <- s[, .(scan = get(scol), interference = get(icol))]
      out <- merge(out, unique(s, by = "scan"), by = "scan", all.x = TRUE)
    }
  }
  out
}
ri <- rbindlist(lapply(ri_files, read_masic), fill = TRUE)

tmt <- merge(spec, ri, by = c("plex", "fraction", "scan"))
cat(sprintf("PSMs with reporter ions: %d of %d\n", nrow(tmt), nrow(spec)))
for (ch in CHANNELS) set(tmt, which(tmt[[ch]] <= 0), ch, NA_real_)
tmt <- tmt[rowSums(!is.na(tmt[, ..CHANNELS])) >= 8]   # at least half the channels quantified
if ("interference" %in% names(tmt)) {
  cat(sprintf("Interference score available: removing %d PSMs below %.2f\n",
              tmt[interference < MIN_INTERFERENCE, .N], MIN_INTERFERENCE))
  tmt <- tmt[is.na(interference) | interference >= MIN_INTERFERENCE]
}
fwrite(tmt, file.path(OUT_DIR, "psm_tmt.tsv.gz"), sep = "\t")
cat(sprintf("Final PSMs: %d; peptides: %d; genes: %d\n",
            nrow(tmt), uniqueN(tmt$sequence), uniqueN(tmt[n_genes == 1, genes])))

# ---- 3. Protein x channel matrix (relative, within plex) ---------------------
# log2; equal loading per channel (column median); each PSM relative to its own
# mean across the 16 channels; protein = median of its unique PSMs.
long <- melt(tmt[n_genes == 1, c("plex", "genes", CHANNELS), with = FALSE],
             id.vars = c("plex", "genes"), variable.name = "channel", value.name = "int")
long[, psm := rowid(plex, channel)]
long <- long[!is.na(int)]
long[, l2 := log2(int)]
long[, l2 := l2 - median(l2), by = .(plex, channel)]
long[, rel := l2 - mean(l2), by = .(plex, psm)]
prot <- long[, .(rel = median(rel)), by = .(genes, plex, channel)]
prot[, col := sprintf("plex%02d_%s", plex, channel)]
M <- as.matrix(dcast(prot, genes ~ col, value.var = "rel"), rownames = "genes")

# ---- 4. Supplementary Data Set 2: published protein values --------------------
meta <- as.data.table(read_excel(file.path(SUPP_DIR, "SupplementaryDataSet1_OmicsSampleMetadata.xlsx"),
                                 sheet = "OmicsSampleMetadata"))
sd2  <- as.data.table(read_excel(file.path(SUPP_DIR, "SupplementaryDataSet2_mRNA_protein_abundances.xlsx"),
                                 sheet = 2, na = "NA"))
pcols <- meta$protein_ColumnHeader
S <- as.matrix(sd2[, ..pcols]); storage.mode(S) <- "numeric"
rownames(S) <- sd2$Cre_GID; colnames(S) <- meta$SampleID
S <- S[rowSums(!is.na(S)) >= 30, ]
S <- sweep(S, 2, apply(S, 2, median, na.rm = TRUE))
S <- S - rowMeans(S, na.rm = TRUE)

common <- intersect(rownames(M), rownames(S))
cat(sprintf("Genes in common with Supplementary Data Set 2: %d\n", length(common)))
R <- cor(M[common, ], S[common, ], use = "pairwise.complete.obs")

# ---- 5. Assignment ------------------------------------------------------------
# Each published sample is assigned to one channel (each channel used at most once,
# Hungarian algorithm); the 3 channels left over are the references.
# The published values still carry the plex effect (a per-protein offset shared by
# the 15 samples of a plex), which our within-plex profiles do not have. So we
# iterate: assign -> centre the published values within the assigned plexes ->
# correlate again -> reassign, until the assignment no longer changes.
assign_once <- function(R) {                     # R: channels x samples
  a <- clue::solve_LSAP(t(R) + 1, maximum = TRUE)
  setNames(rownames(R)[as.integer(a)], colnames(R))
}
centre_within <- function(X, grp) {
  for (g in unique(grp)) X[, grp == g] <- X[, grp == g] - rowMeans(X[, grp == g], na.rm = TRUE)
  X
}
plex_of <- function(cols) sub("_.*", "", cols)

R0  <- R
asg <- assign_once(R)
for (it in 1:10) {
  Sc  <- centre_within(S[common, ], plex_of(asg[colnames(S)]))
  R   <- cor(M[common, ], Sc, use = "pairwise.complete.obs")
  new <- assign_once(R)
  changed <- sum(new != asg[names(new)])
  cat(sprintf("Iteration %d: %d samples reassigned\n", it, changed))
  asg <- new
  if (changed == 0) break
}
ref <- setdiff(rownames(R), asg)
cat("Reference (unassigned) channel per plex:", paste(ref, collapse = ", "), "\n")
if (any(table(plex_of(ref)) != 1)) warning("References are not one per plex: check the assignment")
spread <- apply(M[common, ], 2, sd, na.rm = TRUE)
cat("Lowest-spread channel per plex (a pooled reference should be flat):",
    paste(sapply(split(names(spread), plex_of(names(spread))), function(cc) cc[which.min(spread[cc])]),
          collapse = ", "), "\n")

design <- data.table(SampleID = names(asg), column = unname(asg))
design[, `:=`(plex = as.integer(sub("plex(\\d+)_.*", "\\1", column)),
              channel = sub(".*_", "", column))]
design[, r_initial  := R0[cbind(column, SampleID)]]
design[, r_assigned := R[cbind(column, SampleID)]]
design[, r_second := sapply(seq_len(.N), function(i) max(R[column[i], colnames(R) != SampleID[i]]))]
design[, best_is_assigned := colnames(R)[apply(R, 1, which.max)][match(column, rownames(R))] == SampleID]
design <- rbind(design,
                data.table(column = ref, SampleID = "REFERENCE",
                           plex = as.integer(sub("plex(\\d+)_.*", "\\1", ref)),
                           channel = sub(".*_", "", ref)), fill = TRUE)
design <- merge(design, meta[, .(SampleID, Population, TimeCode, ExperimentalReplicate)],
                by = "SampleID", all.x = TRUE)
design[, channel := factor(channel, CHANNELS)]
setorder(design, plex, channel)
fwrite(design, file.path(OUT_DIR, "channel_design.tsv"), sep = "\t")

cat("\nAssignment quality:\n")
print(design[SampleID != "REFERENCE",
             .(median_r_initial = round(median(r_initial), 3),
               median_r = round(median(r_assigned), 3), min_r = round(min(r_assigned), 3),
               median_second = round(median(r_second), 3),
               argmax_agrees = sum(best_is_assigned), n = .N), by = plex])
cat("\nBalance (each plex should hold one replicate of each light x time state):\n")
print(dcast(design[SampleID != "REFERENCE"], Population + TimeCode ~ plex, fun.aggregate = length))

pheatmap::pheatmap(R[order(rownames(R)), order(colnames(R))], cluster_rows = FALSE, cluster_cols = FALSE,
                   fontsize = 6, main = "Correlation: TMT channels vs published samples",
                   filename = file.path(OUT_DIR, "channel_correlation.png"), width = 11, height = 9)
cat("\nDone. Outputs in", OUT_DIR, "\n")
