# =============================================================================
# Proteomic_Chlamy_TMT.R
# TMTpro 16-plex proteome of Chlamydomonas acclimated to three light intensities
# across the diurnal cycle (Dupuis, Ojeda et al. 2025, Plant Cell 37: koaf086;
# MassIVE MSV000095622): from the public identification and quantification files
# to differential abundance, with two pipelines and two relevance criteria, as in
# the label-free practical (PC3):
#
#   pipeline A  ratio to the plex reference + limma          (step by step)
#   pipeline B  MSstatsTMT                                    (all in one)
#   criterion   FC: |log2FC| >= LFC     Bayesian: P(|log2FC| > LFC) >= P_BAYES
#
# Sections
#   1. Packages
#   2. PSMs: MS-GF+ identifications (mzTab) + MASIC reporter ions        [cached]
#   3. Channel -> sample design (reconstructed if not supplied)          [cached]
#   4. Pipeline A: ratio to reference, QC, limma
#   5. Pipeline B: MSstatsTMT                                            [cached]
#   6. Relevance: FC and Bayesian criteria, both pipelines               [cached]
#   7. Comparison of pipelines and criteria
#
# Heavy steps are cached: if their output exists they are skipped (set the
# FORCE_* switches to TRUE to recompute). Everything that may need editing is in
# section 0.
#
# Folder layout expected in BASE (as downloaded from MassIVE MSV000095622):
#   Result_Files/   72 *.mzTab               MS-GF+ PSMs, one per plex x fraction
#   MASIC_Results/  72 *_ReporterIons.txt    reporter-ion intensities per scan
#   supplementary/  SupplementaryDataSet1_OmicsSampleMetadata.xlsx
#                   SupplementaryDataSet2_mRNA_protein_abundances.xlsx
#                   (from the article; used only to recover the channel design)
# =============================================================================

# ---- 0. Settings (edit here) ---------------------------------------------------
BASE <- Sys.getenv("PROT_BASE", path.expand("~/BMSi/PlantCell2025/TMT"))

PSM_FDR      <- 0.01        # PSM q-value (target-decoy, recomputed over all files)
MIN_CHANNELS <- 12          # pipeline A: a PSM needs >= this many of its 15 sample channels
MIN_PLEX     <- 2           # a protein must be quantified in >= 2 of the 3 plexes
FDR          <- 0.05        # adjusted p-value for the statistical tests
LFC          <- log2(1.5)   # 0.58: relevance threshold (FC 1.5, see note below)
P_BAYES      <- 0.9         # Bayesian criterion: P(|log2FC| > LFC) >= P_BAYES
SEED         <- 2026

# Why FC 1.5 and not 2 (the threshold of PC3): in TMT measured at MS2 level,
# peptides co-isolated with the target add their own reporter ions, which pulls
# every ratio towards 1:1 ("ratio compression"). A true 2-fold change is often
# measured as ~1.5-fold, so FC 2 would be much stricter here than in label-free.

FORCE_PSM      <- FALSE     # re-read the 72 + 72 files (~5 min)
FORCE_DESIGN   <- FALSE     # reconstruct the channel design again
FORCE_MSSTATS  <- FALSE     # rerun MSstatsTMT (long)
FORCE_BAYES    <- FALSE     # refit the Bayesian models (long)

# If the authors provide the true channel design, save it as a tab-separated file
# with columns plex, channel, SampleID (SampleID = REFERENCE for the pool) and set:
DESIGN_FILE  <- NULL        # e.g. file.path(BASE, "design_from_authors.tsv")

# ---- 1. Packages ------------------------------------------------------------------
cran <- c("data.table", "readxl", "clue", "ggplot2", "pheatmap", "rstan")
bioc <- c("limma", "MSstatsTMT")
for (p in cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(limma)
})
theme_set(theme_bw(base_size = 11))
set.seed(SEED)

D  <- function(...) { p <- file.path(BASE, "results_TMT", ...); dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE); p }
PROC <- file.path(BASE, "processed"); dir.create(PROC, showWarnings = FALSE)
CHANNELS  <- c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N",
               "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N")
LIGHTS    <- c("LL", "ML", "HL")
TIMES     <- c("CTm02", "CTp02", "CTp06", "CTp10", "CTm10")     # -2, +2, +6, +10, +14 h after dawn
HOURS     <- c(CTm02 = -2, CTp02 = 2, CTp06 = 6, CTp10 = 10, CTm10 = 14)
LIGHT_COL <- c(LL = "#2A6F97", ML = "#7F8A99", HL = "#B5541B")
CONTRASTS <- c("HL_vs_LL", "HL_vs_ML", "ML_vs_LL")
banner <- function(x) cat("\n==== ", x, " ", strrep("=", max(3, 70 - nchar(x))), "\n", sep = "")

# =================================================================================
# 2. PSMs: identifications + reporter ions
# =================================================================================
banner("2. PSMs")
PSM_FILE <- file.path(PROC, "psm_tmt.tsv.gz")
parse_run <- function(f) {
  m <- regmatches(basename(f), regexec("plex_(\\d+)_f(\\d+)", basename(f)))[[1]]
  list(plex = as.integer(m[2]), fraction = as.integer(m[3]))
}

if (FORCE_PSM || !file.exists(PSM_FILE)) {
  mz_files <- list.files(file.path(BASE, "Result_Files"), "\\.mzTab$", full.names = TRUE, recursive = TRUE)
  mz_files <- mz_files[!duplicated(basename(mz_files))]
  ri_files <- list.files(file.path(BASE, "MASIC_Results"), "_ReporterIons\\.txt$", full.names = TRUE, recursive = TRUE)
  cat("mzTab files:", length(mz_files), "| MASIC ReporterIons files:", length(ri_files), "\n")
  stopifnot(length(mz_files) > 0, length(ri_files) > 0)
  if (length(mz_files) != 72 || length(ri_files) != 72) warning("Expected 72 + 72 files (3 plexes x 24 fractions)")

  # The meaning of search_engine_score[n] is declared in each file's MTD header and
  # its order is NOT the same in every file: find the column by its CV accession.
  score_column <- function(f, cv) {
    mtd <- grep("^MTD\tpsm_search_engine_score\\[", readLines(f, n = 300), value = TRUE)
    hit <- grep(cv, mtd, fixed = TRUE, value = TRUE)
    if (length(hit) != 1) stop("No ", cv, " declared in ", basename(f))
    sub("^MTD\tpsm_(search_engine_score\\[\\d+\\]).*", "\\1", hit)
  }
  read_mztab_psm <- function(f) {
    d <- fread(f, skip = "PSH", sep = "\t", quote = "", fill = TRUE, showProgress = FALSE)[PSH == "PSM"]
    ev <- as.numeric(d[[score_column(f, "MS:1002052")]])     # MS-GF:SpecEValue, lower = better
    if (any(ev < 0 | ev > 1, na.rm = TRUE)) stop("SpecEValue out of [0, 1] in ", basename(f))
    r <- parse_run(f)
    data.table(plex = r$plex, fraction = r$fraction,
               scan = as.integer(sub(".*scan=(\\d+).*", "\\1", d$spectra_ref)),
               sequence = d$sequence, mods = d$modifications, charge = as.integer(d$charge),
               accession = d$accession, spec_evalue = ev,
               decoy = as.integer(d$`opt_global_cv_MS:1002217_decoy_peptide`) == 1L | grepl("^XXX_", d$accession))
  }
  psm <- rbindlist(lapply(mz_files, read_mztab_psm))
  psm[, gene := sub("_4532\\..*$", "", sub("^XXX_", "", accession))]
  # one row per spectrum (a peptide shared by several proteins appears once per protein)
  spec <- psm[, .(sequence = sequence[1], mods = mods[1], charge = charge[1], spec_evalue = min(spec_evalue),
                  decoy = all(decoy), genes = paste(sort(unique(gene[!decoy])), collapse = ";")),
              by = .(plex, fraction, scan)]
  rm(psm)
  # global target-decoy FDR (MS-GF+ q-values are per file)
  setorder(spec, spec_evalue, na.last = TRUE)
  spec[, fdr := cumsum(decoy) / pmax(cumsum(!decoy), 1)]
  spec[, qvalue := rev(cummin(rev(fdr)))]
  fwrite(spec[, .(spectra = .N, decoys = sum(decoy), PSMs_passing = sum(!decoy & qvalue <= PSM_FDR)),
              by = .(plex, fraction)][order(plex, fraction)], D("0_PSMs", "PSMs_per_file.tsv"), sep = "\t")
  cat(sprintf("Identified spectra: %d targets + %d decoys; PSMs at q <= %.2f: %d\n",
              spec[decoy == FALSE, .N], spec[decoy == TRUE, .N], PSM_FDR, spec[decoy == FALSE & qvalue <= PSM_FDR, .N]))
  spec <- spec[decoy == FALSE & qvalue <= PSM_FDR]
  spec[, n_genes := lengths(strsplit(genes, ";"))]

  read_masic <- function(f) {
    d <- fread(f, showProgress = FALSE)
    ion_cols <- grep("^Ion_1[0-9]{2}\\.[0-9]+$", names(d), value = TRUE)   # isotope-corrected intensities
    stopifnot(length(ion_cols) == 16)
    out <- d[, c("ScanNumber", ion_cols), with = FALSE]
    setnames(out, c("scan", CHANNELS))
    r <- parse_run(f); out[, `:=`(plex = r$plex, fraction = r$fraction)]
    out
  }
  ri <- rbindlist(lapply(ri_files, read_masic))
  tmt <- merge(spec, ri, by = c("plex", "fraction", "scan"))
  rm(ri, spec)
  for (ch in CHANNELS) set(tmt, which(tmt[[ch]] <= 0), ch, NA_real_)
  tmt <- tmt[rowSums(!is.na(tmt[, ..CHANNELS])) >= 8]
  fwrite(tmt, PSM_FILE, sep = "\t")
} else cat("Using cached", PSM_FILE, "\n")

tmt <- fread(PSM_FILE)
cat(sprintf("PSMs: %d | peptides: %d | genes (unique peptides): %d\n",
            nrow(tmt), uniqueN(tmt$sequence), uniqueN(tmt[n_genes == 1, genes])))
tmt <- tmt[n_genes == 1 & grepl("^Cre", genes)]     # unique peptides of Chlamydomonas genes (drops contaminants)
cat(sprintf("Kept (peptide unique to one Chlamydomonas gene): %d PSMs\n", nrow(tmt)))

# =================================================================================
# 3. Channel -> sample design
# =================================================================================
banner("3. Channel design")
DESIGN_CACHE <- file.path(PROC, "channel_design.tsv")
meta <- as.data.table(readxl::read_excel(file.path(BASE, "supplementary", "SupplementaryDataSet1_OmicsSampleMetadata.xlsx"),
                                         sheet = "OmicsSampleMetadata"))

if (!is.null(DESIGN_FILE)) {
  design <- fread(DESIGN_FILE, colClasses = list(character = "channel"))
  cat("Design supplied by the user:", DESIGN_FILE, "\n")
} else if (!FORCE_DESIGN && file.exists(DESIGN_CACHE)) {
  design <- fread(DESIGN_CACHE, colClasses = list(character = "channel"))
  cat("Using cached reconstructed design:", DESIGN_CACHE, "\n")
} else {
  # The channel -> sample table was not deposited. It is recovered by correlating each
  # channel (protein profile relative to its plex) with the published protein values.
  long <- melt(tmt[, c("plex", "genes", CHANNELS), with = FALSE], id.vars = c("plex", "genes"),
               variable.name = "channel", value.name = "int", variable.factor = FALSE)
  long[, psm := rowid(plex, channel)]
  long <- long[!is.na(int)][, l2 := log2(int)]
  long[, l2 := l2 - median(l2), by = .(plex, channel)]
  long[, rel := l2 - mean(l2), by = .(plex, psm)]
  pr <- long[, .(rel = median(rel)), by = .(genes, plex, channel)]
  pr[, col := sprintf("plex%02d_%s", plex, channel)]
  Mch <- as.matrix(dcast(pr, genes ~ col, value.var = "rel"), rownames = "genes")
  rm(long, pr)
  sd2 <- as.data.table(readxl::read_excel(file.path(BASE, "supplementary", "SupplementaryDataSet2_mRNA_protein_abundances.xlsx"),
                                          sheet = 2, na = "NA"))
  S <- as.matrix(sd2[, meta$protein_ColumnHeader, with = FALSE]); storage.mode(S) <- "numeric"
  rownames(S) <- sd2$Cre_GID; colnames(S) <- meta$SampleID
  S <- S[rowSums(!is.na(S)) >= 30, ]
  S <- sweep(S, 2, apply(S, 2, median, na.rm = TRUE)); S <- S - rowMeans(S, na.rm = TRUE)
  common <- intersect(rownames(Mch), rownames(S))
  assign_once <- function(R) setNames(rownames(R)[as.integer(clue::solve_LSAP(t(R) + 1, maximum = TRUE))], colnames(R))
  centre_within <- function(X, g) { for (k in unique(g)) X[, g == k] <- X[, g == k] - rowMeans(X[, g == k], na.rm = TRUE); X }
  plex_lab <- function(cols) sub("_.*", "", cols)
  R0 <- R <- cor(Mch[common, ], S[common, ], use = "pairwise.complete.obs")
  asg <- assign_once(R)
  for (it in 1:10) {   # the published values still carry the plex effect: centre within plex and reassign
    R <- cor(Mch[common, ], centre_within(S[common, ], plex_lab(asg[colnames(S)])), use = "pairwise.complete.obs")
    new <- assign_once(R); changed <- sum(new != asg[names(new)]); asg <- new
    cat(sprintf("Iteration %d: %d samples reassigned\n", it, changed))
    if (changed == 0) break
  }
  ref_cols <- setdiff(rownames(R), asg)
  design <- data.table(SampleID = names(asg), column = unname(asg))
  design[, r_initial := R0[cbind(column, SampleID)]][, r_assigned := R[cbind(column, SampleID)]]
  design[, r_second := sapply(seq_len(.N), function(i) max(R[column[i], colnames(R) != SampleID[i]]))]
  design <- rbind(design, data.table(SampleID = "REFERENCE", column = ref_cols), fill = TRUE)
  design[, `:=`(plex = as.integer(sub("plex(\\d+)_.*", "\\1", column)), channel = sub(".*_", "", column))]
  design[, column := NULL]
  setorder(design, plex, channel)
  fwrite(design, DESIGN_CACHE, sep = "\t")
  pheatmap::pheatmap(R[order(rownames(R)), order(colnames(R))], cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 6,
                     main = "Channels vs published samples (after centring within plex)",
                     filename = D("0_PSMs", "channel_correlation.png"), width = 11, height = 9)
  cat(sprintf("Assignment: median r %.2f (initial %.2f), median second best %.2f\n",
              median(design$r_assigned, na.rm = TRUE), median(design$r_initial, na.rm = TRUE), median(design$r_second, na.rm = TRUE)))
}
design[, channel := as.character(channel)]
ref <- design[SampleID == "REFERENCE", .(plex, ref_channel = channel)]
smp <- merge(design[SampleID != "REFERENCE", .(plex, channel, SampleID)],
             meta[, .(SampleID, Light = Population, Time = TimeCode, Replicate = ExperimentalReplicate)], by = "SampleID")
smp[, `:=`(Light = factor(Light, LIGHTS), Time = factor(Time, TIMES))]
smp[, State := factor(paste(Light, Time, sep = "_"), levels = as.vector(outer(LIGHTS, TIMES, paste, sep = "_")))]
setorder(smp, Light, Time, plex)
stopifnot(nrow(smp) == 45, nrow(ref) == 3, all(table(smp$State) == 3))
cat("References:", paste0("plex", ref$plex, " = ", ref$ref_channel, collapse = ", "), "\n")
cat("Balance (one replicate of each state per plex):", all(table(smp$State, smp$plex) == 1), "\n")

# =================================================================================
# 4. Pipeline A: ratio to the plex reference + limma
# =================================================================================
banner("4. Pipeline A: ratio to reference + limma")
tmt[, psm_id := .I]
long <- melt(tmt[, c("psm_id", "plex", "sequence", "genes", CHANNELS), with = FALSE],
             id.vars = c("psm_id", "plex", "sequence", "genes"), measure.vars = CHANNELS,
             variable.name = "channel", value.name = "int", variable.factor = FALSE)
long <- long[!is.na(int) & int > 0][, l2 := log2(int)]
long <- merge(long, ref, by = "plex")
refv <- long[channel == ref_channel, .(psm_id, l2_ref = l2)]
long <- merge(long[channel != ref_channel], refv, by = "psm_id")       # PSMs without reference are dropped
long[, ratio := l2 - l2_ref]            # cancels everything specific to that peptide in that plex
long <- merge(long, smp[, .(plex, channel, SampleID)], by = c("plex", "channel"))
long[, n_ch := .N, by = psm_id]
long <- long[n_ch >= MIN_CHANNELS]
long[, ratio := ratio - median(ratio), by = SampleID]                  # equal loading
protA <- long[, .(ratio = median(ratio), n_psm = .N), by = .(genes, SampleID)]
MA <- as.matrix(dcast(protA, genes ~ SampleID, value.var = "ratio"), rownames = "genes")[, smp$SampleID]
psm_n <- protA[, .(n_psm = sum(n_psm) / 45), by = genes]; psm_n <- setNames(psm_n$n_psm, psm_n$genes)
# without reference (only to show the plex effect)
long[, l2n := l2 - median(l2), by = SampleID]
Mraw <- as.matrix(dcast(long[, .(l2 = median(l2n)), by = .(genes, SampleID)], genes ~ SampleID, value.var = "l2"),
                  rownames = "genes")[, smp$SampleID]
rm(long, refv)

# missing values come in blocks of 15: whole plexes
n_plex <- apply(MA, 1, function(x) sum(tapply(!is.na(x), smp$plex, any)))
tab <- data.table(plexes = 3:1, proteins = as.integer(table(factor(n_plex, 3:1))))
cat("Proteins quantified in 3, 2 and 1 plexes:\n"); print(tab)
fwrite(tab, D("1_QC", "proteins_by_number_of_plexes.tsv"), sep = "\t")
all_states <- apply(MA, 1, function(x) all(tapply(!is.na(x), smp$State, any)))
MA <- MA[n_plex >= MIN_PLEX & all_states, ]
cat(sprintf("Proteins kept (>= %d plexes and every state): %d\n", MIN_PLEX, nrow(MA)))

pca_plot <- function(X, title, file) {
  X <- X[complete.cases(X), ]
  pc <- prcomp(t(X)); ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  d <- cbind(smp, PC1 = pc$x[smp$SampleID, 1], PC2 = pc$x[smp$SampleID, 2])
  g <- ggplot(d, aes(PC1, PC2, colour = Light, shape = factor(plex))) + geom_point(size = 3) +
    scale_colour_manual(values = LIGHT_COL) +
    labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]), shape = "Plex",
         title = title, subtitle = sprintf("%d proteins with no missing value", nrow(X)))
  ggsave(file, g, width = 6.5, height = 5)
  cat(sprintf("PCA %-32s PC1 %.1f%%, PC2 %.1f%%\n", title, ve[1], ve[2]))
}
pca_plot(Mraw[intersect(rownames(Mraw), rownames(MA)), ], "without reference normalization", D("1_QC", "PCA_1_without_reference.png"))
pca_plot(MA, "ratio to the plex reference", D("1_QC", "PCA_2_ratio_to_reference.png"))

# check against centring within plex (valid only because the design is balanced)
centre_plex <- function(X) { for (p in unique(smp$plex)) X[, smp$plex == p] <- X[, smp$plex == p] - rowMeans(X[, smp$plex == p], na.rm = TRUE); X }
Mc <- centre_plex(Mraw[rownames(MA), ]); Mr <- centre_plex(MA)
r_meth <- sapply(seq_len(nrow(MA)), function(i) suppressWarnings(cor(Mr[i, ], Mc[i, ], use = "pairwise.complete.obs")))
cat(sprintf("Ratio to reference vs centring within plex: median per-protein correlation %.2f\n", median(r_meth, na.rm = TRUE)))
fwrite(data.table(Gene = rownames(MA), MA), D("2_A_limma", "A_protein_log2ratio.tsv"), sep = "\t")

# limma: one mean per light x time state; plex as a random block
X <- model.matrix(~ 0 + State, data = smp); colnames(X) <- levels(smp$State)
dc <- duplicateCorrelation(MA, X, block = smp$plex)
cat(sprintf("Within-plex correlation left after normalization: %.3f\n", dc$consensus.correlation))
fitA <- lmFit(MA, X, block = smp$plex, correlation = dc$consensus.correlation)
w <- function(states) { v <- setNames(rep(0, ncol(X)), colnames(X)); v[states] <- 1 / length(states); v }
st <- function(l = LIGHTS, t = TIMES) as.vector(outer(l, t, paste, sep = "_"))
C_light <- cbind(HL_vs_LL = w(st("HL")) - w(st("LL")), HL_vs_ML = w(st("HL")) - w(st("ML")), ML_vs_LL = w(st("ML")) - w(st("LL")))
C_time  <- sapply(TIMES[-1], function(t) w(st(LIGHTS, t)) - w(st(LIGHTS, TIMES[1])))
C_int   <- do.call(cbind, lapply(c("HL", "ML"), function(l) sapply(TIMES[-1], function(t)
  (w(st(l, t)) - w(st("LL", t))) - (w(st(l, TIMES[1])) - w(st("LL", TIMES[1]))))))
eb <- function(C) eBayes(contrasts.fit(fitA, C), trend = TRUE, robust = TRUE)
fl <- eb(C_light); ft <- eb(C_time); fi <- eb(C_int)
png(D("1_QC", "A_limma_mean_variance.png"), 700, 500); plotSA(fl, main = "Pipeline A: mean-variance trend"); dev.off()

resA <- rbindlist(lapply(CONTRASTS, function(k) {
  tt <- topTable(fl, coef = k, number = Inf, sort.by = "none")
  data.table(Protein = rownames(tt), Contrast = k, log2FC = tt$logFC, SE = tt$logFC / tt$t, P = tt$P.Value, adj.P = tt$adj.P.Val)
}))
fwrite(resA, D("2_A_limma", "A_light_contrasts.tsv"), sep = "\t")
ftab <- function(f, name) { tt <- topTable(f, number = Inf, sort.by = "none")
  data.table(Protein = rownames(tt), Test = name, F = tt$F, P = tt$P.Value, adj.P = tt$adj.P.Val) }
timeA <- ftab(ft, "time_of_day"); intA <- ftab(fi, "light_x_time")
fwrite(rbind(timeA, intA), D("2_A_limma", "A_time_and_interaction_Ftests.tsv"), sep = "\t")
cat(sprintf("Time of day (F-test): %d proteins | light x time interaction: %d proteins (adj.P < %.2f)\n",
            sum(timeA$adj.P < FDR), sum(intA$adj.P < FDR), FDR))

top_t <- timeA[order(P)][1:12, Protein]
prof <- melt(data.table(Protein = top_t, MA[top_t, , drop = FALSE]), id.vars = "Protein", variable.name = "SampleID", value.name = "y")
prof <- merge(prof, smp[, .(SampleID, Light, Time)], by = "SampleID")[, Hour := HOURS[as.character(Time)]]
g <- ggplot(prof, aes(Hour, y, colour = Light)) +
  annotate("rect", xmin = 12, xmax = 15, ymin = -Inf, ymax = Inf, fill = "grey88") +
  annotate("rect", xmin = -3, xmax = 0, ymin = -Inf, ymax = Inf, fill = "grey88") +
  geom_point(size = 1, alpha = 0.6) + stat_summary(fun = mean, geom = "line", na.rm = TRUE) +
  facet_wrap(~ Protein, scales = "free_y", ncol = 4) + scale_colour_manual(values = LIGHT_COL) +
  scale_x_continuous(breaks = unname(HOURS), labels = c("-2", "+2", "+6", "+10", "+14")) +
  labs(title = "Pipeline A: strongest time-of-day effects (grey = dark)", x = "Hours after dawn", y = "log2 ratio to reference")
ggsave(D("2_A_limma", "A_profiles_top_time.png"), g, width = 11, height = 7)

# =================================================================================
# 5. Pipeline B: MSstatsTMT
# =================================================================================
banner("5. Pipeline B: MSstatsTMT")
MS_CACHE   <- D("3_B_MSstatsTMT", "B_MSstatsTMT.rds")
QUANT_CACHE <- D("3_B_MSstatsTMT", "B_proteinSummarization.rds")   # the slow step, kept on its own
if (FORCE_MSSTATS || !file.exists(MS_CACHE)) {
  if (FORCE_MSSTATS || !file.exists(QUANT_CACHE)) {
    # MSstatsTMT input: one row per feature (peptide + modifications + charge) and channel.
    # As the MSstatsTMT converters do, a feature seen in several fractions or scans of a
    # plex keeps only its most intense PSM; the plex is the "Mixture" and its "Run".
    feat <- copy(tmt)
    feat[, total := rowSums(.SD, na.rm = TRUE), .SDcols = CHANNELS]
    feat[, feature := paste(sequence, mods, sep = "|")]
    feat <- feat[feat[, .I[which.max(total)], by = .(plex, genes, feature, charge)]$V1]
    inp <- melt(feat[, c("plex", "genes", "feature", "charge", CHANNELS), with = FALSE],
                id.vars = c("plex", "genes", "feature", "charge"), measure.vars = CHANNELS,
                variable.name = "Channel", value.name = "Intensity", variable.factor = FALSE)
    rm(feat)
    ann <- rbind(smp[, .(plex, Channel = channel, BioReplicate = SampleID, Condition = as.character(State))],
                 ref[, .(plex, Channel = ref_channel, BioReplicate = "Norm", Condition = "Norm")])
    inp <- merge(inp, ann, by = c("plex", "Channel"))
    inp[, `:=`(ProteinName = genes, PeptideSequence = feature, Charge = charge,
               PSM = paste(feature, charge, sep = "_"), Mixture = paste0("plex", plex),
               TechRepMixture = 1L, Run = paste0("plex", plex))]
    inp <- inp[, .(ProteinName, PeptideSequence, Charge, PSM, Mixture, TechRepMixture, Run,
                   Channel, Condition, BioReplicate, Intensity)]
    cat(sprintf("MSstatsTMT input: %d rows (%d features x 16 channels)\n", nrow(inp), nrow(inp) / 16))
    t0 <- Sys.time()
    quant <- MSstatsTMT::proteinSummarization(as.data.frame(inp), method = "msstats", global_norm = TRUE,
                                              reference_norm = TRUE, remove_norm_channel = TRUE, MBimpute = TRUE,
                                              use_log_file = FALSE, verbose = FALSE)
    saveRDS(quant, QUANT_CACHE)
    rm(inp)
    cat(sprintf("proteinSummarization took %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  } else quant <- readRDS(QUANT_CACHE)
  cm <- t(C_light)          # rows = contrasts, columns = the 15 conditions (averages over the day)
  comp <- tryCatch(
    MSstatsTMT::groupComparisonTMT(quant, contrast.matrix = cm, moderated = TRUE, use_log_file = FALSE, verbose = FALSE),
    error = function(e) {    # very old limma versions lack an argument MSstatsTMT needs for moderation
      warning("Moderated test failed (", conditionMessage(e), "); using the unmoderated one. Update limma to fix it.")
      MSstatsTMT::groupComparisonTMT(quant, contrast.matrix = cm, moderated = FALSE, use_log_file = FALSE, verbose = FALSE)
    })
  saveRDS(list(protein = as.data.table(quant$ProteinLevelData), comparison = as.data.table(comp$ComparisonResult)), MS_CACHE)
} else cat("Using cached", MS_CACHE, "\n")
msB <- readRDS(MS_CACHE)

pB <- msB$protein[!is.na(Abundance)]
MB <- as.matrix(dcast(pB, Protein ~ BioReplicate, value.var = "Abundance", fun.aggregate = mean), rownames = "Protein")
MB <- MB[, intersect(smp$SampleID, colnames(MB))]
resB <- msB$comparison[, .(Protein = as.character(Protein), Contrast = as.character(Label), log2FC, SE, P = pvalue, adj.P = adj.pvalue, issue)]
resB <- resB[is.finite(log2FC) & !is.na(P)]
# Same rule as pipeline A: proteins quantified in >= MIN_PLEX plexes. A protein seen in a
# single plex has one value per condition: no residual variance, so MSstatsTMT gives it
# the prior variance only ("SingleMeasurePerCondition") and its p-value depends on the
# fold change alone (the curved band at -log10 P ~ 3 in the volcano).
plexB <- pB[, .(n_plex = uniqueN(Mixture)), by = Protein]
keepB <- plexB[n_plex >= MIN_PLEX, as.character(Protein)]
cat(sprintf("MSstatsTMT: %d proteins in >= %d plexes kept; %d single-plex proteins dropped\n",
            length(keepB), MIN_PLEX, uniqueN(resB$Protein) - length(intersect(unique(resB$Protein), keepB))))
resB <- resB[Protein %in% keepB & (is.na(issue) | issue != "SingleMeasurePerCondition")]
MB <- MB[intersect(rownames(MB), keepB), ]
fwrite(data.table(Protein = rownames(MB), MB), D("3_B_MSstatsTMT", "B_protein_abundance.tsv"), sep = "\t")
fwrite(resB, D("3_B_MSstatsTMT", "B_light_contrasts.tsv"), sep = "\t")
cat(sprintf("MSstatsTMT: %d proteins summarized, %d tested\n", nrow(MB), uniqueN(resB$Protein)))
pca_plot(MB, "MSstatsTMT protein abundance", D("1_QC", "PCA_3_MSstatsTMT.png"))

# =================================================================================
# 6. Relevance: fold change and Bayesian criteria
# =================================================================================
banner("6. Relevance")
# Bayesian model, one per protein (Stan), on the proteins significant in at least one
# light contrast. As in PC3, empirical priors are taken from ALL proteins:
#   - sigma: centred on the trend of log(SD) with abundance and number of PSMs
#   - light effects: Student-t(3, 0, s), s = robust spread of all estimated light effects
# Differences from PC3 (label-free, two groups):
#   - no dropout curve: in TMT values are missing in whole plexes because the
#     instrument did not pick the peptide, not because it was low: missing values are
#     simply left out (they carry no information about abundance);
#   - factorial design: y = time-of-day mean + light effect (LL = 0); the three light
#     contrasts are derived from the two light effects. The light x time interaction is
#     not modelled (it goes into the residual: conservative).
stan_code <- "
data {
  int<lower=1> N;  vector[N] y;  array[N] int<lower=1, upper=5> tt;  array[N] int<lower=1, upper=3> ll;
  real mu_loc;  real<lower=0> mu_scale;  real<lower=0> s_beta;  real logsig_loc;  real<lower=0> logsig_scale;
}
parameters { vector[5] mu;  real b_ML;  real b_HL;  real log_sigma; }
transformed parameters { vector[3] gamma = [0, b_ML, b_HL]'; }
model {
  mu ~ normal(mu_loc, mu_scale);
  b_ML ~ student_t(3, 0, s_beta);  b_HL ~ student_t(3, 0, s_beta);
  log_sigma ~ normal(logsig_loc, logsig_scale);
  y ~ normal(mu[tt] + gamma[ll], exp(log_sigma));
}
generated quantities { real HL_vs_LL = b_HL;  real HL_vs_ML = b_HL - b_ML;  real ML_vs_LL = b_ML; }
"
hyper <- function(M) {      # empirical priors from all proteins of one pipeline
  Xa <- model.matrix(~ Time + Light, data = smp[match(colnames(M), SampleID)])
  f <- lmFit(M, Xa)
  dfr <- f$df.residual; ok <- dfr >= 3 & is.finite(f$sigma) & f$sigma > 0
  d <- data.table(Protein = rownames(M)[ok], logsd = log(f$sigma[ok]), n = dfr[ok], mean = f$Amean[ok],
                  psm = log(pmax(1, psm_n[rownames(M)[ok]])))
  d[is.na(psm), psm := 0]
  d[, logsd_c := logsd - (digamma(n / 2) - log(n / 2)) / 2]        # bias of log(sample SD)
  trend <- lm(logsd_c ~ mean + psm, data = d)
  tau <- sqrt(max(var(residuals(trend)) - mean(trigamma(d$n / 2) / 4), 0.05^2))
  eff <- as.vector(f$coefficients[, c("LightML", "LightHL")])
  list(trend = trend, tau = tau, s_beta = max(mad(eff, na.rm = TRUE), 0.05),
       mu_loc = median(M, na.rm = TRUE), mu_scale = 3 * sd(M, na.rm = TRUE))
}
bayes_pipeline <- function(M, res, name) {
  cache <- D("4_Relevance", paste0("Bayes_", name, ".rds"))
  # Only proteins that could pass: significant AND |log2FC| + 3 SE >= LFC in some contrast.
  # If even the upper end of a 99.7% interval stays below LFC, P(|log2FC| > LFC) cannot reach
  # P_BAYES (the shrinkage prior only pulls it lower), so the fit would be wasted time.
  n_sig <- res[adj.P < FDR, uniqueN(Protein)]
  sig <- res[adj.P < FDR & abs(log2FC) + 3 * SE >= LFC, unique(Protein)]
  sig <- intersect(sig, rownames(M))
  if (!FORCE_BAYES && file.exists(cache)) { out <- readRDS(cache); todo <- setdiff(sig, out$Protein) } else { out <- NULL; todo <- sig }
  cat(sprintf("Pipeline %s: %d significant proteins, %d could reach the threshold, %d still to fit\n",
              name, n_sig, length(sig), length(todo)))
  if (length(todo)) {
    if (!exists("stan_model_obj")) stan_model_obj <<- rstan::stan_model(model_code = stan_code)
    h <- hyper(M)
    info <- smp[match(colnames(M), SampleID)]
    fits <- list(); k <- 0
    for (p in todo) {
      y <- M[p, ]; o <- !is.na(y)
      dat <- list(N = sum(o), y = as.array(unname(y[o])), tt = as.array(as.integer(info$Time[o])),
                  ll = as.array(as.integer(info$Light[o])), mu_loc = h$mu_loc, mu_scale = h$mu_scale, s_beta = h$s_beta,
                  logsig_loc = unname(predict(h$trend, data.frame(mean = mean(y, na.rm = TRUE),
                                                                   psm = if (is.na(psm_n[p])) 0 else log(max(1, psm_n[p]))))),
                  logsig_scale = h$tau)
      invisible(capture.output(fit <- rstan::sampling(stan_model_obj, data = dat, chains = 4, iter = 2000,
                                                      seed = SEED, refresh = 0, cores = 1)))
      dr <- rstan::extract(fit, CONTRASTS)
      rh <- max(rstan::summary(fit, pars = c("mu", "b_ML", "b_HL", "log_sigma"))$summary[, "Rhat"])
      fits[[p]] <- rbindlist(lapply(CONTRASTS, function(k2) { b <- dr[[k2]]
        data.table(Protein = p, Contrast = k2, post_mean = mean(b), low95 = quantile(b, 0.025), high95 = quantile(b, 0.975),
                   P_up = mean(b > LFC), P_down = mean(b < -LFC), Rhat = rh) }))
      k <- k + 1
      if (k %% 100 == 0 || k == length(todo)) {
        out <- rbind(out, rbindlist(fits)); fits <- list(); saveRDS(out, cache)
        cat(sprintf("  %d / %d proteins\n", k, length(todo)))
      }
    }
  }
  out[Protein %in% sig]
}
bA <- bayes_pipeline(MA, resA, "A")
bB <- bayes_pipeline(MB, resB, "B")
if (any(c(bA$Rhat, bB$Rhat) > 1.01)) warning("Some Bayesian fits have Rhat > 1.01")

classify <- function(res, bay, pipeline) {
  r <- merge(res, bay[, .(Protein, Contrast, post_mean, low95, high95, P_up, P_down)], by = c("Protein", "Contrast"), all.x = TRUE)
  r[, Pipeline := pipeline]
  r[, FC_relevant := adj.P < FDR & abs(log2FC) >= LFC]
  r[, Bayes_relevant := adj.P < FDR & !is.na(P_up) & pmax(P_up, P_down) >= P_BAYES]
  r[, Direction := fifelse(log2FC > 0, "up", "down")]
  r
}
REL <- rbind(classify(resA, bA, "A_limma"), classify(resB, bB, "B_MSstatsTMT"), fill = TRUE)
fwrite(REL, D("4_Relevance", "relevance_all.tsv"), sep = "\t")
summ <- REL[, .(tested = .N, significant = sum(adj.P < FDR),
                FC_up = sum(FC_relevant & Direction == "up"), FC_down = sum(FC_relevant & Direction == "down"),
                Bayes_up = sum(Bayes_relevant & Direction == "up"), Bayes_down = sum(Bayes_relevant & Direction == "down")),
            by = .(Pipeline, Contrast)]
cat(sprintf("Relevance (adj.P < %.2f; FC: |log2FC| >= %.2f; Bayes: P(|log2FC| > %.2f) >= %.2f):\n", FDR, LFC, LFC, P_BAYES))
print(summ)
fwrite(summ, D("4_Relevance", "relevance_summary.tsv"), sep = "\t")

g <- ggplot(REL, aes(log2FC, -log10(P))) + geom_point(colour = "grey75", size = 0.5) +
  geom_point(data = REL[FC_relevant == TRUE], colour = "#B5541B", size = 0.7) +
  geom_point(data = REL[Bayes_relevant == TRUE], shape = 1, colour = "#14213D", size = 1.2) +
  geom_vline(xintercept = c(-LFC, LFC), linetype = 2, colour = "grey50") +
  facet_grid(Pipeline ~ Contrast) +
  labs(title = "Light acclimation, averaged over the day",
       subtitle = "orange = FC criterion; dark circle = Bayesian criterion", x = "log2 fold change", y = "-log10 P")
ggsave(D("4_Relevance", "volcanoes.png"), g, width = 11, height = 7)

# =================================================================================
# 7. Comparison of pipelines and criteria
# =================================================================================
banner("7. Comparison")
both <- merge(resA[, .(Protein, Contrast, A = log2FC)], resB[, .(Protein, Contrast, B = log2FC)], by = c("Protein", "Contrast"))
corr <- both[, .(proteins = .N, r_log2FC = round(cor(A, B), 3),
                 slope_B_on_A = round(coef(lm(B ~ A))[2], 2)), by = Contrast]
cat("log2FC, pipeline A vs B:\n"); print(corr)
g <- ggplot(both, aes(A, B)) + geom_point(size = 0.4, alpha = 0.4) + geom_abline(linetype = 2, colour = "#B5541B") +
  facet_wrap(~ Contrast) + labs(x = "log2FC, pipeline A (ratio to reference + limma)", y = "log2FC, pipeline B (MSstatsTMT)")
ggsave(D("5_Comparison", "log2FC_A_vs_B.png"), g, width = 11, height = 4)

jac <- rbindlist(lapply(CONTRASTS, function(k) {
  L <- list(`A FC` = REL[Pipeline == "A_limma" & Contrast == k & FC_relevant, paste(Protein, Direction)],
            `A Bayes` = REL[Pipeline == "A_limma" & Contrast == k & Bayes_relevant, paste(Protein, Direction)],
            `B FC` = REL[Pipeline == "B_MSstatsTMT" & Contrast == k & FC_relevant, paste(Protein, Direction)],
            `B Bayes` = REL[Pipeline == "B_MSstatsTMT" & Contrast == k & Bayes_relevant, paste(Protein, Direction)])
  cmb <- combn(names(L), 2)
  data.table(Contrast = k, list1 = cmb[1, ], list2 = cmb[2, ],
             n1 = sapply(cmb[1, ], function(a) length(L[[a]])), n2 = sapply(cmb[2, ], function(a) length(L[[a]])),
             shared = apply(cmb, 2, function(x) length(intersect(L[[x[1]]], L[[x[2]]]))),
             Jaccard = apply(cmb, 2, function(x) { u <- length(union(L[[x[1]]], L[[x[2]]]))
               if (u == 0) NA_real_ else round(length(intersect(L[[x[1]]], L[[x[2]]])) / u, 2) }))
}))
cat("\nOverlap of the relevant lists (same protein, same direction):\n"); print(jac)
fwrite(jac, D("5_Comparison", "list_overlap.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), D("sessionInfo.txt"))
cat("\nDone. Results in", file.path(BASE, "results_TMT"), "\n")
