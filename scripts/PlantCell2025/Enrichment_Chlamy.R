# =============================================================================
# Enrichment_Chlamy.R
# Functional enrichment of the light-acclimation results of Chlamydomonas
# (Dupuis, Ojeda et al. 2025, Plant Cell 37: koaf086), at two levels:
#   RNA     (GEO GSE275433)     pipelines: limma-voom, DESeq2
#   protein (MassIVE MSV000095622) pipelines: A = ratio to reference + limma, B = MSstatsTMT
# each with two relevance criteria (FC and Bayesian) and two annotation sources
# (eggNOG-mapper and the Phytozome v6.1 annotation), as in the enrichment practical (PC4):
#
#   4. ORA of every list and GSEA of every pipeline, same settings
#   5. Descriptive comparison: shared terms (Jaccard), NES of one choice against another
#   6. Differential enrichment: functions where two pipelines disagree systematically
#   7. Resampling (jackknife): is the difference larger than the effect of the samples?
#
# Sections
#   1. Packages
#   2. Data: statistics and lists of the RNA and TMT blocks
#   3. Annotation: eggNOG-mapper and Phytozome, mapped to genes
#   4. ORA and GSEA
#   5. Comparison of pipelines, criteria and annotation sources
#   6. Differential enrichment
#   7. Resampling                                                          [cached]
#   8. Files for Cytoscape / ClueGO
#
# Everything that may need editing is in section 0. The resampling is cached: if its
# output exists it is reused (FORCE_RESAMPLING = TRUE to recompute).
#
# Files expected in BASE/InputData:
#   Annotation/          CreinhardtiiCC_4532_707_v6.1.annotation_info.txt   (Phytozome)
#                        *.emapper.annotations                               (eggNOG-mapper on the
#                                                                             v6.1 primary proteins)
#   TranscriptomicData/  relevance_all.tsv, sample_design.tsv,
#                        GSE275433_counts_reconstructed.tsv(.gz)            (from RNAseq_Chlamy.R)
#   ProteomicData/       4_Relevance/relevance_all.tsv, 2_A_limma/A_protein_log2ratio.tsv,
#                        3_B_MSstatsTMT/B_MSstatsTMT.rds                     (from Proteomic_Chlamy_TMT.R)
# =============================================================================

# ---- 0. Settings (edit here) ---------------------------------------------------
BASE  <- Sys.getenv("ENRICH_BASE", path.expand("~/BMSi/PlantCell2025/Enrichment"))
INPUT <- file.path(BASE, "InputData")

CONTRASTS <- c("HL_vs_LL", "HL_vs_ML", "ML_vs_LL")
FDR       <- 0.05     # adjusted p-value (terms)
MIN_SIZE  <- 10       # smallest gene set tested
MAX_SIZE  <- 500      # largest gene set tested (very general terms are left out)
ORA_MIN   <- 3        # an ORA term needs >= this many genes of the list
DELTA_NES <- 1        # smallest NES difference between two pipelines that counts as relevant
DELTA_ER  <- 1        # ORA: smallest difference of log2 fold enrichment (2-fold) that counts as relevant
LFC       <- c(RNA = 1, Protein = log2(1.5))   # FC criterion of each block, to rebuild the lists in section 7b
SEED      <- 2026

# Resampling (section 7): the unit resampled must be an INDEPENDENT unit.
#   RNA: the 45 libraries are independent (replicate labels are not batches).
#        "sample":    jackknife, leave out one library at a time (45 resamples, 44 df)
#        "replicate": jackknife, leave out BR1, BR2 or BR3 in every state at once (3 resamples, 2 df)
#        "bootstrap": stratified bootstrap, draw 3 libraries WITH replacement within each of the
#                     15 states (BOOT_B resamples; df = 45 libraries - 15 states = 30). The SE is
#                     the SD over resamples x sqrt(3/2), the usual correction for 3 per stratum.
#   TMT: the plex is a batch, so the unit is the plex (jackknife, 3 resamples, 2 df). There is no
#        TMT bootstrap: a plex drawn twice breaks the mixture model, and 3 plexes give 10 distinct draws.
JACK_RNA      <- "bootstrap"  # "sample", "replicate" or "bootstrap"
BOOT_B        <- 50           # bootstrap resamples (each one refits limma-voom and DESeq2: ~2.5 min)
JACK_MODERATE <- TRUE         # also report resampling SEs moderated across terms (limma::squeezeVar)
FORCE_GSEA       <- FALSE     # recompute the GSEA of section 4 (~25 min)
FORCE_RESAMPLING <- FALSE     # recompute the resampling (RNA by sample: ~1.5 h, mostly DESeq2)

# ---- 1. Packages ------------------------------------------------------------------
cran <- c("data.table", "ggplot2", "pheatmap")
bioc <- c("limma", "edgeR", "DESeq2", "fgsea", "GO.db", "AnnotationDbi", "BiocParallel", "MSstatsTMT")
for (p in cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(limma); library(edgeR) })
theme_set(theme_bw(base_size = 11))
set.seed(SEED)

D <- function(...) { p <- file.path(BASE, "results_Enrichment", ...); dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE); p }
banner <- function(x) cat("\n==== ", x, " ", strrep("=", max(3, 70 - nchar(x))), "\n", sep = "")
inp <- function(...) { f <- file.path(INPUT, ...); if (!file.exists(f) && file.exists(paste0(f, ".gz"))) f <- paste0(f, ".gz")
  if (!file.exists(f)) stop("Missing input file: ", f); f }
LIGHTS <- c("LL", "ML", "HL")

# =================================================================================
# 2. Data: statistics and lists of the RNA and TMT blocks
# =================================================================================
banner("2. Data")
# Signed z-score of every gene in every pipeline: the ranking used by GSEA
# (z = sign(log2FC) x the normal quantile of the two-sided p-value).
zscore <- function(lfc, p) sign(lfc) * qnorm(pmax(p, 1e-300) / 2, lower.tail = FALSE)

relR <- fread(inp("TranscriptomicData", "relevance_all.tsv"))
relP <- fread(inp("ProteomicData", "4_Relevance", "relevance_all.tsv"))
relR[, `:=`(Omic = "RNA", id = Gene, Pipeline = fifelse(Pipeline == "A_limma_voom", "RNA_limma", "RNA_DESeq2"))]
relP[, `:=`(Omic = "Protein", id = Protein, Pipeline = fifelse(Pipeline == "A_limma", "TMT_A", "TMT_B"))]
cols <- c("Omic", "Pipeline", "Contrast", "id", "log2FC", "P", "adj.P", "FC_relevant", "Bayes_relevant", "Direction")
rel <- rbind(relR[, ..cols], relP[, ..cols])
rel <- rel[Contrast %in% CONTRASTS & is.finite(log2FC) & is.finite(P)]
rel[, z := zscore(log2FC, P)]
PIPES <- c("RNA_limma", "RNA_DESeq2", "TMT_A", "TMT_B")
omic_of <- function(pipe) if (startsWith(pipe, "RNA")) "RNA" else "Protein"

# Lists: 2 pipelines x 2 criteria per omic, up and down, per contrast. Universe = genes tested.
lists <- list()
for (pp in PIPES) for (crit in c("FC", "Bayes")) for (k in CONTRASTS) {
  r <- rel[Pipeline == pp & Contrast == k]; sel <- r[[paste0(crit, "_relevant")]] %in% TRUE
  lists[[paste(pp, crit, k, sep = "|")]] <- list(pipeline = pp, criterion = crit, contrast = k, universe = unique(r$id),
                                                  up = r[sel & Direction == "up", unique(id)], down = r[sel & Direction == "down", unique(id)])
}
lsum <- rbindlist(lapply(lists, function(l) data.table(Pipeline = l$pipeline, Criterion = l$criterion, Contrast = l$contrast,
                                                     Universe = length(l$universe), Up = length(l$up), Down = length(l$down))))
fwrite(lsum, D("1_Data", "lists_summary.tsv"), sep = "\t")
cat("Lists (genes/proteins):\n"); print(dcast(lsum, Pipeline + Criterion ~ Contrast, value.var = c("Up", "Down")))

# =================================================================================
# 3. Annotation: eggNOG-mapper and Phytozome
# =================================================================================
banner("3. Annotation")
# Two sources, both mapped to the v6.1 gene (Cre01.g000050_4532.1.p -> Cre01.g000050):
#   eggNOG-mapper  transfers annotation from orthologous groups: GO, KEGG pathways, Pfam
#   Phytozome      the genome's own annotation (InterPro2GO, Pfam, PANTHER families)
# GO terms are propagated to their ancestors in both, so that a gene annotated to a specific
# term also counts for the more general ones (and both sources are treated alike).
gene_of <- function(x) sub("_4532.*$", "", x)
emf <- list.files(file.path(INPUT, "Annotation"), pattern = "emapper\\.annotations$", full.names = TRUE)
if (!length(emf)) stop("No *.emapper.annotations file in ", file.path(INPUT, "Annotation"))
em <- fread(emf[1], skip = "#query", sep = "\t", quote = "", na.strings = c("-", ""), fill = TRUE)
setnames(em, "#query", "query"); em <- em[!startsWith(query, "##")][, id := gene_of(query)]
pz <- fread(inp("Annotation", "CreinhardtiiCC_4532_707_v6.1.annotation_info.txt"), sep = "\t", quote = "", na.strings = "")
setnames(pz, 1, "pacId"); pz[, id := gene_of(locusName)]

split_col <- function(d, col, pattern) {   # one row per (gene, item)
  d <- d[!is.na(get(col)), .(id, v = as.character(get(col)))]
  unique(d[, .(item = trimws(unlist(strsplit(v, pattern)))), by = id][item != ""])
}
go_info <- function(ids) suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = unique(ids), columns = c("ONTOLOGY", "TERM"), keytype = "GOID"))
anc <- list(BP = as.list(GO.db::GOBPANCESTOR), MF = as.list(GO.db::GOMFANCESTOR), CC = as.list(GO.db::GOCCANCESTOR))
go_sets <- function(g) {   # g: id, item (GO id) -> GO_BP, GO_MF, GO_CC, propagated to ancestors
  info <- go_info(g$item); g[, Ont := info$ONTOLOGY[match(item, info$GOID)]]; g <- g[!is.na(Ont)]
  out <- list()
  for (o in c("BP", "MF", "CC")) {
    x <- g[Ont == o]
    a <- x[, .(Term = unique(c(item, unlist(anc[[o]][item])))), by = id][Term != "all"]
    nm <- go_info(a$Term); a[, Name := nm$TERM[match(Term, nm$GOID)]]
    out[[paste0("GO_", o)]] <- unique(a[, .(Term, id, Name)])
  }
  out
}
kegg_names <- tryCatch({
  x <- fread("https://rest.kegg.jp/list/pathway", header = FALSE, sep = "\t", col.names = c("Term", "Name"))
  x[, Term := sub("^path:", "", Term)]
}, error = function(e) {   # offline: the KEGG map list shipped with clusterProfiler
  e2 <- new.env(); utils::data("kegg_category", package = "clusterProfiler", envir = e2)
  k <- as.data.table(get("kegg_category", e2)); data.table(Term = paste0("map", k$id), Name = k$name)
})
kegg_names[, Name := trimws(sub("New!$", "", Name))]
sets <- list()
s <- go_sets(split_col(em, "GOs", ","))
kg <- split_col(em, "KEGG_Pathway", ",")
kg[grepl("^[0-9]{5}$", item), item := paste0("map", item)]
s$KEGG <- merge(kg[startsWith(item, "map"), .(Term = item, id)], kegg_names, by = "Term", all.x = TRUE)[, Name := fifelse(is.na(Name), Term, Name)]
s$Pfam <- split_col(em, "PFAMs", ",")[, item := sub("_[0-9]+_[0-9]+$", "", item)][, .(Term = item, id, Name = item)]
sets$eggNOG <- lapply(s, unique)
s <- go_sets(split_col(pz, "GO", "[ ,]"))
s$Pfam    <- split_col(pz, "Pfam", "[ ,]")[, .(Term = item, id, Name = item)]
s$PANTHER <- split_col(pz, "Panther", "[ ,]")[, item := sub(":.*$", "", item)][, .(Term = item, id, Name = item)]
sets$Phytozome <- lapply(s, unique)
rm(anc, em, pz, s, kg); invisible(gc())   # the GO ancestor tables are large

universe_all <- unique(rel$id)
cov <- rbindlist(lapply(names(sets), function(src) rbindlist(lapply(names(sets[[src]]), function(a) {
  x <- sets[[src]][[a]]
  data.table(Source = src, Annotation = a, Terms = uniqueN(x$Term), Genes = uniqueN(x$id),
             RNA_universe_pct = round(100 * mean(unique(rel[Omic == "RNA", id]) %in% x$id), 1),
             Protein_universe_pct = round(100 * mean(unique(rel[Omic == "Protein", id]) %in% x$id), 1))
}))))
fwrite(cov, D("2_Annotation", "annotation_coverage.tsv"), sep = "\t")
cat("Annotation coverage (% of the genes / proteins tested):\n"); print(cov)

# =================================================================================
# 4. ORA and GSEA
# =================================================================================
banner("4. ORA and GSEA")
ora <- function(ids, universe, s) {   # hypergeometric test; universe = annotated genes tested in the pipeline
  u <- intersect(universe, s$id); x <- intersect(ids, u); if (length(x) < ORA_MIN) return(NULL)
  t <- s[id %in% u, .(K = uniqueN(id), k = uniqueN(id[id %in% x]), genes = paste(sort(unique(id[id %in% x])), collapse = "/")), by = .(Term, Name)]
  t <- t[K >= MIN_SIZE & K <= MAX_SIZE]
  if (!nrow(t)) return(NULL)
  t[, `:=`(n = length(x), N = length(u))][, pvalue := phyper(k - 1, K, N - K, n, lower.tail = FALSE)][, padj := p.adjust(pvalue, "BH")]
  t[k >= ORA_MIN & padj < FDR][order(padj)]
}
prep_sets <- function(s, genes) {   # list of gene sets restricted to the ranked genes, size-filtered
  pw <- split(s$id, s$Term); pw <- lapply(pw, function(v) intersect(unique(v), genes))
  pw[lengths(pw) >= MIN_SIZE & lengths(pw) <= MAX_SIZE]
}
gsea <- function(z, s, seed = SEED) {   # fgsea on the signed z-scores
  z <- z[is.finite(z) & !duplicated(names(z))]; z <- sort(z, decreasing = TRUE)
  pw <- prep_sets(s, names(z)); set.seed(seed)
  r <- suppressWarnings(fgsea::fgseaMultilevel(pw, z, minSize = MIN_SIZE, maxSize = MAX_SIZE, eps = 0, BPPARAM = BiocParallel::SerialParam()))
  r <- as.data.table(r)[is.finite(NES), .(Term = pathway, size, ES, NES, pvalue = pval, padj)]
  merge(r, unique(s[, .(Term, Name)]), by = "Term", all.x = TRUE)
}
ANN <- unlist(lapply(names(sets), function(src) paste(src, names(sets[[src]]), sep = ":")))
get_set <- function(a) { p <- strsplit(a, ":")[[1]]; sets[[p[1]]][[p[2]]] }

ora_all <- rbindlist(lapply(names(lists), function(ln) { l <- lists[[ln]]
  rbindlist(lapply(c("up", "down"), function(dirn) rbindlist(lapply(ANN, function(a) {
    r <- ora(l[[dirn]], l$universe, get_set(a))
    if (!is.null(r) && nrow(r)) r[, `:=`(Pipeline = l$pipeline, Criterion = l$criterion, Contrast = l$contrast, Direction = dirn, Annotation = a)]
  }))))
}), fill = TRUE)
fwrite(ora_all[, .(Pipeline, Criterion, Contrast, Direction, Annotation, Term, Name, Ratio = paste0(k, "/", n), Background = paste0(K, "/", N), padj, genes)],
       D("3_Enrichment", "ORA_all.tsv"), sep = "\t")
GSEA_CACHE <- D("3_Enrichment", "GSEA_all.rds")
if (FORCE_GSEA || !file.exists(GSEA_CACHE)) {
  t0 <- Sys.time()
  gsea_all <- rbindlist(lapply(PIPES, function(pp) rbindlist(lapply(CONTRASTS, function(k) {
    r <- rel[Pipeline == pp & Contrast == k]; z <- setNames(r$z, r$id)
    rbindlist(lapply(ANN, function(a) gsea(z, get_set(a))[, `:=`(Pipeline = pp, Contrast = k, Annotation = a)]))
  }))))
  saveRDS(gsea_all, GSEA_CACHE)
  cat(sprintf("GSEA of %d pipeline x contrast x annotation combinations: %.1f min\n",
              uniqueN(gsea_all[, .(Pipeline, Contrast, Annotation)]), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
} else { gsea_all <- readRDS(GSEA_CACHE); cat("Using cached", GSEA_CACHE, "\n") }
gsea_all[, Name := trimws(sub("New!$", "", Name))]
fwrite(gsea_all[order(padj)], D("3_Enrichment", "GSEA_all.tsv"), sep = "\t")
ora_sum <- ora_all[, .N, by = .(Pipeline, Criterion, Contrast, Annotation)]
gsea_sum <- gsea_all[padj < FDR, .N, by = .(Pipeline, Contrast, Annotation)]
fwrite(ora_sum, D("3_Enrichment", "ORA_summary.tsv"), sep = "\t"); fwrite(gsea_sum, D("3_Enrichment", "GSEA_summary.tsv"), sep = "\t")
cat("ORA, significant terms, HL vs LL (up + down):\n")
print(dcast(ora_sum[Contrast == "HL_vs_LL"], Annotation ~ Pipeline + Criterion, value.var = "N", fill = 0))
cat("GSEA, significant terms, HL vs LL:\n")
print(dcast(gsea_sum[Contrast == "HL_vs_LL"], Annotation ~ Pipeline, value.var = "N", fill = 0))

top <- gsea_all[Contrast == "HL_vs_LL" & Annotation %in% c("eggNOG:GO_BP", "eggNOG:KEGG") & padj < FDR]
top <- top[order(-abs(NES)), head(.SD, 12), by = .(Pipeline, Annotation)]
if (nrow(top)) {
  g <- ggplot(top, aes(NES, reorder(substr(Name, 1, 55), NES), fill = NES > 0)) + geom_col() +
    facet_wrap(~ Pipeline + Annotation, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c(`TRUE` = "#B5541B", `FALSE` = "#2A6F97"), labels = c(`TRUE` = "higher in HL", `FALSE` = "lower in HL"), name = NULL) +
    labs(x = "NES", y = NULL, title = "GSEA, high vs low light: strongest enriched functions (eggNOG)")
  ggsave(D("3_Enrichment", "GSEA_top_HL_vs_LL.png"), g, width = 13, height = 14)
}

# =================================================================================
# 5. Comparison of pipelines, criteria and annotation sources
# =================================================================================
banner("5. Comparison")
jac <- function(a, b) { u <- length(union(a, b)); if (u == 0) NA_real_ else length(intersect(a, b)) / u }
# 5.1 ORA: Jaccard of the significant terms of the four lists of each omic (same contrast, same source)
for (ctr in CONTRASTS) for (src in names(sets)) {   # not "k": ora_all has a column k
  ln <- c(outer(c("RNA_limma", "RNA_DESeq2", "TMT_A", "TMT_B"), c("FC", "Bayes"), paste, sep = "_"))
  terms <- lapply(ln, function(x) ora_all[paste(Pipeline, Criterion, sep = "_") == x & Contrast == ctr & startsWith(Annotation, src), paste(Direction, Annotation, Term)])
  names(terms) <- ln
  J <- outer(ln, ln, Vectorize(function(a, b) jac(terms[[a]], terms[[b]]))); dimnames(J) <- list(ln, ln)
  fwrite(data.table(List = ln, round(J, 2)), D("4_Comparison", sprintf("Jaccard_ORA_%s_%s.tsv", ctr, src)), sep = "\t")
  if (ctr == "HL_vs_LL" && sum(is.finite(J)) > 1 && length(unique(J[is.finite(J)])) > 1)
    pheatmap::pheatmap(J, display_numbers = TRUE, number_format = "%.2f", cluster_rows = FALSE, cluster_cols = FALSE, na_col = "grey90",
                       main = paste("Jaccard of significant ORA terms, HL vs LL,", src), filename = D("4_Comparison", sprintf("Jaccard_ORA_HL_vs_LL_%s.png", src)), width = 8, height = 7)
}
# 5.2 GSEA: NES of pairs of choices (all terms tested in both)
PAIRS <- list(
  "RNA: limma-voom vs DESeq2"        = list(a = c("RNA_limma", "eggNOG"), b = c("RNA_DESeq2", "eggNOG")),
  "TMT: A (limma) vs B (MSstatsTMT)" = list(a = c("TMT_A", "eggNOG"),     b = c("TMT_B", "eggNOG")),
  "RNA limma-voom: eggNOG vs Phytozome" = list(a = c("RNA_limma", "eggNOG"), b = c("RNA_limma", "Phytozome")),
  "TMT A: eggNOG vs Phytozome"       = list(a = c("TMT_A", "eggNOG"),     b = c("TMT_A", "Phytozome")))
nes_pairs <- rbindlist(lapply(names(PAIRS), function(pn) { a <- PAIRS[[pn]]$a; b <- PAIRS[[pn]]$b
  ga <- gsea_all[Pipeline == a[1] & startsWith(Annotation, a[2])][, Cat := sub(".*:", "", Annotation)]
  gb <- gsea_all[Pipeline == b[1] & startsWith(Annotation, b[2])][, Cat := sub(".*:", "", Annotation)]
  m <- merge(ga[, .(Contrast, Cat, Term, Name, NES_1 = NES, padj_1 = padj)], gb[, .(Contrast, Cat, Term, NES_2 = NES, padj_2 = padj)], by = c("Contrast", "Cat", "Term"))
  m[, Pair := pn]
}))
nes_pairs[, Significant := fifelse(padj_1 < FDR & padj_2 < FDR, "both", fifelse(padj_1 < FDR, "first only", fifelse(padj_2 < FDR, "second only", "neither")))]
fwrite(nes_pairs, D("4_Comparison", "NES_pairs.tsv"), sep = "\t")
nes_cor <- nes_pairs[, .(Terms = .N, r_NES = round(cor(NES_1, NES_2, use = "complete.obs"), 2), Enriched_any = sum(Significant != "neither"),
                         Both = sum(Significant == "both"), First_only = sum(Significant == "first only"), Second_only = sum(Significant == "second only"),
                         Same_direction_pct = round(100 * mean(sign(NES_1[Significant != "neither"]) == sign(NES_2[Significant != "neither"])), 1)),
                     by = .(Pair, Contrast)]
fwrite(nes_cor, D("4_Comparison", "NES_pairs_summary.tsv"), sep = "\t")
cat("GSEA, pairs of choices (terms tested in both):\n"); print(nes_cor[Contrast == "HL_vs_LL"])
g <- ggplot(nes_pairs[Contrast == "HL_vs_LL"][order(Significant != "neither")], aes(NES_1, NES_2, colour = Significant)) +
  geom_point(size = 0.7, alpha = 0.6) + geom_abline(linetype = 2, colour = "grey50") + facet_wrap(~ Pair, ncol = 2) +
  scale_colour_manual(values = c(both = "#14213D", `first only` = "#B5541B", `second only` = "#2A6F97", neither = "grey80")) +
  labs(x = "NES, first choice", y = "NES, second choice", colour = "Enriched in", title = "High vs low light: NES of the same terms under two choices")
ggsave(D("4_Comparison", "NES_pairs_HL_vs_LL.png"), g, width = 11, height = 9)

# =================================================================================
# 6. Differential enrichment: where do two pipelines disagree?
# =================================================================================
banner("6. Differential enrichment")
# For the genes tested in both pipelines, d = z(first) - z(second); GSEA on d finds functions whose
# genes systematically score higher (NES > 0) or lower (NES < 0) in the first pipeline.
# Caution: genes are permuted, so correlation between the genes of a function makes this test
# optimistic; section 7 answers the same question by resampling the samples.
DPAIRS <- list("RNA: limma-voom vs DESeq2" = c("RNA_limma", "RNA_DESeq2"), "TMT: A vs B" = c("TMT_A", "TMT_B"))
diff_all <- rbindlist(lapply(names(DPAIRS), function(pn) rbindlist(lapply(CONTRASTS, function(k) {
  a <- DPAIRS[[pn]][1]; b <- DPAIRS[[pn]][2]
  m <- merge(rel[Pipeline == a & Contrast == k, .(id, za = z)], rel[Pipeline == b & Contrast == k, .(id, zb = z)], by = "id")
  d <- setNames(m$za - m$zb, m$id)
  rbindlist(lapply(grep("^eggNOG", ANN, value = TRUE), function(an) gsea(d, get_set(an))[, `:=`(Pair = pn, Contrast = k, Annotation = an, Genes = nrow(m))]))
}))))
fwrite(diff_all[order(padj)], D("5_Differential", "Differential_enrichment.tsv"), sep = "\t")
diff_sum <- diff_all[, .(Terms_tested = .N, Discordant = sum(padj < FDR, na.rm = TRUE), Higher_in_first = sum(padj < FDR & NES > 0, na.rm = TRUE),
                         Higher_in_second = sum(padj < FDR & NES < 0, na.rm = TRUE)), by = .(Pair, Contrast)]
fwrite(diff_sum, D("5_Differential", "Differential_summary.tsv"), sep = "\t")
cat("Differential enrichment (eggNOG annotation, gene permutation):\n"); print(diff_sum)

# =================================================================================
# 7. Resampling: jackknife over the independent units
# =================================================================================
banner("7. Resampling")
# Each pipeline is rerun without one independent unit at a time (jackknife). For each function the
# difference of NES between the two pipelines is recomputed every time; its spread is the
# uncertainty that comes only from the choice of samples:
#     SE_jack = sqrt((n - 1) / n * sum((d_i - mean(d))^2)),   t = d_full / SE_jack,   df = n - 1
# A "pipeline effect" needs: adjusted p < FDR, |d| >= DELTA_NES and the function enriched in at
# least one pipeline (significance + relevance, as for genes).
# Two details:
#   - NES of a resample = its ES divided by the normalization constant (ES/NES) of the full data,
#     so that the Monte Carlo noise of fgsea's normalization does not enter the jackknife
#     (the jackknife multiplies every deviation by n - 1).
#   - With few units the SE itself is very uncertain (relative error ~ 1/sqrt(2 df)): 50% with 2 df.
#     JACK_MODERATE squeezes the SEs of all terms towards their common value (empirical Bayes, as
#     limma does for gene variances), which adds "borrowed" degrees of freedom.
counts_file <- inp("TranscriptomicData", "GSE275433_counts_reconstructed.tsv")
cn <- fread(counts_file); counts <- as.matrix(cn[, -(1:5)]); rownames(counts) <- cn$gene
smpR <- fread(inp("TranscriptomicData", "sample_design.tsv"))
smpR[, State := factor(State, levels = as.vector(outer(LIGHTS, c("m02h", "p02h", "p06h", "p10h", "m10h"), paste, sep = "_")))]
MA <- as.matrix(fread(inp("ProteomicData", "2_A_limma", "A_protein_log2ratio.tsv")), rownames = 1)
msB <- readRDS(inp("ProteomicData", "3_B_MSstatsTMT", "B_MSstatsTMT.rds"))
protB <- as.data.table(msB$protein)
smpP <- unique(protB[, .(SampleID = as.character(BioReplicate), Plex = as.character(Mixture), Condition = as.character(Condition))])
smpP <- smpP[SampleID %in% colnames(MA)]
smpP[, Light := c(`0050` = "LL", `0200` = "ML", `1000` = "HL")[sub("^light(\\d+)_.*", "\\1", SampleID)]]
smpP[, State := factor(Condition, levels = as.vector(outer(LIGHTS, c("CTm02", "CTp02", "CTp06", "CTp10", "CTm10"), paste, sep = "_")))]
MA <- MA[, smpP$SampleID]

light_contrasts <- function(states) {   # weights of the 3 light contrasts over the 15 states
  w <- function(l) { v <- setNames(rep(0, length(states)), states); v[startsWith(states, paste0(l, "_"))] <- 1 / 5; v }
  cbind(HL_vs_LL = w("HL") - w("LL"), HL_vs_ML = w("HL") - w("ML"), ML_vs_LL = w("ML") - w("LL"))
}
genesR <- unique(rel[Pipeline == "RNA_limma", id])
fit_rna <- function(pipe, keep) {   # rerun the RNA pipeline on the samples in 'keep'
  sd <- smpR[match(keep, SampleID)]                # keeps duplicates (bootstrap draws)
  sd[, Col := make.unique(SampleID)]                 # a library drawn twice needs two column names
  cm <- counts[genesR, sd$SampleID, drop = FALSE]; colnames(cm) <- sd$Col
  Cm <- light_contrasts(levels(sd$State))
  X <- model.matrix(~ 0 + State, data = sd); colnames(X) <- levels(sd$State)
  if (pipe == "RNA_limma") {
    dge <- calcNormFactors(DGEList(cm))
    v <- voomWithQualityWeights(dge, X); f <- eBayes(contrasts.fit(lmFit(v, X), Cm), robust = TRUE)
    rbindlist(lapply(CONTRASTS, function(k) { t <- topTable(f, coef = k, number = Inf, sort.by = "none")
      data.table(Contrast = k, id = rownames(t), log2FC = t$logFC, P = t$P.Value) }))
  } else {
    dds <- DESeq2::DESeqDataSetFromMatrix(cm, colData = data.frame(sd, row.names = sd$Col), design = ~ 0 + State)
    dds <- DESeq2::DESeq(dds, quiet = TRUE); rn <- DESeq2::resultsNames(dds)
    rbindlist(lapply(CONTRASTS, function(k) { cv <- setNames(Cm[, k], paste0("State", rownames(Cm)))[rn]
      r <- DESeq2::results(dds, contrast = unname(cv)); data.table(Contrast = k, id = rownames(r), log2FC = r$log2FoldChange, P = r$pvalue) }))
  }
}
fit_tmt <- function(pipe, keep) {   # rerun the TMT test without the plex left out (values are kept)
  sd <- smpP[SampleID %in% keep]; Cm <- light_contrasts(levels(sd$State))
  if (pipe == "TMT_A") {
    M <- MA[, sd$SampleID]
    ok <- apply(M, 1, function(x) all(tapply(!is.na(x), sd$State, any)) & uniqueN(sd$Plex[!is.na(x)]) >= 2)
    M <- M[ok, ]
    X <- model.matrix(~ 0 + State, data = sd); colnames(X) <- levels(sd$State)
    dc <- duplicateCorrelation(M, X, block = sd$Plex)
    f <- eBayes(contrasts.fit(lmFit(M, X, block = sd$Plex, correlation = dc$consensus.correlation), Cm), trend = TRUE, robust = TRUE)
    rbindlist(lapply(CONTRASTS, function(k) { t <- topTable(f, coef = k, number = Inf, sort.by = "none")
      data.table(Contrast = k, id = rownames(t), log2FC = t$logFC, P = t$P.Value) }))
  } else {
    pd <- protB[as.character(BioReplicate) %in% sd$SampleID]
    keepP <- pd[!is.na(Abundance), .(n = uniqueN(Mixture)), by = Protein][n >= 2, Protein]
    pd <- pd[Protein %in% keepP]
    cm <- t(Cm); colnames(cm) <- rownames(Cm)
    gc_tmt <- function(mod) { out <- NULL   # capture.output hides MSstatsTMT's progress bars
      invisible(capture.output(out <- MSstatsTMT::groupComparisonTMT(list(ProteinLevelData = as.data.frame(pd)), contrast.matrix = cm,
                                                                     moderated = mod, use_log_file = FALSE, verbose = FALSE)))
      out }
    res <- tryCatch(gc_tmt(TRUE), error = function(e) gc_tmt(FALSE))
    r <- as.data.table(res$ComparisonResult)
    r <- r[is.finite(log2FC) & !is.na(pvalue) & (is.na(issue) | issue != "SingleMeasurePerCondition")]
    r[, .(Contrast = as.character(Label), id = as.character(Protein), log2FC, P = pvalue)]
  }
}
# RNA units: the samples LEFT OUT (jackknife) or the samples DRAWN (bootstrap)
BOOT <- JACK_RNA == "bootstrap"
boot_draws <- function(B) { set.seed(SEED + 1)   # own seed: the same draws whatever ran before
  setNames(lapply(seq_len(B), function(b) unlist(lapply(split(smpR$SampleID, smpR$State), function(x) x[sample.int(length(x), length(x), replace = TRUE)]), use.names = FALSE)),
           sprintf("b%03d", seq_len(B))) }
units <- list(
  RNA = switch(JACK_RNA, sample = setNames(as.list(smpR$SampleID), smpR$SampleID),
               replicate = setNames(lapply(sort(unique(smpR$Rep)), function(r) smpR[Rep == r, SampleID]), sort(unique(smpR$Rep))),
               bootstrap = boot_draws(BOOT_B), stop("JACK_RNA must be sample, replicate or bootstrap")),
  Protein = setNames(lapply(sort(unique(smpP$Plex)), function(p) smpP[Plex == p, SampleID]), sort(unique(smpP$Plex))))
all_samples <- list(RNA = smpR$SampleID, Protein = smpP$SampleID)
# SE and df of a resampling statistic d (one value per resample)
is_boot <- function(om) BOOT && om == "RNA"
rs_se <- function(d, boot) if (boot) sqrt(1.5) * sd(d) else sqrt((length(d) - 1) / length(d) * sum((d - mean(d))^2))
rs_df <- function(om, n) if (is_boot(om)) nrow(smpR) - nlevels(smpR$State) else n - 1

# deterministic NES of a resample: ES / (ES/NES of the full data)
norm_const <- gsea_all[grepl("^eggNOG", Annotation), .(Pipeline, Contrast, Annotation, Term, scale = ES / NES)]
es_only <- function(z, pw) {
  z <- z[is.finite(z) & !duplicated(names(z))]; z <- sort(z, decreasing = TRUE)
  vapply(pw, function(g) { pos <- match(intersect(g, names(z)), names(z))
    if (length(pos) < MIN_SIZE) return(NA_real_) else fgsea::calcGseaStat(z, sort(pos), gseaParam = 1) }, numeric(1))
}
RS_ANN <- grep("^eggNOG", ANN, value = TRUE)
# one cache per pipeline (RNA: per jackknife unit), so that changing JACK_RNA does not redo TMT;
# the gene statistics of every resample are kept too (section 7b rebuilds the lists from them)
cache_of <- function(pp, what = "ES") D("6_Resampling", "cache", if (omic_of(pp) == "RNA") sprintf("%s_%s_%s.rds", what, pp, JACK_RNA) else sprintf("%s_%s.rds", what, pp))
jk <- list(); jst <- list()   # column Left_out = the unit left out (jackknife) or the draw (bootstrap: b001...)
for (pp in PIPES) {
  om <- omic_of(pp); uu <- units[[om]]
  if (!FORCE_RESAMPLING && file.exists(cache_of(pp))) jk[[pp]] <- readRDS(cache_of(pp))
  if (!FORCE_RESAMPLING && file.exists(cache_of(pp, "STATS"))) jst[[pp]] <- readRDS(cache_of(pp, "STATS"))
  # gene statistics are only needed if the ORA resampling (7b) is not cached yet
  need_stats <- FORCE_RESAMPLING || !file.exists(D("6_Resampling", "cache", sprintf("ORA_resampled_RNA-%s.rds", JACK_RNA)))
  done <- if (is.null(jk[[pp]])) character() else unique(jk[[pp]]$Left_out)
  if (need_stats) done <- intersect(done, if (is.null(jst[[pp]])) character() else unique(jst[[pp]]$Left_out))
  todo <- setdiff(names(uu), done)
  if (!length(todo)) next
  cat(sprintf("Resampling %s: %d of %d units to go\n", pp, length(todo), length(uu))); t0 <- Sys.time()
  for (u in todo) {
    keep <- if (is_boot(om)) uu[[u]] else setdiff(all_samples[[om]], uu[[u]])
    st <- if (om == "RNA") fit_rna(pp, keep) else fit_tmt(pp, keep)
    st[, z := zscore(log2FC, P)]
    res <- rbindlist(lapply(CONTRASTS, function(k) { zz <- st[Contrast == k]; z <- setNames(zz$z, zz$id)
      rbindlist(lapply(RS_ANN, function(a) { pw <- prep_sets(get_set(a), names(z))
        data.table(Contrast = k, Annotation = a, Term = names(pw), ES = es_only(z, pw)) })) }))
    if (!is.null(jk[[pp]])) jk[[pp]] <- jk[[pp]][Left_out != u]
    if (!is.null(jst[[pp]])) jst[[pp]] <- jst[[pp]][Left_out != u]
    jk[[pp]] <- rbind(jk[[pp]], res[, `:=`(Pipeline = pp, Left_out = u)])
    jst[[pp]] <- rbind(jst[[pp]], st[, .(Contrast, id, log2FC, P, Pipeline = pp, Left_out = u)])
    saveRDS(jk[[pp]], cache_of(pp)); saveRDS(jst[[pp]], cache_of(pp, "STATS")); rm(st, res); invisible(gc())
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins")); nd <- match(u, todo)
    if (nd %% 5 == 0 || nd == length(todo)) cat(sprintf("  %d / %d (%.0f min, ~%.0f min left)\n", nd, length(todo), el, el / nd * (length(todo) - nd)))
  }
}
jkd <- rbindlist(jk)
jkd <- merge(jkd, norm_const, by = c("Pipeline", "Contrast", "Annotation", "Term"))[, NES := ES / scale]
jk_all <- rbindlist(lapply(names(DPAIRS), function(pn) {
  a <- DPAIRS[[pn]][1]; b <- DPAIRS[[pn]][2]; n <- length(units[[omic_of(a)]])
  full <- merge(gsea_all[Pipeline == a & Annotation %in% RS_ANN, .(Contrast, Annotation, Term, Name, NES_1 = NES, padj_1 = padj)],
                gsea_all[Pipeline == b & Annotation %in% RS_ANN, .(Contrast, Annotation, Term, NES_2 = NES, padj_2 = padj)],
                by = c("Contrast", "Annotation", "Term"))
  full[, `:=`(Delta = NES_1 - NES_2, Enriched = padj_1 < FDR | padj_2 < FDR)]
  m <- merge(jkd[Pipeline == a, .(Contrast, Annotation, Term, Left_out, A = NES)], jkd[Pipeline == b, .(Contrast, Annotation, Term, Left_out, B = NES)],
             by = c("Contrast", "Annotation", "Term", "Left_out"))
  se <- m[is.finite(A) & is.finite(B), .(n_sub = .N, SE = rs_se(A - B, is_boot(omic_of(a)))), by = .(Contrast, Annotation, Term)][n_sub == n]
  r <- merge(full, se, by = c("Contrast", "Annotation", "Term"))
  r[, df := rs_df(omic_of(a), n)][, t := Delta / pmax(SE, 1e-6)][, pvalue := 2 * pt(-abs(t), df)]
  if (JACK_MODERATE) {
    r[, c("SE_mod", "df_mod") := { sq <- squeezeVar(pmax(SE, 1e-6)^2, df = df[1]); list(sqrt(sq$var.post), pmin(df[1] + sq$df.prior, 1000)) },
      by = .(Contrast, Annotation)]
    r[, t_mod := Delta / SE_mod][, pvalue_mod := 2 * pt(-abs(t_mod), df_mod)]
  }
  r[, padj := p.adjust(pvalue, "BH"), by = .(Contrast, Annotation)]
  r[, Pipeline_effect := padj < FDR & abs(Delta) >= DELTA_NES & Enriched]
  if (JACK_MODERATE) { r[, padj_mod := p.adjust(pvalue_mod, "BH"), by = .(Contrast, Annotation)]
    r[, Pipeline_effect_mod := padj_mod < FDR & abs(Delta) >= DELTA_NES & Enriched] }
  r[, Pair := pn]
}), fill = TRUE)
fwrite(jk_all[order(padj)], D("6_Resampling", sprintf("Resampling_NES_difference_RNA-%s.tsv", JACK_RNA)), sep = "\t")
jk_sum <- jk_all[, c(list(Terms = .N, df = df[1], t_crit = round(qt(0.975, df[1]), 2), Median_abs_Delta = round(median(abs(Delta)), 2),
                          Median_SE = round(median(SE), 2), Beyond_noise = sum(padj < FDR), Pipeline_effect = sum(Pipeline_effect)),
                     if (JACK_MODERATE) list(df_moderated = round(median(df_mod), 1), Pipeline_effect_moderated = sum(Pipeline_effect_mod))),
                 by = .(Pair, Contrast)]
fwrite(jk_sum, D("6_Resampling", sprintf("Resampling_summary_RNA-%s.tsv", JACK_RNA)), sep = "\t")
cat("Resampling (eggNOG annotation): is the pipeline effect larger than the effect of the samples?\n"); print(jk_sum)
jk_all[, Class := fifelse(Pipeline_effect, "pipeline effect", fifelse(Enriched, "enriched: within noise or small", "not enriched in either"))]
g <- ggplot(jk_all[Contrast == "HL_vs_LL"][order(Class != "not enriched in either")], aes(Delta, pmin(abs(t), 30), colour = Class)) +
  geom_point(size = 0.8, alpha = 0.7) + facet_wrap(~ Pair) +
  geom_vline(xintercept = c(-DELTA_NES, DELTA_NES), linetype = 2, colour = "grey50") +
  geom_hline(aes(yintercept = qt(0.975, df)), linetype = 3, colour = "grey40") +
  scale_colour_manual(values = c(`pipeline effect` = "#B5541B", `enriched: within noise or small` = "#14213D", `not enriched in either` = "grey80")) +
  labs(x = "NES difference between pipelines (all samples)", y = "|difference| / resampling SE (capped at 30)", colour = NULL,
       title = "High vs low light: pipeline effect vs sampling noise", subtitle = "dotted line: t critical value (95%) for the df of each omic")
ggsave(D("6_Resampling", sprintf("Resampling_volcano_HL_vs_LL_RNA-%s.png", JACK_RNA)), g, width = 11, height = 5)

# ---- 7b. Resampling of ORA: do the LISTS of two pipelines give different functions? ----------
# GSEA ranks all genes; ORA depends on the lists, and lists of ~100 genes can gain or lose a whole
# function with a few genes. Here the FC lists (adj.P < 0.05 and |log2FC| >= LFC) are rebuilt in
# every resample and ORA is rerun. (The Bayesian lists would need the Stan models refitted in
# every resample, too long for a practical, so only the FC criterion is resampled.)
# Two measures per term and direction:
#   - log2 fold enrichment, ER = log2(((k + 0.5)/(n + 1)) / ((K + 0.5)/(N + 1))): the jackknife
#     compares the difference between pipelines with its SE, as for the NES (relevant: |dER| >= DELTA_ER);
#   - stability: in what fraction of the resamples is the term significant in each pipeline.
banner("7b. Resampling of ORA (FC lists)")
ORA_RS_CACHE <- D("6_Resampling", "cache", sprintf("ORA_resampled_RNA-%s.rds", JACK_RNA))
ora_tab <- function(ids, universe, s, terms) {   # ORA statistics of the given terms; BH over ALL terms tested, as in section 4
  u <- intersect(universe, s$id); x <- intersect(ids, u)
  t <- s[id %in% u, .(K = uniqueN(id), k = uniqueN(id[id %in% x])), by = Term][K >= MIN_SIZE & K <= MAX_SIZE]
  if (!nrow(t)) return(NULL)
  t[, `:=`(n = length(x), N = length(u))]
  t[, pvalue := phyper(k - 1, K, N - K, n, lower.tail = FALSE)][, padj := p.adjust(pvalue, "BH")]
  t[, ER := log2(((k + 0.5) / (n + 1)) / ((K + 0.5) / (N + 1)))]
  t[Term %in% terms]
}
fc_lists <- function(st, om) {   # rebuild the FC lists of one resample
  st <- copy(st)[, adj.P := p.adjust(P, "BH"), by = Contrast]
  st[, .(up = list(id[adj.P < FDR & log2FC >= LFC[om]]), down = list(id[adj.P < FDR & log2FC <= -LFC[om]]), universe = list(unique(id))), by = Contrast]
}
ORA_RS_ANN <- grep("^eggNOG", ANN, value = TRUE)
# the terms followed: significant in ORA of the FC lists of either pipeline of a pair, all samples
tracked <- ora_all[Criterion == "FC" & Annotation %in% ORA_RS_ANN, .(Pipeline, Contrast, Direction, Annotation, Term)]
if (!FORCE_RESAMPLING && file.exists(ORA_RS_CACHE)) ora_rs <- readRDS(ORA_RS_CACHE) else {
  t0 <- Sys.time()
  ora_rs <- rbindlist(lapply(names(DPAIRS), function(pn) { pr <- DPAIRS[[pn]]; om <- omic_of(pr[1])
    tr <- unique(tracked[Pipeline %in% pr, .(Contrast, Direction, Annotation, Term)])
    rbindlist(lapply(pr, function(pp) rbindlist(lapply(c("ALL", names(units[[om]])), function(u) {
      st <- if (u == "ALL") rel[Pipeline == pp, .(Contrast, id, log2FC, P)] else jst[[pp]][Left_out == u]
      L <- fc_lists(st, om)
      rbindlist(lapply(seq_len(nrow(L)), function(i) { ctr <- L$Contrast[i]   # not "k": ORA tables have a column k
        rbindlist(lapply(c("up", "down"), function(dirn) rbindlist(lapply(ORA_RS_ANN, function(a) {
          terms <- tr[Contrast == ctr & Direction == dirn & Annotation == a, Term]
          if (!length(terms)) return(NULL)
          r <- ora_tab(L[[dirn]][[i]], L$universe[[i]], get_set(a), terms)
          if (!is.null(r)) r[, `:=`(Contrast = ctr, Direction = dirn, Annotation = a)] }))))
      }))[, `:=`(Pair = pn, Pipeline = pp, Left_out = u)]
    }))))
  }))
  saveRDS(ora_rs, ORA_RS_CACHE)
  cat(sprintf("ORA of the resampled FC lists: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
ora_jk <- rbindlist(lapply(names(DPAIRS), function(pn) { pr <- DPAIRS[[pn]]; n <- length(units[[omic_of(pr[1])]])
  x <- ora_rs[Pair == pn]
  w <- dcast(x, Contrast + Direction + Annotation + Term + Left_out ~ Pipeline, value.var = c("ER", "padj", "k"))
  setnames(w, c(paste0("ER_", pr), paste0("padj_", pr), paste0("k_", pr)), c("ER_1", "ER_2", "padj_1", "padj_2", "k_1", "k_2"))
  full <- w[Left_out == "ALL"]; rs <- w[Left_out != "ALL" & is.finite(ER_1) & is.finite(ER_2)]
  se <- rs[, .(n_sub = .N, SE = rs_se(ER_1 - ER_2, is_boot(omic_of(pr[1]))),
               Stable_1 = mean(padj_1 < FDR & k_1 >= ORA_MIN, na.rm = TRUE), Stable_2 = mean(padj_2 < FDR & k_2 >= ORA_MIN, na.rm = TRUE)),
           by = .(Contrast, Direction, Annotation, Term)][n_sub == n]
  r <- merge(full[, .(Contrast, Direction, Annotation, Term, ER_1, ER_2, padj_1, padj_2, k_1, k_2)], se, by = c("Contrast", "Direction", "Annotation", "Term"))
  r[, `:=`(Sig_1 = padj_1 < FDR & k_1 >= ORA_MIN, Sig_2 = padj_2 < FDR & k_2 >= ORA_MIN)]
  r[, Delta := ER_1 - ER_2][, df := rs_df(omic_of(pr[1]), n)][, t := Delta / pmax(SE, 1e-6)][, pvalue := 2 * pt(-abs(t), df)]
  r[, padj := p.adjust(pvalue, "BH"), by = .(Contrast, Annotation)]
  r[, Pipeline_effect := padj < FDR & abs(Delta) >= DELTA_ER & (Sig_1 | Sig_2)]
  r[, Pair := pn]
}))
ora_jk <- merge(ora_jk, unique(rbindlist(lapply(ORA_RS_ANN, function(a) unique(get_set(a)[, .(Term, Name)])[, Annotation := a]))),
                by = c("Annotation", "Term"), all.x = TRUE)
fwrite(ora_jk[order(padj)], D("6_Resampling", sprintf("Resampling_ORA_difference_RNA-%s.tsv", JACK_RNA)), sep = "\t")
ora_jk_sum <- ora_jk[, .(Terms = .N, df = df[1], Only_one_pipeline = sum(xor(Sig_1, Sig_2)),
                         Beyond_noise = sum(padj < FDR), Pipeline_effect = sum(Pipeline_effect),
                         Median_abs_dER = round(median(abs(Delta)), 2), Median_SE = round(median(SE), 2),
                         Stable_when_significant_1 = round(median(Stable_1[Sig_1]), 2), Stable_when_significant_2 = round(median(Stable_2[Sig_2]), 2)),
                     by = .(Pair, Contrast)]
fwrite(ora_jk_sum, D("6_Resampling", sprintf("Resampling_ORA_summary_RNA-%s.tsv", JACK_RNA)), sep = "\t")
cat("ORA of the FC lists: does the pipeline change the enriched functions beyond sampling noise?\n"); print(ora_jk_sum)

# =================================================================================
# 8. Files for Cytoscape / ClueGO
# =================================================================================
banner("8. ClueGO files")
for (ln in names(lists)) for (dirn in c("up", "down")) { l <- lists[[ln]]
  writeLines(l[[dirn]], D("7_ClueGO", sprintf("List_%s_%s_%s_%s.txt", l$pipeline, l$criterion, l$contrast, dirn))) }
for (a in ANN) { g <- get_set(a)[, .(genes = paste(unique(id), collapse = "\t")), by = .(Term, Name)]
  writeLines(paste(g$Term, g$Name, g$genes, sep = "\t"), D("7_ClueGO", paste0("Annotation_", sub(":", "_", a), ".gmt"))) }

writeLines(capture.output(sessionInfo()), D("sessionInfo.txt"))
cat("\nDone. Results in", file.path(BASE, "results_Enrichment"), "\n")
