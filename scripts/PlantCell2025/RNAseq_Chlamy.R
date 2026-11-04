# =============================================================================
# RNAseq_Chlamy.R
# Transcriptome of Chlamydomonas acclimated to three light intensities across
# the diurnal cycle (Dupuis, Ojeda et al. 2025, Plant Cell 37: koaf086;
# GEO GSE275433): from the published FPKM table to differential expression,
# with two pipelines and two relevance criteria, as in the RNA-seq practical (PC2)
# and the TMT block:
#
#   pipeline A  limma-voom (with sample quality weights)    (step by step)
#   pipeline B  DESeq2                                       (all in one)
#   criterion   FC: |log2FC| >= LFC     Bayesian: P(|log2FC| > LFC) >= P_BAYES
#
# Sections
#   1. Packages
#   2. Counts recovered from the FPKM table                           [cached]
#   3. Design and quality control
#   4. Pipeline A: limma-voom
#   5. Pipeline B: DESeq2                                              [cached]
#   6. Relevance: FC and Bayesian criteria, both pipelines            [cached]
#   7. Comparison of pipelines and criteria
#   8. What if we had used the FPKM? (limma-trend on log2 FPKM)
#
# Heavy steps are cached: if their output exists they are skipped (set the
# FORCE_* switches to TRUE to recompute). Everything that may need editing is in
# section 0.
#
# Files expected in BASE:
#   GSE275433_Cr073.light_intensity.FPKMs.replicates.v6.1.tsv.gz   (GEO, supplementary)
#   CreinhardtiiCC_4532_707_v6.1.gene_exons.gff3.gz                (Phytozome, annotation)
# =============================================================================

# ---- 0. Settings (edit here) ---------------------------------------------------
BASE <- Sys.getenv("RNA_BASE", path.expand("~/BMSi/PlantCell2025/RNAseq"))
FPKM_FILE <- file.path(BASE, "GSE275433_Cr073.light_intensity.FPKMs.replicates.v6.1.tsv.gz")
GFF_FILE  <- file.path(BASE, "CreinhardtiiCC_4532_707_v6.1.gene_exons.gff3.gz")

FDR      <- 0.05        # adjusted p-value for the statistical tests
LFC      <- 1           # relevance threshold: FC 2 (log2FC 1), as in PC2/PC3
P_BAYES  <- 0.9         # Bayesian criterion: P(|log2FC| > LFC) >= P_BAYES
SEED     <- 2026

# Why FC 2 here and FC 1.5 in TMT: read counts have no "ratio compression";
# a 2-fold change in the cell is measured as ~2-fold.

FORCE_COUNTS <- FALSE   # rebuild the count matrix from the FPKM (~1 min)
FORCE_DESEQ  <- FALSE   # rerun DESeq2 (~5 min)
FORCE_BAYES  <- FALSE   # refit the Bayesian models (long)

# ---- 1. Packages ------------------------------------------------------------------
cran <- c("data.table", "ggplot2", "pheatmap", "rstan")
bioc <- c("limma", "edgeR", "DESeq2")
for (p in cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(limma); library(edgeR)
})
theme_set(theme_bw(base_size = 11))
set.seed(SEED)

D  <- function(...) { p <- file.path(BASE, "results_RNAseq", ...); dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE); p }
PROC <- file.path(BASE, "processed"); dir.create(PROC, showWarnings = FALSE)
LIGHTS    <- c("LL", "ML", "HL")
TIMES     <- c("m02h", "p02h", "p06h", "p10h", "m10h")      # -2, +2, +6, +10, +14 h after dawn
HOURS     <- c(m02h = -2, p02h = 2, p06h = 6, p10h = 10, m10h = 14)
LIGHT_COL <- c(LL = "#2A6F97", ML = "#7F8A99", HL = "#B5541B")
CONTRASTS <- c("HL_vs_LL", "HL_vs_ML", "ML_vs_LL")
banner <- function(x) cat("\n==== ", x, " ", strrep("=", max(3, 70 - nchar(x))), "\n", sep = "")

# =================================================================================
# 2. Counts recovered from the FPKM table
# =================================================================================
# GEO only provides FPKM (Cufflinks):  FPKM = count * 1e9 / (Leff * N)
#   Leff = effective length = transcript length - mean fragment length (mu)
#   N    = fragments of the sample (library size)
# Counts are whole numbers, so for each sample we look for the N (and a common mu)
# that turn the FPKM of low-expressed genes into integers. Exact for single-isoform
# genes; a close approximation for multi-isoform genes (Cufflinks adds up isoforms,
# each with its own length) and very short genes.
banner("2. Counts from FPKM")
COUNTS_FILE <- file.path(PROC, "GSE275433_counts_reconstructed.tsv.gz")
fp <- fread(FPKM_FILE)
fpkm <- as.matrix(fp[, -(1:3)]); rownames(fpkm) <- fp$Cre_GID
if (FORCE_COUNTS || !file.exists(COUNTS_FILE)) {
  # gene length = exons of the primary transcript ('longest=1') in the gff3
  g <- fread(GFF_FILE, header = FALSE, sep = "\t", quote = "", skip = "chromosome",
             col.names = c("chr", "src", "type", "start", "end", "score", "strand", "phase", "attr"))
  mrna <- g[type == "mRNA"]
  mrna[, tx := sub(".*ID=([^;]+).*", "\\1", attr)]
  mrna[, gene := sub("_4532\\..*", "", tx)]
  mrna[, primary := grepl("longest=1", attr)]
  ex <- g[type == "exon"][, tx := sub(".*Parent=([^;]+).*", "\\1", attr)]
  mrna <- merge(mrna[, .(tx, gene, primary)], ex[, .(len = sum(end - start + 1)), by = tx], by = "tx")
  genes <- mrna[, .(length = len[primary][1], n_isoforms = .N), by = gene]
  L    <- setNames(genes$length, genes$gene)[rownames(fpkm)]
  niso <- setNames(genes$n_isoforms, genes$gene)[rownames(fpkm)]
  stopifnot(!anyNA(L)); rm(g, ex, mrna)

  # integrality: mean distance to the nearest integer of the low counts (0.5-30) of
  # single-isoform genes > 1 kb; ~0.25 when N is wrong, ~0.02-0.05 when right
  frac_err <- function(f, Lg, N, mu) {
    cc <- f * (Lg - mu) * N / 1e9; cc <- cc[cc > 0.5 & cc < 30]
    if (length(cc) < 50) return(NA_real_)
    mean(abs(cc - round(cc)))
  }
  fit_N <- function(f, Lg, mu) {
    grid <- seq(1e6, 60e6, by = 2e4)
    err  <- vapply(grid, function(n) frac_err(f, Lg, n, mu), 0)
    # 2N, 3N... also give whole numbers (every count doubled): take the SMALLEST
    # N whose error is close to the best one, then the minimum of that basin
    i0  <- which(err < min(err, na.rm = TRUE) + 0.05)[1]
    win <- which(abs(grid / grid[i0] - 1) < 0.02)
    N <- grid[win][which.min(err[win])]
    for (s in c(2e-3, 4e-4, 8e-5)) {
      gr <- N * (1 + seq(-s, s, length.out = 81))
      N  <- gr[which.min(vapply(gr, function(n) frac_err(f, Lg, n, mu), 0))]
    }
    N
  }
  use <- niso == 1 & L > 1000
  f1 <- fpkm[use, 1]; L1 <- L[use]; k1 <- f1 > 0
  mu_err <- function(m) { N <- fit_N(f1[k1], L1[k1], m); frac_err(f1[k1], L1[k1], N, m) }
  mu <- seq(150, 260, by = 10)[which.min(vapply(seq(150, 260, by = 10), mu_err, 0))]
  mu <- seq(mu - 10, mu + 10, by = 0.5)[which.min(vapply(seq(mu - 10, mu + 10, by = 0.5), mu_err, 0))]
  cat(sprintf("Mean fragment length (mu): %.1f bp\n", mu))
  qc <- rbindlist(lapply(seq_len(ncol(fpkm)), function(s) {
    f <- fpkm[use, s]; k <- f > 0; N <- fit_N(f[k], L[use][k], mu)
    data.table(SampleID = colnames(fpkm)[s], N_fragments = N, int_error = frac_err(f[k], L[use][k], N, mu))
  }))
  q <- fpkm * pmax(L - mu, 1) * rep(qc$N_fragments, each = nrow(fpkm)) / 1e9
  counts <- round(q); storage.mode(counts) <- "integer"
  dev <- abs(q - round(q)); low <- q > 0.5 & q < 50
  cat(sprintf("Low counts (<50) that are whole numbers: single-isoform genes %.0f%%, multi-isoform genes %.0f%%\n",
              100 * mean(dev[niso == 1 & L > 500, ][low[niso == 1 & L > 500, ]] < 0.05),
              100 * mean(dev[niso > 1, ][low[niso > 1, ]] < 0.05)))
  fwrite(data.table(gene = rownames(counts), symbol = fp$Cre_symbol, defline = fp$Cre_defline,
                    length = L, n_isoforms = niso, counts), COUNTS_FILE, sep = "\t")
  fwrite(qc, file.path(PROC, "GSE275433_counts_reconstruction_QC.tsv"), sep = "\t")
  rm(q, dev, low)
} else cat("Using cached", COUNTS_FILE, "\n")
cn <- fread(COUNTS_FILE)
counts <- as.matrix(cn[, -(1:5)]); rownames(counts) <- cn$gene
ann <- cn[, .(gene, symbol, defline, length, n_isoforms)]
cat(sprintf("Counts: %d genes x %d samples\n", nrow(counts), ncol(counts)))

# =================================================================================
# 3. Design and quality control
# =================================================================================
banner("3. Design and QC")
# sample names: Cr073_<letter>_<light>_<Dk|Lt>_<time>___BR<replicate>
smp <- data.table(SampleID = colnames(counts))
smp[, c("Light", "Phase", "Time", "Rep") := tstrsplit(sub("^Cr073_[A-O]_", "", SampleID), "_+")]
smp[, `:=`(Light = factor(Light, LIGHTS), Time = factor(Time, TIMES), Rep = factor(Rep))]
smp[, State := factor(paste(Light, Time, sep = "_"), levels = as.vector(outer(LIGHTS, TIMES, paste, sep = "_")))]
smp[, Hour := HOURS[as.character(Time)]]
stopifnot(all(table(smp$State) == 3))
fwrite(smp, D("1_QC", "sample_design.tsv"), sep = "\t")

lib <- data.table(smp, M_reads = colSums(counts) / 1e6)
cat("Library size (million fragments):\n"); print(summary(lib$M_reads))
print(lib[M_reads < 0.5 * median(M_reads), .(SampleID, M_reads = round(M_reads, 2))])
g <- ggplot(lib, aes(reorder(SampleID, M_reads), M_reads, fill = Light)) + geom_col() + coord_flip() +
  scale_fill_manual(values = LIGHT_COL) + labs(x = NULL, y = "Million fragments", title = "Library sizes")
ggsave(D("1_QC", "library_sizes.png"), g, width = 6.5, height = 8)

dge <- DGEList(counts, samples = as.data.frame(smp))
keep <- filterByExpr(dge, group = smp$State)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
cat(sprintf("Genes expressed (filterByExpr): %d of %d\n", nrow(dge), nrow(counts)))
logcpm <- cpm(dge, log = TRUE, prior.count = 2)

pc <- prcomp(t(logcpm)); ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
pcd <- cbind(smp, pc$x[, 1:3])
g1 <- ggplot(pcd, aes(PC1, PC2, colour = Light, shape = Time)) + geom_point(size = 3) +
  scale_colour_manual(values = LIGHT_COL) + labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
  title = "PCA, log CPM", subtitle = "PC1 and PC2: time of day")
g2 <- ggplot(pcd, aes(PC1, PC3, colour = Light, shape = Time)) + geom_point(size = 3) +
  scale_colour_manual(values = LIGHT_COL) + labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC3 (%.1f%%)", ve[3]),
  title = "PCA, log CPM", subtitle = "PC3: light intensity")
ggsave(D("1_QC", "PCA_PC1_PC2.png"), g1, width = 6.5, height = 5)
ggsave(D("1_QC", "PCA_PC1_PC3.png"), g2, width = 6.5, height = 5)
cat(sprintf("PCA: PC1 %.1f%%, PC2 %.1f%%, PC3 %.1f%%\n", ve[1], ve[2], ve[3]))
print(pcd[, .(PC1 = round(mean(PC1)), PC2 = round(mean(PC2)), PC3 = round(mean(PC3))), by = Time])
print(pcd[, .(PC1 = round(mean(PC1)), PC2 = round(mean(PC2)), PC3 = round(mean(PC3))), by = Light])

ann_col <- data.frame(Light = smp$Light, Time = smp$Time, row.names = smp$SampleID)
pheatmap::pheatmap(cor(logcpm), annotation_col = ann_col, show_colnames = FALSE, fontsize_row = 6,
                   annotation_colors = list(Light = LIGHT_COL), main = "Sample correlation (log CPM)",
                   filename = D("1_QC", "sample_correlation_heatmap.png"), width = 9, height = 8)
# each sample against the mean of its own state: an outlier stands out
st_mean <- sapply(levels(smp$State), function(s) rowMeans(logcpm[, smp$State == s, drop = FALSE]))
r_own <- sapply(seq_len(ncol(logcpm)), function(i) cor(logcpm[, i], st_mean[, as.character(smp$State[i])]))
lib[, r_own_state := round(r_own, 4)]
cat("Lowest correlation with the mean of their own state:\n"); print(lib[order(r_own_state)][1:4, .(SampleID, M_reads = round(M_reads, 2), r_own_state)])
fwrite(lib, D("1_QC", "sample_QC.tsv"), sep = "\t")

# Is "BR1/BR2/BR3" a batch (e.g. one culture sampled all day)? If so, it should be blocked.
X  <- model.matrix(~ 0 + State, data = smp); colnames(X) <- levels(smp$State)
dc <- duplicateCorrelation(voom(dge, X), X, block = smp$Rep)
cat(sprintf("Correlation within replicate label: %.3f (close to 0: replicates are not batches, no blocking needed)\n",
            dc$consensus.correlation))

# =================================================================================
# 4. Pipeline A: limma-voom
# =================================================================================
banner("4. Pipeline A: limma-voom")
# voom turns counts into log CPM with a precision weight per observation (from the
# mean-variance trend). voomWithQualityWeights also gives each SAMPLE a weight, so
# a noisier sample (e.g. the shallow LL -2 h BR1 library) counts less, without removing it.
png(D("1_QC", "A_voom_mean_variance.png"), 1100, 500)
vA <- voomWithQualityWeights(dge, X, plot = TRUE)
dev.off()
sw <- data.table(SampleID = colnames(vA), weight = round(vA$targets$sample.weights, 2))
cat("Sample quality weights (lowest):\n"); print(sw[order(weight)][1:4])
fwrite(sw, D("2_A_limma_voom", "A_sample_weights.tsv"), sep = "\t")
fitA <- lmFit(vA, X)

w  <- function(states) { v <- setNames(rep(0, ncol(X)), colnames(X)); v[states] <- 1 / length(states); v }
st <- function(l = LIGHTS, t = TIMES) as.vector(outer(l, t, paste, sep = "_"))
C_light <- cbind(HL_vs_LL = w(st("HL")) - w(st("LL")), HL_vs_ML = w(st("HL")) - w(st("ML")), ML_vs_LL = w(st("ML")) - w(st("LL")))
C_time  <- sapply(TIMES[-1], function(t) w(st(LIGHTS, t)) - w(st(LIGHTS, TIMES[1])))
C_int   <- do.call(cbind, lapply(c("HL", "ML"), function(l) sapply(TIMES[-1], function(t)
  (w(st(l, t)) - w(st("LL", t))) - (w(st(l, TIMES[1])) - w(st("LL", TIMES[1]))))))
eb <- function(C) eBayes(contrasts.fit(fitA, C), robust = TRUE)
fl <- eb(C_light); ft <- eb(C_time); fi <- eb(C_int)

resA <- rbindlist(lapply(CONTRASTS, function(k) {
  tt <- topTable(fl, coef = k, number = Inf, sort.by = "none")
  data.table(Gene = rownames(tt), Contrast = k, log2FC = tt$logFC, SE = tt$logFC / tt$t,
             AveLogCPM = tt$AveExpr, P = tt$P.Value, adj.P = tt$adj.P.Val)
}))
fwrite(merge(resA, ann, by.x = "Gene", by.y = "gene"), D("2_A_limma_voom", "A_light_contrasts.tsv"), sep = "\t")
ftab <- function(f, name) { tt <- topTable(f, number = Inf, sort.by = "none")
  data.table(Gene = rownames(tt), Test = name, F = tt$F, P = tt$P.Value, adj.P = tt$adj.P.Val) }
timeA <- ftab(ft, "time_of_day"); intA <- ftab(fi, "light_x_time")
fwrite(rbind(timeA, intA), D("2_A_limma_voom", "A_time_and_interaction_Ftests.tsv"), sep = "\t")
cat(sprintf("Time of day (F-test): %d genes | light x time interaction: %d genes (adj.P < %.2f) | of %d\n",
            sum(timeA$adj.P < FDR), sum(intA$adj.P < FDR), FDR, nrow(timeA)))

top_t <- timeA[order(P)][1:12, Gene]
lab <- ann[match(top_t, gene), fifelse(is.na(symbol) | symbol == "NA", gene, paste0(symbol, " (", gene, ")"))]
prof <- melt(data.table(Gene = lab, logcpm[top_t, , drop = FALSE]), id.vars = "Gene", variable.name = "SampleID", value.name = "y")
prof <- merge(prof, smp[, .(SampleID, Light, Hour)], by = "SampleID")
g <- ggplot(prof, aes(Hour, y, colour = Light)) +
  annotate("rect", xmin = 12, xmax = 15, ymin = -Inf, ymax = Inf, fill = "grey88") +
  annotate("rect", xmin = -3, xmax = 0, ymin = -Inf, ymax = Inf, fill = "grey88") +
  geom_point(size = 1, alpha = 0.6) + stat_summary(fun = mean, geom = "line") +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) + scale_colour_manual(values = LIGHT_COL) +
  scale_x_continuous(breaks = unname(HOURS), labels = c("-2", "+2", "+6", "+10", "+14")) +
  labs(title = "Pipeline A: strongest time-of-day effects (grey = dark)", x = "Hours after dawn", y = "log2 CPM")
ggsave(D("2_A_limma_voom", "A_profiles_top_time.png"), g, width = 12, height = 7)

# =================================================================================
# 5. Pipeline B: DESeq2
# =================================================================================
banner("5. Pipeline B: DESeq2")
DS_CACHE <- D("3_B_DESeq2", "B_DESeq2.rds")
if (FORCE_DESEQ || !file.exists(DS_CACHE)) {
  suppressPackageStartupMessages(library(DESeq2))
  t0 <- Sys.time()
  # DESeq2 uses raw counts (its own normalization: median of ratios) and the same genes as A
  dds <- DESeqDataSetFromMatrix(counts[rownames(dge), ], colData = as.data.frame(smp), design = ~ 0 + State)
  dds <- DESeq(dds, quiet = TRUE)
  rn <- resultsNames(dds)                       # "StateLL_m02h", ...
  cvec <- function(col) { v <- setNames(C_light[, col], paste0("State", rownames(C_light))); unname(v[rn]) }
  light <- rbindlist(lapply(CONTRASTS, function(k) {
    r <- results(dds, contrast = cvec(k), alpha = FDR)
    data.table(Gene = rownames(r), Contrast = k, log2FC = r$log2FoldChange, SE = r$lfcSE,
               baseMean = r$baseMean, P = r$pvalue, adj.P = r$padj)
  }))
  # time of day and interaction: likelihood-ratio tests between nested models
  lrt <- function(full, reduced, name) {
    d2 <- DESeqDataSetFromMatrix(counts[rownames(dge), ], colData = as.data.frame(smp), design = full)
    d2 <- DESeq(d2, test = "LRT", reduced = reduced, quiet = TRUE)
    r <- results(d2); data.table(Gene = rownames(r), Test = name, stat = r$stat, P = r$pvalue, adj.P = r$padj)
  }
  ftest <- rbind(lrt(~ Light + Time, ~ Light, "time_of_day"), lrt(~ Light * Time, ~ Light + Time, "light_x_time"))
  disp <- data.table(Gene = rownames(dds), baseMean = rowData(dds)$baseMean, disp_fit = rowData(dds)$dispFit,
                     disp_gene = rowData(dds)$dispGeneEst, disp_final = dispersions(dds))
  vst_m <- assay(vst(dds, blind = FALSE))
  saveRDS(list(light = light, ftest = ftest, disp = disp, sf = sizeFactors(dds),
               norm = counts(dds, normalized = TRUE), vst = vst_m,
               disp_prior_var = attr(dispersionFunction(dds), "dispPriorVar"),
               disp_fun = dispersionFunction(dds)), DS_CACHE)
  png(D("1_QC", "B_DESeq2_dispersion.png"), 700, 550); plotDispEsts(dds, main = "Pipeline B: dispersion estimates"); dev.off()
  cat(sprintf("DESeq2 took %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(dds)
} else cat("Using cached", DS_CACHE, "\n")
dsB <- readRDS(DS_CACHE)
resB <- dsB$light[!is.na(adj.P)]        # DESeq2 sets adj.P = NA for low-count genes (independent filtering)
fwrite(merge(resB, ann, by.x = "Gene", by.y = "gene"), D("3_B_DESeq2", "B_light_contrasts.tsv"), sep = "\t")
fwrite(dsB$ftest, D("3_B_DESeq2", "B_time_and_interaction_LRT.tsv"), sep = "\t")
cat(sprintf("Time of day (LRT): %d genes | light x time interaction: %d genes (adj.P < %.2f)\n",
            dsB$ftest[Test == "time_of_day", sum(adj.P < FDR, na.rm = TRUE)],
            dsB$ftest[Test == "light_x_time", sum(adj.P < FDR, na.rm = TRUE)], FDR))
cat(sprintf("Size factors (DESeq2): %.2f - %.2f\n", min(dsB$sf), max(dsB$sf)))

# =================================================================================
# 6. Relevance: fold change and Bayesian criteria
# =================================================================================
banner("6. Relevance")
# One Bayesian model per gene (Stan), only for genes significant in a light contrast
# whose |log2FC| + 3 SE reaches LFC (the others cannot reach P_BAYES). Same factorial
# structure as in TMT: y = time-of-day mean + light effect (LL = 0); the light x time
# interaction is not modelled (it goes into the noise: conservative). Empirical priors
# come from ALL genes of each pipeline.
#   A (log CPM, voom weights): normal likelihood; each observation's SD = sigma / sqrt(voom
#     weight); prior on sigma centred on the trend of log(residual SD) with mean log CPM.
#   B (counts): negative binomial likelihood with DESeq2 size factors; prior on the
#     dispersion centred on DESeq2's fitted trend, with DESeq2's own prior width.
#   Both: light effects ~ Student-t(3, 0, s), s = robust spread of all light effects.
stan_A <- "
data { int<lower=1> N; vector[N] y; vector<lower=0>[N] wt; array[N] int<lower=1,upper=5> tt; array[N] int<lower=1,upper=3> ll;
       real mu_loc; real<lower=0> mu_scale; real<lower=0> s_beta; real logsig_loc; real<lower=0> logsig_scale; }
parameters { vector[5] mu; real b_ML; real b_HL; real log_sigma; }
transformed parameters { vector[3] gamma = [0, b_ML, b_HL]'; }
model { mu ~ normal(mu_loc, mu_scale); b_ML ~ student_t(3, 0, s_beta); b_HL ~ student_t(3, 0, s_beta);
        log_sigma ~ normal(logsig_loc, logsig_scale);
        y ~ normal(mu[tt] + gamma[ll], exp(log_sigma) ./ sqrt(wt)); }
generated quantities { real HL_vs_LL = b_HL; real HL_vs_ML = b_HL - b_ML; real ML_vs_LL = b_ML; }
"
stan_B <- "
data { int<lower=1> N; array[N] int<lower=0> y; vector[N] log_sf; array[N] int<lower=1,upper=5> tt; array[N] int<lower=1,upper=3> ll;
       real mu_loc; real<lower=0> mu_scale; real<lower=0> s_beta; real logdisp_loc; real<lower=0> logdisp_scale; }
parameters { vector[5] mu; real b_ML; real b_HL; real log_disp; }
transformed parameters { vector[3] gamma = [0, b_ML, b_HL]'; }
model { mu ~ normal(mu_loc, mu_scale); b_ML ~ student_t(3, 0, s_beta); b_HL ~ student_t(3, 0, s_beta);
        log_disp ~ normal(logdisp_loc, logdisp_scale);
        y ~ neg_binomial_2_log(log_sf + mu[tt] + gamma[ll], inv(exp(log_disp))); }
generated quantities {            // natural log -> log2
  real HL_vs_LL = b_HL / log(2); real HL_vs_ML = (b_HL - b_ML) / log(2); real ML_vs_LL = b_ML / log(2); }
"
info <- smp
Xa <- model.matrix(~ Time + Light, data = smp)
hyperA <- function() {
  f <- lmFit(vA, Xa)
  d <- data.table(logsd = log(f$sigma), n = f$df.residual, mean = f$Amean)
  d[, logsd_c := logsd - (digamma(n / 2) - log(n / 2)) / 2]          # bias of log(sample SD)
  trend <- loess(logsd_c ~ mean, data = d, span = 0.5)
  tau <- sqrt(max(var(residuals(trend)) - mean(trigamma(d$n / 2) / 4), 0.05^2))
  list(trend = trend, tau = tau, s_beta = max(mad(as.vector(f$coefficients[, c("LightML", "LightHL")])), 0.05),
       mu_loc = median(vA$E), mu_scale = 3 * sd(vA$E))
}
hyperB <- function() {
  lnorm <- log(dsB$norm + 0.5)
  f <- lmFit(lnorm, Xa)
  list(s_beta = max(mad(as.vector(f$coefficients[, c("LightML", "LightHL")])), 0.05),
       mu_loc = median(lnorm), mu_scale = 3 * sd(lnorm),
       logdisp_scale = sqrt(max(dsB$disp_prior_var, 0.05^2)))
}
bayes_pipeline <- function(res, name) {
  cache <- D("4_Relevance", paste0("Bayes_", name, ".rds"))
  n_sig <- res[adj.P < FDR, uniqueN(Gene)]
  sig <- res[adj.P < FDR & abs(log2FC) + 3 * SE >= LFC, unique(Gene)]
  if (!FORCE_BAYES && file.exists(cache)) { out <- readRDS(cache); todo <- setdiff(sig, out$Gene) } else { out <- NULL; todo <- sig }
  cat(sprintf("Pipeline %s: %d significant genes, %d could reach the threshold, %d still to fit\n",
              name, n_sig, length(sig), length(todo)))
  if (length(todo)) {
    mod <- rstan::stan_model(model_code = if (name == "A") stan_A else stan_B)
    h <- if (name == "A") hyperA() else hyperB()
    fits <- list(); k <- 0; t0 <- Sys.time()
    for (p in todo) {
      if (name == "A") {
        y <- vA$E[p, ]
        dat <- list(N = length(y), y = unname(y), wt = unname(vA$weights[match(p, rownames(vA)), ]),
                    tt = as.integer(info$Time), ll = as.integer(info$Light), mu_loc = h$mu_loc, mu_scale = h$mu_scale,
                    s_beta = h$s_beta, logsig_loc = unname(predict(h$trend, data.frame(mean = mean(y)))), logsig_scale = h$tau)
      } else {
        y <- counts[p, ]
        dat <- list(N = length(y), y = as.integer(unname(y)), log_sf = unname(log(dsB$sf)),
                    tt = as.integer(info$Time), ll = as.integer(info$Light), mu_loc = h$mu_loc, mu_scale = h$mu_scale,
                    s_beta = h$s_beta, logdisp_loc = log(dsB$disp[Gene == p, disp_fit]), logdisp_scale = h$logdisp_scale)
      }
      invisible(capture.output(fit <- rstan::sampling(mod, data = dat, chains = 4, iter = 1500,
                                                      seed = SEED, refresh = 0, cores = 1)))
      dr <- rstan::extract(fit, CONTRASTS)
      rh <- max(rstan::summary(fit, pars = c("mu", "b_ML", "b_HL"))$summary[, "Rhat"])
      fits[[p]] <- rbindlist(lapply(CONTRASTS, function(k2) { b <- dr[[k2]]
        data.table(Gene = p, Contrast = k2, post_mean = mean(b), low95 = quantile(b, 0.025), high95 = quantile(b, 0.975),
                   P_up = mean(b > LFC), P_down = mean(b < -LFC), Rhat = rh) }))
      k <- k + 1
      if (k %% 100 == 0 || k == length(todo)) {
        out <- rbind(out, rbindlist(fits)); fits <- list(); saveRDS(out, cache)
        el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
        cat(sprintf("  %d / %d genes (%.0f min, ~%.0f min left)\n", k, length(todo), el, el / k * (length(todo) - k)))
      }
    }
  }
  out[Gene %in% sig]
}
bA <- bayes_pipeline(resA, "A")
bB <- bayes_pipeline(resB, "B")
if (any(c(bA$Rhat, bB$Rhat) > 1.01)) warning("Some Bayesian fits have Rhat > 1.01")

classify <- function(res, bay, pipeline) {
  r <- merge(res[, .(Gene, Contrast, log2FC, SE, P, adj.P)],
             bay[, .(Gene, Contrast, post_mean, low95, high95, P_up, P_down)], by = c("Gene", "Contrast"), all.x = TRUE)
  r[, Pipeline := pipeline]
  r[, FC_relevant := adj.P < FDR & abs(log2FC) >= LFC]
  r[, Bayes_relevant := adj.P < FDR & !is.na(P_up) & pmax(P_up, P_down) >= P_BAYES]
  r[, Direction := fifelse(log2FC > 0, "up", "down")]
  r
}
REL <- rbind(classify(resA, bA, "A_limma_voom"), classify(resB, bB, "B_DESeq2"))
REL <- merge(REL, ann[, .(Gene = gene, symbol, defline)], by = "Gene", all.x = TRUE)
fwrite(REL, D("4_Relevance", "relevance_all.tsv"), sep = "\t")
summ <- REL[, .(tested = .N, significant = sum(adj.P < FDR),
                FC_up = sum(FC_relevant & Direction == "up"), FC_down = sum(FC_relevant & Direction == "down"),
                Bayes_up = sum(Bayes_relevant & Direction == "up"), Bayes_down = sum(Bayes_relevant & Direction == "down")),
            by = .(Pipeline, Contrast)][order(Contrast, Pipeline)]
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
both <- merge(resA[, .(Gene, Contrast, A = log2FC)], resB[, .(Gene, Contrast, B = log2FC)], by = c("Gene", "Contrast"))
corr <- both[, .(genes = .N, r_log2FC = round(cor(A, B), 3), slope_B_on_A = round(coef(lm(B ~ A))[2], 2)), by = Contrast]
cat("log2FC, pipeline A vs B:\n"); print(corr)
g <- ggplot(both, aes(A, B)) + geom_point(size = 0.4, alpha = 0.4) + geom_abline(linetype = 2, colour = "#B5541B") +
  facet_wrap(~ Contrast) + labs(x = "log2FC, pipeline A (limma-voom)", y = "log2FC, pipeline B (DESeq2)")
ggsave(D("5_Comparison", "log2FC_A_vs_B.png"), g, width = 11, height = 4)

lists <- function(k) list(`A FC` = REL[Pipeline == "A_limma_voom" & Contrast == k & FC_relevant, paste(Gene, Direction)],
                          `A Bayes` = REL[Pipeline == "A_limma_voom" & Contrast == k & Bayes_relevant, paste(Gene, Direction)],
                          `B FC` = REL[Pipeline == "B_DESeq2" & Contrast == k & FC_relevant, paste(Gene, Direction)],
                          `B Bayes` = REL[Pipeline == "B_DESeq2" & Contrast == k & Bayes_relevant, paste(Gene, Direction)])
jac <- rbindlist(lapply(CONTRASTS, function(k) {
  L <- lists(k); cmb <- combn(names(L), 2)
  data.table(Contrast = k, list1 = cmb[1, ], list2 = cmb[2, ],
             n1 = sapply(cmb[1, ], function(a) length(L[[a]])), n2 = sapply(cmb[2, ], function(a) length(L[[a]])),
             shared = apply(cmb, 2, function(x) length(intersect(L[[x[1]]], L[[x[2]]]))),
             Jaccard = apply(cmb, 2, function(x) { u <- length(union(L[[x[1]]], L[[x[2]]]))
               if (u == 0) NA_real_ else round(length(intersect(L[[x[1]]], L[[x[2]]])) / u, 2) }))
}))
cat("\nOverlap of the relevant lists (same gene, same direction):\n"); print(jac)
fwrite(jac, D("5_Comparison", "list_overlap.tsv"), sep = "\t")
core <- rbindlist(lapply(CONTRASTS, function(k) { x <- Reduce(intersect, lists(k))
  if (!length(x)) return(NULL)
  data.table(Contrast = k, Gene = sub(" .*", "", x), Direction = sub(".* ", "", x)) }))
core <- merge(core, ann[, .(Gene = gene, symbol, defline)], by = "Gene", all.x = TRUE)
core <- merge(core, resA[, .(Gene, Contrast, log2FC_A = round(log2FC, 2))], by = c("Gene", "Contrast"))
fwrite(core[order(Contrast, -abs(log2FC_A))], D("5_Comparison", "core_genes_all_four_lists.tsv"), sep = "\t")
cat("Core genes (in all four lists):\n"); print(core[, .N, by = .(Contrast, Direction)])

# =================================================================================
# 8. What if we had used the FPKM?
# =================================================================================
banner("8. FPKM instead of counts")
# The usual shortcut when only FPKM is available: limma-trend on log2(FPKM + 1).
# It runs, but it ignores that a gene with few reads is much noisier than its FPKM says.
lf <- log2(fpkm[rownames(dge), ] + 1)
fitF <- eBayes(contrasts.fit(lmFit(lf, X), C_light), trend = TRUE, robust = TRUE)
resF <- rbindlist(lapply(CONTRASTS, function(k) { tt <- topTable(fitF, coef = k, number = Inf, sort.by = "none")
  data.table(Gene = rownames(tt), Contrast = k, log2FC = tt$logFC, adj.P = tt$adj.P.Val) }))
resF[, FC_relevant := adj.P < FDR & abs(log2FC) >= LFC][, Direction := fifelse(log2FC > 0, "up", "down")]
cmpF <- rbindlist(lapply(CONTRASTS, function(k) {
  a <- REL[Pipeline == "A_limma_voom" & Contrast == k & FC_relevant, paste(Gene, Direction)]
  f <- resF[Contrast == k & FC_relevant == TRUE, paste(Gene, Direction)]
  m <- merge(resA[Contrast == k, .(Gene, A = log2FC)], resF[Contrast == k, .(Gene, F = log2FC)], by = "Gene")
  data.table(Contrast = k, r_log2FC = round(cor(m$A, m$F), 3), slope_F_on_A = round(coef(lm(F ~ A, m))[2], 2),
             FC_list_counts = length(a), FC_list_FPKM = length(f), shared = length(intersect(a, f)),
             only_FPKM = length(setdiff(f, a)), only_counts = length(setdiff(a, f)))
}))
print(cmpF)
fwrite(cmpF, D("5_Comparison", "FPKM_vs_counts.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), D("sessionInfo.txt"))
cat("\nDone. Results in", file.path(BASE, "results_RNAseq"), "\n")
