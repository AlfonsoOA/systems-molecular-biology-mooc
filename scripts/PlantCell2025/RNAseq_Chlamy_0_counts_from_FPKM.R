# =============================================================================
# RNA-seq Chlamydomonas (Dupuis, Ojeda et al. 2025, Plant Cell; GEO GSE275433)
# Step 0 - Recover read counts from the published FPKM table
# -----------------------------------------------------------------------------
# GEO only provides FPKM (Cufflinks). FPKM = count * 1e9 / (Leff * N), where
#   Leff = effective length = transcript length - mean fragment length (mu)
#   N    = mapped fragments of the sample (library size)
# Counts are integers, so for each sample we look for the N (and a common mu)
# that turn the FPKM of low-expressed genes into whole numbers.
# Exact for single-isoform genes; close approximation for multi-isoform genes
# (Cufflinks sums isoforms with their own lengths) and very short genes.
# Inputs (in BASE):
#   GSE275433_Cr073.light_intensity.FPKMs.replicates.v6.1.tsv.gz  (GEO)
#   CreinhardtiiCC_4532_707_v6.1.gene_exons.gff3.gz               (Phytozome)
# Output: processed/GSE275433_counts_reconstructed.tsv.gz (+ QC table)
# =============================================================================
BASE <- Sys.getenv("RNA_BASE", path.expand("~/BMSi/PlantCell2025/RNAseq"))
FPKM_FILE <- file.path(BASE, "GSE275433_Cr073.light_intensity.FPKMs.replicates.v6.1.tsv.gz")
GFF_FILE  <- file.path(BASE, "CreinhardtiiCC_4532_707_v6.1.gene_exons.gff3.gz")
OUT_DIR   <- file.path(BASE, "processed"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(data.table))

## 1. FPKM table -----------------------------------------------------------
d <- fread(FPKM_FILE)
fpkm <- as.matrix(d[, -(1:3)]); rownames(fpkm) <- d$Cre_GID
cat("FPKM:", nrow(fpkm), "genes x", ncol(fpkm), "samples\n")

## 2. Gene lengths from the gff3 (primary transcript = 'longest=1') --------
g <- fread(GFF_FILE, header = FALSE, sep = "\t", quote = "", skip = "chromosome",
           col.names = c("chr","src","type","start","end","score","strand","phase","attr"))
mrna <- g[type == "mRNA"]
mrna[, tx := sub(".*ID=([^;]+).*", "\\1", attr)]
mrna[, gene := sub("_4532\\..*", "", tx)]
mrna[, primary := grepl("longest=1", attr)]
ex <- g[type == "exon"]
ex[, tx := sub(".*Parent=([^;]+).*", "\\1", attr)]
txlen <- ex[, .(len = sum(end - start + 1)), by = tx]
mrna <- merge(mrna[, .(tx, gene, primary)], txlen, by = "tx")
genes <- mrna[, .(length = len[primary][1], n_isoforms = .N), by = gene]
L    <- setNames(genes$length, genes$gene)[rownames(fpkm)]
niso <- setNames(genes$n_isoforms, genes$gene)[rownames(fpkm)]
stopifnot(!anyNA(L))

## 3. Integrality score ------------------------------------------------------
# mean distance to the nearest integer of the low counts (0.5-30) of
# single-isoform genes longer than 1 kb: ~0.25 if wrong, ~0 if right.
frac_err <- function(f, Lg, N, mu) {
  cc <- f * (Lg - mu) * N / 1e9
  cc <- cc[cc > 0.5 & cc < 30]
  if (length(cc) < 50) return(NA_real_)
  mean(abs(cc - round(cc)))
}
fit_N <- function(f, Lg, mu) {
  grid <- seq(1e6, 60e6, by = 2e4)                  # coarse grid
  err  <- vapply(grid, function(n) frac_err(f, Lg, n, mu), 0)
  # 2N, 3N... also give whole numbers (all counts doubled): take the SMALLEST
  # N whose error is close to the best one, not the global minimum.
  i0 <- which(err < min(err, na.rm = TRUE) + 0.05)[1]
  win <- which(abs(grid / grid[i0] - 1) < 0.02)       # minimum of that basin
  N <- grid[win][which.min(err[win])]
  for (s in c(2e-3, 4e-4, 8e-5)) {                   # refine around the minimum
    gr <- N * (1 + seq(-s, s, length.out = 81))
    N  <- gr[which.min(vapply(gr, function(n) frac_err(f, Lg, n, mu), 0))]
  }
  N
}
use <- niso == 1 & L > 1000

## 4. Common mean fragment length (mu), estimated on sample 1 -------------
f1 <- fpkm[use, 1]; L1 <- L[use]; keep <- f1 > 0
mu_grid <- seq(150, 260, by = 10)
mu_err  <- vapply(mu_grid, function(m) { N <- fit_N(f1[keep], L1[keep], m)
                                         frac_err(f1[keep], L1[keep], N, m) }, 0)
mu <- mu_grid[which.min(mu_err)]
mu_f <- seq(mu - 10, mu + 10, by = 0.5)
mu <- mu_f[which.min(vapply(mu_f, function(m) { N <- fit_N(f1[keep], L1[keep], m)
                                                frac_err(f1[keep], L1[keep], N, m) }, 0))]
cat(sprintf("Mean fragment length (mu): %.1f bp\n", mu))

## 5. Library size per sample -------------------------------------------------
qc <- rbindlist(lapply(seq_len(ncol(fpkm)), function(s) {
  f <- fpkm[use, s]; k <- f > 0
  N <- fit_N(f[k], L[use][k], mu)
  data.table(sample = colnames(fpkm)[s], N_fragments = N,
             int_error = frac_err(f[k], L[use][k], N, mu))
}))

## 6. Counts -------------------------------------------------------------------
Leff <- pmax(L - mu, 1)
q <- fpkm * Leff * rep(qc$N_fragments, each = nrow(fpkm)) / 1e9
counts <- round(q); storage.mode(counts) <- "integer"
dev <- abs(q - round(q)); low <- q > 0.5 & q < 50
qc[, counted_M := colSums(counts) / 1e6]
cat(sprintf("Low counts (<50) that are whole numbers (|x - round(x)| < 0.05): single-isoform %.0f%%, multi-isoform %.0f%%\n",
            100 * mean(dev[niso == 1 & L > 500, ][low[niso == 1 & L > 500, ]] < 0.05),
            100 * mean(dev[niso > 1, ][low[niso > 1, ]] < 0.05)))
print(qc[, .(sample, N_M = round(N_fragments / 1e6, 2), int_error = round(int_error, 3))], nrows = 50)

out <- data.table(gene = rownames(counts), symbol = d$Cre_symbol, defline = d$Cre_defline,
                  length = L, n_isoforms = niso, counts)
fwrite(out, file.path(OUT_DIR, "GSE275433_counts_reconstructed.tsv.gz"), sep = "\t")
fwrite(qc,  file.path(OUT_DIR, "GSE275433_counts_reconstruction_QC.tsv"), sep = "\t")
cat("Written:", file.path(OUT_DIR, "GSE275433_counts_reconstructed.tsv.gz"), "\n")
