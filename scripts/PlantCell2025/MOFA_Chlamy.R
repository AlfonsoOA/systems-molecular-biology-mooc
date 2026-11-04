# =============================================================================
# MOFA_Chlamy.R
# Integration of RNA and protein in the light-acclimation experiment of
# Chlamydomonas (Dupuis, Ojeda et al. 2025, Plant Cell 37: koaf086) with MOFA+.
#
# Three questions, none of them about single genes:
#   Q1. Without being told the design, does MOFA find the time of day and the light?
#   Q2. Which variation is shared by RNA and protein, and which belongs to one layer?
#   Q3. Where do RNA and protein disagree (genes and functions)? The article reports
#       that only ~11% of the light-responsive changes occur at both levels.
#
# Inputs (folder InputData, outputs of the previous blocks):
#   TranscriptomicData/ GSE275433_counts_reconstructed.tsv, sample_design.tsv, relevance_all.tsv
#   ProteomicData/      A_protein_log2ratio.tsv (TMT pipeline A), B_MSstatsTMT.rds (plex of each
#                       sample), relevance_all.tsv
#   Annotation/         Annotation_eggNOG_KEGG.gmt, Annotation_eggNOG_GO_BP.gmt (enrichment block)
# One input per layer, chosen in the benchmark block: limma-voom-style log CPM for RNA and
# TMT pipeline A log2 ratios, plex removed, for protein.
# =============================================================================

# ---- 0. Settings ------------------------------------------------------------------
BASE  <- Sys.getenv("MOFA_BASE", path.expand("~/BMSi/PlantCell2025/MOFA"))
INPUT <- file.path(BASE, "InputData")

N_FACTORS    <- 10      # factors to fit (factors explaining < 1% in both layers are reported, not used)
MAX_NA_PROT  <- 0.2     # a protein may miss at most 20% of the 45 samples
SCALE_VIEWS  <- TRUE    # give RNA and protein the same total variance (else RNA dominates)
SEED         <- 2026
USE_BASILISK <- as.logical(Sys.getenv("MOFA_BASILISK", "TRUE"))  # TRUE: MOFA2 installs its own Python (first run ~10 min). FALSE: use the
                        # Python that reticulate finds (needs 'pip install mofapy2')
FORCE_TRAIN  <- FALSE   # retrain even if the model file exists (~2-5 min)
FDR          <- 0.05

# ---- 1. Packages ------------------------------------------------------------------
cran <- c("data.table", "ggplot2", "ggrepel")
bioc <- c("limma", "edgeR", "MOFA2", "fgsea")
for (p in cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(limma); library(edgeR); library(MOFA2) })
theme_set(theme_bw(base_size = 11))
set.seed(SEED)
D <- function(...) { p <- file.path(BASE, "results_MOFA", ...); dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE); p }
banner <- function(x) cat("\n==== ", x, " ", strrep("=", max(3, 70 - nchar(x))), "\n", sep = "")
inp <- function(...) { f <- file.path(INPUT, ...); if (!file.exists(f)) stop("Missing input file: ", f); f }
TIMES <- c("m02h", "p02h", "p06h", "p10h", "m10h"); LIGHTS <- c("LL", "ML", "HL")

# =================================================================================
# 2. Data and sample pairing
# =================================================================================
banner("2. Data and sample pairing")
cn <- fread(inp("TranscriptomicData", "GSE275433_counts_reconstructed.tsv"))
counts <- as.matrix(cn[, -(1:5)]); rownames(counts) <- cn$gene
smp <- fread(inp("TranscriptomicData", "sample_design.tsv"))
smp[, `:=`(Light = factor(Light, LIGHTS), Time = factor(Time, TIMES))]
counts <- counts[, smp$SampleID]
keep <- filterByExpr(counts, group = smp$State)
rna <- cpm(calcNormFactors(DGEList(counts[keep, ])), log = TRUE, prior.count = 2)

prot <- as.matrix(fread(inp("ProteomicData", "A_protein_log2ratio.tsv")), rownames = 1)
msB <- readRDS(inp("ProteomicData", "B_MSstatsTMT.rds"))
pm <- unique(as.data.table(msB$protein)[, .(Prot_ID = as.character(BioReplicate), Plex = as.character(Mixture))])
rm(msB)
pm <- pm[Prot_ID %in% colnames(prot)]
pm[, Light := c(`0050` = "LL", `0200` = "ML", `1000` = "HL")[sub("^light(\\d+)_.*", "\\1", Prot_ID)]]
pm[, Time := c(CTm02 = "m02h", CTp02 = "p02h", CTp06 = "p06h", CTp10 = "p10h", CTm10 = "m10h")[sub(".*_(CT...)_\\d$", "\\1", Prot_ID)]]
pm[, Rep := paste0("BR", as.integer(sub(".*_(\\d)$", "\\1", Prot_ID)) + 1)]   # TMT replicate 0,1,2 -> BR1,BR2,BR3
pm[, State := paste(Light, Time, sep = "_")]

# Are the RNA and protein samples of a state the same cultures? If so, the deviation of
# each culture from its state mean should correlate between the two layers for the right
# pairing. Test the 6 possible pairings in each of the 15 states.
shared0 <- intersect(rownames(rna), rownames(prot)); shared0 <- shared0[rowSums(is.na(prot[shared0, ])) == 0]
resid <- function(X, st) { for (s in unique(st)) X[, st == s] <- X[, st == s] - rowMeans(X[, st == s, drop = FALSE]); X }
Rr <- resid(rna[shared0, smp$SampleID], smp$State); Pr <- resid(prot[shared0, pm$Prot_ID], pm$State)
perms <- list(c(1, 2, 3), c(1, 3, 2), c(2, 1, 3), c(2, 3, 1), c(3, 1, 2), c(3, 2, 1))
pairing <- rbindlist(lapply(unique(smp$State), function(s) {
  a <- smp[State == s][order(Rep)]; b <- pm[State == s][order(Rep)]
  cc <- cor(Rr[, a$SampleID], Pr[, b$Prot_ID])
  sc <- sapply(perms, function(p) sum(cc[cbind(1:3, p)]))
  data.table(State = s, identity_r = round(mean(diag(cc)), 3), best = paste(perms[[which.max(sc)]], collapse = ""),
             identity_is_best = which.max(sc) == 1) }))
fwrite(pairing, D("1_Data", "pairing_check.tsv"), sep = "\t")
n_id <- sum(pairing$identity_is_best)
cat(sprintf("Pairing BRk (RNA) = replicate k-1 (TMT) is the best of 6 in %d of 15 states (by chance ~2.5; binomial p = %.1e)\n",
            n_id, pbinom(n_id - 1, 15, 1 / 6, lower.tail = FALSE)))
if (n_id < 6) warning("Weak support for the pairing: consider integrating state means instead")

# paired matrices with the RNA sample names; plex removed from protein, keeping the state means
pm <- pm[match(paste(smp$State, smp$Rep), paste(State, Rep))]
P <- prot[, pm$Prot_ID]; colnames(P) <- smp$SampleID
P <- removeBatchEffect(P, batch = pm$Plex, design = model.matrix(~ 0 + State, data = smp))
P <- P[rowMeans(is.na(P)) <= MAX_NA_PROT, ]
genes <- intersect(rownames(rna), rownames(P))    # the same genes in both layers: weights comparable gene by gene
R <- rna[genes, ]; P <- P[genes, ]
cat(sprintf("Samples: %d paired | genes in both layers: %d (RNA expressed: %d, proteins with <= %d%% missing: %d)\n",
            ncol(R), length(genes), nrow(rna), round(100 * MAX_NA_PROT), sum(rowMeans(is.na(prot)) <= MAX_NA_PROT)))

# =================================================================================
# 3. MOFA+ model
# =================================================================================
banner("3. MOFA+")
MODEL_FILE <- D("2_Model", "MOFA_Chlamy.hdf5")
if (FORCE_TRAIN || !file.exists(MODEL_FILE)) {
  views <- list(RNA = R, Protein = P)
  for (v in names(views)) rownames(views[[v]]) <- paste0(rownames(views[[v]]), "_", v)   # feature names unique across views
  mo <- create_mofa(views)
  dopt <- get_default_data_options(mo); dopt$scale_views <- SCALE_VIEWS
  mopt <- get_default_model_options(mo); mopt$num_factors <- N_FACTORS
  topt <- get_default_training_options(mo); topt$seed <- SEED; topt$convergence_mode <- "medium"
  topt$drop_factor_threshold <- -1; topt$verbose <- FALSE
  mo <- prepare_mofa(mo, data_options = dopt, model_options = mopt, training_options = topt)
  t0 <- Sys.time()
  model <- run_mofa(mo, outfile = MODEL_FILE, use_basilisk = USE_BASILISK, save_data = TRUE)
  cat(sprintf("Trained in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
} else model <- load_model(MODEL_FILE)
meta <- copy(smp)[, sample := SampleID][, Plex := pm$Plex]
samples_metadata(model) <- as.data.frame(meta[match(samples_names(model)[[1]], sample)])

ve <- calculate_variance_explained(model)
r2f <- as.data.table(ve$r2_per_factor[[1]], keep.rownames = "Factor")
r2t <- ve$r2_total[[1]]
cat("Variance explained in total (%):\n"); print(round(r2t, 1))
use <- r2f[RNA >= 1 | Protein >= 1, Factor]
cat(sprintf("Factors explaining >= 1%% in at least one layer: %d of %d\n", length(use), nrow(r2f)))

# =================================================================================
# 4. Q1: what do the factors capture?
# =================================================================================
banner("4. Q1: factors vs time, light and technical variables")
Fm <- get_factors(model)[[1]]
md <- meta[match(rownames(Fm), sample)]
assoc <- rbindlist(lapply(colnames(Fm), function(f) {
  y <- Fm[, f]
  a <- anova(lm(y ~ Time * Light + Rep + Plex, data = cbind(md, y = y)))
  ss <- setNames(a$`Sum Sq`, trimws(rownames(a))); tot <- sum(ss)
  data.table(Factor = f, Time = ss["Time"] / tot, Light = ss["Light"] / tot, `Time:Light` = ss["Time:Light"] / tot,
             Replicate = ss["Rep"] / tot, Plex = ss["Plex"] / tot, Residual = ss["Residuals"] / tot) }))
assoc <- merge(r2f, assoc, by = "Factor")[order(-(RNA + Protein))]
assoc[, Label := fifelse(Time >= 0.5, "time of day", fifelse(Light >= 0.5, "light", fifelse(`Time:Light` >= 0.3, "light x time",
                 fifelse(Time + Light + `Time:Light` >= 0.5, "time + light", "other"))))]
num <- names(assoc)[sapply(assoc, is.numeric)]; assoc[, (num) := lapply(.SD, round, 3), .SDcols = num]
fwrite(assoc, D("3_Factors", "factor_variance_and_design.tsv"), sep = "\t")
cat("Variance explained per layer (%) and share of each factor's values explained by the design:\n"); print(assoc)

LIGHT_F <- assoc[Factor %in% use][which.max(Light), Factor]                      # the factor most explained by light
TIME_F  <- assoc[Factor %in% use & Label == "time of day"][order(-(RNA + Protein)), head(Factor, 2)]   # the two largest time factors
cat(sprintf("Time factors: %s | light factor: %s\n", paste(TIME_F, collapse = ", "), LIGHT_F))
fd <- cbind(md, as.data.table(Fm))
tcol <- c(m02h = "#14213D", p02h = "#E8A87C", p06h = "#B5541B", p10h = "#7A3A16", m10h = "#2A6F97")
g <- ggplot(fd, aes(.data[[TIME_F[1]]], .data[[TIME_F[2]]], colour = Time, shape = Light)) + geom_point(size = 3) +
  scale_colour_manual(values = tcol) + labs(title = "The two largest factors draw the clock", x = TIME_F[1], y = TIME_F[2])
ggsave(D("3_Factors", "factors_clock.png"), g, width = 7, height = 5.5)
g <- ggplot(fd, aes(Time, .data[[LIGHT_F]], colour = Light)) + geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c(LL = "#8A96A3", ML = "#E8A87C", HL = "#B5541B")) +
  labs(title = sprintf("The light factor (%s) separates the light levels at every time", LIGHT_F), x = "time", y = LIGHT_F)
ggsave(D("3_Factors", "factor_light.png"), g, width = 7, height = 5)
cat(sprintf("Mean of the light factor by light level: %s\n", paste(fd[, .(m = round(mean(get(LIGHT_F)), 2)), by = Light][, paste(Light, m)], collapse = ", ")))
r2l <- melt(r2f, id.vars = "Factor", variable.name = "Layer", value.name = "R2")
r2l[, Factor := factor(Factor, levels = rev(r2f$Factor))]
g <- ggplot(r2l, aes(Layer, Factor, fill = R2)) + geom_tile(colour = "white") + geom_text(aes(label = round(R2, 1)), size = 3.5) +
  scale_fill_gradient(low = "#F1F4F7", high = "#2A6F97", name = "% variance") + labs(x = NULL, y = NULL, title = "Variance explained per factor and layer")
ggsave(D("3_Factors", "variance_explained.png"), g, width = 4.5, height = 5)

# =================================================================================
# 5. Q2-Q3: shared and layer-specific signals, genes and functions
# =================================================================================
banner("5. Q2-Q3: weights in RNA and protein")
W <- get_weights(model, scale = TRUE)
for (v in names(W)) rownames(W[[v]]) <- sub(paste0("_", v, "$"), "", rownames(W[[v]]))
wr <- rbindlist(lapply(use, function(f) data.table(Factor = f, Gene = genes, RNA = W$RNA[genes, f], Protein = W$Protein[genes, f])))
wcor <- wr[, .(r_weights = round(cor(RNA, Protein), 2)), by = Factor]
print(merge(assoc[Factor %in% use, .(Factor, RNA, Protein, Label)], wcor, by = "Factor")[order(-(RNA + Protein))])
fwrite(wr, D("4_Weights", "weights_RNA_protein.tsv"), sep = "\t")

# genes that load on the light factor in one layer only (top 2% of |weight| in one, bottom half in the other)
sym <- unique(fread(inp("TranscriptomicData", "relevance_all.tsv"), select = c("Gene", "symbol", "defline")))
lw <- merge(wr[Factor == LIGHT_F], sym, by = "Gene", all.x = TRUE)
lw[, `:=`(qR = frank(abs(RNA)) / .N, qP = frank(abs(Protein)) / .N)]
lw[, Class := fifelse(qR >= 0.98 & qP >= 0.98, "both", fifelse(qP >= 0.98 & qR < 0.5, "protein only", fifelse(qR >= 0.98 & qP < 0.5, "RNA only", "")))]
fwrite(lw[Class != ""][order(Class, -abs(Protein) - abs(RNA))], D("4_Weights", sprintf("light_factor_%s_discordant_genes.tsv", LIGHT_F)), sep = "\t")
cat("Light factor, top-2% genes by layer:\n"); print(lw[Class != "", .N, by = Class])
lab <- lw[Class != "" & !is.na(symbol) & symbol != ""][order(-pmax(abs(RNA), abs(Protein)))][, head(.SD, 6), by = Class]
g <- ggplot(lw, aes(RNA, Protein)) + geom_point(colour = "grey75", size = 0.6) +
  geom_point(data = lw[Class != ""], aes(colour = Class), size = 1.2) +
  ggrepel::geom_text_repel(data = lab, aes(label = symbol), size = 3, max.overlaps = 30) +
  scale_colour_manual(values = c(both = "#14213D", `protein only` = "#B5541B", `RNA only` = "#2A6F97")) +
  labs(title = sprintf("Light factor (%s): weight of each gene in RNA and in protein", LIGHT_F), x = "RNA weight (scaled)", y = "protein weight (scaled)", colour = NULL)
ggsave(D("4_Weights", sprintf("light_factor_%s_weights.png", LIGHT_F)), g, width = 7.5, height = 6)

# functions: GSEA on the weights of each layer, for each factor in use
read_gmt <- function(f) { l <- strsplit(readLines(f), "\t"); setNames(lapply(l, `[`, -(1:2)), paste0(sapply(l, `[`, 1), "|", sapply(l, `[`, 2))) }
sets <- c(lapply(read_gmt(inp("Annotation", "Annotation_eggNOG_KEGG.gmt")), identity), read_gmt(inp("Annotation", "Annotation_eggNOG_GO_BP.gmt")))
sets <- lapply(sets, intersect, genes); sets <- sets[lengths(sets) >= 10 & lengths(sets) <= 500]
gs <- rbindlist(lapply(use, function(f) rbindlist(lapply(c("RNA", "Protein"), function(v) {
  st <- sort(setNames(W[[v]][genes, f], genes), decreasing = TRUE)
  r <- suppressWarnings(fgsea::fgseaMultilevel(sets, st, minSize = 10, maxSize = 500, eps = 0, nPermSimple = 2000))
  data.table(Factor = f, Layer = v, Term = r$pathway, NES = r$NES, padj = r$padj) }))))
gs[, c("ID", "Name") := tstrsplit(Term, "|", fixed = TRUE)]
gw <- dcast(gs, Factor + ID + Name ~ Layer, value.var = c("NES", "padj"))
gw[, Class := fifelse(padj_RNA < FDR & padj_Protein < FDR & sign(NES_RNA) == sign(NES_Protein), "both, same sign",
              fifelse(padj_RNA < FDR & padj_Protein < FDR, "both, opposite sign",
              fifelse(padj_Protein < FDR, "protein only", fifelse(padj_RNA < FDR, "RNA only", "neither"))))]
fwrite(gw[order(Factor, Class, padj_Protein)], D("5_Functions", "GSEA_weights_by_layer.tsv"), sep = "\t")
gsum <- gw[Class != "neither", .N, by = .(Factor, Class)]
cat("Functions loading on each factor (GSEA on weights, eggNOG KEGG + GO_BP):\n"); print(dcast(gsum, Factor ~ Class, value.var = "N", fill = 0))
kg <- gw[Factor == LIGHT_F & grepl("^map", ID) & Class != "neither"]
cat(sprintf("KEGG pathways on the light factor (%s):\n", LIGHT_F)); print(kg[order(Class, -abs(NES_Protein)), .(Name, NES_RNA = round(NES_RNA, 2), NES_Protein = round(NES_Protein, 2), Class)])

# =================================================================================
# 6. The article's figure: how many light responses occur at both levels?
# =================================================================================
banner("6. RNA vs protein light effects (pipelines A)")
rr <- fread(inp("TranscriptomicData", "relevance_all.tsv"))[grepl("^A_", Pipeline), .(Gene, Contrast, lfc_R = log2FC, rel_R = FC_relevant)]
pp <- fread(inp("ProteomicData", "relevance_all.tsv"))[grepl("^A_", Pipeline), .(Gene = Protein, Contrast, lfc_P = log2FC, rel_P = FC_relevant)]
rp <- merge(rr, pp, by = c("Gene", "Contrast"))
both <- rp[, .(genes = .N, r_log2FC = round(cor(lfc_R, lfc_P, use = "complete.obs"), 2),
               relevant_either = sum(rel_R | rel_P), relevant_both = sum(rel_R & rel_P),
               pct_both = round(100 * sum(rel_R & rel_P) / sum(rel_R | rel_P), 1),
               same_sign_when_both = round(100 * mean(sign(lfc_R[rel_R & rel_P]) == sign(lfc_P[rel_R & rel_P])), 1)), by = Contrast]
fwrite(both, D("6_Article", "RNA_vs_protein_light_effects.tsv"), sep = "\t")
print(both)
cat("\nDone. Results in", D(""), "\n")
writeLines(capture.output(sessionInfo()), D("sessionInfo.txt"))
