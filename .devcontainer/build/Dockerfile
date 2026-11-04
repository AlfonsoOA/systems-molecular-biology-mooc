# Image of the MOOC: the image of the course (same R, packages, conda environments and
# RStudio Server, so the Core sessions run exactly as in the course) plus the packages
# that only the Plant Cell 2025 re-analysis needs.
FROM ghcr.io/alfonsooa/bmsi-practicals:2026-27
USER root
RUN Rscript -e 'BiocManager::install(c("MSstatsTMT", "MOFA2", "clue"), ask = FALSE, update = FALSE); \
      m <- setdiff(c("MSstatsTMT", "MOFA2", "clue"), rownames(installed.packages())); \
      if (length(m)) stop("Not installed: ", paste(m, collapse = ", ")); \
      for (p in c("MSstatsTMT", "MOFA2", "clue")) suppressPackageStartupMessages(library(p, character.only = TRUE))' \
 && rm -rf /tmp/* /root/.cache
