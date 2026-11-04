#!/bin/bash
# Instructor side, run once in WSL: builds the four data packages of the Plant Cell 2025
# re-analysis from the MOOC folder, to upload to the release "data-2026-27" of this
# repository. Copies only; the MOOC folder is not modified.
#   bash make_MOOC_data.sh
set -e
M=/mnt/c/Users/Alfonso/Desktop/AlfonsoOA_MSI/2Docencia/15_MOOC
R=$M/PC2_RNAseq/PlantCell2025/Transcriptomic_Chlamy
P=$M/PC3_Proteomics/PlantCell2025/Proteomic_Chlamy
E=$M/PC4_Enrichment_MOFA/PlantCell2025/Enrichment
F=$M/PC4_Enrichment_MOFA/PlantCell2025/MOFA
OUT=$M/release_assets; mkdir -p "$OUT"
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT

echo "1/4 RNA-seq"
d=$T/RNAseq; mkdir -p $d/processed
cp "$R"/GSE275433_Cr073.light_intensity.FPKMs.replicates.v6.1.tsv.gz "$R"/CreinhardtiiCC_4532_707_v6.1.gene_exons.gff3.gz $d/
cp "$R"/processed/GSE275433_counts_reconstructed.tsv.gz "$R"/processed/GSE275433_counts_reconstruction_QC.tsv $d/processed/
python3 -c "import zipfile,sys; zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2])" "$R/RNAseq_results_and_caches.zip" $d

echo "2/4 TMT"
d=$T/TMT; mkdir -p $d/processed $d/supplementary $d/results_TMT/3_B_MSstatsTMT $d/results_TMT/4_Relevance
cp "$P"/processed/psm_tmt.tsv.gz "$P"/processed/channel_design.tsv "$P"/processed/qc_files.tsv $d/processed/
cp "$P"/supplementary/*.xlsx $d/supplementary/
cp "$E"/InputData/ProteomicData/3_B_MSstatsTMT/B_MSstatsTMT.rds $d/results_TMT/3_B_MSstatsTMT/
cp "$E"/InputData/ProteomicData/4_Relevance/Bayes_A.rds "$E"/InputData/ProteomicData/4_Relevance/Bayes_B.rds $d/results_TMT/4_Relevance/

echo "3/4 Enrichment"
d=$T/Enrichment; mkdir -p $d/results_Enrichment/6_Resampling/cache
cp -r "$E"/InputData $d/
rm -f $d/InputData/ProteomicData/3_B_MSstatsTMT/B_proteinSummarization.rds   # 171 MB, not read by the enrichment
python3 -c "import zipfile,sys; zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2])" "$E/Enrichment_bootstrap_caches.zip" $d/results_Enrichment/6_Resampling/cache

echo "4/4 MOFA"
d=$T/MOFA; mkdir -p $d/results_MOFA/2_Model
cp "$F"/MOFA_Chlamy.hdf5 $d/results_MOFA/2_Model/

for m in RNAseq TMT Enrichment MOFA; do
  python3 - "$T/$m" "$OUT/MOOC_PlantCell_$m.zip" <<'PY'
import os, sys, zipfile
src, out = sys.argv[1], sys.argv[2]
with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, allowZip64=True) as z:
    for root, _, files in os.walk(src):
        for f in files:
            p = os.path.join(root, f); z.write(p, os.path.relpath(p, src))
PY
  echo "  $(ls -l "$OUT/MOOC_PlantCell_$m.zip" | awk '{printf "%8.1f MB", $5/1048576}')  MOOC_PlantCell_$m.zip"
done
echo "Done: upload the four zips in $OUT to the release data-2026-27 of AlfonsoOA/systems-molecular-biology-mooc."
