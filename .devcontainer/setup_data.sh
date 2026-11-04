#!/bin/bash
# Runs once, when a learner's codespace is created.
#  1. Core: the data and scripts of the course sessions, taken from the course repository
#     itself (AlfonsoOA/bmsi-practicals), so the MOOC always has the same materials.
#  2. Plant Cell 2025: the data of the re-analysis, from the release "data-2026-27" of
#     this repository, and its scripts.
# Everything goes to ~/BMSi. Files already there are never overwritten.
REPO=$(cd "$(dirname "$0")/.." && pwd)
COURSE=https://github.com/AlfonsoOA/bmsi-practicals
RELEASE=https://github.com/AlfonsoOA/systems-molecular-biology-mooc/releases/download/data-2026-27

# 1. Core
rm -rf /tmp/bmsi-practicals
if git clone --depth 1 -q "$COURSE" /tmp/bmsi-practicals; then
  bash /tmp/bmsi-practicals/.devcontainer/setup_data.sh
else
  echo "WARNING: the course repository could not be reached; run this script again later: bash $0"
fi

# 2. Plant Cell 2025
PC=~/BMSi/PlantCell2025; mkdir -p $PC
for m in RNAseq TMT Enrichment MOFA; do
  if [ -d "$PC/$m" ]; then echo "PlantCell2025/$m already there"; continue; fi
  z=MOOC_PlantCell_$m.zip
  echo "Downloading $z..."
  if curl -fsSL --retry 3 -o "$PC/$z" "$RELEASE/$z"; then
    mkdir -p "$PC/$m" && unzip -n -q "$PC/$z" -d "$PC/$m" && rm -f "$PC/$z"
  else
    rm -f "$PC/$z"; echo "WARNING: $z could not be downloaded; run this script again later: bash $0"
  fi
done
# MOFA reads the outputs of the RNA-seq, TMT and enrichment blocks: its InputData points to them
I=$PC/MOFA/InputData; mkdir -p $I/TranscriptomicData $I/ProteomicData
E=$PC/Enrichment/InputData
for f in GSE275433_counts_reconstructed.tsv sample_design.tsv relevance_all.tsv; do
  ln -sfn $E/TranscriptomicData/$f $I/TranscriptomicData/$f; done
ln -sfn $E/ProteomicData/2_A_limma/A_protein_log2ratio.tsv $I/ProteomicData/A_protein_log2ratio.tsv
ln -sfn $E/ProteomicData/3_B_MSstatsTMT/B_MSstatsTMT.rds   $I/ProteomicData/B_MSstatsTMT.rds
ln -sfn $E/ProteomicData/4_Relevance/relevance_all.tsv     $I/ProteomicData/relevance_all.tsv
ln -sfn $PC/Enrichment/results_Enrichment/7_ClueGO         $I/Annotation   # written by Enrichment_Chlamy.R

mkdir -p ~/BMSi/scripts/PlantCell2025
cp -n "$REPO"/scripts/PlantCell2025/*.R ~/BMSi/scripts/PlantCell2025/ 2>/dev/null || true
echo "Data and scripts are in ~/BMSi:"; ls ~/BMSi ~/BMSi/PlantCell2025
