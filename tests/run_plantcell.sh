#!/bin/bash
# Test run of the Plant Cell 2025 scripts inside the MOOC image, as a learner would run
# them (heavy steps taken from the caches of the release). Prints a markdown table with
# exit status, time and peak memory. The Core sessions are tested in the course repository.
REPO=$(pwd); LOGS=$REPO/test_logs; mkdir -p $LOGS
bash .devcontainer/setup_data.sh > $LOGS/setup.log 2>&1
run() {   # name, script
  /usr/bin/time -v -o $LOGS/$1.time Rscript "$2" > $LOGS/$1.log 2>&1
  local st=$?
  local mem=$(awk -F': ' '/Maximum resident/ {printf "%.2f", $2/1048576}' $LOGS/$1.time)
  local t=$(awk -F': ' '/Elapsed/ {print $2}' $LOGS/$1.time)
  echo "| $1 | $([ $st -eq 0 ] && echo OK || echo "FAILED ($st)") | $t | $mem GB |"
}
S=~/BMSi/scripts/PlantCell2025
echo "| Script | Result | Time (h:mm:ss) | Peak memory |"
echo "| --- | --- | --- | --- |"
run RNAseq_Chlamy     $S/RNAseq_Chlamy.R
run Proteomic_TMT     $S/Proteomic_Chlamy_TMT.R
run Enrichment_Chlamy $S/Enrichment_Chlamy.R
run MOFA_Chlamy       $S/MOFA_Chlamy.R
for d in RNAseq TMT Enrichment MOFA; do
  if [ -d ~/BMSi/PlantCell2025/$d ]; then echo "| data PlantCell2025/$d | OK | | $(du -shL ~/BMSi/PlantCell2025/$d | cut -f1) |"; else echo "| data PlantCell2025/$d | MISSING | | |"; fi
done
for d in PC1_assembly PC2_RNAseq PC3_Proteomics PC4_Enrichment; do
  [ -d ~/BMSi/$d ] && echo "| Core data $d | OK | | |" || echo "| Core data $d | MISSING | | |"
done

# Caches of the slow steps of the enrichment (GSEA, resampling), to add to
# MOOC_PlantCell_Enrichment.zip so that learners do not recompute them (~35 min)
mkdir -p $LOGS/enrichment_caches
( cd ~/BMSi/PlantCell2025/Enrichment && find results_Enrichment -name "*.rds" -path "*3_Enrichment*" -o -name "*.rds" -path "*6_Resampling/cache*" | \
  zip -q $LOGS/enrichment_caches/Enrichment_caches.zip -@ ) 2>/dev/null
