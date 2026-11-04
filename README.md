# Same Data, Different Answers: Robust Omics Data Analysis from Reads to Biology

*Every number depends on a choice.* A MOOC based on the computer practicals of
*Systems Molecular Biology* (University of Córdoba), hosted by DTU.

The analysis of omics data, from raw reads and mass-spectrometry reports to
biology: genome assembly, RNA-seq, proteomics, functional enrichment and multi-omics
integration. Each session has two parts:

- **Core**: the computer practicals of the course *Systems Molecular Biology* (University of
  Córdoba), the same materials the students of the degree use.
- **Plant Cell 2025**: a re-analysis of a published experiment, *Chlamydomonas reinhardtii*
  under three light intensities along the day (Dupuis, Ojeda et al. 2025, *The Plant Cell*
  37: koaf086), with transcriptomics, TMT proteomics, enrichment and MOFA+.

Everything runs in your web browser: R, RStudio, all the packages, the command-line tools
and the data. Nothing to install.

## Start (first session, about 10 minutes)

1. Create a free account at [github.com](https://github.com) if you do not have one.
2. **Before creating the codespace:** go to [github.com/settings/codespaces](https://github.com/settings/codespaces)
   → **Default idle timeout** → write **240** minutes → **Save**. Working in RStudio does
   not count as activity for GitHub; with the default (30 minutes) the codespace stops
   in the middle of a session. The setting only applies to codespaces created after it.
3. On this page, click the green **Code** button → **Codespaces** → **Create codespace on main**.
4. Wait while the codespace is created (5–10 minutes the first time: it downloads the data
   of both parts into `~/BMSi`, in the cloud, not on your computer).
5. When the **Terminal** at the bottom shows "RStudio is ready", Ctrl + click on the address
   (`https://<your-codespace>-8787.app.github.dev`). Never use `localhost:8787`.

## Every following session

Go to [github.com/codespaces](https://github.com/codespaces) and open **the same codespace**.
Do not create a new one: your results are kept in the one you already have.

## Where things are

| What | Where |
| --- | --- |
| Core data and scripts | `~/BMSi/PC1_assembly` … `~/BMSi/PC4_Enrichment`, `~/BMSi/scripts` |
| Plant Cell 2025 data | `~/BMSi/PlantCell2025/RNAseq`, `TMT`, `Enrichment`, `MOFA` |
| Plant Cell 2025 scripts | `~/BMSi/scripts/PlantCell2025` |
| Keep your Core results | `bash ~/BMSi/scripts/pack_results.sh PC1` (PC2, PC3, PC4), then Files → More → Export |

Run the Plant Cell scripts in this order: `RNAseq_Chlamy.R`, `Proteomic_Chlamy_TMT.R`,
`Enrichment_Chlamy.R`, `MOFA_Chlamy.R` (MOFA uses the annotation files the enrichment
script writes). Before each script: **Session → Restart R**.

## Free usage

A GitHub account includes 120 core-hours a month; this codespace has 4 cores, so that is
**30 hours of use a month**. Stop the codespace when you finish (github.com/codespaces → **…**
→ **Stop codespace**) and delete it when the course is over.

---

## For the instructor

- **Core is not copied here.** At codespace creation, `.devcontainer/setup_data.sh` clones
  the course repository ([AlfonsoOA/bmsi-practicals](https://github.com/AlfonsoOA/bmsi-practicals))
  and runs its own setup, so both groups always have the same materials.
- **Image:** `.devcontainer/build/Dockerfile` starts from the course image and adds only
  MSstatsTMT, MOFA2 and clue. `.github/workflows/build-image.yml` builds it, publishes it to
  GHCR (`ghcr.io/alfonsooa/bmsi-mooc:2026-27`) and test-runs the Plant Cell scripts at
  codespace size; the results table is in the summary of each run. After the first build,
  make the package public once: GitHub profile → Packages → `bmsi-mooc` → Package settings →
  Change visibility → Public.
- **Data:** release `data-2026-27`, four zips built with `tools/make_MOOC_data.sh` from the
  MOOC folder (raw mzTab and MASIC files are not distributed: the processed PSM table and the
  caches of the long steps are).
- **Scripts:** `scripts/PlantCell2025`. Their data folder is `~/BMSi/PlantCell2025/<block>`
  unless the variables `RNA_BASE`, `PROT_BASE`, `ENRICH_BASE` or `MOFA_BASE` say otherwise.

## Licence and citation

Teaching materials: CC BY 4.0 (LICENSE-CONTENT.md). Code: MIT (LICENSE). How to cite: CITATION.cff.
