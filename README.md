# BAFFR CAR-T Single-Cell Pipeline

Snakemake + R pipeline for processing 10x single-cell (CITE-seq / multimodal) data
from a BAFFR CAR-T construct study, from Cell Ranger output through per-sample QC
to multi-sample integration and downstream analysis.

Part of a personal collection of Snakemake-based bioinformatics pipelines run on
an HPC cluster (Slurm scheduler, `module load` environment modules).

## Pipeline flow

```
cellrange_multi.smk              Cell Ranger multi (GEX + Antibody Capture + CAR reference)
        |
BAFFR_CART.fa / .gtf              custom CAR transgene reference (CAR_BAFFR_1/2, WPRE)
feature_reference_mod.csv         CITE-seq antibody capture feature reference (Anti_EGFR)
        |
        v
BAFFR_singlecell.smk  --include-->  Seurat.smk
        |                              |- createSeuratobj  -> Rscript/CreatSeuratObj_multimodal.R
        |                              '- SeuratQC          -> Rscript/SeuratQC.R (doublet removal via DoubletFinder)
        v
submit_sample_BAFFR.smk           generates one per-sample sbatch script (`smk_code/{sample}.sh`)
        |
        v
Rscript/BAFFR_Data_Integration.R  multi-sample merge/integration (TCR, RPCA/harmony), QC plots
        |
        v
Rscript/BAFFR_post_integration.R  downstream analysis on the integrated object
                                   (cluster markers, pseudobulk DE, score/dotplots)
```

## Files

| File | Purpose |
|---|---|
| `cellrange_multi.smk` | Snakemake rule running `cellranger multi` per sample (GEX + Antibody Capture + VDJ/custom reference). |
| `BAFFR_CART.fa`, `BAFFR_CART.gtf` | Custom transgene sequences (`CAR_BAFFR_1`, `CAR_BAFFR_2`, `WPRE`) appended to the reference genome so Cell Ranger can quantify CAR/WPRE expression. |
| `feature_reference_mod.csv` | CITE-seq Antibody Capture feature reference (Anti_EGFR) passed to `cellranger multi`. |
| `BAFFR_singlecell.smk` | Main per-sample Snakemake pipeline; includes `Seurat.smk` and targets `SeuratQC.bmk` per sample. |
| `Seurat.smk` | Rules `createSeuratobj` (build Seurat object from the Cell Ranger multimodal matrix) and `SeuratQC` (doublet removal + QC). |
| `Rscript/CreatSeuratObj_multimodal.R` | Builds the per-sample Seurat object from the `sample_filtered_feature_bc_matrix`, separating Gene Expression / Antibody Capture modalities. |
| `Rscript/SeuratQC.R` | Removes doublets (DoubletFinder classifications) and recomputes QC metrics. |
| `submit_sample_BAFFR.smk` | Per-sample job dispatcher: generates one sbatch script per sample under `{dir_out}/smk_code/{sample}.sh`, each invoking `BAFFR_singlecell.smk` for that sample. |
| `Rscript/BAFFR_Data_Integration.R` | Merges/integrates per-sample Seurat objects across the cohort (TCR annotation via scRepertoire, RPCA/harmony integration), with QC boxplots. |
| `Rscript/BAFFR_post_integration.R` | Downstream analysis on the integrated object: cluster markers, pseudobulk DESeq2 comparisons, signature score dot plots. |
| `Rscript/F_sc.R` | General single-cell helper library: Seurat object creation, QC metrics and cutoffs, DoubletFinder wrapper, SCTransform clustering, SingleR annotation (single/double reference), marker-density cutoff helpers, and plotting/TSV I/O utilities. |
| `Rscript/CART_singlecell_functions.R` | CAR-T analysis helper library: T-cell signature gene sets, ProjecTILs reference mapping, integration and resolution selection, CAR-positive and group/sample proportion plots, TCR clonotype (CDR3, top/shared clone) comparisons, pseudobulk DE, pseudotime analysis, signature scoring, heatmaps, and GSEA. |

## Helper libraries

`F_sc.R` and `CART_singlecell_functions.R` define functions only and are loaded with `source()`:

- `F_sc.R` is sourced by `CreatSeuratObj_multimodal.R`, `SeuratQC.R`, and `BAFFR_Data_Integration.R`.
- `CART_singlecell_functions.R` is sourced by `BAFFR_Data_Integration.R`.

The `source()` calls and some reference files inside these libraries (for example, the housekeeping
gene list in `F_sc.R` and the CD8 ProjecTILs reference in `CART_singlecell_functions.R`) use absolute
cluster paths. Update these paths before running the scripts in a different environment.

## Running

```bash
# 1. Generate one sbatch script per sample
snakemake -s submit_sample_BAFFR.smk -j1

# 2. Submit the generated per-sample scripts
sbatch <dir_out>/smk_code/<sample>.sh
```

`Rscript/BAFFR_Data_Integration.R` and `Rscript/BAFFR_post_integration.R` are run manually (not via Snakemake)
once all per-sample QC objects are ready, and expect per-project working directories/paths
to be set near the top of each script.

## Requirements

- Snakemake (run on HPC via `module load`; no per-rule conda envs in this pipeline)
- Cell Ranger (`module load cellranger`)
- R (`module load R/RStudio_R-4.4.1`) with Seurat, DoubletFinder, scRepertoire, DESeq2, and
  the other packages loaded at the top of each script
