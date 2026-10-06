# Spatial-SCB

This repository contains the custom R scripts used for the spatial transcriptomics analysis of subchondral bone in this study. The workflow includes Visium HD/Visium preprocessing, single-cell reference re-annotation, spatial cell-state identification, RCTD-based cell type deconvolution, core-halo region identification, regional functional analysis, and spatial cell-cell communication analysis.

## Repository structure

```text
Spatial-SCB/
├── README.md
├── LICENSE
└── src/
    ├── 01_data_processing.R
    ├── 02_sc_reference_reannotated.R
    ├── 03_extract_oc_ad_bins_16um.R
    ├── 04_RCTD.R
    ├── 05_Core_identification.R
    ├── 06_Area_description.R
    ├── 07_Spatial_ccc.R
    ├── 08_Inter-ecotype_communication.R
    └── utils.R
```

All scripts are stored in the `src/` folder. The scripts are designed to be run in numerical order, with shared helper functions stored in `utils.R`.

## Workflow overview

| Step | Script | Purpose |
|---|---|---|
| 1 | `01_data_processing.R` | Load and preprocess Visium HD/Visium spatial transcriptomics data, add sample metadata, filter spatial bins/spots, normalize data, perform dimensionality reduction and integration, and generate initial spatial visualizations. |
| 2 | `02_sc_reference_reannotated.R` | Re-annotate the single-cell bone marrow reference, assign broad and refined cell type labels, validate marker expression, and save the refined reference object. |
| 3 | `03_extract_oc_ad_bins_16um.R` | Identify osteoclast- and adipocyte-enriched spatial bins using marker-based module scores, generate spatial plots, differential expression results, and GO enrichment summaries. |
| 4 | `04_RCTD.R` | Build the RCTD reference and spatial query objects, run RCTD deconvolution, normalize cell type weights, add the RCTD assay to the spatial Seurat object, and visualize cell type proportions. |
| 5 | `05_Core_identification.R` | Define osteoblast-enriched core regions based on RCTD-estimated osteoblast proportions, remove low-connectivity regions, define halo regions around core areas, and save the core-halo annotated object. |
| 6 | `06_Area_description.R` | Compare Core, Halo, and Other regions using differential expression, GO enrichment, GSVA/limma pathway analysis, and cell type composition visualization. |
| 7 | `07_Spatial_ccc.R` | Perform spatial ligand-receptor co-localization analysis and infer sender-receiver cell type pairs for spatial cell-cell communication. |
| 8 | `08_Inter-ecotype_communication.R` | Analyze inter-ecotype communication using spatial ecotype labels and CellChat-based communication inference. |
| Helper | `utils.R` | Contains shared functions for data loading, plotting, RCTD preparation, core-halo visualization, enrichment analysis, ligand-receptor co-localization, and cell type attribution. |

## Software requirements

The analysis was performed in R. The following R packages are required or used by at least one script:

```r
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Matrix)
library(UCell)
library(spacexr)
library(msigdbr)
library(GSVA)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(randomcoloR)
library(scRNAtoolVis)
library(gghalves)
library(ggbreak)
library(ggprism)
library(ggpubr)
library(RColorBrewer)
library(cols4all)
library(FNN)
library(spdep)
library(igraph)
library(arrow)
library(furrr)
library(future.apply)
library(metap)
library(UpSetR)
library(glmnet)
library(CellChat)
library(liana)
library(CrossTalkeR)
library(jsonlite)
library(tibble)
library(purrr)
```

Some packages are available from CRAN, some from Bioconductor, and some from GitHub. Install them before running the scripts.

Example installation commands:

```r
install.packages(c(
  "Seurat", "dplyr", "tidyverse", "ggplot2", "ggrepel", "patchwork",
  "Matrix", "UCell", "msigdbr", "GSVA", "limma", "randomcoloR",
  "gghalves", "ggbreak", "ggprism", "ggpubr", "RColorBrewer",
  "cols4all", "FNN", "spdep", "igraph", "arrow", "furrr",
  "future.apply", "metap", "UpSetR", "glmnet", "jsonlite", "tibble", "purrr"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "clusterProfiler", "org.Hs.eg.db", "ComplexHeatmap", "circlize"
))
```

For `spacexr`, `CellChat`, `liana`, `CrossTalkeR`, and `scRNAtoolVis`, please follow the installation instructions from the corresponding package repositories.

---

## How to run the workflow

Because the scripts call `source("utils.R")`, run them from inside the `src/` folder:

```bash
cd src
```

Then run the scripts in order:

```bash
Rscript 01_data_processing.R
Rscript 02_sc_reference_reannotated.R
Rscript 03_extract_oc_ad_bins_16um.R
Rscript 04_RCTD.R
Rscript 05_Core_identification.R
Rscript 06_Area_description.R
Rscript 07_Spatial_ccc.R
Rscript 08_Inter-ecotype_communication.R
```

The scripts can also be run interactively in RStudio by opening each script and running the sections sequentially.

## Parameters and paths to modify before running

Several scripts contain project-specific path settings near the beginning of the file. Users should update these paths before running the workflow.

Common path variables include:

```r
DATA_DIR <- "path/to/local/project/data"
OUT_DIR  <- "path/to/local/project/results"
```

Other commonly modified parameters include:

| Parameter | Used in | Description |
|---|---|---|
| `bin.size` | `01_data_processing.R` | Spatial bin size used when loading Visium HD data. |
| `sampleID` | Multiple scripts | Sample identifiers used to subset, merge, and plot spatial objects. |
| `BoneLocation` | Multiple scripts | Metadata variable used to distinguish subchondral and distant/control regions. |
| `group` | Multiple scripts | Case/control group label used for downstream comparison. |


## Troubleshooting

### `source("utils.R")` cannot find the file

Run the scripts from inside the `src/` folder:

```bash
cd src
Rscript 01_data_processing.R
```

Alternatively, update the source line in each script to the full path of `utils.R`.

### Object or file not found

Check the path variables at the beginning of each script, especially:

```r
DATA_DIR
OUT_DIR
obj_file
```

Also confirm that the output from the previous step was successfully generated before running the next step.

### Package not found

Install the missing package from CRAN, Bioconductor, or GitHub, depending on the package. For example:

```r
install.packages("Seurat")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
```

### Memory or parallelization error

Reduce the number of parallel workers in scripts that use `future`, `furrr`, or `future.apply`. For example:

```r
future::plan("multisession", workers = 4)
```

Also adjust:

```r
options(future.globals.maxSize = 100 * 1024^3)
```

according to the available memory.

## Citation

If you use or adapt this code, please cite the associated study.
