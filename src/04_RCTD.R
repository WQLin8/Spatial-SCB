############################################################
# Spatial Deconvolution using RCTD
# Author: Weiqiang Lin
############################################################

# ===================== 1. Libraries =====================
library(spacexr)
library(Seurat)
library(dplyr)
library(Matrix)
library(patchwork)

source("utils.R")

# ===================== 2. Paths =====================
DATA_DIR <- "data/"
OUT_DIR  <- "results/"

# ===================== 3. Load Data =====================
sc   <- readRDS(file.path(DATA_DIR, "sc_reference_refined.rds"))
obj1 <- readRDS(file.path(DATA_DIR, "hd_marrow_processed_56um.rds"))
obj2 <- readRDS(file.path(DATA_DIR, "hd_marrow_cellstates.rds"))


reference <- build_reference(sc, obj2)

query <- build_query(obj1)

RCTD_res <- run_rctd_pipeline(query, reference)

saveRDS(RCTD_res, file.path(DATA_DIR, "rctd_full.rds"))

weights <- RCTD_res@results$weights
norm_weights <- normalize_weights(weights)

obj1[["RCTD"]] <- CreateAssayObject(data = t(as.matrix(norm_weights)))
DefaultAssay(obj1) <- "RCTD"

plot_rctd <- function(obj, cell_types) {
  
  for (ct in cell_types) {
    plot_spatial_feature_grid(
      obj,
      feature = ct,
      path = file.path(OUT_DIR,
                       paste0("RCTD_", gsub(" ", "_", ct), ".pdf"))
    )
  }
}

cell_types <- colnames(obj1[["RCTD"]])

plot_rctd(obj1, cell_types)

saveRDS(obj1, file.path(DATA_DIR, "integrated_RCTD_obj.rds"))

