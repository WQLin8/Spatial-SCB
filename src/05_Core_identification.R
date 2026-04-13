############################################################
# Spatial Core-Halo Architecture Analysis
# Author: Weiqiang Lin
############################################################

library(Seurat)
library(dplyr)
library(patchwork)

source("utils.R")

DATA_DIR <- "data/"
OUT_DIR  <- "results/"

obj <- readRDS(file.path(DATA_DIR, "integrated_RCTD_obj.rds"))

obj <- identify_Core(obj, "Osteoblast", 0.85)

plot_Core(obj, "Osteoblast_Core", "osteoblast_Core.pdf")

filter_all_samples <- function(obj) {
  sample_paths <- list(
    "hd_68_1" = "Human_HD_68_batch1/outs/binned_outputs/square_056um/spatial/tissue_positions.parquet",
    "hd_68_2" = "Human_HD_68_batch2/outs/binned_outputs/square_056um/spatial/tissue_positions.parquet",
    "hd_82_1" = "Human_HD_82_batch1/outs/binned_outputs/square_056um/spatial/tissue_positions.parquet",
    "hd_82_2" = "Human_HD_82_batch2/outs/binned_outputs/square_056um/spatial/tissue_positions.parquet",
    "visium_12" = "Human_Visium_12/spatial/tissue_positions.csv"
  )
  
  obj_list <- lapply(names(sample_paths), function(sid) {
    
    filter_low_connectivity_areas(
      sample_id = sid,
      location_path = file.path(DATA_DIR, sample_paths[[sid]]),
      obj = obj,
      size = 4
    )
    
  })
  
  names(obj_list) <- names(sample_paths)
  
  obj_new <- merge(obj_list[[1]], y = obj_list[-1])
  
  DefaultAssay(obj_new) <- "Spatial.056um"
  obj_new <- JoinLayers(obj_new)
  
  return(obj_new)
}

obj <- filter_all_samples(obj)

plot_ring(obj, "core_halo_structure.pdf")

saveRDS(obj, file.path(DATA_DIR, "corehalo_obj.rds"))

