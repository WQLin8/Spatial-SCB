############################################################
# Human Bone Spatial Transcriptomics Analysis Pipeline
# Author: Weiqiang Lin
# Description: Integration and analysis of Visium HD + Visium
############################################################

library(Seurat)
library(dplyr)
library(msigdbr)
library(UCell)
library(clustree)
library(patchwork)
library(RColorBrewer)

options(future.globals.maxSize = 1e9)

source("utils.R")

DATA_DIR <- "~/spatial_bone_HD/Processed_data_via_SpaceRanger"
OUT_DIR  <- "~/spatial_bone_HD"

hd_68_1 <- load_hd_sample("Human_HD_68", 1, "Subchondral", "Case", TRUE)
hd_68_2 <- load_hd_sample("Human_HD_68", 2, "Subchondral", "Case", FALSE, "trab_segmentation")

hd_82_1 <- load_hd_sample("Human_HD_82", 1, "Distant", "Control", TRUE)
hd_82_2 <- load_hd_sample("Human_HD_82", 2, "Subchondral", "Case", FALSE, "trab.segmantation")


visium_12 <- Load10X_Spatial(file.path(DATA_DIR, "Human_Visium_12"))

DefaultAssay(visium_12) <- "Spatial"
visium_12[["Spatial.056um"]] <- visium_12[["Spatial"]]
visium_12@assays$Spatial <- NULL

DefaultAssay(visium_12) <- "Spatial.056um"

# Fix image slot
visium_12@images[["slice1.056um"]] <- visium_12@images[["slice1"]]
visium_12@images[["slice1"]] <- NULL
visium_12@images[["slice1.056um"]]@assay <- "Spatial.056um"

# Metadata
trab <- read.csv(file.path(DATA_DIR, "trab_location/trab-12.csv"))
rownames(trab) <- trab$Barcode
trab$trab[trab$trab == ""] <- "marrow"

visium_12$RegionType <- trab[colnames(visium_12), "trab"]
visium_12 <- RenameCells(visium_12, add.cell.id = "visium_12")
visium_12$sampleID <- "visium_12"
visium_12$BoneLocation <- "Distant"
visium_12$group <- "Control"

hd <- merge(
  hd_68_1,
  y = list(hd_68_2, hd_82_1, hd_82_2, visium_12),
  project = "Human_Bone"
)

hd <- JoinLayers(hd)

obj <- hd

obj[["Spatial.056um"]] <- split(obj[["Spatial.056um"]], f = obj$sampleID)

obj <- subset(obj, nCount_Spatial.056um >= 50 & RegionType == "marrow")

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, percent.mt <= 20)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 3000)
obj <- ScaleData(obj)

obj <- RunPCA(obj)

obj <- IntegrateLayers(
  object = obj,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony"
)

obj <- JoinLayers(obj)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.3)
obj <- RunUMAP(obj, dims = 1:15, reduction = "harmony")

osteoblast <- c("BGLAP", "SPP1", "IBSP", "COL1A1", "MEPE", "DMP1")
endothelial <- c("PECAM1", "VWF", "CDH5", "EMCN", "CLDN5", "ENG")
vsmc <- c("ACTA2", "TAGLN", "MYH11", "CNN1")

obj <- AddModuleScore(obj, features = list(osteoblast), name = "OB_score")
obj <- AddModuleScore(obj, features = list(endothelial), name = "endo_score")
obj <- AddModuleScore(obj, features = list(vsmc), name = "vsmc_score")

saveRDS(obj, file.path(OUT_DIR, "processed/hd_marrow_processed.rds"))

plot_spatial_feature_grid(obj, feature = "BGLAP",
                          path = file.path(OUT_DIR, "results/BGLAP.pdf"))

plot_spatial_dim_grid(
  hd,
  path = file.path(OUT_DIR, "results/region.pdf"),
  cols = c("marrow" = "#85B8D5", "trab" = "#DC9C9B")
)
