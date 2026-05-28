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

hd_03 <- load_hd_sample("Human_HD_03", 1, "Subchondral", "Case", TRUE)
hd_01 <- load_hd_sample("Human_HD_01", 2, "Subchondral", "Case", FALSE, "trab_segmentation")

hd_04 <- load_hd_sample("Human_HD_04", 1, "Distant", "Control", TRUE)
hd_01 <- load_hd_sample("Human_HD_01", 2, "Subchondral", "Case", FALSE, "trab.segmantation")


v_01 <- Load10X_Spatial(file.path(DATA_DIR, "Human_Visium_12"))

DefaultAssay(v_01) <- "Spatial"
v_01[["Spatial.056um"]] <- v_01[["Spatial"]]
v_01@assays$Spatial <- NULL

DefaultAssay(visium_12) <- "Spatial.056um"

# Fix image slot
v_01@images[["slice1.056um"]] <- v_01@images[["slice1"]]
v_01@images[["slice1"]] <- NULL
v_01@images[["slice1.056um"]]@assay <- "Spatial.056um"

# Metadata
trab <- read.csv(file.path(DATA_DIR, "trab_location/trab-12.csv"))
rownames(trab) <- trab$Barcode
trab$trab[trab$trab == ""] <- "marrow"

v_01$RegionType <- trab[colnames(v_01), "trab"]
v_01 <- RenameCells(v_01, add.cell.id = "v_01")
v_01$sampleID <- "v_01"
v_01$BoneLocation <- "Distant"
v_01$group <- "Control"

hd <- merge(
  hd_03,
  y = list(hd_01, hd_04, hd_02, v_01),
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
