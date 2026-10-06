############################################################
# Spatial Cell Identification (Osteoclast & Adipocyte)
# Author: Weiqiang Lin
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(UCell)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cols4all)
library(RColorBrewer)

source("utils.R")

DATA_DIR <- "data/"
OUT_DIR  <- "results/"

obj <- readRDS(file.path(DATA_DIR, "hd_marrow_combined_16um.rds"))
DefaultAssay(obj) <- "Spatial.016um"


# ===================== Osteoclast =====================
osteoclast_markers <- c(
  "ACP5","CTSK","CALCR","TCIRG1","CLCN7",
  "TNFRSF11A","ATP6V1C1","ATP6V0D2",
  "DCSTAMP","OCSTAMP","TM7SF4"
)

obj <- identify_cell_state(obj, osteoclast_markers, "Osteoclast")

plot_spatial_binary(obj, "Osteoclast", "osteoclast_spatial.pdf")

de_oc <- run_volcano(obj, "Osteoclast", "volcano_osteoclast.pdf")

run_go(de_oc, "go_osteoclast.pdf")

# ===================== Adipocyte =====================

adipocyte_markers <- c(
  "PLIN1","PLIN4","ADIPOQ","CIDEC",
  "G0S2","THRSP","LEP","PNPLA2","LIPE","CIDEA"
)

obj <- identify_cell_state(obj, adipocyte_markers, "Adipocyte", 0.9997)

plot_spatial_binary(obj, "Adipocyte", "adipocyte_spatial.pdf")

de_ad <- run_volcano(obj, "Adipocyte", "volcano_adipocyte.pdf")

run_go(de_ad, "go_adipocyte.pdf")

saveRDS(obj, file.path(DATA_DIR, "hd_marrow_cellstates.rds"))

