############################################################
# Single-cell Bone Marrow Reference Processing
# Author: Weiqiang Lin
# Description: Cell type re-annotation and validation
############################################################

# ===================== 1. Libraries =====================
library(Seurat)
library(dplyr)
library(ggplot2)
library(randomcoloR)
library(ComplexHeatmap)

source("utils.R")

# ===================== 2. Paths =====================
DATA_DIR <- "data/"
OUT_DIR  <- "results/"

# ===================== 3. Load Reference =====================
reference <- readRDS(file.path(DATA_DIR, "bone_marrow_reference.rds"))

Idents(reference) <- "cluster_anno_l2"

# ===================== 4. Cell Type Mapping =====================

# ---- Level 1 (broad lineage) ----
mapping_lvl1 <- c(
  "Plasma Cell" = "B Cell",
  "Pre-B" = "B Cell",
  "Pro-B" = "B Cell",
  "Pre-Pro B" = "B Cell",
  "Mature B" = "B Cell",
  "CD8+ T-Cell" = "T Cell",
  "CD4+ T-Cell" = "T Cell",
  "Neutrophil" = "Granulocytic Lineage",
  "Ba/Eo/Ma" = "Granulocytic Lineage",
  "Monocyte" = "Monocyte",
  "Macrophages" = "Macrophages",
  "GMP" = "Myeloid Progenitor",
  "Late Myeloid" = "Myeloid Progenitor",
  "Early Myeloid Progenitor" = "Myeloid Progenitor",
  "pDC" = "Dendritic Cell Lineage",
  "Cycling DCs" = "Dendritic Cell Lineage",
  "Adipo-MSC" = "Mesenchymal Stem Cell",
  "Osteo-MSC" = "Mesenchymal Stem Cell",
  "APOD+ MSC" = "Mesenchymal Stem Cell",
  "THY1+ MSC" = "Mesenchymal Stem Cell",
  "RNAlo MSC" = "Mesenchymal Stem Cell",
  "Fibro-MSC" = "Mesenchymal Stem Cell",
  "Osteoblast" = "Osteoblast",
  "Late Erythroid" = "Megakaryocyte-Erythroid Lineage",
  "RBC" = "Megakaryocyte-Erythroid Lineage",
  "Erythroblast" = "Megakaryocyte-Erythroid Lineage",
  "MEP" = "Megakaryocyte-Erythroid Lineage",
  "Megakaryocyte" = "Megakaryocyte-Erythroid Lineage",
  "MPP" = "HSPC",
  "CLP" = "HSPC",
  "HSC" = "HSPC",
  "Cycling HSPC" = "HSPC",
  "VSMC" = "VSMC",
  "AEC" = "Endothelial Cell",
  "SEC" = "Endothelial Cell"
)

# ---- Level 2 (refined) ----
mapping_lvl2 <- c(
  "Plasma Cell" = "B Lineage",
  "Pre-B" = "B Lineage",
  "Pro-B" = "B Lineage",
  "Pre-Pro B" = "B Lineage",
  "Mature B" = "B Lineage",
  "CD8+ T-Cell" = "T Cell",
  "CD4+ T-Cell" = "T Cell",
  "Neutrophil" = "Granulocytic Lineage",
  "Ba/Eo/Ma" = "Granulocytic Lineage",
  "Monocyte" = "Monocyte",
  "Macrophages" = "Macrophages",
  "GMP" = "MP",
  "Late Myeloid" = "MP",
  "Early Myeloid Progenitor" = "MP",
  "pDC" = "CLP",
  "Cycling DCs" = "DC",
  "Adipo-MSC" = "Adipo-MSC",
  "Osteo-MSC" = "Osteo-MSC",
  "APOD+ MSC" = "Adipo-MSC",
  "THY1+ MSC" = "THY1+ MSC",
  "RNAlo MSC" = "RNAlo MSC",
  "Fibro-MSC" = "Fibro-MSC",
  "Osteoblast" = "Osteoblast",
  "Late Erythroid" = "Erythroid",
  "RBC" = "Erythroid",
  "Erythroblast" = "Erythroid",
  "MEP" = "MEP",
  "Megakaryocyte" = "Megakaryocyte",
  "MPP" = "HSPC",
  "HSC" = "HSPC",
  "Cycling HSPC" = "HSPC",
  "VSMC" = "VSMC",
  "AEC" = "AEC",
  "SEC" = "SEC"
)

# ===================== 5. Apply Mapping =====================

reference$cell_type_lvl1 <- mapping_lvl1[as.character(reference$cluster_anno_l2)]
reference$cell_type_lvl2 <- mapping_lvl2[as.character(reference$cluster_anno_l2)]

# Remove unwanted cells
reference <- subset(reference, cell_type_lvl2 != "RNAlo MSC")

Idents(reference) <- "cell_type_lvl2"

# ===================== 6. UMAP Visualization =====================

palette <- distinctColorPalette(length(unique(reference$cell_type_lvl2)))

umap_df <- reference@reductions$UMAP_dim30@cell.embeddings %>%
  as.data.frame() %>%
  cbind(cell_type = reference$cell_type_lvl2)

p_umap <- ggplot(umap_df, aes(UMAP_1, UMAP_2, fill = cell_type)) +
  geom_point(size = 1.8, shape = 21, colour = "grey30", stroke = 0.05) +
  scale_fill_manual(values = palette) +
  theme_classic() +
  labs(fill = "Cell type")

ggsave(file.path(OUT_DIR, "sc_reference_umap.pdf"),
       p_umap, width = 8, height = 6)

# ===================== 7. Save =====================

saveRDS(reference,
        file.path(DATA_DIR, "sc_reference_refined.rds"))

# ===================== 8. Marker Gene Validation =====================

markers <- c(
  "MS4A1","CD79A","CD19",
  "PPARG","LPL","PLIN1",
  "CSF3R","S100A8",
  "HBB","HBA1","KLF1",
  "BGLAP","SPP1","COL1A1",
  "RUNX2","SP7",
  "THY1","PDGFRA",
  "ACTA2","TAGLN","MYH11",
  "GATA1","ITGA2B",
  "KDR","PECAM1","VWF",
  "CD3D","CD4",
  "GJA5","DLL4",
  "CD68","CD163",
  "PROM1","GATA2",
  "LYZ","CTSS",
  "HLA-DRA","FCER1A",
  "DCN","LUM",
  "IL7R","DNTT",
  "PF4"
)

avg_expr <- AverageExpression(reference, features = markers)$RNA
avg_expr <- t(scale(t(avg_expr)))

# Color function
col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

ht <- Heatmap(
  avg_expr,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  rect_gp = gpar(col = "black", lwd = 0.3),
  name = "Z-score"
)

pdf(file.path(OUT_DIR, "sc_marker_heatmap.pdf"), width = 6, height = 10)
draw(ht)
dev.off()

