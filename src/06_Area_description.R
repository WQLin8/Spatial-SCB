###############################################################################
# Spatial Core-Halo description analysis in subchondral bone
# Author: Weiqiang Lin
###############################################################################
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggprism)
library(ggpubr)
library(ComplexHeatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(limma)
library(GSVA)
library(msigdbr)
library(scRNAtoolVis)
library(randomcoloR)
library(gghalves)
library(ggbreak)
library(stringr)
library(grid)

source("utils.R")
obj_file <- "~/spatial_bone_HD/data_for_downstream_analysis_include_visium/corehalo_obj.rds"
out_dir  <- "~/spatial_bone_HD/results"

###############################################################################
# Global settings
###############################################################################
SPATIAL_ASSAY <- "Spatial.056um"
RCTD_ASSAY    <- "rctd_full"
RING_LEVELS   <- c("Core", "Halo", "Other")
RING_COLORS   <- c(Core = "#F07B70", Halo = "#32B93D", Other = "#6F9AFD")

CORE_GENES <- c("BGLAP", "IBSP", "SP7", "PHOSPHO1", "IFITM5",
                "ALPL", "COL1A1", "DMP1", "RUNX2", "SATB2")
HALO_GENES <- c("POSTN", "COL3A1", "TEK", "MCTP1", "PECAM1",
                "CDH5", "THBS2", "PLXDC2", "OGN")
OTHER_GENES <- c("IGKC", "CD74", "FABP4", "ADIPOQ", "PLIN1",
                 "CD14", "TYROBP", "APOE", "LYZ")
MARK_GENES <- c(CORE_GENES, HALO_GENES, OTHER_GENES)

obj <- readRDS(obj_file)
DefaultAssay(obj) <- SPATIAL_ASSAY
obj$ring <- factor(obj$ring, levels = RING_LEVELS)
Idents(obj) <- "ring"

markers_all <- FindAllMarkers(
  obj,
  assays          = "data",
  only.pos        = FALSE,
  min.pct         = 0.01,
  logfc.threshold = 0.01
)

top_markers <- markers_all %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 500) %>%
  ungroup()

gene_keep <- top_markers$gene
obj_sub   <- subset(obj, features = gene_keep)

st_data <- prepareDataFromscRNA(
  object        = obj_sub,
  assays        = SPATIAL_ASSAY,
  diffData      = top_markers,
  showAverage   = TRUE,
  keep.uniqGene = FALSE
)

enrich <- enrichCluster(
  object        = st_data,
  OrgDb         = org.Hs.eg.db,
  type          = "BP",
  organism      = "hsa",
  pvalueCutoff  = 0.05,
  topn          = 15,
  seed          = 5201314
)

ha <- ComplexHeatmap::HeatmapAnnotation(
  region = c("Core", "Halo", "Other"),
  col = list(region = RING_COLORS),
  gp = grid::gpar(col = "white"),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    region = list(
      title     = "Region",
      title_gp  = grid::gpar(fontsize = 14),
      labels_gp = grid::gpar(fontsize = 12)
    )
  )
)

pdf(file.path(out_dir, "06_hotspot-deg-path.pdf"), height = 10, width = 14, onefile = FALSE)
ht <- visCluster(
  object               = st_data,
  plot.type            = "both",
  column_names_rot     = 0,
  column_names_gp      = grid::gpar(fontsize = 16, just = "left"),
  show_row_dend        = TRUE,
  add.mline            = TRUE,
  markGenes            = MARK_GENES,
  markGenes.side       = "right",
  annoTerm.data        = enrich,
  line.side            = "right",
  panel.arg            = c(4, 0.25, 6, "grey90", NA),
  textbox.size         = 10,
  add.line             = TRUE,
  add.box              = TRUE,
  add.point            = FALSE,
  textbox.pos          = c(0.5, 0.1),
  sample.col           = RING_COLORS,
  annnoblock.gp        = c("white", 12),
  bar.width            = 8,
  cluster.order        = 1:3,
  add.sampleanno       = TRUE,
  annnoblock.text      = FALSE,
  go.size              = "pval",
  genes.gp             = c("italic", 10, NA),
  go.col               = rep(unname(RING_COLORS), each = 15),
  ctAnno.col           = unname(RING_COLORS),
  HeatmapAnnotation    = ha,
  heatmap_legend_param = list(
    title         = "Z-score",
    title_gp      = grid::gpar(fontsize = 14),
    labels_gp     = grid::gpar(fontsize = 12),
    legend_height = unit(3, "cm"),
    grid_width    = unit(0.6, "cm")
  ),
  add.bar              = FALSE
)

draw(ht, annotation_legend_side = "right", heatmap_legend_side = "right")
grid::grid.text(
  "GO Biological Process",
  x    = unit(0.88, "npc"),
  y    = unit(0.96, "npc"),
  just = "right",
  gp   = grid::gpar(fontsize = 16, fontface = "bold", col = "firebrick")
)
dev.off()

###############################################################################
# 2. Subchondral vs Distant: overall hotspot DE + GO
###############################################################################

# Subchondral-focused samples
obj_subchondral <- subset(obj, sampleID != "visium_12" & sampleID != "hd_82_1")
res_sub <- run_region_specific_go(obj_subchondral,
                                  logfc_threshold_findall = 0.1,
                                  sig_logfc_cutoff = 0.5)

p_sub <- jjVolcano(
  diffData         = res_sub$markers,
  log2FC.cutoff    = 0.1,
  fontface         = "italic",
  topGeneN         = 10,
  tile.col         = distinctColorPalette(3),
  legend.position  = c(0.8, 0.95),
  flip             = FALSE,
  pSize            = 1,
  height           = 0.4,
  celltypeSize     = 5
)

ggsave(file.path(out_dir, "volcano-Core-halo-DE_Subchondral.pdf"), p_sub, height = 10, width = 6)
plot_go_dotplot(res_sub$go_df, file.path(out_dir, "De-path_Subchondral.pdf"))

# Distant-focused samples
obj_distant <- subset(obj, sampleID %in% c("visium_12", "hd_82_1"))
res_dis <- run_region_specific_go(obj_distant,
                                  logfc_threshold_findall = 0.01,
                                  sig_logfc_cutoff = 0.5)

p_dis <- jjVolcano(
  diffData         = res_dis$markers,
  log2FC.cutoff    = 0.01,
  fontface         = "italic",
  topGeneN         = 10,
  tile.col         = distinctColorPalette(3),
  legend.position  = c(0.8, 0.95),
  flip             = FALSE,
  pSize            = 1,
  height           = 0.4,
  celltypeSize     = 5
)

ggsave(file.path(out_dir, "volcano-Core-halo-DE_Distant.pdf"), p_dis, height = 10, width = 6)
plot_go_dotplot(res_dis$go_df, file.path(out_dir, "Dde-path_Distant.pdf"))

###############################################################################
# 3. Core area: DEG + GSVA/limma
###############################################################################
obj_core <- subset(obj, ring == "Core")

Idents(obj_core) <- "ring"
deg_core <- FindMarkers(
  obj_core,
  ident.1  = "Subchondral",
  ident.2  = "Distant",
  group.by = "BoneLocation"
)

plot_volcano_with_labels(
  deg_df        = deg_core,
  title_text    = "Core",
  outfile       = file.path(out_dir, "volcano-Core-DE.pdf"),
  fdr_cutoff    = 0.05,
  logfc_cutoff  = 0.1,
  xlim_range    = c(-15, 10),
  n_label       = 25
)

# Original code ultimately overwrote GO-BP with Hallmark; keep Hallmark here
core_gsva_diff <- run_gsva_limma(obj_core, assay = SPATIAL_ASSAY, gene_set_type = "hallmark")
plot_gsva_bar(core_gsva_diff, file.path(out_dir, "gsva_bar_core.pdf"))

###############################################################################
# 4. Halo area: DEG + ranked plot + GSVA/limma
###############################################################################
obj_halo <- subset(obj, ring == "Halo")

Idents(obj_halo) <- "ring"
deg_halo <- FindMarkers(
  obj_halo,
  ident.1  = "Subchondral",
  ident.2  = "Distant",
  group.by = "BoneLocation"
)
write.csv(deg_halo, file.path(out_dir, "Halo-DE-list.csv"))

deg_halo_tbl <- make_deg_table(deg_halo, fdr_cutoff = 0.05, logfc_cutoff = 0.1)
up_genes <- deg_halo_tbl %>% filter(Sig == "Up") %>% arrange(desc(avg_log2FC))
write.csv(up_genes, file.path(out_dir, "Halo-DE-Up_list.csv"), row.names = FALSE)

# Ranked DEG plot
rank_df <- deg_halo_tbl %>%
  rownames_to_column(var = "gene_name") %>%
  mutate(rank = -1 * rank(avg_log2FC, ties.method = "max"))

genes_to_label <- c(
  "POSTN", "C1QTNF8", "GJA5", "ENPP1", "MYH11", "ACTN2",
  "GPIHBP1", "EGFL7", "CD276", "PTH1R", "MMP13", "DLX5",
  "SPARC", "IGHD", "MZB1", "FABP5", "CEBPA", "PLIN4", "ADIPOQ"
)

rank_label_df <- rank_df %>% filter(gene_name %in% genes_to_label)

p_rank <- rank_df %>%
  ggplot(aes(x = rank, y = avg_log2FC, color = p_val_adj, size = abs(avg_log2FC))) +
  geom_point(alpha = 0.85) +
  scale_x_continuous(
    breaks = c(-15000, -10000, -5000, 0),
    labels = c(0, 5000, 10000, 15000)
  ) +
  scale_color_gradient2(
    low      = "#f88072",
    high     = "#8bb1d3",
    mid      = "#ffffff",
    midpoint = 0.5,
    name     = "adj P-value"
  ) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_text_repel(
    data = rank_label_df,
    aes(label = gene_name),
    box.padding       = 0.4,
    segment.curvature = -0.1,
    segment.ncp       = 3,
    segment.angle     = 20,
    max.overlaps      = Inf
  ) +
  scale_size(range = c(1, 8), name = "abs(log2FC)") +
  labs(
    x = "Rank of differentially expressed genes",
    y = "log2(Fold Change)"
  ) +
  theme_bw() +
  theme(
    panel.border      = element_rect(linewidth = 1.25),
    legend.background = element_rect(color = "#888888", linetype = 1),
    axis.text         = element_text(color = "#000000", size = 12.5),
    axis.title        = element_text(color = "#000000", size = 15)
  )

ggsave(file.path(out_dir_manu, "halo_DEG_plot.pdf"), p_rank, height = 7, width = 9)

# Volcano
# Prevent infinite -log10 values in plotting
safe_deg_halo <- deg_halo
safe_deg_halo$p_val_adj[safe_deg_halo$p_val_adj < 1e-300] <- 1e-300
plot_volcano_with_labels(
  deg_df        = safe_deg_halo,
  title_text    = "Halo",
  outfile       = file.path(out_dir, "volcano-Halo-DE.pdf"),
  fdr_cutoff    = 0.05,
  logfc_cutoff  = 0.1,
  xlim_range    = c(-15, 10),
  n_label       = 25
)

# GSVA/limma
halo_gsva_diff <- run_gsva_limma(obj_halo, assay = SPATIAL_ASSAY, gene_set_type = "hallmark")
plot_gsva_bar(halo_gsva_diff, file.path(out_dir, "gsva_bar_halo.pdf"))

###############################################################################
# 5. Osteogenic area: RCTD cell proportion comparison across hotspot regions
###############################################################################
obj_rctd <- readRDS(obj_file)
df_long <- prepare_rctd_long(obj_rctd, assay = RCTD_ASSAY)

# Example: single cell-type plot
celltype_name <- "Hematopoietic Stem-Progenitor Cell"
plot_single_celltype_violin(
  df_long,
  celltype_name,
  file.path(out_dir, paste0("ring_", celltype_name, "_cell_proportion.pdf"))
)

# All cell types in one plot with y-axis breaks
plot_all_celltype_boxbreak(
  df_long,
  file.path(out_dir, "area_cell_proportion.pdf")
)

