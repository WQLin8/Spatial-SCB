#############################################################################
# Inter-ecotype communication analysis
# Author: Weiqiang Lin
#############################################################################

library(CellChat)
library(patchwork)
library(Seurat)
library(jsonlite)
library(arrow)
library(tibble)
library(purrr)
library(liana)
library(tidyverse)
library(CrossTalkeR)
library(dplyr)

options(future.globals.maxSize = 250 * 1024^3)

source("utils.R")

obj <- readRDS("integrated_obj_RCTD_hotspot_ISCHIA.rds")

Idents(obj) <- "cc_10"


##################################################
#### 1. Spatial visualization function ###########
##################################################

paletteMartin <- c(
  '#e6194b','#3cb44b','#ffe119','#4363d8','#f0995bff',
  '#911eb4','#09f309ff','#46f0f0','#f032e6','#bcf60c',
  '#fabebe','#008080','#e6beff','#f35f09ff'
)

all_ccs <- unique(obj$CompositionCluster_CC)
color_mapping <- setNames(paletteMartin[seq_along(all_ccs)], all_ccs)

plot_spatial <- function(obj, sample_ids, color_mapping){

  plots <- list()

  for(i in seq_along(sample_ids)){

    sample_id <- sample_ids[i]
    sample_obj <- subset(obj, sampleID == sample_id)

    aspect_ratio <- get_image_aspect_ratio(sample_obj)

    pt_size <- c(5,3.2,3.2,3.2,2.5)[i]
    shape <- ifelse(sample_id == "visium_12", 21, 22)

    plots[[i]] <- SpatialDimPlot(
      sample_obj,
      group.by = "CompositionCluster_CC",
      pt.size.factor = pt_size,
      shape = shape,
      cols = color_mapping
    ) +
      theme(
        aspect.ratio = aspect_ratio,
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.position = "right"
      )
  }

  wrap_plots(plots, ncol = 2)
}

pdf("01_spatial_map.pdf", 10, 8)
print(plot_spatial(obj,
                    c("hd_68_1","hd_68_2","hd_82_1","hd_82_2","visium_12"),
                    color_mapping))
dev.off()


##################################################
#### 2. CellChat builder function ################
##################################################

build_cellchat <- function(obj, samples, location_paths, assay="Spatial.056um"){

  data.input <- NULL
  meta <- NULL
  spatial.locs <- NULL
  spatial.factors <- data.frame()

  for(i in seq_along(samples)){

    samp <- samples[i]
    location_path <- location_paths[i]

    spatial.loc <- if (grepl("\\.parquet$", location_path)) {
      read_parquet(location_path) %>%
        as.data.frame() %>%
        select(barcode, imagerow = 3, imagecol = 4)
    } else {
      read.csv(location_path) %>%
        as.data.frame() %>%
        select(barcode, imagerow = 3, imagecol = 4)
    }

    rownames(spatial.loc) <- paste0(samp, "_", spatial.loc$barcode)
    spatial.loc <- spatial.loc[, -1]

    obj.sub <- subset(obj, sampleID == samp)

    data.input <- cbind(
      data.input,
      GetAssayData(obj.sub, slot="data", assay=assay)
    )

    meta <- rbind(meta, obj.sub@meta.data)

    spatial.loc <- spatial.loc[colnames(obj.sub), ]
    spatial.locs <- rbind(spatial.locs, spatial.loc)

    spatial.factors <- rbind(spatial.factors,
                             data.frame(ratio=1, tol=0.75, row.names=samp))
  }

  meta$samples <- as.factor(meta$sampleID)
  meta$cc_10 <- as.factor(meta$cc_10)

  colnames(spatial.locs) <- c("x","y")

  cellchat <- createCellChat(
    object = data.input,
    meta = meta[colnames(data.input), ],
    group.by = "cc_10",
    datatype = "spatial",
    coordinates = spatial.locs[colnames(data.input), ],
    spatial.factors = spatial.factors
  )

  return(cellchat)
}


##################################################
#### 3. Case analysis ###########################
##################################################

samples_case <- c("hd_68_1","hd_68_2","hd_82_2")

location_case <- c(
  "/Human_HD_68_batch1/...parquet",
  "/Human_HD_68_batch2/...parquet",
  "/Human_HD_82_batch2/...parquet"
)

cellchat.case <- build_cellchat(obj, samples_case, location_case)

CellChatDB <- CellChatDB.human
cellchat.case@DB <- subsetDB(CellChatDB)

cellchat.case <- subsetData(cellchat.case)

future::plan("multisession", workers=4)

cellchat.case <- identifyOverExpressedGenes(cellchat.case)
cellchat.case <- identifyOverExpressedInteractions(cellchat.case)

cellchat.case <- computeCommunProb(
  cellchat.case,
  trim=0.1,
  k.min=5,
  type="truncatedMean",
  distance.use=TRUE,
  interaction.range=6,
  scale.distance=1,
  contact.range=2,
  contact.dependent=FALSE
)

cellchat.case <- filterCommunication(cellchat.case, min.cells=0)
cellchat.case <- computeCommunProbPathway(cellchat.case)
cellchat.case <- aggregateNet(cellchat.case)


pdf("02_case_circle.pdf",6,6)
netVisual_circle(cellchat.case@net$count,
                 vertex.weight=rowSums(cellchat.case@net$count),
                 weight.scale=TRUE,
                 label.edge=TRUE)
dev.off()


##################################################
#### 4. Control analysis #########################
##################################################

samples_control <- c("hd_82_1","visium_12")

location_control <- c(
  "/82_batch1/...parquet",
  "/visium12.csv"
)

cellchat.control <- build_cellchat(obj, samples_control, location_control)

cellchat.control@DB <- subsetDB(CellChatDB)
cellchat.control <- subsetData(cellchat.control)

cellchat.control <- identifyOverExpressedGenes(cellchat.control)
cellchat.control <- identifyOverExpressedInteractions(cellchat.control)

cellchat.control <- computeCommunProb(
  cellchat.control,
  trim=0.1,
  k.min=5,
  type="truncatedMean",
  distance.use=TRUE,
  interaction.range=6,
  scale.distance=1,
  contact.range=2,
  contact.dependent=FALSE
)

cellchat.control <- filterCommunication(cellchat.control, min.cells=0)
cellchat.control <- computeCommunProbPathway(cellchat.control)
cellchat.control <- aggregateNet(cellchat.control)


pdf("03_control_circle.pdf",6,6)
netVisual_circle(cellchat.control@net$count,
                 vertex.weight=rowSums(cellchat.control@net$count),
                 weight.scale=TRUE,
                 label.edge=FALSE)
dev.off()


##################################################
#### 5. LR table + chord plot ###################
##################################################

get_lr_table <- function(cellchat_obj){

  subsetCommunication(cellchat_obj) %>%
    filter(pval <= 0.05) %>%
    select(source, target, ligand, receptor, prob, pval) %>%
    as_tibble()
}

df.case <- get_lr_table(cellchat.case)
df.control <- get_lr_table(cellchat.control)


pdf("04_case_chord.pdf",8,8)
chord_freq(df.case, facing="bending.outside", cex=1.2, adj=c(0.5,2.3))
dev.off()

pdf("05_control_chord.pdf",8,8)
chord_freq(df.control, facing="bending.outside", cex=1.2, adj=c(0.5,2.3))
dev.off()


##################################################
#### 6. Centrality analysis ######################
##################################################

cellchat.case <- netAnalysis_computeCentrality(cellchat.case)

pdf("06_case_outgoing.pdf",8,8)
netAnalysis_signalingRole_scatter(cellchat.case)
dev.off()

pdf("07_case_heatmap.pdf",8,8)
netAnalysis_signalingRole_heatmap(cellchat.case, pattern="outgoing")
dev.off()


##################################################
#### 7. Merge Case + Control #####################
##################################################

object.list <- list(Case=cellchat.case, Control=cellchat.control)

cellchat.merge <- mergeCellChat(object.list, add.names=names(object.list))


##################################################
#### 8. Differential interaction #################
##################################################

pdf("08_diff_interaction.pdf",10,6)
netVisual_diffInteraction(
  cellchat.merge,
  weight.scale=TRUE,
  comparison=c(2,1)
)
dev.off()


##################################################
#### 9. Bubble comparison ########################
##################################################

pdf("09_bubble_compare.pdf",10,5)
netVisual_bubble(
  cellchat.merge,
  sources.use=c("CC10","CC2","CC11","CC5","CC9"),
  targets.use=c("CC10","CC2","CC11","CC5","CC9"),
  comparison=c(1,2),
  angle.x=45
)
dev.off()