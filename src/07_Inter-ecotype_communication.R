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

obj <- readRDS("data/integrated_obj_RCTD_hotspot_ISCHIA.rds")

Idents(obj) <- "cc_13"


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
    shape <- ifelse(sample_id == "V_01", 21, 22)

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

pdf("spatial_map.pdf", 10, 8)
print(plot_spatial(obj,
                    c("HD_03","HD_01","HD_04","HD_02","V_01"),
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
  meta$cc_13 <- as.factor(meta$cc_13)

  colnames(spatial.locs) <- c("x","y")

  cellchat <- createCellChat(
    object = data.input,
    meta = meta[colnames(data.input), ],
    group.by = "cc_13",
    datatype = "spatial",
    coordinates = spatial.locs[colnames(data.input), ],
    spatial.factors = spatial.factors
  )

  return(cellchat)
}


##################################################
#### 3. Case analysis ###########################
##################################################

samples_case <- c("HD_03","HD_01","HD_02")

location_case <- c(
  "Human_HD_03/tissue_positions.parquet",
  "Human_HD_01/tissue_positions.parquet",
  "Human_HD_04/tissue_positions.parquet"
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

##################################################
#### 4. Control analysis #########################
##################################################

samples_control <- c("HD_04","V_01")

location_control <- c(
  "Human_HD_04/tissue_positions.parquet",
  "Human_V_01/tissue_positions.csv"
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


object.list <- list(Case=cellchat.case, Control=cellchat.control)

cellchat.merge <- mergeCellChat(object.list, add.names=names(object.list))


##################################################
#### 5. Differential interaction #################
##################################################

pdf("diff_interaction.pdf",10,6)
netVisual_diffInteraction(
  cellchat.merge,
  weight.scale=TRUE,
  comparison=c(2,1)
)
dev.off()


##################################################
#### 6. Bubble comparison ########################
##################################################

pdf("bubble_compare.pdf",10,5)
netVisual_bubble(
  cellchat.merge,
  sources.use=c("CC10","CC2","CC11","CC5","CC9"),
  targets.use=c("CC10","CC2","CC11","CC5","CC9"),
  comparison=c(1,2),
  angle.x=45
)
dev.off()
