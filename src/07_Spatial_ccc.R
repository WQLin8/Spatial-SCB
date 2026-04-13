#############################################################################
# Spatial Cell-Cell Communication Analysis
# Author: Weiqiang Lin
#############################################################################
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(furrr)
library(future.apply)
library(metap)
library(UpSetR)
library(arrow)
library(dplyr)

plan(multisession, workers = 20)
options(future.globals.maxSize = 100 * 1024^3)

source("utils.R")

LRpair <- read.csv("data/LRpair_human.csv", row.names = 1)
obj <- readRDS("data/05_integrated_obj_RCTD_hotspot.rds")

DefaultAssay(obj) <- "Spatial.056um"
Idents(obj) <- "ring"
obj_sub <- subset(obj, ring != "Other")

#############################################################################
########################### 1. LR Colocalization ###########################
#############################################################################

sample_info <- tibble(
  sample_id = c("hd_68_1","hd_68_2","hd_82_1","hd_82_2","visium_12"),
  location_path = c(
    "spatial/Human_HD_68_batch1/tissue_positions.parquet",
    "spatial/Human_HD_68_batch2/tissue_positions.parquet",
    "spatial/Human_HD_82_batch1/tissue_positions.parquet",
    "spatial/Human_HD_82_batch2/tissue_positions.parquet",
    "spatial/Human_Visium_12/tissue_positions.csv"
  )
)

run_LR_coloc <- function(sample_id, location_path) {
  LR_Colco(
    obj = obj_sub,
    sample_id = sample_id,
    location_path = location_path,
    distance = 6,
    LRpair = LRpair
  )
}

LR_results <- future_pmap(
  list(sample_info$sample_id, sample_info$location_path),
  run_LR_coloc
)

names(LR_results) <- sample_info$sample_id

#############################################################################
########################### 2. Extract LR pairs ############################
#############################################################################

lr_pairs_list <- lapply(names(LR_results), function(sid) {
  res <- LR_results[[sid]]
  if (is.null(res) || nrow(res) == 0) return(NULL)
  rownames(res)
})

names(lr_pairs_list) <- names(LR_results)

lr_long <- bind_rows(lapply(names(lr_pairs_list), function(sid) {
  data.frame(
    LR_pair = lr_pairs_list[[sid]],
    sample_id = sid
  )
}))

write.csv(lr_long, "results/LR_colocalization_all_samples.csv", row.names = FALSE)

#############################################################################
########################### 3. UpSet plot ##################################
#############################################################################

all_pairs <- unique(unlist(lr_pairs_list))

upset_mat <- matrix(
  0,
  nrow = length(all_pairs),
  ncol = length(lr_pairs_list),
  dimnames = list(all_pairs, names(lr_pairs_list))
)

for (sid in names(lr_pairs_list)) {
  upset_mat[lr_pairs_list[[sid]], sid] <- 1
}

upset_df <- as.data.frame(upset_mat)

pdf("results/LR_upset.pdf", width = 8, height = 5)
upset(
  upset_df,
  sets = names(lr_pairs_list),
  order.by = "freq",
  mb.ratio = c(0.6, 0.4),
  text.scale = c(1.3, 1.2, 1.2, 1.2, 1.5, 1.2),
  mainbar.y.label = "Intersection Size",
  sets.x.label = "LR Pairs per Sample"
)
dev.off()

#############################################################################
########################### 4. Cell-cell communication #####################
#############################################################################

run_CCC <- function(sample_ids, location_paths) {
  
  LR_results <- future_pmap(
    list(sample_ids, location_paths),
    function(sample_id, location_path) {
      LR2CellComm(
        obj = obj_sub,
        sample_id = sample_id,
        location_path = location_path,
        distance = 6,
        LRpair = LRpair
      )
    }
  )
  
  names(LR_results) <- sample_ids
  
  lr_all <- bind_rows(lapply(names(LR_results), function(sid) {
    df <- LR_results[[sid]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df %>% mutate(sample_id = sid)
  }))
  
  lr_summary <- lr_all %>%
    group_by(ligand, receptor, sender, receiver) %>%
    summarise(
      n_samples = n_distinct(sample_id),
      sampleID = paste(unique(sample_id), collapse = ";"),
      mean_strength = mean(interaction_strength, na.rm = TRUE),
      fisher_p = tryCatch(
        metap::sumlog(adjp[!is.na(adjp)])$p,
        error = function(e) NA
      ),
      .groups = "drop"
    ) %>%
    filter(n_samples >= 2)
  
  return(lr_summary)
}

#############################################################################
########################### 5. Subchondral CCC #############################
#############################################################################

sub_samples <- c("hd_68_1","hd_68_2","hd_82_2")

sub_summary <- run_CCC(
  sub_samples,
  sample_info$location_path[sample_info$sample_id %in% sub_samples]
)

sub_plot <- sub_summary %>%
  mutate(
    interaction = paste0(sender, "->", receiver),
    LR_pair = paste0(ligand, "-", receptor),
    log_p = -log10(fisher_p)
  )

pdf("results/CCC_Subchondral.pdf", width = 6, height = 3)
ggplot(sub_plot, aes(interaction, LR_pair)) +
  geom_point(aes(size = mean_strength, fill = log_p),
             shape = 21, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()

#############################################################################
########################### 6. Distant CCC #################################
#############################################################################

dis_samples <- c("hd_82_1","visium_12")

dis_summary <- run_CCC(
  dis_samples,
  sample_info$location_path[sample_info$sample_id %in% dis_samples]
)

dis_plot <- dis_summary %>%
  mutate(
    interaction = paste0(sender, "->", receiver),
    LR_pair = paste0(ligand, "-", receptor),
    log_p = -log10(fisher_p)
  )

pdf("results/CCC_Distant.pdf", width = 6, height = 4)
ggplot(dis_plot, aes(interaction, LR_pair)) +
  geom_point(aes(size = mean_strength, fill = log_p),
             shape = 21, color = "black") +
  scale_fill_gradient2(low="blue", mid="white", high="red") +
  theme_minimal()
dev.off()
