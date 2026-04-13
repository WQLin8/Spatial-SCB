load_hd_sample <- function(sample_name, batch, bone_loc, group,
                           has_capture = TRUE,
                           trab_col = "trab.segmentation") {
  
  message("Loading: ", sample_name, " batch", batch)
  
  # Load data
  obj <- Load10X_Spatial(
    file.path(DATA_DIR, paste0(sample_name, "_batch", batch, "/outs")),
    bin.size = 56
  )
  
  # Load trab annotation
  trab_file <- file.path(DATA_DIR, "trab_location",
                         paste0("trab-", gsub("Human_HD_", "", sample_name),
                                "-batch", batch, ".csv"))
  
  trab <- read.csv(trab_file)
  rownames(trab) <- trab$Barcode
  
  # Fix empty annotation
  trab[[trab_col]][trab[[trab_col]] == ""] <- "marrow"
  
  # Optional: filter captured area
  if (has_capture) {
    cap_file <- file.path(DATA_DIR, "trab_location",
                          paste0("capture_area-",
                                 gsub("Human_HD_", "", sample_name),
                                 "-", batch, ".csv"))
    
    cap <- read.csv(cap_file)
    rownames(cap) <- cap$Barcode
    
    captured_cells <- rownames(cap)[cap$capture_area == "captured"]
    obj <- subset(obj, cells = captured_cells)
  }
  
  # Add metadata
  obj$RegionType <- trab[colnames(obj), trab_col]
  obj <- RenameCells(obj, add.cell.id = paste0(sample_name, "_", batch))
  obj$sampleID <- paste0(sample_name, "_", batch)
  obj$BoneLocation <- bone_loc
  obj$group <- group
  
  return(obj)
}


## ---- Cell state identification ----
identify_cell_state <- function(obj, markers, name, quantile_cutoff = 0.9998) {
  
  markers <- intersect(markers, rownames(obj))
  
  # UCell
  obj_ucell <- AddModuleScore_UCell(obj,
                                    features = list(tmp = markers),
                                    name = paste0(name, "_UCell"))
  
  # Seurat score
  obj_seurat <- AddModuleScore(obj,
                               features = list(markers),
                               name = paste0(name, "_Score"))
  
  # Threshold
  thr_ucell  <- quantile(obj_ucell[[paste0(name, "_UCell")]], quantile_cutoff)
  thr_seurat <- quantile(obj_seurat[[paste0(name, "_Score1")]], quantile_cutoff)
  
  ucell_pos  <- colnames(obj_ucell)[obj_ucell[[paste0(name, "_UCell")]] > thr_ucell]
  seurat_pos <- colnames(obj_seurat)[obj_seurat[[paste0(name, "_Score1")]] > thr_seurat]
  
  common <- intersect(ucell_pos, seurat_pos)
  
  obj[[name]] <- colnames(obj) %in% common
  obj[[name]] <- factor(obj[[name]], levels = c(FALSE, TRUE))
  
  message(name, " spots: ", sum(obj[[name]] == TRUE))
  
  return(obj)
}

## ---- Spatial plotting ----
plot_spatial_binary <- function(obj, group, filename) {
  
  sample_ids <- unique(obj$sampleID)
  plots <- list()
  
  for (i in seq_along(sample_ids)) {
    sample_obj <- subset(obj, sampleID == sample_ids[i])
    aspect_ratio <- get_image_aspect_ratio(sample_obj)
    
    plots[[i]] <- SpatialDimPlot(
      sample_obj,
      group.by = group,
      pt.size.factor = ifelse(i == 1, 5, 3.2),
      shape = 22
    ) + theme(aspect.ratio = aspect_ratio)
  }
  
  p <- wrap_plots(plots, ncol = 2)
  
  ggsave(file.path(OUT_DIR, filename), p, width = 8, height = 8)
}

## ---- Differential analysis + volcano ----
run_volcano <- function(obj, group, filename) {
  
  de <- FindMarkers(obj,
                    ident.1 = "TRUE",
                    ident.2 = "FALSE",
                    group.by = group,
                    min.pct = 0.1)
  
  de$p_val_adj[de$p_val_adj < 1e-150] <- 1e-150
  
  up <- subset(de, avg_log2FC > 0 & p_val_adj < 0.05)
  down <- subset(de, avg_log2FC < 0 & p_val_adj < 0.05)
  
  labeldata <- rbind(head(up, 20), head(down, 10))
  
  p <- ggplot(de, aes(avg_log2FC, -log10(p_val_adj))) +
    geom_point(aes(color = -log10(p_val_adj), size = -log10(p_val_adj))) +
    scale_color_gradientn(colours = rev(brewer.pal(11, 'RdYlBu'))) +
    geom_text_repel(data = labeldata,
                    aes(label = rownames(labeldata)),
                    size = 2.5) +
    theme_bw()
  
  ggsave(file.path(OUT_DIR, filename), p, width = 7, height = 5)
  
  return(de)
}

## ---- GO enrichment ----
run_go <- function(de, filename) {
  
  sig <- rownames(subset(de, avg_log2FC > 0.5 & p_val_adj < 0.05))
  
  ego <- enrichGO(
    gene = sig,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP"
  )
  
  df <- as.data.frame(ego) %>% slice_min(pvalue, n = 20)
  
  p <- ggplot(df, aes(-log10(pvalue), reorder(Description, -log10(pvalue)))) +
    geom_col(aes(fill = Count)) +
    theme_classic()
  
  ggsave(file.path(OUT_DIR, filename), p, width = 8, height = 6)
}

# ===================== 4. Build Reference =====================
build_reference <- function(sc, obj_spatial, label = "cell_type_lvl2") {
  
  common_genes <- intersect(rownames(sc), rownames(obj_spatial))
  
  sc_counts  <- GetAssayData(sc, slot = "counts")[common_genes, ]
  sp_counts  <- GetAssayData(obj_spatial, slot = "counts")[common_genes, ]
  
  # ---- Extract special cell states ----
  adipo_idx <- which(obj_spatial$Adipocyte == "TRUE")
  oc_idx    <- which(obj_spatial$Osteoclast == "TRUE")
  
  adipo_counts <- sp_counts[, adipo_idx]
  oc_counts    <- sp_counts[, oc_idx]
  
  # ---- Merge ----
  merged_counts <- cbind(sc_counts, adipo_counts, oc_counts)
  
  # ---- Labels ----
  Idents(sc) <- label
  cluster_sc <- sc[[label]][,1]
  
  cluster_ad <- rep("Adipocyte", ncol(adipo_counts))
  cluster_oc <- rep("Osteoclast", ncol(oc_counts))
  
  names(cluster_sc) <- colnames(sc_counts)
  names(cluster_ad) <- colnames(adipo_counts)
  names(cluster_oc) <- colnames(oc_counts)
  
  cluster_all <- c(cluster_sc, cluster_ad, cluster_oc)
  cluster_all <- factor(cluster_all[colnames(merged_counts)])
  
  # ---- nUMI ----
  nUMI <- Matrix::colSums(merged_counts)
  
  reference <- Reference(
    counts = merged_counts,
    cell_types = cluster_all,
    nUMI = nUMI,
    min_UMI = 50
  )
  
  return(reference)
}

run_rctd_pipeline <- function(query, reference, cores = 20) {
  
  rctd <- create.RCTD(
    query,
    reference,
    max_cores = cores,
    UMI_min = 50,
    CELL_MIN_INSTANCE = 1
  )
  
  rctd <- run.RCTD(rctd, doublet_mode = "full")
  
  return(rctd)
}

build_query <- function(obj) {
  
  counts <- GetAssayData(obj, slot = "counts", assay = "Spatial.056um")
  
  image_names <- names(obj@images)
  
  coords_list <- lapply(image_names, function(img) {
    pos <- GetTissueCoordinates(obj, image = img)
    colnames(pos) <- c("x", "y")
    return(pos)
  })
  
  coords <- do.call(rbind, coords_list)
  coords <- coords[, c("x", "y")]
  
  query <- SpatialRNA(coords, counts, colSums(counts))
  
  return(query)
}

# ===================== 4. Core Identification =====================

identify_Core <- function(obj, celltype = "Osteoblast", quantile_cutoff = 0.85) {
  
  DefaultAssay(obj) <- "RCTD"
  
  mat <- obj[["RCTD"]]@data
  
  prop <- mat[celltype, ]
  obj[[paste0(celltype, "_prop")]] <- prop
  
  thr <- quantile(prop, quantile_cutoff)
  
  flag_name <- paste0(celltype, "_Core")
  
  obj[[flag_name]] <- prop >= thr
  obj[[flag_name]] <- factor(obj[[flag_name]], levels = c(FALSE, TRUE))
  
  message(celltype, " Core spots: ", sum(obj[[flag_name]] == TRUE))
  
  return(obj)
}

# ===================== 5. Spatial Visualization =====================

plot_Core <- function(obj, group, filename) {
  
  sample_ids <- unique(obj$sampleID)
  plots <- list()
  
  for (i in seq_along(sample_ids)) {
    
    sample_obj <- subset(obj, sampleID == sample_ids[i])
    aspect_ratio <- get_image_aspect_ratio(sample_obj)
    
    pt_size <- ifelse(sample_ids[i] == "visium_12", 2.5,
                      ifelse(i == 1, 5, 3.2))
    
    shape_type <- ifelse(sample_ids[i] == "visium_12", 21, 22)
    
    plots[[i]] <- SpatialDimPlot(
      sample_obj,
      group.by = group,
      pt.size.factor = pt_size,
      shape = shape_type
    ) + theme(aspect.ratio = aspect_ratio)
  }
  
  p <- wrap_plots(plots, ncol = 2)
  
  ggsave(file.path(OUT_DIR, filename), p, width = 8, height = 10)
}

# ===================== 7. Core-Halo Visualization =====================

plot_ring <- function(obj, filename) {
  
  sample_ids <- unique(obj$sampleID)
  plots <- list()
  
  for (i in seq_along(sample_ids)) {
    
    sample_obj <- subset(obj, sampleID == sample_ids[i])
    aspect_ratio <- get_image_aspect_ratio(sample_obj)
    
    pt_size <- ifelse(sample_ids[i] == "visium_12", 2.5,
                      ifelse(i == 1, 5, 3.2))
    
    shape_type <- ifelse(sample_ids[i] == "visium_12", 21, 22)
    
    plots[[i]] <- SpatialDimPlot(
      sample_obj,
      group.by = "ring",
      pt.size.factor = pt_size,
      shape = shape_type
    ) +
      labs(fill = "Region") +
      theme(
        aspect.ratio = aspect_ratio,
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)
      )
  }
  
  p <- wrap_plots(plots, ncol = 2)
  
  ggsave(file.path(OUT_DIR, filename), p, width = 8, height = 10)
}


get_hallmark_gene_sets <- function(expr_mat) {
  hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
  gset <- split(x = hallmark_sets$gene_symbol, f = hallmark_sets$gs_name)
  gset <- lapply(gset, intersect, rownames(expr_mat))
  gset[sapply(gset, length) >= 10 & sapply(gset, length) <= 5000]
}

get_go_bp_gene_sets <- function(expr_mat) {
  human_go_bp <- msigdbr(species = "Homo sapiens",
                         collection = "C5",
                         subcategory = "BP") %>%
    dplyr::select(gs_name, gene_symbol)
  gset <- split(x = human_go_bp$gene_symbol, f = human_go_bp$gs_name)
  gset <- lapply(gset, intersect, rownames(expr_mat))
  gset[sapply(gset, length) >= 10 & sapply(gset, length) <= 5000]
}

run_go_enrichment <- function(genes, group_name, top_n = 10) {
  ego <- enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )

  ego_df <- data.frame(ego)
  if (nrow(ego_df) == 0) return(NULL)

  ego_df$Group <- group_name

  if (group_name == "Other") {
    ego_df <- ego_df %>% slice_min(order_by = p.adjust, n = top_n)
  } else {
    ego_df <- ego_df %>%
      arrange(p.adjust, BgRatio) %>%
      filter(row_number(p.adjust) <= top_n)
  }

  ego_df %>%
    separate(GeneRatio, into = c("num", "denom"), sep = "/", convert = TRUE) %>%
    mutate(
      GeneRatio   = num / denom,
      Description = factor(Description, levels = Description[order(-log10(p.adjust))])
    ) %>%
    dplyr::select(-num, -denom) %>%
    arrange(-log10(p.adjust)) %>%
    mutate(Description = factor(Description, levels = Description))
}

plot_go_dotplot <- function(go_df, outfile, width = 7, height = 6) {
  if (is.null(go_df) || nrow(go_df) == 0) return(invisible(NULL))

  go_df$Description <- as.factor(go_df$Description)
  go_df$Description <- forcats::fct_inorder(go_df$Description)

  p <- ggplot(go_df, aes(Group, Description)) +
    geom_point(aes(color = p.adjust, size = GeneRatio)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.6, vjust = 0.5),
      legend.title = element_text(margin = margin(b = 8))
    ) +
    scale_color_gradient(
      high   = "#6699CC",
      low    = "#CC3333",
      name   = "p.adjust",
      breaks = c(min(go_df$p.adjust), max(go_df$p.adjust)),
      labels = c("low", "high")
    ) +
    labs(x = NULL, y = NULL) +
    guides(size = guide_legend(order = 1))

  ggsave(outfile, p, width = width, height = height)
}

make_deg_table <- function(deg_df, fdr_cutoff = 0.05, logfc_cutoff = 0.1) {
  deg_df <- data.frame(deg_df)
  deg_df$gene_name <- rownames(deg_df)
  deg_df$Sig <- ifelse(
    deg_df$p_val_adj < fdr_cutoff & abs(deg_df$avg_log2FC) >= logfc_cutoff,
    ifelse(deg_df$avg_log2FC > 0, "Up", "Down"),
    "noSig"
  )
  deg_df
}

plot_volcano_with_labels <- function(deg_df, title_text, outfile,
                                     fdr_cutoff = 0.05,
                                     logfc_cutoff = 0.1,
                                     xlim_range = c(-15, 10),
                                     n_label = 25) {
  plot_df <- make_deg_table(deg_df, fdr_cutoff, logfc_cutoff)

  p1 <- ggplot(plot_df, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = Sig)) +
    geom_point(alpha = 0.65, size = 2) +
    scale_color_manual(values = c("Down" = "#546de5", "noSig" = "#d2dae2", "Up" = "#ff4757")) +
    xlim(xlim_range) +
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), lty = 4, col = "black", lwd = 0.8) +
    geom_hline(yintercept = -log10(fdr_cutoff), lty = 4, col = "black", lwd = 0.8) +
    labs(x = "log2FC", y = "-log10FDR") +
    ggtitle(title_text) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position = "right",
      legend.title = element_blank(),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 12, color = "black")
    )

  p2 <- p1 + geom_text_repel(
    data = plot_df %>%
      filter(p_val_adj < fdr_cutoff & abs(avg_log2FC) >= logfc_cutoff) %>%
      mutate(direction = ifelse(avg_log2FC > 0, "Up", "Down")) %>%
      group_by(direction) %>%
      arrange(p_val_adj, .by_group = TRUE) %>%
      slice_head(n = n_label) %>%
      ungroup(),
    aes(label = gene_name),
    size = 3.5,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"),
    segment.color = "black",
    show.legend = FALSE
  )

  ggsave(outfile, p2, height = 6, width = 8)
  invisible(plot_df)
}

run_gsva_limma <- function(seu_obj, assay = SPATIAL_ASSAY,
                           gene_set_type = c("hallmark", "go_bp")) {
  gene_set_type <- match.arg(gene_set_type)

  DefaultAssay(seu_obj) <- assay
  Idents(seu_obj) <- "BoneLocation"

  expr <- GetAssayData(seu_obj, assay = assay, slot = "data")
  expr <- as.matrix(expr)
  expr <- expr[rowSums(expr) > 0, ]

  if (any(duplicated(rownames(expr)))) {
    expr <- rowsum(expr, group = rownames(expr), reorder = FALSE)
  }

  gset <- switch(
    gene_set_type,
    hallmark = get_hallmark_gene_sets(expr),
    go_bp    = get_go_bp_gene_sets(expr)
  )

  gsva_res <- gsva(
    expr,
    gset.idx.list = gset,
    kcdf          = "Gaussian",
    method        = "gsva",
    mx.diff       = TRUE,
    parallel.sz   = 20
  )

  bone_loc <- Idents(seu_obj)
  group_list <- data.frame(
    celltype = colnames(gsva_res),
    group    = ifelse(bone_loc[colnames(gsva_res)] == "Subchondral", "case", "control")
  )

  design <- model.matrix(~ 0 + factor(group_list$group))
  colnames(design) <- levels(factor(group_list$group))
  rownames(design) <- colnames(gsva_res)

  contrast.matrix <- makeContrasts(case - control, levels = design)

  fit  <- lmFit(gsva_res, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)

  topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
}

plot_gsva_bar <- function(diff_df, outfile) {
  dat_plot <- data.frame(id = row.names(diff_df), t = diff_df$t)
  dat_plot$id <- str_replace(dat_plot$id, "GOBP_", "")
  dat_plot$threshold <- factor(
    ifelse(dat_plot$t > -2, ifelse(dat_plot$t >= 2, "Up", "NoSig"), "Down"),
    levels = c("Up", "Down", "NoSig")
  )

  dat_plot <- dat_plot %>% arrange(t)
  dat_plot$id <- factor(dat_plot$id, levels = dat_plot$id)

  p <- ggplot(data = dat_plot, aes(x = id, y = t, fill = threshold)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("Up" = "#36638a", "NoSig" = "#cccccc", "Down" = "#7bcd7b")) +
    geom_hline(yintercept = c(-2, 2), color = "white", size = 0.5, lty = "dashed") +
    xlab("") +
    ylab("t value of GSVA score, Case versus Control") +
    guides(fill = FALSE) +
    theme_prism(border = TRUE) +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )

  low1  <- dat_plot %>% filter(t < -2) %>% nrow()
  low0  <- dat_plot %>% filter(t < 0)  %>% nrow()
  high0 <- dat_plot %>% filter(t < 2)  %>% nrow()
  high1 <- nrow(dat_plot)

  if (low1 > 0) {
    p <- p + geom_text(data = dat_plot[1:low1, ], aes(x = id, y = 0.1, label = id),
                       hjust = 0, color = "black")
  }
  if (low0 > low1) {
    p <- p + geom_text(data = dat_plot[(low1 + 1):low0, ], aes(x = id, y = 0.1, label = id),
                       hjust = 0, color = "grey")
  }
  if (high0 > low0) {
    p <- p + geom_text(data = dat_plot[(low0 + 1):high0, ], aes(x = id, y = -0.1, label = id),
                       hjust = 1, color = "grey")
  }
  if (high1 > high0) {
    p <- p + geom_text(data = dat_plot[(high0 + 1):high1, ], aes(x = id, y = -0.1, label = id),
                       hjust = 1, color = "black")
  }

  ggsave(outfile, p, width = 16, height = 10)
}

run_region_specific_go <- function(seu_obj, logfc_threshold_findall = 0.1,
                                   sig_logfc_cutoff = 0.5) {
  markers <- FindAllMarkers(
    seu_obj,
    assays          = "data",
    only.pos        = FALSE,
    min.pct         = 0.01,
    logfc.threshold = logfc_threshold_findall
  )

  sig_dge_all <- subset(markers, p_val_adj < 0.05 & abs(avg_log2FC) > sig_logfc_cutoff)

  df_core  <- run_go_enrichment(subset(sig_dge_all, avg_log2FC > 0 & cluster == "Core")$gene,  "Core")
  df_halo  <- run_go_enrichment(subset(sig_dge_all, avg_log2FC > 0 & cluster == "Halo")$gene,  "Halo")
  df_other <- run_go_enrichment(subset(sig_dge_all, avg_log2FC > 0 & cluster == "Other")$gene, "Other")

  list(markers = markers, go_df = bind_rows(df_core, df_halo, df_other))
}

prepare_rctd_long <- function(seu_obj, assay = RCTD_ASSAY) {
  DefaultAssay(seu_obj) <- assay

  core  <- subset(seu_obj, ring == "Core")
  halo  <- subset(seu_obj, ring == "Halo")
  other <- subset(seu_obj, ring == "Other")

  n <- length(rownames(core@assays[[assay]]@data))

  bind_region_block <- function(x) {
    df <- rbind(
      x@assays[[assay]]@data,
      x$sampleID[colnames(x@assays[[assay]]@data)],
      x$ring[colnames(x@assays[[assay]]@data)]
    )
    rownames(df)[(n + 1):(n + 2)] <- c("sampleID", "ring")
    df
  }

  df_all <- cbind(bind_region_block(core), bind_region_block(halo), bind_region_block(other))

  expr_data <- df_all[1:n, ]
  expr_data <- tibble::rownames_to_column(as.data.frame(expr_data), var = "celltype")

  cell_data <- expr_data %>%
    pivot_longer(cols = -celltype, names_to = "spot", values_to = "proportion")

  meta_info <- as.data.frame(df_all)[(n + 1):(n + 2), ] %>%
    tibble::rownames_to_column("meta") %>%
    pivot_longer(cols = -meta, names_to = "spot", values_to = "value") %>%
    pivot_wider(names_from = meta, values_from = value)

  df_long <- left_join(cell_data, meta_info, by = "spot")
  df_long$proportion <- as.numeric(df_long$proportion)
  df_long$ring <- factor(df_long$ring, levels = RING_LEVELS)
  df_long
}

plot_single_celltype_violin <- function(df_long, celltype_name, outfile) {
  data <- df_long[df_long$celltype == celltype_name, ]
  max_y <- max(data$proportion, na.rm = TRUE)
  y_limit <- max_y * 1.5
  colours <- c("#1b9e77", "#d95f02", "#7570b3")

  p <- ggplot(data, aes(x = ring, y = proportion, fill = ring, color = ring)) +
    scale_color_manual(values = rev(colours), name = "ring") +
    scale_fill_manual(values = rev(colours), name = "ring") +
    geom_half_violin(
      position = position_nudge(x = 0.15, y = 0),
      side     = "R",
      adjust   = 2.5,
      trim     = TRUE,
      color    = NA,
      na.rm    = TRUE,
      alpha    = 0.8
    ) +
    geom_boxplot(outlier.shape = NA, width = 0.2, alpha = 0.7) +
    xlab("") +
    ylab("Proportion") +
    scale_y_continuous(limits = c(0, y_limit)) +
    ggpubr::stat_compare_means(
      comparisons = list(c("Core", "Halo"), c("Core", "Other"), c("Halo", "Other")),
      method      = "wilcox.test",
      label       = "p.signif"
    ) +
    theme_classic() +
    theme(
      axis.text.x   = element_blank(),
      axis.text     = element_text(size = 12),
      axis.ticks.x  = element_blank(),
      plot.margin   = unit(c(1, 1, 1, 1), "cm"),
      plot.title    = element_text(hjust = 0.5)
    ) +
    ggtitle(celltype_name)

  ggsave(outfile, p, width = 6, height = 5)
}

plot_all_celltype_boxbreak <- function(df_long, outfile) {
  colours <- c("#F07B70", "#32B93D", "#6F9AFD")

  p <- ggplot(df_long, aes(x = celltype, y = proportion, fill = ring)) +
    geom_boxplot(
      width         = 0.4,
      alpha         = 0.6,
      position      = position_dodge(0.75),
      show.legend   = TRUE,
      outlier.shape = NA
    ) +
    scale_fill_manual(values = colours, name = "Area") +
    scale_y_break(c(0.0003, 0.005), scales = "free") +
    scale_y_break(c(0.05, 0.1), scales = "free") +
    labs(x = "Cell type", y = "proportion") +
    theme_minimal() +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 1),
      panel.border  = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    )

  ggsave(outfile, p, width = 12, height = 4, device = cairo_pdf)
}
