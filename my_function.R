res_dir <- "/public/home/GENE_proc/wth/res/"
options(future.globals.maxSize = 128*1024^3)
set.seed(520)
my_print <- function(..., sep = ""){  # 输出时间和信息
  formatted_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  str <- paste(as.vector(list(...)), collapse = sep)
  message(paste0(formatted_time, "  ", str))
}
identity_seu <- function(seu, target_col = "celltype"){  # 判断该Seurat对象是哪种细胞类型
  celltype_name <- unique(as.character(seu[[target_col]][, 1]))
  celltype_name <- ifelse(length(celltype_name) == 1, celltype_name[1], "All")
  return(celltype_name)
}
p_to_star <- function(p){
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}
fmt_p <- function(p){
  dplyr::case_when(
    is.na(p) ~ "p=NA",
    p < 0.001 ~ paste0("p=", format(p, scientific = TRUE, digits = 2)),
    TRUE ~ paste0("p=", sprintf("%.3f", p))
  )
}
align_to_cells <- function(mat, cells){
  mat <- as(mat, "dgCMatrix")
  extra_cells <- setdiff(colnames(mat), cells)
  if (length(extra_cells) > 0) {  # 如果herv计数结果中有多余细胞：删掉
    mat <- mat[, setdiff(colnames(mat), extra_cells), drop = FALSE]
  }
  miss_cells <- setdiff(cells, colnames(mat))
  if (length(miss_cells) > 0) {  # 如果herv计数结果中缺少某些细胞：补0
    zero <- Matrix(
      0,
      nrow = nrow(mat),
      ncol = length(miss_cells),
      sparse = TRUE,
      dimnames = list(rownames(mat), miss_cells)
    )
    mat <- cbind(mat, zero)
  }
  mat[, cells, drop = FALSE]
}
replace_assay <- function(seu, assay_name, counts){
  counts <- as(counts, "dgCMatrix")
  seu[[assay_name]] <- CreateAssayObject(counts = counts)
  seu
}
normalize_herv_by_rna_depth <- function(seu, assay, scale.factor = 10000) {
  rna_counts <- GetAssayData(seu, assay = "RNA",  layer = "counts")
  herv_counts <- GetAssayData(seu, assay = assay, layer = "counts")
  lib_rna <- Matrix::colSums(rna_counts)
  herv_norm <- t(t(herv_counts) / lib_rna) * scale.factor
  herv_norm@x <- log1p(herv_norm@x)
  seu <- SetAssayData(seu, assay = assay, layer = "data", new.data = herv_norm)
  return(seu)
}
run_go_one <- function(gene_vec, universe_genes = NULL, ont = "BP") {
  gene_vec <- unique(na.omit(gene_vec))
  my_print(paste0("  - n_gene:", length(gene_vec), " n_universe_gene:", length(universe_genes)))
  ego <- enrichGO(
    gene          = gene_vec,
    universe      = universe_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE
  )
  if (is.null(ego)){
    my_print("  - no res")
    return(NULL)
  }
  ego
}
mat_to_df <- function(mat, value_name) {
  df <- as.data.frame(as.table(as.matrix(mat)), stringsAsFactors = FALSE)
  colnames(df) <- c("trait", "module", value_name)
  df
}
residualize_vec <- function(y, meta) {
  fit <- lm(y ~ log10_nCount_RNA + percent_mito, data = meta)
  resid(fit)
}
herv_gene_module_corr <- function(seu, herv_assay, target_herv, is_plot = T){
  celltype_name <- identity_seu(seu)
  DefaultAssay(seu) <- "RNA"
  seu$AD_binary <- ifelse(seu$group == "AD", 1, 0)
  mat <- LayerData(seu, assay = herv_assay, layer = "data")[target_herv, , drop = FALSE]
  mat <- t(as.matrix(mat))
  colnames(mat) <- paste0("herv_", make.names(target_herv))
  seu <- AddMetaData(seu, metadata = as.data.frame(mat))
  trait_cols <- c(
    paste0("herv_", make.names(target_herv)),
    "AD_binary"
  )
  seu <- ModuleTraitCorrelation(
    seu,
    traits = trait_cols,
    group.by = "celltype",
    wgcna_name = "gene_net"
  )
  mt <- GetModuleTraitCorrelation(seu)
  # pdf(file.path(res_dir, paste(celltype_name, herv_assay, "trait_corr.pdf", sep = "_")), width = 16, height = 16)
  # p <- PlotModuleTraitCorrelation(seu, label = "fdr", label_symbol = "stars")
  # print(p)
  # dev.off()
  if(is_plot) return(PlotModuleTraitCorrelation(seu, label = "fdr", label_symbol = "stars", combine = FALSE)[[2]])
  else return(seu)
}
get_plot <- function(net_obj_name_list, target_herv_assay, target_herv){
  plot_list <- lapply(net_obj_name_list, function(net_obj_name){
    my_print(net_obj_name)
    my_print("- reading RData")
    load(file.path(res_dir, net_obj_name))
    my_print("- calc corr")
    herv_gene_module_corr(seu, target_herv_assay, target_herv, T)
  })
  pdf(file.path(res_dir, paste0(target_herv_assay, "_subcluster_network_corr.pdf")), width = 32, height = 24)
  print(wrap_plots(plot_list, ncol = 3))
  dev.off()
  return(plot_list)
}
run_herv_geneset_coexpr_one_celltype <- function(seu_ct, herv_assay, herv_features, gene_features, agg = "mean", is_donor = FALSE, donor_var = "orig.ident") {
  ct <- identity_seu(seu_ct)
  level_name <- ifelse(is_donor, "donor", "cell")
  res_name <- paste(ct, herv_features[1], gene_features[1], agg, level_name, sep = " ")
  my_print(res_name, ":")
  herv_use <- intersect(rownames(seu_ct[[herv_assay]]), herv_features)
  gene_use <- intersect(rownames(seu_ct[["RNA"]]), gene_features)
  my_print("- hERV features used:", length(herv_use), " genes used:", length(gene_use))
  if(length(herv_use)==0 | length(gene_use)==0) return(NULL)
  my_print("- calculating residualize_vec")
  # 构建模型，设置协变量
  meta <- seu_ct@meta.data %>%
    mutate(log10_nCount_RNA = log10(nCount_RNA + 1))
  herv_mat <- GetAssayData(seu_ct, assay = herv_assay, layer = "data")[herv_use, , drop = FALSE]
  rna_mat  <- GetAssayData(seu_ct, assay = "RNA",  layer = "data")[gene_use, , drop = FALSE]
  herv_counts <- GetAssayData(seu_ct, assay = herv_assay, layer = "counts")[herv_use, , drop = FALSE]
  rna_counts  <- GetAssayData(seu_ct, assay = "RNA", layer = "counts")[gene_use, , drop = FALSE]
  herv_pct <- Matrix::rowMeans(herv_counts > 0)
  gene_set_pct <- mean(Matrix::colSums(rna_counts > 0) > 0)
  if (agg == "mean") {
    gene_score <- Matrix::colMeans(rna_mat)
  } else {
    gene_score <- Matrix::colSums(rna_mat)
  }
  gene_set_pct <- mean(  # 该组基因中任意一个表达>0即记为检测到
    Matrix::colSums(GetAssayData(seu_ct, assay = "RNA", layer = "counts")[gene_use, , drop = FALSE] > 0) > 0
  )
  # 残差化
  if (!is_donor) {
    my_print("- cell-level residualize")
    herv_resid <- t(apply(as.matrix(herv_mat), 1, residualize_vec, meta = meta))
    gene_resid <- residualize_vec(as.numeric(gene_score), meta = meta)
    n_used <- ncol(herv_resid)
  } else {
    my_print("- donor-level aggregation")
    donors <- as.character(meta[[donor_var]])
    donor_levels <- unique(donors)
    meta_donor <- meta %>%  # 协变量：对每个样本取平均
      mutate(donor = donors) %>%
      group_by(donor) %>%
      summarise(
        log10_nCount_RNA = mean(log10_nCount_RNA, na.rm = TRUE),
        percent_mito = mean(percent_mito, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      as.data.frame()
    rownames(meta_donor) <- meta_donor$donor
    herv_donor <- sapply(donor_levels, function(d) {  # hERV表达：每个样本内按细胞平均
      idx <- which(donors == d)
      Matrix::rowMeans(herv_mat[, idx, drop = FALSE])
    })
    if (is.vector(herv_donor)) {
      herv_donor <- matrix(herv_donor, nrow = length(herv_use), dimnames = list(herv_use, donor_levels))
    }
    colnames(herv_donor) <- donor_levels
    gene_donor <- sapply(donor_levels, function(d) {  # gene set表达：每个样本内按细胞平均
      idx <- which(donors == d)
      mean(gene_score[idx], na.rm = TRUE)
    })
    gene_donor <- as.numeric(gene_donor)
    names(gene_donor) <- donor_levels
    meta_donor <- meta_donor[colnames(herv_donor), , drop = FALSE]
    my_print("- donor-level residualize")
    herv_resid <- t(apply(as.matrix(herv_donor), 1, residualize_vec, meta = meta_donor))
    gene_resid <- residualize_vec(as.numeric(gene_donor), meta = meta_donor)
    n_used <- ncol(herv_resid)
  }
  my_print("- calculating corr")
  # 计算rho+p值
  res_list <- vector("list", length = nrow(herv_resid))
  for (i in seq_len(nrow(herv_resid))) {
    x <- as.numeric(herv_resid[i, ])
    y <- as.numeric(gene_resid)
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    if (length(x) < 10 || sd(x) == 0 || sd(y) == 0) {
      res_list[[i]] <- data.frame(
        feature = rownames(herv_resid)[i],
        gene_set = paste(gene_use, collapse = ";"),
        n_gene = length(gene_use),
        agg = agg,
        level = level_name,
        rho = NA_real_,
        p_val = NA_real_,
        n = length(x),
        herv_pct = herv_pct[rownames(herv_resid)[i]],
        gene_set_pct = gene_set_pct
      )
    } else {
      ct_res <- suppressWarnings(
        cor.test(x, y, method = "spearman", exact = FALSE)
      )
      res_list[[i]] <- data.frame(
        feature = rownames(herv_resid)[i],
        gene_set = paste(gene_use, collapse = ";"),
        n_gene = length(gene_use),
        agg = agg,
        level = level_name,
        rho = unname(ct_res$estimate),
        p_val = ct_res$p.value,
        n = length(x),
        herv_pct = herv_pct[rownames(herv_resid)[i]],
        gene_set_pct = gene_set_pct
      )
    }
  }
  bind_rows(res_list) %>%
    mutate(
      p_adj = p.adjust(p_val, method = "BH")
    ) %>%
    arrange(p_adj, p_val, desc(abs(rho)))
}
get_gene_herv_corr_res <- function(seu, herv_assay, target_herv, gene_panel, res_name = NA, agg = "mean", is_donor = FALSE, donor_var = "orig.ident", write_csv = T){
  if(is.na(res_name)) res_name <- identity_seu(seu)
  loop_tag <- names(gene_panel)
  if(is.null(loop_tag)) loop_tag <- 1:length(gene_panel)
  res <- lapply(loop_tag, function(sig) {
    row <- run_herv_geneset_coexpr_one_celltype(
      seu_ct = seu,
      herv_assay = herv_assay,
      herv_features = target_herv,
      gene_features = gene_panel[[sig]],
      agg = agg,
      is_donor = is_donor,
      donor_var = donor_var
    )
    if(is.null(row)) return(NULL)
    else return(row %>% mutate(signature = sig, .before = 1))
  }) %>% bind_rows()
  if(write_csv) write.csv(res, file = file.path(res_dir, paste0(res_name, "_geneset_", herv_assay, ifelse(is_donor, "_donor", ""), "_corr.csv")), row.names = F)
  else return(res)
}
get_gene_herv_corr_res_plus <- function(seu, herv_assay, gene_panel_list, res_name = NA, agg = "mean", is_donor = FALSE, donor_var = "orig.ident", write_csv = T){
  if(is.na(res_name)) res_name <- identity_seu(seu)
  res <- lapply(names(gene_panel_list), function(herv){
    get_gene_herv_corr_res(
      seu,
      herv_assay = herv_assay,
      target_herv = herv,
      gene_panel = gene_panel_list[[herv]],
      res_name = res_name,
      agg = agg,
      is_donor = is_donor,
      donor_var = donor_var, 
      write_csv = F
    )
  }) %>% bind_rows()
  if(write_csv) write.csv(res, file = file.path(res_dir, paste0(res_name, "_geneset_", herv_assay, ifelse(is_donor, "_donor", ""), "_corr_plus.csv")), row.names = F)
  else return(res)
}
plot_herv_gene_corr <- function(seu_ct, herv_assay, herv_feature, gene, res_name = NA, agg = "mean", is_donor = FALSE, donor_var = "orig.ident", is_plot = T) {
  if(is.na(res_name)) res_name <- identity_seu(seu_ct)
  my_print("calc corr: ", herv_feature[1], "-", gene[1])
  meta <- seu_ct@meta.data %>%
    mutate(log10_nCount_RNA = log10(nCount_RNA + 1))
  # hERV
  herv_mat <- GetAssayData(seu_ct, assay = herv_assay, layer = "data")
  if (!herv_feature %in% rownames(herv_mat)) stop("hERV feature not found in assay.")
  x_raw <- as.numeric(herv_mat[herv_feature, ])
  # gene/gene set
  rna_mat <- GetAssayData(seu_ct, assay = "RNA", layer = "data")
  if (length(gene) == 1 && gene %in% colnames(meta)) {
    y_raw <- meta[[gene]]
    gene_name <- gene
  } else {
    gene_use <- intersect(gene, rownames(rna_mat))
    if (length(gene_use) == 0) stop("gene not found in RNA assay.")
    y_mat <- rna_mat[gene_use, , drop = FALSE]
    y_raw <- if (agg == "mean") Matrix::colMeans(y_mat) else Matrix::colSums(y_mat)
    gene_name <- if (length(gene_use) == 1) gene_use else paste(gene_use, collapse = "+")
  }
  # cell-level/donor-level
  if (!is_donor) {
    x <- residualize_vec(x_raw, meta)
    y <- residualize_vec(as.numeric(y_raw), meta)
    df_plot <- data.frame(x = x, y = y)
  } else {
    donors <- as.character(meta[[donor_var]])
    meta_donor <- meta %>%
      mutate(donor = donors) %>%
      group_by(donor) %>%
      summarise(
        log10_nCount_RNA = mean(log10_nCount_RNA, na.rm = TRUE),
        percent_mito = mean(percent_mito, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      as.data.frame()
    rownames(meta_donor) <- meta_donor$donor
    donor_levels <- meta_donor$donor
    x_donor <- sapply(donor_levels, function(d) mean(x_raw[donors == d], na.rm = TRUE))
    y_donor <- sapply(donor_levels, function(d) mean(y_raw[donors == d], na.rm = TRUE))
    x <- residualize_vec(as.numeric(x_donor), meta_donor)
    y <- residualize_vec(as.numeric(y_donor), meta_donor)
    df_plot <- data.frame(x = x, y = y, donor = donor_levels)
  }
  ct_res <- suppressWarnings(cor.test(df_plot$x, df_plot$y, method = "spearman", exact = FALSE))
  my_print("- drawing")
  rho_txt <- round(unname(ct_res$estimate), 3)
  p_txt <- signif(ct_res$p.value, 3)
  if(!is_plot) return(list(df_plot, rho_txt, p_txt))
  ggplot(df_plot, aes(x = x, y = y)) +
    geom_point(size = 0.5, alpha = 0.35) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.7) +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = paste0("rho=", rho_txt, "\np=", p_txt),
      hjust = -0.1, vjust = 1.2,
      size = 4.5
    ) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    labs(
      title = paste0(res_name, ":    ", gene_name, " ~ ", herv_feature),
      x = herv_feature,
      y = gene_name
    )
}
plot_herv_gene_corr_df <- function(seu_list, plot_specs, ncol = 4, res_name = "", title = "", is_plot = F){
  plots <- vector("list", nrow(plot_specs))
  for(i in seq_len(nrow(plot_specs))) {
    ct_i <- plot_specs$ct[i]
    herv_i <- plot_specs$herv[i]
    gene_i <- plot_specs$gene[i]
    p <- plot_herv_gene_corr(
      seu_ct = seu_list[[ct_i]],
      herv_assay = "HERV",
      herv_feature = herv_i,
      gene = gene_i,
      res_name = ifelse(ct_i=="Mic2", "Mic cluster2", NA),
      is_donor = TRUE,
    ) +
      geom_point(size = 1, alpha = 0.85) +
      labs(x = NULL, y = NULL) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(size = 13.5, face = "bold"),
        axis.title = element_blank(),
        axis.text = element_text(size = 9)
      )
    plots[[i]] <- p
  }
  res <- wrap_plots(plots, ncol = ncol)
  if(title != "") res <- res + plot_annotation(title = title)
  if(!is_plot) return(res)
  pdf(file.path(res_dir, paste0("herv_gene_donor_corr", res_name, ".pdf")), width = 24, height = ((length(plots)-1)%/%4+1)*6)
  print(res)
  dev.off()
}
