library(ProjecTILs)
ref_CD8 <- readRDS("/net/isi-dcnl/ifs/user_data/lkwak/Data/Jibin/reference/SingleR_Ref/CD8T_human_ref_v1.rds")
list_of_ref.maps <- get.reference.maps(reference = "CD4")
ref_CD4 <- list_of_ref.maps[["human"]][["CD4"]]
ref_combined <- merge(ref_CD4, ref_CD8)

gdT_sig <- c("TRDC","TRGC1","TRGC2","KLRD1","GZMB","NKG7")
prolif_sig <- c("MKI67","TOP2A","STMN1","TYMS","PCNA","BIRC5")
Treg_sig <- c("FOXP3","IL2RA","CTLA4","TNFRSF4","TNFRSF18","IKZF2","LAG3","CD27")
Tex_sig <- c("PDCD1","LAG3","TIGIT","CTLA4", "HAVCR2","TOX","CXCL13","ENTPD1","LAYN")
CTL_sig <- c("GZMB","GZMH","GZMA", "GNLY","NKG7","PRF1","KLRD1","KLRG1","FGFBP2","CCL4","CCL5")
Tem_sig <- c("GZMK","CCL5","IFNG","TNF","IL2","CXCR3","KLRB1","RORA")
Temra_sig <- c("FGFBP2","NKG7","GZMB","GNLY","PRF1","KLRD1","KLRG1","CX3CR1")
Tpex_sig <- c("TCF7","LEF1","IL7R","CCR7","SLAMF6","CXCR5")
Tcm_sig <- c("IL7R","CCR7", "SELL","TCF7","LTB","BCL2","PASK")
Tnaive_sig <- c("CCR7","IL7R","LRRN3", "NPM1", "TCF7","LEF1","LTB","MAL","SELL","TRAC")
glyco_sig <- c("HK2","GPI","PFKP","ALDOA","TPI1","GAPDH","PGK1","PGAM1","ENO1","PKM","LDHA")
OXPHOS_sig <-c("COX4I1","COX5A","CYCS","ATP5MC1","UQCRC1","UQCRH","SDHA","SDHB")
Cytokine1_sig <- c("IL2","IFNG","TNF","CSF2","TBX21")
Cytokine2_sig <- c("IL4","IL5","IL6","IL10","IL13","GATA3")

CD8_signatures <- list(
  Naive          = c("CCR7","SELL","TCF7","LEF1","IL7R","KLF2","CD27","CD28","SATB1","BACH2"),
  Tcm            = c("CCR7","IL7R","TCF7","CD27","CD28","SELL","BCL2","EOMES","CXCR3","IL2RB"),
  Tem            = c("GZMK","GZMA","CCL5","CXCR3","ITGA1","CD44","EOMES","TBX21","PRF1","NKG7","CX3CR1"),
  Temra          = c("GZMB","GZMH","PRF1","GNLY","CX3CR1","KLRG1","TBX21","B3GAT1","FGFBP2","ZEB2","KLRD1"),
  Effector       = c("GZMB","PRF1","IFNG","TNF","TBX21","KLRG1","CX3CR1","GNLY","NKG7","CCL4"),
  Tpex           = c("TCF7","SLAMF6","CXCR5","PDCD1","TOX","EOMES","ID3","BCL6","SELL","CD28","TIGIT"),
  Tex            = c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","TOX","TOX2","NR4A1","ENTPD1","CXCL13","CD38","LAYN"),
  Ttex           = c("TOX","PDCD1","EOMES","TBX21","GZMK","GZMB","LAG3","TIGIT","TCF7","CXCR3","CD28"),
  Trm            = c("ITGAE","ITGA1","CXCR6","CD69","ZNF683","ZNF683","BHLHE40","RORA","RUNX3","GZMK"),
  NKlike         = c("NCAM1","TYROBP","FCER1G","KIR2DL3","KLRC1","KLRC2","NCR1","GNLY","B3GAT1","FCGR3A"),
  EarlyActivated = c("CD69","CD38","HLA-DRA","ICOS","NR4A1","FOS","JUN","IRF4","BATF","TNFRSF9"),
  Proliferating  = c("MKI67","TOP2A","PCNA","CDK1","CCNB1","CCNA2","MCM2","MCM4","STMN1","UBE2C")
)

CD4_signatures <- list(
  Naive     = c("CCR7","SELL","TCF7","LEF1","IL7R","KLF2","CD27","CD28","SATB1","BACH2","CD40LG"),
  Tcm       = c("CCR7","IL7R","TCF7","CD27","CD28","SELL","BCL2","IL2RB","CXCR3","CD40LG"),
  Tem       = c("CXCR3","CCR5","ITGA1","CD44","IFNG","TNF","GZMK","NKG7","CCL5","EOMES","TBX21"),
  Th1       = c("TBX21","IFNG","TNF","CXCR3","CCR5","STAT1","STAT4","GZMB","PRF1","CCL3","CCL4"),
  Th2       = c("GATA3","IL4","IL5","IL13","CCR4","CCR3","STAT6","IRF4","PTGDR2","IL1RL1"),
  Th17      = c("RORC","RORA","IL17A","IL17F","IL22","CCR6","STAT3","IRF4","BATF","CXCR6"),
  Th22      = c("IL22","AHR","CCR10","CCR4","CCR6","FPR1","S100A8","S100A9"),
  Tfh       = c("CXCR5","BCL6","PDCD1","ICOS","IL21","CXCL13","MAF","SH2D1A","TIGIT","TOX2"),
  Treg      = c("FOXP3","IL2RA","CTLA4","IKZF2","RTKN2","FCRL3","TIGIT","LAYN","ENTPD1","PLCL1"),
  eTreg     = c("FOXP3","IL2RA","TNFRSF9","TNFRSF18","CCR8","LAYN","ENTPD1","TIGIT","BATF","CCR4","MAGEH1"),
  Trm       = c("CD69","ITGAE","CXCR6","ITGA1","BHLHE40","ZNF683","CXCL13","PDCD1","GZMK"),
  CD4_CTL   = c("GZMB","GZMH","PRF1","GNLY","NKG7","CX3CR1","EOMES","TBX21","SLAMF7","CRTAM","FGFBP2"),
  Tex       = c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","TOX","TOX2","NR4A1","CXCL13","ENTPD1","LAYN"),
  Tr1       = c("IL10","IFNG","EOMES","GZMB","PRF1","LAG3","CTLA4","LILRB4","CCL5","CXCR3","PRDM1"),
  Tscm      = c("TCF7","LEF1","CCR7","SELL","IL7R","CD27","CD28","BACH2","ID3","KLF2","BCL2"),
  Prolif    = c("MKI67","TOP2A","PCNA","CDK1","CCNB1","CCNA2","MCM2","MCM4","STMN1","UBE2C")
)

CDR3_stat <- function(seu,combine.CR,celltype) {
  #select samples in this project
  Idents(seu) <- "orig.ident"
  #select CART cells for downstream
  seu <- subset(seu, subset = CAR_positive=="TRUE")
  cells_keep <- colnames(seu) 
  combined.CR.filtered <- lapply(combined.CR, function(df) {df[df$barcode %in% cells_keep, ]})
  cdr3_stats <- map_df(names(combined.CR.filtered),
                       ~ combined.CR.filtered[[.x]] %>% mutate(TCR=paste0(TCR1,"_",TCR2), cdr3_aa=paste0(cdr3_aa1,"_",cdr3_aa2),cdr3_nt=paste0(cdr3_nt1,"_",cdr3_nt2)) %>%
                         summarise(sample = .x,cdr3a_NAs     = sum(is.na(cdr3_aa1)),cdr3a_nonNA  = sum(!is.na(cdr3_aa1)),
                                   cdr3b_NA     = sum(is.na(cdr3_aa2)),cdr3b_nonNA  = sum(!is.na(cdr3_aa2)),
                                   cdr3_aa1_n = n_distinct(cdr3_aa1, na.rm = TRUE),cdr3_aa2_n = n_distinct(cdr3_aa2, na.rm = TRUE),
                                   cdr3_nt1_n = n_distinct(cdr3_nt1, na.rm = TRUE),cdr3_nt2_n = n_distinct(cdr3_nt2, na.rm = TRUE),
                                   TCR1_n = n_distinct(TCR1, na.rm = TRUE),TCR2_n = n_distinct(TCR2, na.rm = TRUE),
                                   cdr3_aa_n = n_distinct(cdr3_aa, na.rm = TRUE),cdr3_nt_n = n_distinct(cdr3_nt, na.rm = TRUE),TCR_n = n_distinct(TCR, na.rm = TRUE)))
  write.table(cdr3_stats,paste0(celltype,"number of clones in each sample.csv"),sep=",",row.names = F,quote=F)
  return(combined.CR.filtered)
}

find_dims <- function(combined,celltype) {
  Tcells <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
  Tcells <- FindVariableFeatures(Tcells, selection.method = "vst")
  VariableFeatures(Tcells) <- grep("TRAV|TRBV|TRDV|TRGV|TRGC|TRAJ|TRAC|TRDC|ENSG|LINC|-AS1|^MT-|^RP[SL]",VariableFeatures(Tcells),value = TRUE, invert = TRUE)
  Tcells <- ScaleData(Tcells, vars.to.regress = c("nCount_RNA","percent.mt"))
  Tcells <- RunPCA(Tcells)
  ElbowPlot(object = Tcells,ndims=50)
  ggsave(paste0(celltype," Elbowplot of 50 dims.png"), dpi=300, width=6, height=5)
  pct_var <- Tcells[["pca"]]@stdev^2 / sum(Tcells[["pca"]]@stdev^2) * 100
  cumvar <- cumsum(pct_var)
  # Plot cumulative variance
  plot(cumvar, xlab = "PC", ylab = "Cumulative variance (%)", type = "b")
  abline(h = 80, col = "red", lty = 2)  # common 80% threshold
  # Find PC where you cross 80% explained
  pc80 <- which(cumvar > 80)[1]
  message(celltype, ": PC explaining >80% cumvar = ", pc80)
  return(Tcells)
}

merge_seurat <- function(Tcells,dims,celltype) {
  Tcells <- FindNeighbors(Tcells, dims = dims, reduction = "pca")
  Tcells <- FindClusters(Tcells, resolution = 0.3, cluster.name = "unintegrated_clusters")
  Tcells <- RunUMAP(Tcells,reduction = "pca", dims = dims,reduction.name = "umap.unintegrated")
  DimPlot(Tcells, reduction = "umap.unintegrated", group.by = c("group", "seurat_clusters"))
  ggsave(paste0(celltype," merged all samples UMAP group.png"), dpi=300, width=10, height=5)
  #plot to see cluster distribution amont group and samples
  DimPlot(Tcells, reduction = "umap.unintegrated",pt.size = 1, ncol = 2, label = T, split.by = "group")
  ggsave(paste0(celltype," Dimplot split by group after merge.png"), dpi=300, width=9, height=4)
  DimPlot(Tcells, reduction = "umap.unintegrated",pt.size = 1, ncol = 2, label = F, split.by = "orig.ident")
  ggsave(paste0(celltype," Dimplot split by sample after merge.png"), dpi=300, width=7, height=10)
  return(Tcells)
}

plot_sample_portion <- function(pt, sample, name, reduction,cluster_colors) {
  table <- pt %>% filter(grepl(sample,Var3))
  ggplot(table, aes(x = Var3, y = Freq, fill = Var1)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.5) +
    xlab("Sample") +
    ylab("Proportion") +
    scale_fill_manual(values = cluster_colors) +
    theme(legend.title = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  ggsave(paste("percentage of clusters in", sample, "after",name, reduction, "integration plan B.png"), dpi=300,width=5, height=5)
}

plot_group_portion <- function(pt, group, name, reduction,cluster_colors) {
  table <- pt %>% filter(grepl(group,Var2))
  ggplot(table, aes(x = Var2, y = Freq, fill = Var1)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.5) +
    xlab("group") + ylab("Proportion") +
    scale_fill_manual(values = cluster_colors) +
    theme(legend.title = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(paste("percentage of clusters in",group,"group after",name, reduction, "integration plan B.png"), dpi=300,width=4, height=5)
}

CARpos_barplot <- function(feature, group, name, reduction, feature_name, tag,
                           width = 8, height = 5, dpi = 300,
                           obj = NULL, variable = NULL) {
  
  tab <- as.data.frame(table(feature, group))
  colnames(tab) <- c("Var1", "Var2", "Freq")
  tab <- tab %>%
    group_by(Var2) %>%
    mutate(Prop = Freq / sum(Freq)) %>%
    ungroup()
  # Use UMAP colors if obj and variable are provided
  if (!is.null(obj) && !is.null(variable)) {
    p_umap <- DimPlot(obj, reduction = reduction, group.by = variable)
    umap_colors <- unique(ggplot_build(p_umap)$data[[1]]$colour)
    cluster_ids <- unique(p_umap[[1]][["data"]][[variable]])
    cluster_colors <- setNames(umap_colors, cluster_ids)
    cluster_colors <- cluster_colors[sort(names(cluster_colors))]
    my_cols <- cluster_colors
  } else {
    my_cols <- setNames(viridis(length(unique(tab$Var1))), sort(unique(tab$Var1)))
  }
  p <- ggplot(tab, aes(x = Var2, y = Freq, fill = Var1)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.5) +
    geom_text(aes(label = ifelse(Prop > 0.01, Freq, "")),
              position = position_fill(vjust = 0.5),
              size = 2.5) +
    xlab(tag) +
    ylab("Proportion") +
    scale_fill_manual(values = my_cols) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
  
  filename <- paste0("percentage of ", feature_name, " in each ", tag,
                     " after ", name, reduction, " integration.png")
  ggsave(filename, plot = p, dpi = dpi, width = width, height = height)
  return(p)
}

sample_heatmap <- function(seu_object, group,ident, width, height) {
  Idents(seu_object) <- ident
  cluster.markers <- FindAllMarkers(object = seu_object,only.pos = F, min.pct = 0, logfc.threshold = 0,return.thresh = 1)
  write.table(cluster.markers,paste("marker genes list of each",group, "sample.csv"),sep=",",quote=F,row.names=F)
  top30 <- cluster.markers %>% filter(pct.1>0.3, p_val_adj<0.05, !grepl("TRAV|TRBV|TRDV|TRGV|TRGC|TRAJ|TRAC|TRDC|ENSG|LINC|-AS1",gene)) 
  scaled_genes <- rownames(LayerData(seu_object, assay="RNA", layer="scale.data"))
  top30 <- top30 %>% filter(gene %in% scaled_genes) %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
  cells_use <- WhichCells(seu_object, downsample = 1000)
  DoHeatmap(seu_object, cells = cells_use, features = top30$gene,group.by=ident,disp.max = 4,raster = TRUE,slot = "scale.data") + scale_fill_gradientn(colors = c("#89cff0", "white", "red"))
  ggsave(paste("Heatmap of top 30 variable genes in each", group, "sample.png"), dpi=300,width=9, height=height)
  # ── heatmap 2: pseudobulk heatmap ───────────────────────────────────────
  genes_pb <- unique(top30$gene)
  # aggregate raw counts per sample (pseudobulk)
  samples    <- unique(seu_object@meta.data[[ident]])
  pb_matrix  <- sapply(samples, function(s) {
    cells_s <- colnames(seu_object)[seu_object@meta.data[[ident]] == s]
    # sum raw counts across cells for each gene
    Matrix::rowSums(seu_object[["RNA"]]$counts[genes_pb, cells_s, drop=FALSE])
  })
  
  # normalize: log CPM per sample
  pb_cpm <- log1p(sweep(pb_matrix, 2, colSums(pb_matrix), "/") * 1e6)
  
  # z-score across samples per gene for visualization
  pb_scaled <- t(scale(t(pb_cpm)))
  pb_scaled  <- pmin(pmax(pb_scaled, -4), 4)  # clip at ±4 like disp.max
  
  # order genes by cluster (same order as top30)
  gene_order <- top30 %>% ungroup() %>%
    arrange(cluster, desc(avg_log2FC)) %>%
    pull(gene) %>% unique()
  gene_order <- gene_order[gene_order %in% rownames(pb_scaled)]
  pb_scaled  <- pb_scaled[gene_order, , drop=FALSE]
  
  pheatmap::pheatmap(pb_scaled,
                     color            = colorRampPalette(c("#89cff0","white","red"))(100),
                     breaks           = seq(-4, 4, length.out=101),
                     cluster_rows     = FALSE,   # keep gene order by cluster
                     cluster_cols     = FALSE,   # keep sample order
                     show_rownames    = TRUE,
                     show_colnames    = TRUE,
                     fontsize_row     = 6,
                     fontsize_col     = 10,
                     main             = paste("Pseudobulk heatmap —", group),
                     filename         = paste("Pseudobulk heatmap of top 30 genes in each", group, "sample.png"),
                     width  = width,
                     height = height,
                     dpi    = 300)
  return(cluster.markers)
}

DE_analyse <- function(seu, name, substring, cells, group, contrast) {
  if (cells %in% c("CD4","CD8")) {
    sub <- subset(seu, subset=Cell_type==cells)
  }  else if (cells =="") { sub = seu } else {
    sub <- subset(seu, subset=annotation == cells)}
  pb <- make_pseudobulk(sub, group_var = "orig.ident") 
  pb <- pb[-c(1:3),grepl(substring,colnames(pb))]
  meta <- sub@meta.data %>% select(orig.ident, disease) %>% distinct()
  meta <- meta[match(colnames(pb), meta$orig.ident), ] 
  rownames(meta) <- meta$orig.ident
  dds <- DESeqDataSetFromMatrix(countData = pb,colData = meta, design = ~ disease )
  dds <- DESeq(dds)
  res <- results(dds, contrast = contrast)
  res <- as.data.frame(res[order(res$padj),])
  write.table(res,paste(name,"DEseq2 DEGs with pseudo bulk of",group, cells, "cells plan B.csv"),sep=",",quote=F)
  vsd <- vst(dds, blind = F)
  write.table(assay(vsd),paste(name,"vsd normalized counts of",group, cells,"cells plan B.csv"), sep=",",quote=F)
  res$diffexpressed <- "NO"
  res$diffexpressed[res$log2FoldChange > 0.5 & res$padj < 0.05] <- "UP"
  res$diffexpressed[res$log2FoldChange < -0.5 & res$padj < 0.05] <- "DOWN"
  res$delabel <- NA
  res$delabel[res$diffexpressed != "NO"] <- rownames(res[res$diffexpressed != "NO",])
  res <- res %>% filter(!grepl("TRAV|TRBV|TRDV|TRGV|TRGC|TRAJ|TRAC|TRDC|ENSG|LINC|-AS1|orf|HLA|MIR|-DT",delabel))
  ggplot(data=res, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
    geom_point(size=0.01) + 
    theme(legend.position = "none") +
    ggtitle(paste("ICAN vs noICANS DEGs", group, cells, name)) +
    geom_text_repel(aes(label=delabel),size=3,min.segment.size = 0.3, seed=40, segment.alpha = 0.5, segment.color = "grey50", force = 10, min.segment.length = 0.1) +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-0.5, 0.5), col="red", size=0.3, alpha=0.5) +
    geom_hline(yintercept=-log10(0.05), col="red", size=0.3, alpha=0.5)
  ggsave(paste(name, "volcano plot of DEGs with pseudo bulk of", group, cells, "T cells plan B.png"), width=6, height=5, dpi=300)
}

p_to_stars <- function(p) {if (is.na(p)) return("")
  if (p < 1e-4) "****"
  else if (p < 1e-3) "***"
  else if (p < 1e-2) "**"
  else if (p < 0.05) "*"
  else ""}

barplot_group <- function(obj, reduction, name, group, variable, subgroup, width, height) {
  p <- DimPlot(obj, reduction = reduction, group.by = variable)
  umap_colors <- unique(ggplot_build(p)$data[[1]]$colour)
  cluster_ids <- unique(p[[1]][["data"]][[variable]])
  cluster_colors <- setNames(umap_colors, cluster_ids)
  cluster_colors <- cluster_colors[sort(names(cluster_colors))]
  
  pt <- as.data.frame(table(obj@meta.data[[variable]], obj@meta.data[[subgroup]]))
  pt$Var1 <- as.character(pt$Var1)
  
  df <- pt %>%filter(grepl(group, Var2)) %>% group_by(Var2) %>% mutate(pct = 100 * Freq / sum(Freq), lab = sprintf("%.1f%%", pct)) %>% ungroup()
  
  # ensure exactly two bars
  df$Var2 <- droplevels(factor(df$Var2))
  x_levels <- levels(df$Var2)
  
  # composition test: cluster (Var1) x subgroup (Var2)
  ct <- xtabs(Freq ~ Var1 + Var2, data = df)
  stat_test <- chisq.test(ct)
  pval <- stat_test$p.value
  stars <- p_to_stars(pval)
  x_pos <- seq_along(x_levels)
  # plot
  g <- ggplot(df, aes(x = Var2, y = Freq, fill = Var1)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.5) +
    geom_text(aes(label = ifelse(pct >= 2, lab, "")),position = position_fill(vjust = 0.5), size = 4) +
    xlab("group") + ylab("Proportion") +
    scale_fill_manual(values = cluster_colors) +
    theme(legend.title = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    coord_cartesian(ylim = c(0, 1.12)) + 
    labs(caption = if (!is.na(pval))
      sprintf("Chi-square test: p = %.3g", pval)
      else NULL)
  
  # add bridge only if significant
  if (stars != "") {
    g <- g + geom_segment(x = x_pos[1], xend = x_pos[1],y = 1.02, yend = 1.05,linewidth = 0.8) + 
      geom_segment(x = x_pos[2], xend = x_pos[2],y = 1.02, yend = 1.05,linewidth = 0.8) + 
      geom_segment(x = x_pos[1], xend = x_pos[2],y = 1.05, yend = 1.05,linewidth = 0.8) + 
      annotate("text", x = mean(x_pos), y = 1.08,label = stars,size = 4)}
  outfile <- paste(name, "percentage of", variable, "in", group, "after harmony integration.pdf")
  ggsave(outfile, plot = g, dpi = 300, width = width, height = height)
}

find_resolution <- function(Tcells,dims,method,reduction,celltype) {
  integrated <- IntegrateLayers(object = Tcells, method = method, orig.reduction = "pca", new.reduction = reduction, verbose = FALSE)
  integrated <- FindNeighbors(integrated, reduction = reduction, dims = dims, k.param=50)
  res <- seq(0.1, 1.0, by = 0.1)
  integrated <- FindClusters(integrated, resolution = res)
  clustree(integrated, prefix = "RNA_snn_res.")
  ggsave(paste(celltype,reduction,"cluster tree plot.png"), dpi=300, width=10, height=12)
  return(integrated)
}

integration <- function(seu, res,reduction, reduction_name, clustername, celltype, dims, tcr, CAR_features) {
  
  #integrated <- FindNeighbors(seu, reduction = reduction, dims = dims, k.param=50)
  
  # Then cluster using the graph just built
  integrated <- FindClusters(seu, resolution = res)
  
  # UMAP from the correct reduction
  integrated <- RunUMAP(integrated, reduction = reduction, dims = dims, reduction.name = reduction_name,  n.neighbors = 50L,min.dist = 0.5, seed.use = 42)
  
  integrated <- JoinLayers(integrated)
  nrows <- ceiling(length(unique(integrated$orig.ident)) / 3)
  var_genes <- Trex::quietTCRgenes(VariableFeatures(integrated,assay = "RNA"))
  var_genes <- var_genes[!grepl("TRAV|TRBV|TRDV|TRGV|TRGC|TRAJ|TRAC|TRDC|ENSG|LINC|-AS1",var_genes)]
  VariableFeatures(integrated) <- grep("TRAV|TRBV|TRDV|TRGV|TRGC|TRAJ|TRAC|TRDC|ENSG|LINC|-AS1|^MT-|^RP[SL]",VariableFeatures(integrated),value = TRUE, invert = TRUE)
  plan("sequential")
  cluster.markers <- FindAllMarkers(object = integrated, features = var_genes, only.pos = T, min.pct = 0.2, logfc.threshold = 0.1)
  write.table(cluster.markers,paste0(celltype," marker genes list of ", reduction," integration ", res, ".csv"),sep=",",quote=F,row.names=F)
  #check marker genes of each cluster
  top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  cells_use <- WhichCells(integrated, downsample = 300)
  DoHeatmap(integrated, cells = cells_use, features = top10$gene,group.by=clustername,disp.max = 4,raster = TRUE,slot = "scale.data") + scale_fill_gradientn(colors = c("#89cff0", "white", "red"))
  ggsave(paste0(celltype,reduction, " Heatmap of top 10 variable genes in each Harmony integrated cluster resolution ", res, ".png"), dpi=300,width=12, height=18)
  # Save original UMAP embeddings
  combined <- combineExpression(tcr, integrated, cloneCall="aa", proportion=F,  cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
  combined <- combineExpression(tcr, combined, cloneCall="aa", group.by="sample", proportion=T)
  DimPlot(combined, reduction = reduction_name,group.by = "cloneSize",split.by = "orig.ident",ncol=2)
  ggsave(paste0(celltype,reduction," integrated UMAP TCR.png"), dpi=300, width=7, height=9)
  FeaturePlot(combined, features = CAR_features, reduction = reduction_name, ncol=3)
  ggsave(paste(celltype,"CAR features after",reduction, "integration.png"),dpi=300,width=3*3,height=3*nrows)
  DotPlot(combined, features = CAR_features) + RotatedAxis() + theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(paste(celltype,"dot plot after",reduction, "integration.png"), dpi=300,width=6, height=9,bg="white")
  combined = CellCycleScoring(object = combined, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes)
  DimPlot(combined,group.by='Phase',reduction=reduction_name)
  ggsave(paste(celltype,reduction,"CellCycle_FeaturePlot_all_samples.png"),dpi=300,width=5,height=5)
  DimPlot(combined, reduction = reduction_name, group.by = c("orig.ident", clustername,"Phase","patient","group"))
  ggsave(paste(celltype,reduction,res,"integration cluster umap.png"), dpi=300, width=10, height=6)
  #CARpos_barplot(combined$cloneSize,combined$orig.ident, celltype,"harmony", "Clonesize","sample")
  return(combined)
}

projectil_analysis <- function(seu,ref,celltype) {
  seu <- ProjecTILs.classifier(query = seu,ref = ref, filter.cells=F)
  DimPlot(seu, group.by = "functional.cluster",reduction = "umap.harmony")
  ggsave(paste(celltype,"UMAP of ProjectTILs annotaiton with reference.png"),dpi=300,width=10, height=9)
  return(seu)
}

top_clonal_compare <- function(combined,cr_list, patient, celltype) {
  sampP <- paste0(patient, "-P")
  sampB <- paste0(patient, "-B")
  Idents(combined) <- "patient"
  sub <- subset(combined, idents = patient)
  meta <- sub@meta.data %>% filter(!is.na(CTaa))
  countP <- meta %>% filter(orig.ident == sampP) %>% dplyr::count(CTaa, name = "nP") %>% arrange(desc(nP)) %>% slice_head(n = 10)
  countB <- meta %>% filter(orig.ident == sampB) %>% dplyr::count(CTaa, name = "nB") %>% arrange(desc(nB)) %>% slice_head(n = 10)
  clonalCompare(cr_list, top.clones=NULL, clones = unique(c(countP$CTaa,countB$CTaa)), order.by = c(sampP, sampB), cloneCall = "aa",graph = "alluvial")
  ggsave(paste(celltype,"Top10 CDR3 aa clonal comparison between",patient,"Product and Blood.png"), dpi=300, width=6, height=6)
  countP_seu <- meta %>% filter(orig.ident == sampP) %>% dplyr::count(CTaa, name = "nP") %>% arrange(desc(nP)) %>% slice_head(n = 3)
  countB_seu <- meta %>% filter(orig.ident == sampB) %>% dplyr::count(CTaa, name = "nB") %>% arrange(desc(nB)) %>% slice_head(n = 3)
  top3_seurat <- unique(c(countP_seu$CTaa, countB_seu$CTaa))
  sub$orig.ident <- factor(sub$orig.ident, levels = c(sampP, sampB))
  sub <- highlightClones(sub, cloneCall = "aa", sequence = top3_seurat)
  Seurat::DimPlot(sub, group.by = "highlight",reduction="umap.harmony",split.by = "orig.ident") +
    guides(color = guide_legend(
      ncol          = 2,          # number of columns (change to 2, 3 etc)
      byrow         = TRUE,
      override.aes  = list(size = 3)  # size of legend dots
    ))  +
    ggplot2::theme(
      plot.title      = element_blank(),
      legend.position = "bottom",
      legend.text     = element_text(size = 8),   # font size of labels
      legend.title    = element_text(size = 8),   # font size of "highlight" title
      legend.key.size = unit(0.4, "cm")           # size of legend key boxes
    )
  ggsave(paste(patient,celltype,"top 3 clones UMAP.png"),dpi=300,width=9,height=5)
}

shared_clonal_compare <- function(combined,cr_list, patient, celltype) {
  sampP <- paste0(patient, "-P")
  sampB <- paste0(patient, "-B")
  clones <- intersect(unique(cr_list[[sampP]][["CTaa"]]),unique(cr_list[[sampB]][["CTaa"]]))
  countP <- cr_list[[sampP]] %>% filter(CTaa %in% clones) %>% dplyr::count(CTaa, name = "nP")
  countB <- cr_list[[sampB]] %>% filter(CTaa %in% clones) %>% dplyr::count(CTaa, name = "nB")
  # Merge and rank by total cells (nP + nB)
  top10_tbl <- full_join(countP, countB, by = "CTaa") %>%
    mutate(nP = coalesce(nP, 0L), nB = coalesce(nB, 0L), nTotal = nP + nB) %>%
    arrange(desc(nTotal)) %>% slice_head(n = 10)
  clonalCompare(cr_list, clones = top10_tbl$CTaa, order.by = c(sampP, sampB), cloneCall = "aa",graph = "alluvial")
  ggsave(paste(celltype,"Top10 shared CDR3 aa clonal comparison between",patient,"Product and Blood.png"), dpi=300, width=6, height=6)
  top_tbl <- top10_tbl %>% slice_head(n = 5)
  Idents(combined) <- "patient"
  sub <- subset(combined, idents = patient)
  sub$orig.ident <- factor(sub$orig.ident, levels = c(sampP,sampB))
  sub <- highlightClones(sub,cloneCall = "aa",sequence  = top_tbl$CTaa)
  Seurat::DimPlot(sub, group.by = "highlight",reduction="umap.harmony",split.by = "orig.ident") + guides(color=guide_legend(nrow=3,byrow=TRUE)) +
    ggplot2::theme(plot.title = element_blank(), legend.position = "bottom")
  ggsave(paste(patient,celltype,"shared top 5 clones UMAP.png"),dpi=300,width=9,height=6)
}

signature_scores <- function(seu, signatures,celltype, cluster) {
  for (sig_name in names(signatures)) {
    seu <- AddModuleScore(seu,features = list(signatures[[sig_name]]),
                          name = paste0(sig_name, "_score") )}
  score_cols <- paste0(names(signatures), "_score1")
  
  FeaturePlot(seu,features  = score_cols,
              reduction = "umap.harmony",
              ncol      = 4) & theme(plot.title = element_text(size = 8))
  ggsave(paste(celltype,"module signature scores featureplot.pdf"),width=12,height=9)
  # Dot plot summary per cluster
  DotPlot(seu,features   = score_cols,
          group.by   = cluster) + RotatedAxis() +
    scale_color_gradientn(colors = c("blue", "white", "red"))
  ggsave(paste(celltype,"module signature scores dotplot.pdf"),width=7,height=5)
  return(seu)
}

add_annotation <- function(seu,new.cluster.ids,cluster.col) {
  annotation <- new.cluster.ids[as.character(seu@meta.data[[cluster.col]])]
  names(annotation) <- colnames(seu)
  seu <- AddMetaData(object = seu,metadata = annotation,col.name = 'annotation')
  return(seu)
}

Peudotime_analysis <- function(seu,celltype, root,genes_of_interest) {
  #monocle analysis
  cds <- as.cell_data_set(seu,reductions="umap.harmony")
  reducedDims(cds)[["UMAP"]] <- reducedDims(cds)[["UMAP.HARMONY"]]
  cds@clusters[["UMAP"]]$clusters <- seu$seurat_clusters
  # Now cluster_cells should work
  cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = 1e-3)
  cds <- learn_graph(cds, use_partition = TRUE)
  plot_cells(cds,
             color_cells_by = "cluster",
             label_groups_by_cluster = TRUE,
             label_leaves = TRUE,
             label_principal_points = TRUE)
  cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == root]))
  plot_cells(cds, color_cells_by = "pseudotime",
             group_cells_by = "cluster",
             label_cell_groups = FALSE,
             label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE,
             label_roots = FALSE,
             trajectory_graph_color = "grey60")
  ggsave(paste(celltype,"Monocle Trajectory on UMAP.pdf"), width=10, height=10)
  # Add pseudotime and cluster info back to Seurat metadata
  seu$pseudotime <- pseudotime(cds)
  
  # Basic pseudotime distribution per patient
  ggplot(seu@meta.data %>% filter(!is.infinite(pseudotime)),
         aes(x = pseudotime, fill = patient, color = patient)) +
    geom_density(alpha = 0.4) +
    facet_wrap(~ group, ncol = 1) +  # CR, CD19+relapse, CD19-relapse
    theme_classic() +
    labs(title = "Pseudotime Distribution by Patient Outcome",
         x = "Pseudotime", y = "Density") +
    theme(strip.text = element_text(face = "bold"))
  ggsave(paste(celltype,"Peudotime comparison between patients.pdf"), width=10, height=10)
  
  # Bin pseudotime into quartiles
  seu$pt_bin <- cut(seu$pseudotime,
                    breaks = quantile(seu$pseudotime,
                                      probs = seq(0, 1, 0.25),
                                      na.rm = TRUE),
                    labels = c("Early", "Early-mid", "Late-mid", "Late"),
                    include.lowest = TRUE)
  
  # Cluster composition per pseudotime bin per patient
  seu@meta.data %>%
    filter(!is.na(pt_bin), !is.infinite(pseudotime)) %>%
    group_by(group, patient, pt_bin, annotation) %>%
    summarise(n = n(), .groups = "drop") %>%
    # Key fix — normalize within each patient+group combination
    group_by(patient, group, pt_bin) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup() %>%
    # Set factor levels for correct ordering
    mutate(patient = factor(patient,
                       levels = c("UPN788", "UPN774", "UPN764")),
      group   = factor(group,
                       levels = c("Product", "Blood"))) %>%
    ggplot(aes(x = pt_bin, y = pct, fill = annotation)) +
    geom_bar(stat = "identity") +
    # 3 rows (patient) x 2 cols (Product | Blood)
    facet_grid(patient ~ group) +
    theme_classic() +
    labs(title = paste(celltype,"Cluster Composition Along Pseudotime"),
         x     = "Pseudotime Bin",
         y     = "% Cells",
         fill  = "Cluster") +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
          strip.text       = element_text(face = "bold", size = 10),
          strip.background = element_rect(fill = "grey90"),
          legend.position  = "bottom") +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE))
  ggsave(paste0(celltype,"_cluster_pseudotime_bins_6panels.pdf"), dpi = 300, width = 10, height = 12)
  # Where does the top clone sit in pseudotime?
  # This is unique to your CART clonotype data
  top_clone_per_patient <- seu@meta.data %>%
    filter(!is.na(CTaa)) %>%
    group_by(orig.ident, CTaa) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(orig.ident) %>%
    slice_max(n, n = 10) %>%
    select(orig.ident, top_clone = CTaa)
  
  seu$is_top_clone <- FALSE
  
  for (p in unique(seu$orig.ident)) {
    top <- top_clone_per_patient %>%
      filter(orig.ident == p) %>%
      pull(top_clone)
    
    cells <- rownames(seu@meta.data)[
      seu$orig.ident == p &
        !is.na(seu$CTaa) &
        seu$CTaa %in% top  # fixed: %in% not ==
    ]
    seu$is_top_clone[
      rownames(seu@meta.data) %in% cells
    ] <- TRUE
  }
  
  # Compare pseudotime of top clone vs other clones
  seu@meta.data %>% filter(!is.infinite(pseudotime)) %>%
    ggplot(aes(x = orig.ident, y = pseudotime,
               fill = is_top_clone)) +
    geom_violin(alpha = 0.7, position = position_dodge(0.8)) +
    scale_fill_manual(values = c("TRUE"  = "firebrick",
                                 "FALSE" = "grey70"),
                      labels = c("TRUE"  = "Top 10 Clone",
                                 "FALSE" = "Other Clones")) +
    theme_classic() +
    labs(title = "Top Clone Pseudotime Position by sample",
         x = NULL, y = "Pseudotime", fill = NULL) +
    theme(legend.position = "bottom")
  ggsave(paste(celltype,"TopClone_pseudotime_by_sample.pdf"), width = 7, height = 4)
  
  # Show how top clone cells are distributed across clusters
  # comparing product vs blood per patient
  alluvial_data <- seu@meta.data %>%
    filter(!is.na(CTaa), is_top_clone == TRUE) %>%
    mutate(compartment = ifelse(grepl("-P", orig.ident),
                                "Product", "Blood")) %>%
    group_by(patient, compartment, annotation) %>%
    summarise(n = n(), .groups = "drop")
  
  ggplot(alluvial_data,
         aes(x      = compartment,
             y      = n,
             alluvium = annotation,
             stratum  = annotation,
             fill     = annotation,
             label    = annotation)) +
    geom_stratum(alpha = 0.8) +
    geom_flow(alpha = 0.5) +
    geom_text(stat = "stratum", size = 2.5, fontface = "bold") +
    facet_wrap(~ patient, ncol = 3) +
    theme_classic() +
    labs(title = "Top Clone Distribution: Product → Blood",
         x     = NULL,
         y     = "Number of Cells") +
    theme(legend.position = "none",
          strip.text       = element_text(face = "bold"))
  
  ggsave(paste(celltype,"Top10_clone_cluster_flow.pdf"),width = 12, height = 6)
  
  prop_data <- seu@meta.data %>%
    filter(!is.na(CTaa)) %>%                    # TCR-positive cells only
    group_by(orig.ident, patient, annotation,   # annotation = your cluster names
             is_top_clone) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(orig.ident, annotation) %>%
    mutate(pct = n / sum(n) * 100) %>%
    filter(is_top_clone == TRUE)                # keep only top clone rows
  
  ggplot(prop_data,aes(x    = annotation,
                       y    = pct,
                       fill = orig.ident)) +
    geom_bar(stat     = "identity",
             position = position_dodge(0.8),
             width    = 0.7) +
    facet_wrap(~ patient, ncol = 3) +
    scale_fill_brewer(palette = "Set2") +
    theme_classic() +
    labs(title = "Top 10 Clone Proportion per Cluster per Sample",
         x     = "Cluster",
         y     = "% Top Clone Cells within Cluster",
         fill  = "Sample") +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1,
                                         size = 8),
          strip.text       = element_text(face = "bold"),
          legend.position  = "bottom")
  ggsave(paste(celltype,"Top10_clone_proportion_per_cluster.pdf"),width = 14, height = 6)
  # Compare pseudotime between product and blood per patient
  # Higher pseudotime in blood = more differentiation occurred in vivo
  
  seu@meta.data %>%
    filter(!is.infinite(pseudotime)) %>%
    ggplot(aes(x = group, y = pseudotime, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
    stat_compare_means(method = "wilcox.test",
                       label  = "p.signif") +
    facet_wrap(~ patient, ncol = 3) +
    scale_fill_manual(values = c("Product" = "steelblue",
                                 "Blood"   = "firebrick")) +
    theme_classic() +
    labs(title = "Pseudotime Shift: Product → Blood by Outcome",
         x = NULL, y = "Pseudotime") +
    theme(strip.text = element_text(face = "bold"))
  ggsave(paste(celltype,"Pseudotime_product_vs_blood_by_outcome.pdf"),width = 9, height = 5)
  
  # Extract expression + pseudotime
  expr_pt <- FetchData(seu,vars = c(genes_of_interest,
                                    "pseudotime", "group",
                                    "patient", "orig.ident")) %>%
    filter(!is.infinite(pseudotime))
  
  # Plot each gene along pseudotime per outcome group
  plot_list <- lapply(genes_of_interest, function(gene) {
    ggplot(expr_pt, aes(x    = pseudotime, y    = .data[[gene]], color = orig.ident)) +
      geom_smooth(method = "loess", se = TRUE, span = 0.75) +
      theme_classic() +
      labs(title = gene, x = "Pseudotime", y = "Expression") +
      theme(plot.title    = element_text(face = "bold", size = 10),
            legend.position = "none")
  })
  combined <- wrap_plots(plot_list, ncol = 4) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  print(combined)
  ggsave(paste0(celltype,"Gene_expression_pseudotime_by_outcome.pdf"), plot = combined, width = 9, height = 4)
  return(seu)
}

compare_top_clone <- function(cr_list,patient, top_number, name) {
  sampP <- paste0(patient, "-P")
  sampB <- paste0(patient, "-B")
  clonalCompare(cr_list, top.clones = top_number,samples = c(sampP, sampB), order.by = c(sampP, sampB), cloneCall="aa", graph = "alluvial",palette = "viridis") +
    theme(axis.text = element_text(face = "bold", color = "black",size = 14),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title = element_text(face = "bold", color = "black", size = 12),
          legend.text = element_text(color = "black", size = 14),
          legend.title = element_text(face = "bold", color = "black", size = 14)) 
  ggsave(paste(name,"Top10 CDR3 aa clonal comparison between",patient,"Product and Blood plan B.pdf"), dpi=300, width=10, height=6)
  clonalScatter(cr_list, cloneCall ="gene", x.axis = sampP, y.axis = sampB,dot.size = "total",graph = "proportion")
  ggsave(paste(name,"aa clonal scatter plot for",patient,"Product and Blood plan B.png"), dpi=300, width=6, height=6)
  vizGenes(combined.CR.filtered[c(sampP, sampB)],x.axis = "TRBV",y.axis = NULL,plot = "barplot") 
  ggsave(paste(name,"TRBV vizGene plot for",patient,"Product and Blood plan B.png"), dpi=300, width=10, height=6)
}

make_dimplot_with_pecent <- function(obj, reduction, name, split_variable, variable, width, height) {
  md <- obj@meta.data %>% rownames_to_column("cell") %>% mutate(split = as.character(.data[[split_variable]]),cluster = as.factor(.data[[variable]])) 
  meta_pct <-  md %>% group_by(split, cluster) %>% summarise(n = n(), .groups = "drop") %>% group_by(split) %>% mutate(percent = 100 * n / sum(n))
  umap_df <- as.data.frame(Embeddings(obj, reduction)) %>% rownames_to_column("cell") %>% inner_join(md, by = "cell")
  label_pos <- umap_df %>% group_by(split, cluster) %>% summarise(umap_1 = median(umap_1), umap_2 = median(umap_2),.groups = "drop") %>%
    left_join(meta_pct, by = c("split", "cluster")) %>% mutate(label = sprintf("%.1f%%", percent)) %>% filter(split != "NA")
  label_pos[[split_variable]] <- label_pos$split
  DimPlot(obj, reduction = reduction,pt.size = 1, ncol = 2, label = F, group.by = variable, split.by = split_variable) +
    geom_text(data = label_pos,aes(x = umap_1, y = umap_2, label = label),size = 3,color = "black", fontface = "bold")
  ggsave(paste(name,"harmony integration umap split by", split_variable, "with percentage of",variable,"plan B.svg"), dpi=300, width=width, height=height)
}

make_pseudobulk <- function(seu, group_var = "orig.ident") {
  counts <- GetAssayData(seu, layer = "counts")
  meta <- seu@meta.data
  # split cells by sample
  sample_cells <- split(colnames(seu), meta[[group_var]])
  # sum counts per sample (pseudobulk)
  pseudobulks <- lapply(sample_cells, function(cells) {
    Matrix::rowSums(counts[, cells, drop = FALSE])
  })
  pseudobulks <- do.call(cbind, pseudobulks)
  return(pseudobulks)
}

# Statistical test of pseudotime differences between outcomes
pseudotime_test <- function(seu,tissue,celltype) {
  max_pt <- max(seu@meta.data$pseudotime[
    !is.infinite(seu@meta.data$pseudotime)],
    na.rm = TRUE)
  seu@meta.data %>%
    filter(!is.infinite(pseudotime),group==tissue) %>%
    ggplot(aes(x = patient, y = pseudotime, fill = patient)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
    stat_compare_means(
      comparisons = list(
        c("UPN788", "UPN774"),
        c("UPN788", "UPN764"),
        c("UPN764", "UPN774")),
      method = "wilcox.test",
      label  = "p.signif",hide.ns = TRUE,label.y  = max_pt * 0.85, step.increase = 0.06,vjust = 0.01
    ) + facet_wrap(~ annotation, ncol = 4) +
    theme_classic() +
    labs(title = "Pseudotime Comparison by Patient",
         x = NULL, y = "Pseudotime") +
    theme(strip.text      = element_text(face = "bold"),
          legend.position = "bottom",
          axis.text.x     = element_text(angle = 45, hjust = 1))
  ggsave(paste(tissue,celltype,"pseudotime_stats_by_patient.pdf"),width = 8, height = 10)
}

# Heatmap
CART_metrix_heatmap <- function(summary_df,celltype) {
  mat <- summary_df %>% tibble::column_to_rownames("orig.ident") %>% as.matrix() %>% scale()
  pheatmap(mat,cluster_rows = FALSE,
           cluster_cols = TRUE,
           color        = colorRampPalette(c("blue","white","red"))(100),
           main         = paste(celltype,"CART Metrics by Patient"),
           fontsize      = 10,
           filename = paste0(celltype,"_CART_metrics_heatmap.pdf"),  # add this
           width    = 8,    # inches
           height   = 6)}

make_dimplot_with_percent_by_sample <- function(obj, name, reduction, sample_var = "orig.ident",cluster_var = "annotation", width = 7, height = 7,outdir = ".") {
  p_full <- DimPlot(obj, group.by = cluster_var, reduction = reduction)
  cluster_colors <- p_full$data %>% select(all_of(cluster_var)) %>% bind_cols(color = ggplot_build(p_full)$data[[1]]$colour) %>% distinct() %>% deframe() 
  samples <- obj@meta.data[[sample_var]] %>% as.character() %>% unique() 
  for (s in samples) {
    cells_use <- rownames(obj@meta.data)[obj@meta.data[[sample_var]] == s]
    obj_s <- subset(obj, cells = cells_use)
    md <- obj_s@meta.data %>% rownames_to_column("cell") %>% mutate(cluster = as.factor(.data[[cluster_var]]))
    # % of each cluster within this sample
    meta_pct <- md %>% group_by(cluster) %>% summarise(n = n(), .groups = "drop") %>% mutate(percent = 100 * n / sum(n),label = sprintf("%.1f%%", percent))
    umap_mat <- Embeddings(obj_s, reduction)
    # keep only first 2 dimensions, regardless of how many are present
    umap_df <- as.data.frame(umap_mat[, 1:2, drop = FALSE]) %>% rownames_to_column("cell") %>% setNames(c("cell", "umap_1", "umap_2")) %>% inner_join(md, by = "cell")
    # label position per cluster (median)
    label_pos <- umap_df %>% group_by(cluster) %>%summarise(umap_1 = median(umap_1),umap_2 = median(umap_2),.groups = "drop") %>% left_join(meta_pct, by = "cluster")
    if (reduction=="umap") {
      p <- DimPlot(obj_s, reduction = reduction,pt.size = 1, label = FALSE,group.by = cluster_var) + ggtitle(s) +
        scale_color_manual(values = cluster_colors, drop = TRUE) +
        ggrepel::geom_label_repel(
          data         = label_pos,
          aes(x        = umap_1,
              y        = umap_2,
              label    = meta_pct$label),
          size         = 3,
          fontface     = "bold",
          box.padding  = 0.3,
          label.padding = 0.2,
          label.size   = 0.3,          # border thickness
          fill         = "white",
          alpha        = 0.85,
          max.overlaps = 30,
          inherit.aes  = FALSE         # critical — don't inherit DimPlot aes
        ) + theme(plot.title = element_text(face = "bold"))
    } else {
      p <- DimPlot(obj_s, reduction = reduction,pt.size = 1, label = FALSE,group.by = cluster_var) + 
        geom_text(data = label_pos,aes(x = umap_1, y = umap_2, label = label),size = 5,color = "black", fontface = "bold") +
        theme(axis.text = element_text(face = "bold", color = "black",size = 14),
              axis.title = element_text(face = "bold", color = "black", size = 14),
              legend.text = element_text(face = "bold", color = "black", size = 14),
              legend.title = element_text(face = "bold", color = "black", size = 14)) 
    }
    outfile <- file.path(outdir,sprintf("%s_%s_%s_%s.svg",reduction, name, sample_var, gsub("[^A-Za-z0-9_.-]", "_", s))) 
    ggsave(outfile, plot = p, dpi = 300, width = width, height = height)
  }
}

sample_heatmap <- function(seu_object, group,celltype, width, height) {
  Idents(seu_object) <- celltype
  plan("sequential") 
  cluster.markers <- FindAllMarkers(object = seu_object,only.pos = F, min.pct = 0.1, logfc.threshold = 0,return.thresh = 1)
  write.table(cluster.markers,paste(celltype,"marker genes list of each",group, "sample.csv"),sep=",",quote=F,row.names=F)
  top30 <- cluster.markers %>% filter(pct.1>0.3, p_val_adj<0.05, !grepl("TRAV|TRBV|TRDV|TRGV|TRGC|TRAJ|TRAC|TRDC|ENSG|LINC|-AS1",gene)) 
  scaled_genes <- rownames(LayerData(seu_object, assay="RNA", layer="scale.data"))
  top30 <- top30 %>% filter(gene %in% scaled_genes) %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
  cells_use <- WhichCells(seu_object, downsample = 1000)
  DoHeatmap(seu_object, cells = cells_use, features = top30$gene,group.by=celltype,disp.max = 4,raster = TRUE,slot = "scale.data") + scale_fill_gradientn(colors = c("#89cff0", "white", "red"))
  ggsave(paste(celltype,"Heatmap of top 30 variable genes in each", group, "sample.png"), dpi=300,width=9, height=12)
  # ── heatmap 2: pseudobulk heatmap ───────────────────────────────────────
  genes_pb <- unique(top30$gene)
  # aggregate raw counts per sample (pseudobulk)
  samples    <- unique(seu_object@meta.data[[celltype]])
  pb_matrix  <- sapply(samples, function(s) {
    cells_s <- colnames(seu_object)[seu_object@meta.data[[celltype]] == s]
    # sum raw counts across cells for each gene
    Matrix::rowSums(seu_object[["RNA"]]$counts[genes_pb, cells_s, drop=FALSE])
  })
  
  # normalize: log CPM per sample
  pb_cpm <- log1p(sweep(pb_matrix, 2, colSums(pb_matrix), "/") * 1e6)
  
  # z-score across samples per gene for visualization
  pb_scaled <- t(scale(t(pb_cpm)))
  pb_scaled  <- pmin(pmax(pb_scaled, -4), 4)  # clip at ±4 like disp.max
  
  # order genes by cluster (same order as top30)
  gene_order <- top30 %>% ungroup() %>%
    arrange(cluster, desc(avg_log2FC)) %>%
    pull(gene) %>% unique()
  gene_order <- gene_order[gene_order %in% rownames(pb_scaled)]
  pb_scaled  <- pb_scaled[gene_order, , drop=FALSE]
  
  pheatmap::pheatmap(pb_scaled,
                     color            = colorRampPalette(c("#89cff0","white","red"))(100),
                     breaks           = seq(-4, 4, length.out=101),
                     cluster_rows     = FALSE,   # keep gene order by cluster
                     cluster_cols     = FALSE,   # keep sample order
                     show_rownames    = TRUE,
                     show_colnames    = TRUE,
                     fontsize_row     = 6,
                     fontsize_col     = 10,
                     main             = paste(celltype,"Pseudobulk heatmap —", group),
                     filename         = paste(celltype,"Pseudobulk heatmap of top 30 genes in each", group, "sample.pdf"),
                     width  = width,
                     height = height,
                     dpi    = 300)
  return(cluster.markers)
}

gsea_analysis <- function(markerlist, celltype, group, category, gs, qcutoff, width, height) {
  sig_GSEA <- list()
  clusters <- unique(markerlist$cluster) 
  for (i in seq_along(clusters)) {
    DEG <- markerlist %>% filter(cluster==clusters[i]) %>% select(gene,avg_log2FC)
    DEG_entrez <- merge(DEG,entrez, by.x="gene",by.y="Gene.name") %>% dplyr::select(NCBI.gene..formerly.Entrezgene..ID,avg_log2FC) %>% 
      group_by(NCBI.gene..formerly.Entrezgene..ID) %>% slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>% arrange(desc(avg_log2FC)) %>% ungroup()
    geneList <- setNames(DEG_entrez$avg_log2FC,  as.integer(DEG_entrez$NCBI.gene..formerly.Entrezgene..ID))
    gse <- genGSEA(genelist = geneList, geneset = gs,p_cutoff = 1,q_cutoff = 1,min_gset_size =10,max_gset_size = 2000)
    sig_GSEA[[i]] <- gse[["gsea_df"]]
    sig_GSEA[[i]]$sample <- clusters[i]
  }
  GSEA_results <- do.call(rbind,sig_GSEA)
  signif_results <- GSEA_results %>% filter(qvalue<qcutoff)
  GSEA_results <- GSEA_results %>% filter(ID %in% unique(signif_results$ID))
  ggplot(GSEA_results, aes(x= sample, y=reorder(ID,NES), size=-log10(qvalue), color=NES, group=sample)) + 
    geom_point(alpha = 0.8) + theme_classic() +  geom_point(data = subset(GSEA_results, qvalue < 0.05),aes(size = -log10(qvalue)),
                                                            shape = 21, fill = NA, color = "black", stroke = 0.8) + labs(y = "Pathways", x = "Sample") +
    scale_color_gradient(high = "red2",  low = "mediumblue", space = "Lab") +
    scale_size(range = c(1, 7)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggsave(paste(celltype,group,"dot plot of GSEA",category,"pathways.png"),dpi=300,width=width, height = height)
}
