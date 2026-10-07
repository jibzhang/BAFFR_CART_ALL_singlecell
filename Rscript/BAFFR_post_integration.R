library(pheatmap)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(svglite)
library(tidyverse)
library(genekitr)
library(DESeq2)
library(ggrepel)
workdir="/net/isi-dcnl/ifs/user_data/lkwak/Data/Jibin/project/BAFFR_singlecell/"
setwd(paste0(workdir,"T_cell_integration"))
harmony_inte <- readRDS("CART harmony integrated CART cells planB.rds")
Idents(harmony_inte) <- "orig.ident"
product <- subset(harmony_inte, idents = c("UPN787-P","UPN775-P","UPN685-P","UPN749-P","UPN746-P"))
blood <- subset(harmony_inte, idents = c("UPN787-B","UPN775-B","UPN685-B","UPN749-B","UPN746-B"))

scores=c("CTL_score","Tem_score","Temra_score","Tex_score","Prolif_score","Tcm_score","Tn_score","Tpex_score")
meta_scores=c("Glyco_score","OXPHOS_score")
cytokine_scores=c("Cytokine1_socre","Cytokine2_score")
make_dotplot <- function(seu, desired_order, scores, index, group, width, height) {
  seu$orig.ident <- factor(seu$orig.ident, levels = desired_order)
  DotPlot(seu, features = scores,group = "orig.ident") + RotatedAxis() 
  ggsave(paste0("dot plot of ", index, " index scores for ", group,".svg"), dpi=300,width=width, height=height,bg="white")
  DotPlot(seu, features = scores,group = "disease") + RotatedAxis() +  theme(axis.text = element_text(size = 14,face = "bold", color = "black"),
                                                                          axis.title = element_text(face = "bold", color = "black", size = 14),
                                                                           legend.text = element_text(face = "bold", color = "black", size = 14),
                                                                           legend.title = element_text(face = "bold", color = "black", size = 14)) 
  ggsave(paste0("dot plot of ", index, " index scores for pooled ", group,".svg"), dpi=300,width=width, height=height,bg="white")}

desired_order <- c("UPN685-P", "UPN775-P", "UPN787-P", "UPN749-P", "UPN746-P")
make_dotplot(product, desired_order, meta_scores, "metabolism", "product", 4,4)
make_dotplot(product, desired_order, cytokine_scores, "cytokine", "product", 4,4)
make_dotplot(product, desired_order, scores, "cell type", "product", 6,4)
desired_order <- c("UPN685-B", "UPN775-B", "UPN787-B", "UPN749-B", "UPN746-B")
make_dotplot(blood, desired_order,meta_scores, "metabolism", "blood", 5,5)
make_dotplot(blood, desired_order,cytokine_scores, "cytokine", "blood", 5,5)
make_dotplot(blood, desired_order,scores, "cell type", "blood", 6,5)

genes <- c("CD4","IL2RA","GATA3","ICOS","MX1","IFIT3","SELL","JUN","BACH2","GZMB","CD8A","CD8B","LEF1","PRF1","FOS","FOSB","CXCR3","BATF3","CCR7","IL7R","TCF7","FASLG",
           "ITGA1", "MCM2","MCM4","MCM7","PCNA","MYB","MCM10","TYMS","STMN1","MKI67","TOP2A","CDK1","HMGB2","HMGB3","HMGN2","BIRC5","EZH2","CCNB2","CCNB1",
           "PLK1","TNF","TOX2", "IFI6","GZMA","CXCR4","CD27","RUNX3","CCL5","KLF2","CD44","GZMK","NKG7","HAVCR2","LAG3","SLAMF6","ZEB2","CCL4","FAS","GNLY","GZMH","TIGIT",
           "FGFBP2","CTLA4","IFNG","KLRG1","CX3CR1","CD83","ADGRG1","TBX21","KLRD1","FOXO1","BCL11B","TOX","RORA","EOMES","FOXP3")
dp <- DotPlot(harmony_inte, features = genes, group.by = "annotation")
# Step 2: find which cluster each gene is highest in
gene_order <- dp$data %>%
  group_by(features.plot) %>%
  slice_max(avg.exp.scaled, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(id = factor(id, levels = c(
    "CD4 Tn","Tem-like","Temra","Tnk",
    "Tpex","Tprol-M","Tprol-S","Treg"))) %>%
  arrange(id, desc(avg.exp.scaled)) %>%
  pull(features.plot) %>%
  as.character()
# Step 3: replot with ordered genes
DotPlot(harmony_inte, features = gene_order, group.by = "annotation") + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text = element_text(face = "bold", color = "black",size = 14),
        axis.title = element_text(face = "bold", color = "black", size = 14),
        legend.text = element_text(face = "bold", color = "black", size = 14),
        legend.title = element_text(face = "bold", color = "black", size = 14)) +
  scale_color_gradient2(low = "blue", mid = "lightgrey", high = "red") 
ggsave("signature genes dotplot for each cluster.pdf",dpi = 300, width = 5, height = 15, bg = "white")

gene_sets <- list( CTL = c("GZMB","GZMH","GZMA","GZMK","IFNG","PRF1","CCL3","CCL4","CCL5"),
  Temra = c("NKG7","GNLY","KLRD1","KLRG1","CX3CR1"),  
  Tem = c("TNF","IL2","CXCR3","KLRB1","RORA"),
  Tcm = c("CCR7", "IL7R", "TCF7", "SELL"),
  Proliferation = c("MKI67","TOP2A","STMN1","TYMS","PCNA","BIRC5","MCM5", "CENPF", "CDC20", "CCNB1", "CCNB2", "CDK1"))

make_dotplot_gene_with_bar <- function(seu, gene_sets, index, group, width, height) {
  suppressPackageStartupMessages({library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(cowplot)})
  genes <- unlist(gene_sets, use.names = FALSE)
  gene_group <- stack(gene_sets)
  colnames(gene_group) <- c("gene", "set")
  set_colors <- c(CTL = "#D73027", Temra = "#4575B4", Tem = "#1A9850", Tcm="yellow", Proliferation = "#984EA3")
  # --- 1) Seurat DotPlot (keep Seurat scaling) ---
  p_dot <- DotPlot(seu, features = genes, group.by = "disease") +
    coord_flip() + RotatedAxis() +
    theme(plot.margin = margin(5.5, 2, 5.5, 5.5),axis.text = element_text(face = "bold", color = "black",size = 14),
                axis.title = element_text(face = "bold", color = "black", size = 14),
                legend.text = element_text(face = "bold", color = "black", size = 14),
                legend.title = element_text(face = "bold", color = "black", size = 14)) # small right margin
  # extract gene order used by Seurat
  dpdat <- p_dot$data
  gene_order <- unique(as.character(dpdat$features.plot))
  # --- 2) Gene-set annotation bar (no legends here) ---
  ann <- gene_group %>% filter(gene %in% gene_order) %>% mutate(features.plot = factor(gene, levels = gene_order),x = 1)
  p_bar <- ggplot(ann, aes(x = x, y = features.plot, fill = set)) +
    geom_tile() + scale_fill_manual(values = set_colors) +
    scale_y_discrete(limits = gene_order) + scale_x_continuous(expand = c(0, 0)) +
    theme_void() + theme(plot.margin = margin(5.5, 0, 5.5, 0),legend.position = "none") 
  # --- 3) Legend extraction: take DotPlot legend OUT ---
  legend_grob <- cowplot::get_legend( p_dot + theme(legend.position = "right"))
  legend_bar <- cowplot::get_legend(p_bar + theme(legend.position = "right") + guides(fill = guide_legend(title = "Gene set")))
  legend_all <- cowplot::plot_grid(legend_grob,legend_bar,ncol = 1,align = "v",rel_heights = c(1, 0.6))
  # remove legend from the DotPlot itself so it cannot overlap
  p_dot_noleg <- p_dot + theme(legend.position = "none", plot.margin = margin(5.5, 0, 5.5, 5.5))
  # --- 4) Assemble: DotPlot | bar | legend ---
  p_all <- cowplot::plot_grid(p_dot_noleg, p_bar,NULL,legend_all, nrow = 1,
    rel_widths = c(0.4, 0.06, 0.03, 0.45),align = "h", axis = "tb")
  ggsave(paste0("dot plot of ", index, " index genes for pooled ", group, ".pdf"),
         plot = p_all, dpi = 300, width = width, height = height, bg = "white")
  ggsave(paste0("dot plot of ", index, " index genes for pooled ", group, ".svg"),
         plot = p_all, dpi = 300, width = width, height = height, bg = "white")
}
make_dotplot_gene_with_bar(product, gene_sets, "index genes", "product", 4, 7)

pseudobulk_DEGs <- read.csv("CART DEseq2 DEGs with pseudo bulk of product  cells plan B.csv",sep=",") %>% filter(abs(log2FoldChange)>1,padj<0.05) 
Idents(product) <- "disease"
group.markers <- FindAllMarkers(object = product, only.pos = T)
DEGs <- group.markers %>% filter(gene %in% rownames(pseudobulk_DEGs)) %>% filter(gene !="XIST") %>%
  filter( p_val_adj<0.1, !grepl("TRAV|TRBV|TRDV|TRGV|TRGC|TRAJ|TRAC|TRDC|ENSG|LINC|-AS1|orf|HLA|MIR|-DT",gene)) %>% 
  mutate(abs_log2FC = abs(avg_log2FC)) %>% slice_max(order_by = abs_log2FC, n = 100, with_ties = FALSE)
signatures = list(CTL = c("GZMB","GZMH","GZMA", "GNLY","NKG7","PRF1","KLRD1","KLRG1","FGFBP2","CCL4","CCL5"),
Proliferation = c("MKI67","TOP2A","STMN1","TYMS","PCNA","BIRC5"),
Tpex=c("TCF7","LEF1","IL7R","CCR7","SLAMF6","CXCR5"),
Tex=c("PDCD1","LAG3","TIGIT","CTLA4", "HAVCR2","TOX","CXCL13","ENTPD1","LAYN"))
genes <- unlist(signatures, use.names = FALSE)
gene_group <- stack(signatures)
colnames(gene_group) <- c("gene", "signature")
annotation_col <- gene_group %>% column_to_rownames(var="gene")
ann_colors <- list(set = c(CTL = "#D73027",Proliferation = "#984EA3", Tpex = "#1A9850",Tex = "#4575B4" ))
count <- read.csv("CART vsd normalized counts of product  cells plan B.csv",sep=",",header=T) %>% select(UPN749.P,UPN746.P,UPN787.P,UPN775.P,UPN685.P)
colnames(count) <- sub("..$", "",colnames(count))
sigcount <- scale(t(count[gene_group$gene,]))
cols <- colorRampPalette(colors = brewer.pal(9,"RdBu") )
svglite("product signature heatmap.svg", width = 9, height = 2)
pheatmap(sigcount,color=rev(cols(100)),fontsize = 8,cluster_rows=F,cluster_cols=F,clustering_distance_rows = "correlation",
         annotation_col = annotation_col,annotation_colors = ann_colors,clustering_distance_cols = "correlation", show_rownames = T, clustering_method = "ward.D2")
dev.off()
DEGcount <- scale(t(count[DEGs$gene,]))
svglite("product DEGs heatmap.svg", width = 18, height = 3)
pheatmap(DEGcount,color=rev(cols(100)),fontsize = 8,cluster_rows=F,cluster_cols=T,clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation", show_rownames = T, clustering_method = "ward.D2")
dev.off()

df <- FetchData(harmony_inte, vars = c("orig.ident","TCR_positive","Cell_type"))
count_table <- df %>%  group_by(Cell_type, orig.ident) %>% summarise(n = n())
write.table(count_table, "cell counts for the CD4 and CD8 plan A.csv",sep=",",row.names=F,quote=F)
label_genes <- c(
  ## ---- CR-enriched ----
  # metabolism / mTORC1 & nutrient sensing
  "RRAGD", "RHEBL1", "TBC1D7", "PASK", "SLC7A3", "SLC28A3","H1-1",
  "DGKA", "ENPP2", "CHROMR", "INPP5A", "PIK3R6", "PIK3IP1","STAG3",
  # cytokine receptors & inflammatory signaling
  "IFNGR2", "IL23R", "IL1RAP", "IL1R2", "IL31RA", "IL1RN",
  "TRAF3IP2", "PELI2", "TIFA", "ALPK1",
  # transcriptional regulators
  "MYB", "IRF1", "STAT1", "ATF3", "ARID5A", "PLAG1", "NPAS2", "DACH1",
  # TCR-proximal / GPCR signaling
  "DAPP1", "TIAM2", "DENND6B", "CAMK2D", "RAPGEF4", "GRK5", "RGS2","TMIGD2",
  
  ## ---- PD-enriched ----
  # TEMRA / terminal effector differentiation
  "CX3CR1", "S1PR5", "KIR3DL1", "GZMH", "GZMK", "EOMES", "TYROBP",
  # senescence / growth arrest / apoptosis
  "CDKN2A", "TP73", "FOXO3", "PHLDA1", "HIPK2", "CASP1",
  # growth & metabolic signaling
  "HRAS", "SOS1", "PTPN12", "ESRRA",
  # costimulation / TNF superfamily
  "CD40LG", "TNFSF9", "TNFSF14", "TNFRSF18", "CD80", "CD86",
  # type 2 module
  "IL4", "IL5", "HPGDS", "PTGER2", "IL9R",
  # inhibitory / regulatory
  "SIRPG", "TGFBR3", "BST2",
  # chemokine receptors
  "CCR1", "CCR2", "CCR5", "CCR6", "CCR8", "GPR15",
  # chromatin / transcription
  "KAT2B", "ZBTB7A", "RARG"
)
label_genes <- unique(label_genes)
DE_analyse <- function(seu, name, substring, cells, group,label_sig_only = TRUE ) {
  if (cells %in% c("CD4","CD8")) {
    sub <- subset(seu, subset=Cell_type==cells)
  }  else if (cells =="") { sub = seu } else {
    sub <- subset(seu, subset=annotation == cells)}
  pb <- make_pseudobulk(sub, group_var = "orig.ident") 
  pb <- pb[-c(1:3),grepl(substring,colnames(pb))]
  exclude_pattern <- "^MT-|^TRAV|^TRBV|^TRDV|^TRGV|^TRGC|^TRAJ|^TRAC|^TRDC|^ENSG|^LINC|^MIR|-AS1|-DT|orf|^HLA|^RPS|^RPL"
  pb <- pb[!grepl(exclude_pattern, rownames(pb)), ]
  meta <- sub@meta.data %>% dplyr::select(orig.ident, disease) %>% distinct()
  meta <- meta[match(colnames(pb), meta$orig.ident), ] 
  rownames(meta) <- meta$orig.ident
  dds <- DESeqDataSetFromMatrix(countData = pb,colData = meta, design = ~ disease )
  dds <- DESeq(dds)
  # Check coefficient name first
  print(resultsNames(dds))
  
  # LFC shrinkage for small n
  res <- results(dds, contrast = c("disease", "CR", "PD"))
  res <- lfcShrink(dds, contrast = c("disease", "CR", "PD"), type = "ashr") 
  
  res <- as.data.frame(res[order(res$padj),])
  write.table(res,paste(name,"DEseq2 DEGs with pseudo bulk of",group, cells, "cells plan B.csv"),sep=",",quote=F)
  rld <- rlog(dds, blind = FALSE)
  write.table(assay(rld),  paste(name, "rlog normalized counts of", group, cells, "cells plan B.csv"), sep=",", quote=F)
  res$diffexpressed <- "NO"
  res$diffexpressed[res$log2FoldChange > 0 & res$padj < 0.05] <- "UP"
  res$diffexpressed[res$log2FoldChange < 0 & res$padj < 0.05] <- "DOWN"
  ## ---- CHANGED: label only genes of interest --------------------
  res$gene <- rownames(res)
  res$to_label <- res$gene %in% label_genes
  if (label_sig_only) res$to_label <- res$to_label & res$diffexpressed != "NO"
  
  res$delabel <- NA
  res$delabel[res$to_label] <- res$gene[res$to_label]
  res_filtered <- res %>% filter(abs(log2FoldChange)<10,padj>10^-12)
  ggplot(data=res_filtered, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) +
    geom_point(size=0.5, alpha=0.3) +labs(x = "Log2 fold change", y = "-log10(FDR)") +
    geom_point(data = subset(res_filtered, !is.na(delabel)), aes(x=log2FoldChange, y=-log10(padj)), size=1.2) +
    theme_classic()+
    theme(legend.position = "none",axis.text = element_text(size = 14,face = "bold", color = "black"),
          axis.title = element_text(face = "bold", color = "black", size = 14),
          legend.text = element_text(face = "bold", color = "black", size = 14),
          legend.title = element_text(face = "bold", color = "black", size = 14)) +
    geom_text_repel(data = subset(res_filtered, !is.na(delabel) & log2FoldChange>0),
                    aes(label=delabel),color="black",
                    size=3, fontface="bold",
                    max.overlaps = 10, seed = 1, box.padding = 0.1, point.padding = 0,
                    force = 10, force_pull = 0.05, max.iter = 50000, max.time = 5,
                    nudge_y = max(-log10(res_filtered$padj[!is.na(res_filtered$delabel)]), na.rm=TRUE) * 0.05,
                    min.segment.length = 0, segment.size = 0.2, segment.color = "grey50",
                    xlim = c(0, NA),ylim = c(0, 12)) +
    geom_text_repel(data = subset(res_filtered, !is.na(delabel) & log2FoldChange<0),
                    aes(label=delabel),color="black",
                    size=3, fontface="bold",
                    max.overlaps = 10, seed = 1, box.padding = 0.1, point.padding = 0,
                    force = 10, force_pull = 0.05, max.iter = 50000, max.time = 5,
                    nudge_y = max(-log10(res_filtered$padj[!is.na(res_filtered$delabel)]), na.rm=TRUE) * 0.05,
                    min.segment.length = 0, segment.size = 0.2, segment.color = "grey50",
                    xlim = c(NA, 0),ylim = c(0, 12)) +
    theme(legend.position = "right", legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 3, shape = 16, label=""))) +
    scale_color_manual(values=c("DOWN"="blue", "NO"="black", "UP"="red"), breaks=c("DOWN", "UP"),labels=c("Down", "Up")) +
    geom_hline(yintercept=-log10(0.05), col="red", size=0.3, alpha=0.5) +
    coord_cartesian(xlim = c(-10, 10), ylim = c(0, 12))
  ggsave(paste(name, "volcano plot of DEGs with pseudo bulk of", group, cells, "T cells.svg"), width=6, height=5, dpi=300)
  return(res_filtered)
}
product_DEG <- DE_analyse(harmony_inte,"CART","-P","","product")
anno <- unique(harmony_inte$annotation)
for (i in c(1:length(anno))) {
  DE_analyse(harmony_inte,"CART","-P",anno[i],"product")
  #DE_analyse(harmony_inte,"CART","-B",anno[i],"blood")
}

gs_hallmark <- geneset::getMsigdb(org = "human",category = "H")
gs_kegg <- geneset::getMsigdb(org = "human",category = "C2-CP-KEGG")
entrez=read.table("/net/nfs-irwrsrchnas01/labs/sforman/Seq/CD19CART_singlecell/T_cell_integration/human genes Ensembl115.txt",header = T,sep="\t",stringsAsFactors = F) %>% dplyr::select(Gene.name,NCBI.gene..formerly.Entrezgene..ID) %>% filter(NCBI.gene..formerly.Entrezgene..ID!="NA") %>% distinct()
gsea_analysis <- function(markerlist, group, category, gs, qcutoff, width, height) {
  DEG <- markerlist %>% dplyr::select(log2FoldChange) %>% rownames_to_column(var="gene")  %>% filter(!is.na(log2FoldChange))       
  DEG_entrez <- merge(DEG,entrez, by.x="gene",by.y="Gene.name") %>% dplyr::select(NCBI.gene..formerly.Entrezgene..ID,log2FoldChange) %>% 
      group_by(NCBI.gene..formerly.Entrezgene..ID) %>% slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>% arrange(desc(log2FoldChange)) %>% ungroup()
  lapply(sets, function(s) intersect(s, label_genes)) 
  geneList <- setNames(DEG_entrez$log2FoldChange,  as.integer(DEG_entrez$NCBI.gene..formerly.Entrezgene..ID))
  gse <- genGSEA(genelist = geneList, geneset = gs,p_cutoff = 1,q_cutoff = 1,min_gset_size =10,max_gset_size = 2000,set_seed = TRUE)
  gse_df <- gse[["gsea_df"]]
  signif_results <- gse_df %>% filter(pvalue<qcutoff) %>% filter(!ID %in% c("KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION","KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",
                                                                            "KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS","KEGG_ASTHMA","KEGG_PRION_DISEASES",
                                                                           "KEGG_TYPE_I_DIABETES_MELLITUS","KEGG_NON_SMALL_CELL_LUNG_CANCER","KEGG_VIRAL_MYOCARDITIS",
                                                                           "KEGG_HYPERTROPHIC_CARDIOMYOPATHY_HCM") )
  ggplot(signif_results, aes(x= -log10(pvalue), y=reorder(ID,-log10(pvalue)), size=-log10(pvalue), color=NES)) + 
    geom_point(alpha = 0.8) + theme_classic(base_size = 12) +  geom_point(data = subset(signif_results, pvalue < 0.05),aes(size = -log10(pvalue)),
                                                            shape = 21, fill = NA, color = "black", stroke = 0.8) + labs(y = "Pathways", x = "-log10(pvalue)") +
    scale_color_gradient(high = "red2",  low = "mediumblue", space = "Lab") +
    scale_size(range = c(1, 7)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, face = "bold", color = "black"),
      axis.text.y  = element_text(size = 12, face = "bold", color = "black"),
      axis.title   = element_text(size = 14, face = "bold", color = "black"),
      legend.text  = element_text(size = 12, color = "black"),
      legend.title = element_text(size = 14, face = "bold", color = "black"))
  ggsave(paste(group,"CRvsPD dot plot of GSEA",category,"pathways.svg"),dpi=300,width=width, height = height)
  return(gse_df)
}
product_pathways <- gsea_analysis(product_DEG,"product","Hallmark",gs_hallmark,0.05,9.5,3.6)
product_kegg_pathways <- gsea_analysis(product_DEG,"product","KEGG",gs_kegg,0.05,9.5,4.5)

md <- harmony_inte@meta.data %>% rownames_to_column("cell") %>% mutate(cluster = as.factor(.data[["annotation"]])) 
meta_pct <-  md %>% group_by(cluster) %>% summarise(n = n(), .groups = "drop") %>% mutate(percent = 100 * n / sum(n))

#make split umap with percentage label
make_dimplot_with_pecent <- function(obj, name, split_variable, variable, width, height, ncols) {
  md <- obj@meta.data %>% rownames_to_column("cell") %>% mutate(split = as.character(.data[[split_variable]]),cluster = as.factor(.data[[variable]])) 
  meta_pct <-  md %>% group_by(split, cluster) %>% summarise(n = n(), .groups = "drop") %>% group_by(split) %>% mutate(percent = 100 * n / sum(n))
  umap_df <- as.data.frame(Embeddings(obj, "umap.harmony")) %>% rownames_to_column("cell") %>% inner_join(md, by = "cell")
  label_pos <- umap_df %>% group_by(split, cluster) %>% summarise(umap_1 = median(umapharmony_1), umap_2 = median(umapharmony_2),.groups = "drop") %>%
    left_join(meta_pct, by = c("split", "cluster")) %>% mutate(label = sprintf("%.1f%%", percent)) %>% filter(split != "NA")
  label_pos[[split_variable]] <- label_pos$split
  DimPlot(obj, reduction = "umap.harmony",pt.size = 1, ncol = ncols, label = F, group.by = variable, split.by = split_variable) +
    geom_text(data = label_pos,aes(x = umap_1, y = umap_2, label = label),size = 3,color = "black", fontface = "bold")
  ggsave(paste(name,"harmony integration umap split by", split_variable, "with percentage of",variable,"plan B.svg"), dpi=300, width=width, height=height)
}
make_dimplot_with_pecent(harmony_inte,"CART","Cell_type","annotation", width=8,height=4)
make_dimplot_with_pecent(harmony_inte,"CART","disease","annotation", width=8,height=4)
make_dimplot_with_pecent(harmony_inte,"CART","Cell_type","Phase", width=8,height=4)
make_dimplot_with_pecent(harmony_inte,"CART","disease","Phase", width=8,height=4)
make_dimplot_with_pecent(product,"product","Cell_type","annotation", width=8,height=4)
make_dimplot_with_pecent(blood,"blood","Cell_type","annotation", width=8,height=4)
make_dimplot_with_pecent(product,"product","disease","annotation", width=8,height=4)
make_dimplot_with_pecent(blood,"blood","disease","annotation", width=8,height=4)
make_dimplot_with_pecent(blood,"blood","Cell_type","Phase", width=8,height=4)
make_dimplot_with_pecent(blood,"blood","disease","Phase", width=8,height=4)
make_dimplot_with_pecent(product,"product","CAR_positive","annotation", width=6,height=3,ncols=2)
make_dimplot_with_pecent(product,"product","Cell_type","Phase", width=6,height=3,ncols=2)
make_dimplot_with_pecent(product,"product","disease","Phase", width=6,height=3,ncols=2)

make_dimplot_with_percent_by_sample(harmony_inte, name = "CART",reduction="umap.harmony",sample_var = "orig.ident",cluster_var = "annotation",width = 6, height = 4, outdir = ".")
make_dimplot_with_percent_by_sample(harmony_inte, name = "CART",reduction="umap.harmony",sample_var = "group",cluster_var = "annotation",width = 5, height = 4, outdir = ".")

score_product <- c("Prolif_score","CTL_score","Tex_score","Tpex_score")
score_blood <- c("Tem_score","Temra_score", "Tex_score","Tpex_score")
for (i in 1:length(score_product)) {
  FeaturePlot(product, features = score_product[i],reduction = "umap.harmony", split.by = "disease")+ theme(legend.position = "left") & NoAxes() & 
    scale_color_gradientn(colors = inferno(n = 100, direction = -1) )
  ggsave(paste(name, reduction, "integrated UMAP", score_product[i], "scores in product.svg"), dpi=300, width=6, height=3)
  max_value <- max(product[[score_product[i]]], na.rm = TRUE)
  VlnPlot(product, score_product[i],group.by = "disease",pt.size = 0) + theme(legend.position = "none",plot.title = element_blank()) +
    stat_compare_means(comparisons = list(c("CR", "PD")), method = "wilcox.test",label = "p.signif", hide.ns = F,size = 8,fontface = "bold") + 
    ylim(c(0, max_value * 1.2)) + xlab("") + ylab(score_product[i]) +
    theme(axis.text = element_text(face = "bold", color = "black",size = 14),
          axis.title = element_text(face = "bold", color = "black", size = 14),
          legend.text = element_text(face = "bold", color = "black", size = 14),
          legend.title = element_text(face = "bold", color = "black", size = 14)) 
  ggsave(paste(score_product[i],"score in each product cluster after",name, reduction, "integration.svg"),dpi=300,width=3,height=4)}

for (i in 1:length(score_blood)) {
  FeaturePlot(blood, features = score_blood[i],reduction ="umap.harmony", split.by = "disease" ) + theme(legend.position = "left") & NoAxes() &
    scale_color_gradientn(colors = inferno(n = 100, direction = -1) )
  ggsave(paste(name, reduction, "integrated UMAP", score_blood[i], "scores in blood plan B.svg"), dpi=300, width=6, height=3)
  max_value <- max(blood[[score_blood[i]]], na.rm = TRUE)
  VlnPlot(blood, score_blood[i],group.by = "disease",pt.size = 0 )+ theme(legend.position = "none",plot.title = element_blank()) +
    stat_compare_means(comparisons = list(c("CR", "PD")), method = "wilcox.test", label = "p.signif",hide.ns = F,size = 8,fontface = "bold") + 
    ylim(c(0, max_value * 1.2)) + xlab("") + ylab(score_blood[i]) +
    theme(axis.text = element_text(face = "bold", color = "black",size = 14),
          axis.title = element_text(face = "bold", color = "black", size = 14),
          legend.text = element_text(face = "bold", color = "black", size = 14),
          legend.title = element_text(face = "bold", color = "black", size = 14)) 
  ggsave(paste(score_blood[i],"score in each blood cluster after",name, reduction, "integration plan B.svg"),dpi=300,width=3,height=4)}
#make barplot for groups and samples
p_to_stars <- function(p) {if (is.na(p)) return("")
  if (p < 1e-4) "****"
  else if (p < 1e-3) "***"
  else if (p < 1e-2) "**"
  else if (p < 0.05) "*"
  else ""}
barplot_group <- function(obj, name, group, variable, subgroup, width, height) {
  p <- DimPlot(obj, reduction = "umap.harmony", group.by = variable)
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
    geom_text(aes(label = ifelse(pct >= 2, lab, "")),position = position_fill(vjust = 0.5), size = 6) +
    xlab("group") + ylab("Proportion") +
    scale_fill_manual(values = cluster_colors) +
    theme(legend.title = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text = element_text(size = 14,face = "bold", color = "black"),
          axis.title = element_text(face = "bold", color = "black", size = 14),
          legend.text = element_text(face = "bold", color = "black", size = 14)) +
    coord_cartesian(ylim = c(0, 1.12)) + 
    labs(caption = if (!is.na(pval))
        sprintf("Chi-square test: p = %.3g", pval)
      else NULL)
  
  # add bridge only if significant
  if (stars != "") {
    g <- g + geom_segment(x = x_pos[1], xend = x_pos[1],y = 1.02, yend = 1.05,linewidth = 0.8) + 
      geom_segment(x = x_pos[2], xend = x_pos[2],y = 1.02, yend = 1.05,linewidth = 0.8) + 
      geom_segment(x = x_pos[1], xend = x_pos[2],y = 1.05, yend = 1.05,linewidth = 0.8) + 
      annotate("text", x = mean(x_pos), y = 1.08,label = stars,size = 6)}
  outfile <- paste(name, "percentage of", variable, "in", group, "after harmony integration plan B.svg")
  ggsave(outfile, plot = g, dpi = 300, width = width, height = height)
}

groups <- c("Product","Blood")
for (i in c(1:length(groups))) {
  barplot_group(harmony_inte,"CART",groups[i],"annotation","subgroup",5,5)
  barplot_group(harmony_inte,"CART",groups[i],"Phase","subgroup",5,5)
}
groups <- c("685","749","775","787", "746")
for (i in c(1:length(groups))) {
  barplot_group(harmony_inte,"CART",groups[i],"annotation","orig.ident",4,5)
  barplot_group(harmony_inte,"CART",groups[i],"Phase","orig.ident",4,5)
}

#make barplot for groups and samples
DEG_lis <- list.files(pattern="DEseq2 DEGs with pseudo bulk of.*cells plan B.csv")
DEG_lis <- DEG_lis[c(1,2,4,6,11,12,14,16)]
count_lis <- list.files(pattern="vsd normalized counts of.*cells plan B.csv")
count_lis <- count_lis[c(1:4,6:9)]
DEGs <- list()
DEGs_count <- list()
names <- c("blood","blood CD4","blood CD8","blood Cycling T","product","product CD4","product CD8","product Cycling T")
for (i in c(1:length(DEG_lis))) {
  DEGs[[i]] <- read.table(DEG_lis[i],sep=",",header = T,row.names = 1) %>% filter(abs(log2FoldChange)>1,padj<0.05)
  genes <- rownames(DEGs[[i]])[!grepl("ENSG|LINC|AS-1|AS1",rownames(DEGs[[i]]))]
  DEGs_count[[i]] <- read.table(count_lis[i],sep=",",header = T,row.names = 1)
  DEGs_count[[i]] <- DEGs_count[[i]][genes,]
  counts <- t(scale(t(DEGs_count[[i]])))
  colnames(counts) <- c("UPN685","PD","UPN775","UPN787")
  counts <- as.data.frame(counts) %>% mutate(CR=(UPN685+UPN775+UPN787)/3) %>% select(CR, PD)
  cols <- colorRampPalette(colors = brewer.pal(9,"RdBu") )
  svglite(paste(names[i],"DEGs heatmap.svg"), width = 8, height = 20)
  #annotation_col <- data.frame(row.names=colnames(DEGs_count[[i]]),group=c("CR","PD","CR","CR"))
  pheatmap(counts,color=rev(cols(100)),fontsize = 6,cluster_rows=TRUE,cluster_cols=F,clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation", show_rownames = T, clustering_method = "ward.D2")
  dev.off()
}

#separate CD4 and CD8 cells

celltypist_count <- function(seu,name) {
  subset <- subset(seu,subset=Cell_type==name)
  raw_mat <- GetAssayData(subset, assay = "RNA", slot = "counts")
  celltypist_mat <- t(as.matrix(raw_mat))
  write.csv(celltypist_mat, paste0(name,"_celltypist_counts.csv"))
}
celltypist_count(harmony_inte,"CD8")

CD4 <- subset(harmony_inte,subset=Cell_type=="CD4")
CD8 <- subset(harmony_inte, subset=Cell_type=="CD8")

data_integration <- function(seu, name, method, reduction, resolution) {
  celltypist_pred <- read.table(paste0(name,"_celltypist/predicted_labels.csv"),sep=",",header = T,row.names = 1)
  seu$celltypist <- celltypist_pred$majority_voting
  options(future.globals.maxSize = 100 * 1024^3)
  #integrated <- IntegrateLayers(object = seu, method = method, orig.reduction = "pca", new.reduction = reduction,verbose = FALSE)
  #integrated <- JoinLayers(integrated)
  integrated <- FindNeighbors(seu, reduction = reduction, dims = 1:30)
  res <- seq(0.1, 1.0, by = 0.1)
  integrated <- FindClusters(integrated, resolution = res)
  clustree(integrated, prefix = "RNA_snn_res.")
  ggsave(paste(name,"cluster tree plot plan B.svg"), dpi=300, width=10, height=12)
  integrated <- FindClusters(integrated, resolution = resolution, cluster.name = "integrated.cluster")
  integrated <- RunUMAP(integrated, dims = 1:30, reduction = reduction)
  DimPlot(integrated, reduction = "umap", group.by = c("orig.ident", "integrated.cluster","subgroup","disease","celltypist","Phase"))
  ggsave(paste(name,reduction,"integration cluster umap resolution",resolution, "plan B.svg"), dpi=300, width=15, height=10)
  DimPlot(integrated, reduction = "umap",pt.size = 1, ncol = 2, label = F, group.by = "integrated.cluster", split.by = "orig.ident")
  ggsave(paste(name,reduction,"integration sample split umap", resolution, "plan B.svg"), dpi=300, width=7, height=14)
  var_genes <- VariableFeatures(integrated)
  var_genes <- var_genes[!grepl("TRAV|TRBV|TRDV|TRGV|TRGC|TRAJ|TRAC|TRDC|ENSG|LINC|-AS1",var_genes)]
  plan("sequential")
  cluster.markers <- FindAllMarkers(object = integrated, features = var_genes, only.pos = T, min.pct = 0.25, logfc.threshold = 0.15)
  write.table(cluster.markers,paste(name,"marker genes list of Harmony integration resolution", resolution, "plan B.csv"),sep=",",quote=F,row.names=F)
  top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  cells_use <- WhichCells(integrated, downsample = 300)
  DoHeatmap(integrated, cells = cells_use, features = top10$gene,group.by="integrated.cluster",disp.max = 4,raster = TRUE,slot = "scale.data") + scale_fill_gradientn(colors = c("#89cff0", "white", "red"))
  ggsave(paste("Heatmap of top 10 variable genes in each Harmony integrated cluster", name, "resolution", resolution, "plan B.png"), dpi=300,width=12, height=18)
  saveRDS(integrated,file = paste(name, reduction, "integrated CART cells planB.rds"))
  return(integrated)
}

CD4_integrated <- data_integration(CD4,"CD4", HarmonyIntegration,"harmony",0.7)
CD8_integrated <- data_integration(CD8,"CD8", HarmonyIntegration,"harmony",0.7)
celltypist_pred <- read.table(paste0("CD8","_celltypist/predicted_labels.csv"),sep=",",header = T,row.names = 1)
CD8_integrated$celltypist <- celltypist_pred$majority_voting
DimPlot(CD8_integrated, reduction = "umap", group.by = c("orig.ident", "integrated.cluster","subgroup","disease","celltypist","Phase"))
ggsave(paste("CD8","harmony","integration cluster umap resolution",0.5, "plan B.svg"), dpi=300, width=15, height=10)
saveRDS(CD8_integrated,file = paste("CD8","harmony", "integrated CART cells planB.rds"))

ref_ <- readRDS("/net/isi-dcnl/ifs/user_data/lkwak/Data/Jibin/reference/SingleR_Ref/ref_humanstimCD4T_17types.rds")
pred <- SingleR(test = as.SingleCellExperiment(CD4_integrated),ref = ref_,labels = ref_$cluster.id)
cluster_pred <- SingleR(test = as.SingleCellExperiment(CD4_integrated),ref = ref_,labels = ref_$cluster.id, clusters = CD4_integrated$integrated.cluster)
CD4_integrated$labels <- cluster_pred$labels
CD4_integrated$pruned.labels<- cluster_pred$pruned.labels
DimPlot(CD4_integrated, reduction = "umap", group.by = "labels")

CD4_integrated <- readRDS("CD4 harmony integrated CART cells planB.rds")
score_cols <- c("Prolif_score","Treg_score","Tex_score","CTL_score","Tem_score","Tcm_score","Tn_score","Tfh_score","Th1_score","Th2_score","Th17_score")
CD4_integrated@meta.data <- CD4_integrated@meta.data %>% dplyr::select(-any_of(score_cols))
Tn_sig   <- c("CCR7","IL7R","TCF7","LEF1","LTB","NOSIP")
Tcm_sig  <- c("CCR7","IL7R","TCF7","SELL","MAL")
Tem_sig  <- c("GZMK","CXCR3","IFNG","CCL5")
CTL_sig  <- c("GZMB","PRF1","NKG7","GNLY")
Treg_sig <- c("FOXP3","IL2RA","CTLA4","IKZF2")
Tfh_sig  <- c("CXCR5","PDCD1","ICOS","BCL6")
Tex_sig  <- c("PDCD1","TOX","LAG3","TIGIT","HAVCR2")
Prolif_sig   <- c("MKI67","TOP2A","HMGB2","TYMS","STMN1")
Th1_sig <- c("TBX21","IL12RB2","CXCR3","IFNG")
Th2_sig <- c("GATA3","IL4","IL5","IL13","CCR4")
Th17_sig <- c("RORC","IL17A","IL17F","CCR6","KLRB1")
CD4_integrated <- AddModuleScore(object=CD4_integrated,seed.use = 123, features = list(Prolif_sig,Treg_sig,Tex_sig,CTL_sig,Tem_sig,Tcm_sig,Tn_sig,Tfh_sig,Th1_sig,Th2_sig,Th17_sig),name=score_cols)
CD4_integrated@meta.data <- CD4_integrated@meta.data %>% dplyr::rename(Prolif_score = Prolif_score1,Treg_score = Treg_score2,Tex_score = Tex_score3, CTL_score = CTL_score4,Tem_score = Tem_score5, 
                                                               Tcm_score = Tcm_score6, Tn_score = Tn_score7, Tfh_score = Tfh_score8, Th1_score = Th1_score9, Th2_score = Th2_score10, Th17_score=Th17_score11)

genes <- c("CCR7","IL7R","SELL","TCF7","LEF1","LTB","MAL","GZMK","CXCR3","IFNG","CCL5","GZMB","PRF1","NKG7","GNLY","CTSW","KLRD1","FOXP3","IL2RA","CTLA4",
           "IKZF2","TNFRSF18","TIGIT","CXCR5","PDCD1","ICOS","BCL6","SH2D1A","IL21","TBX21","IL12RB2","GATA3","IL4","IL5","IL13","CCR4","RORC","IL17A",
           "IL17F","CCR6","KLRB1","TOX","LAG3","HAVCR2","MKI67","TOP2A","HMGB2","TYMS","STMN1")
DoHeatmap(CD4_integrated, features = genes,group.by="integrated.cluster",disp.max = 4,raster = TRUE,slot = "scale.data") + scale_fill_gradientn(colors = c("#89cff0", "white", "red"))
ggsave(paste("Heatmap of signature genes in each Harmony integrated cluster", "CD4", "resolution", 0.7, "plan B.png"), dpi=300,width=12, height=18)
FeaturePlot(CD4_integrated, features = genes, reduction = "umap",ncol=7)
ggsave(paste("Featureplot of signature genes in each Harmony integrated cluster", "CD4", "resolution", 0.7, "plan B.png"), dpi=300,width=14, height=14)
new.cluster.ids <- c("0" = "Tn","1" = "Cycling late S histone high","2" = "Cycling Early T","3" = "Activated Cycling T","4" = "Treg","5" = "Tcm","6" = "Cycling G2M phase", "7" ="Tem_Tex")
annotation <- new.cluster.ids[as.character(CD4_integrated$integrated.cluster)]
names(annotation) <- colnames(CD4_integrated)
CD4_integrated <- AddMetaData(object = CD4_integrated,metadata = annotation,col.name = 'annotation')
make_dimplot_with_pecent(CD4_integrated,"CD4","orig.ident","annotation", width=7,height=14)
make_dimplot_with_pecent(CD4_integrated,"CD4","Cell_type","annotation", width=8,height=4)
make_dimplot_with_pecent(CD4_integrated,"CD4","Cell_type","Phase", width=8,height=4)
groups <- c("Product","Blood")
for (i in c(1:length(groups))) {
  barplot_group(CD4_integrated,"CD4",groups[i],"annotation","subgroup")
  barplot_group(CD4_integrated,"CD4",groups[i],"Phase","subgroup")
}
groups <- c("685","749","775","787")
for (i in c(1:length(groups))) {
  barplot_group(CD4_integrated,"CD4",groups[i],"annotation","orig.ident")
  barplot_group(CD4_integrated,"CD4",groups[i],"Phase","orig.ident")
}

anno <- unique(CD4_integrated$annotation)
for (i in c(1:length(anno))) {
  DE_analyse(CD4_integrated,"CD4","-P",anno[i],"product")
}
saveRDS(CD4_integrated,file = paste("CD4","harmony", "integrated CART cells planB.rds"))

CD8_integrated <- readRDS("CD8 harmony integrated CART cells planB.rds")
score_cols <- c("Prolif_score","Tpex_score","Tex_score","CTL_score","Tem_score","Tcm_score","Tn_score","Trm_score","Temra_score","ISG_score")
CD8_integrated@meta.data <- CD8_integrated@meta.data %>% dplyr::select(-any_of(score_cols))
Tn_sig   <- c("CCR7", "IL7R", "LTB", "LST1", "TCF7", "LEF1", "MALAT1", "JUN", "FOS")
Tcm_sig  <- c("CCR7", "IL7R","TCF7", "LEF1", "MALAT1")
Tem_sig  <- c("GZMK", "GZMH","IFNG", "CCL4", "CCL5", "CXCR3")
CTL_sig  <- c("GZMB", "PRF1", "NKG7", "GNLY", "IFNG", "TNF")
Temra_sig <- c("GZMB", "PRF1", "FGFBP2", "NKG7", "KLRD1", "KLRG1", "GNLY")
Tpex_sig  <- c("TCF7", "IL7R", "CXCR5", "SLAMF6", "TOX")
Tex_sig  <- c("PDCD1","TOX","LAG3","TIGIT","HAVCR2","CTLA4","TOX2","CXCL13")
Prolif_sig   <- c("MKI67","TOP2A","HMGB2","TYMS","STMN1", "PLK1", "CCNB1", "CDC20")
Trm_sig <- c("ITGAE","CD69", "CXCR6", "ZNF683", "ITGA1")
ISG_sig <- c("ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "OAS1", "OAS2", "OAS3")
CD8_integrated <- AddModuleScore(object=CD8_integrated,seed.use = 123, features = list(Prolif_sig,Tpex_sig,Tex_sig,CTL_sig,Tem_sig,Tcm_sig,Tn_sig,Trm_sig,Temra_sig,ISG_sig),name=score_cols)
CD8_integrated@meta.data <- CD8_integrated@meta.data %>% dplyr::rename(Prolif_score = Prolif_score1,Tpex_score = Tpex_score2,Tex_score = Tex_score3, CTL_score = CTL_score4,Tem_score = Tem_score5, 
                                                                       Tcm_score = Tcm_score6, Tn_score = Tn_score7, Trm_score = Trm_score8, Temra_score = Temra_score9,
                                                                       ISG_score = ISG_score10)
DotPlot(CD8_integrated, features = c(score_cols,"Cytokine1_socre","Cytokine2_score"),group = "integrated.cluster") + RotatedAxis() 
ggsave("CD8 dot plot of index scores for each cluster plan B.svg", dpi=300,width=8, height=7,bg="white")

genes <- c("IFNG", "TNF", "GZMH", "CCL4", "CCL5", "CXCR3","CCR7", "IL7R", "LTB", "LST1",  "LEF1", "MALAT1", "JUN", "FOS", "MALAT1",
           "TCF7", "IL7R", "CXCR5", "SLAMF6", "TOX","GZMB", "PRF1", "FGFBP2", "NKG7", "KLRD1", "KLRG1", "GNLY",
           "PDCD1","TOX","LAG3","TIGIT","HAVCR2","CTLA4","TOX2","CXCL13","MKI67","TOP2A","HMGB2","TYMS","STMN1", 
           "PLK1", "CCNB1", "CDC20","ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "OAS1", "OAS2", "OAS3")
DoHeatmap(CD8_integrated, features = genes,group.by="integrated.cluster",disp.max = 4,raster = TRUE,slot = "scale.data") + scale_fill_gradientn(colors = c("#89cff0", "white", "red"))
ggsave(paste("Heatmap of signature genes in each Harmony integrated cluster", "CD8", "resolution", 0.7, "plan B.svg"), dpi=300,width=12, height=14)
FeaturePlot(CD8_integrated, features = genes, reduction = "umap",ncol=7)
ggsave(paste("Featureplot of signature genes in each Harmony integrated cluster", "CD8", "resolution", 0.7, "plan B.png"), dpi=300,width=14, height=14)
new.cluster.ids <- c("0" = "Tem","1" = "Tn_Tcm","2" = "Cycling Early S","3" = "Cycling CTL","4" = "Temra","5" = "Tpex","6" = "Cycling S", "7" ="Cycling G2M", "8"="Stress activated")
annotation <- new.cluster.ids[as.character(CD8_integrated$integrated.cluster)]
names(annotation) <- colnames(CD8_integrated)
CD8_integrated <- AddMetaData(object = CD8_integrated,metadata = annotation,col.name = 'annotation')
make_dimplot_with_pecent(CD8_integrated,"CD8","orig.ident","annotation", width=7,height=14)
make_dimplot_with_pecent(CD8_integrated,"CD8","Cell_type","annotation", width=6,height=5)
make_dimplot_with_pecent(CD8_integrated,"CD8","Cell_type","Phase", width=6,height=5)
groups <- c("Product","Blood")
for (i in c(1:length(groups))) {
  barplot_group(CD8_integrated,"CD8",groups[i],"annotation","subgroup")
  barplot_group(CD8_integrated,"CD8",groups[i],"Phase","subgroup")
}
groups <- c("685","749","775","787")
for (i in c(1:length(groups))) {
  barplot_group(CD8_integrated,"CD8",groups[i],"annotation","orig.ident")
  barplot_group(CD8_integrated,"CD8",groups[i],"Phase","orig.ident")
}

anno <- unique(CD8_integrated$annotation)
for (i in c(1:length(anno))) {
  DE_analyse(CD8_integrated,"CD8","-P",anno[i],"product")
}
for (i in c(1,3,5,7,9)) {
  DE_analyse(CD8_integrated,"CD8","-B",anno[i],"blood")
}
saveRDS(CD8_integrated,file = paste("CD8","harmony", "integrated CART cells planB.rds"))

colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)
DimPlot(integrated, reduction = "umap",group.by = "cloneSize") + scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))
ggsave(paste(name, reduction, "integrated all samples UMAP CR plan B.png"), dpi=300, width=8, height=5)
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)
cells_keep <- colnames(harmony_inte) 
combined.CR.filtered <- lapply(combined.CR, function(df) {df[df$barcode %in% cells_keep, ]})

plot_clonaloverlay <- function(seu,group, name, reduction) {
  p_umap <- DimPlot(seu, reduction = "umap", group.by = "annotation")
  umap_colors <- unique(ggplot_build(p_umap)$data[[1]]$colour)
  cluster_ids <- unique(p_umap[[1]][["data"]][["annotation"]])
  cluster_colors <- setNames(umap_colors, cluster_ids)
  Idents(seu) <- "annotation"
  clonalOverlay(seu, reduction = "umap",  cutpoint = 5, bins = 10, facet.by = "orig.ident")+ scale_color_manual(values = cluster_colors,breaks = sort(names(cluster_colors))) +
    guides(color = guide_legend(override.aes = list(size = 5) )) + theme(axis.title = element_text(size = 16),
  axis.text  = element_text(size = 14),
  strip.text = element_text(size = 14, face = "bold"),
  legend.title = element_text(size = 16),
  legend.text  = element_text(size = 14),
  plot.title = element_text(size = 18, face = "bold"))
  ggsave(paste("clonalfrequency overlay in", group, "samples in", name, reduction, "integration plan B.svg"),dpi=300,width=6,height=6)
}
clonal_plots <- function(integrated, name, reduction, combined.CR) {
  Idents(integrated) <- "orig.ident"
  product <- subset(integrated, idents = c("UPN787-P","UPN775-P","UPN749-P","UPN685-P"))
  blood <- subset(integrated, idents = c("UPN787-B","UPN775-B","UPN749-B","UPN685-B"))
  Idents(product) <- "annotation"
  Idents(blood) <- "annotation"
  plot_clonaloverlay(product,"product",name,reduction)
  plot_clonaloverlay(blood,"blood",name,reduction)
  DimPlot(product, reduction = "umap",pt.size = 1, ncol = 2, label = F, split.by = "subgroup")
  ggsave(paste(name,reduction,"integration product subgroup split umap plan B.png"), dpi=300, width=8, height=4)
  DimPlot(blood, reduction = "umap",pt.size = 1, ncol = 2, label = F, split.by = "subgroup")
  ggsave(paste(name,reduction,"integration blood subgroup split umap plan B.png"), dpi=300, width=8, height=4)
  VlnPlot(integrated, c("clonalFrequency"),group.by = "orig.ident")
  ggsave(paste("clonal frequency in each sample after",name, reduction,"integration plan B.png"),dpi=300,width=7,height=7)
  VlnPlot(integrated, c("percent.mt"),group.by = "orig.ident")
  ggsave(paste("mitochondrial gene expression in each sample after",name, reduction, "integration.png"),dpi=300,width=7,height=7)
  VlnPlot(integrated, c("percent.mt"),group.by = "integrated.cluster")
  ggsave(paste("mitochondrial gene expression in each cluster after",name, reduction, "integration.png"),dpi=300,width=7,height=7)
  cells_keep <- colnames(integrated) 
  combined.CR.filtered <- lapply(combined.CR, function(df) {df[df$barcode %in% cells_keep, ]})
  clonalProportion(combined.CR.filtered, cloneCall = "aa") 
  ggsave(paste(name, reduction,"clonal proportion Product and Blood.png"), dpi=300, width=10, height=6)
  clonalHomeostasis(combined.CR.filtered, cloneCall = "aa")
  ggsave(paste(name, reduction,"clonal homeostasis Product and Blood.png"), dpi=300, width=10, height=6)
  percentKmer(combined.CR.filtered, cloneCall = "aa", chain = "TRB", motif.length = 5, top.motifs = 25)
  ggsave(paste(name, reduction,"percentKmer TRB all samples.png"), dpi=300, width=10, height=6)
  percentKmer(combined.CR.filtered, cloneCall = "aa", chain = "TRA", motif.length = 5, top.motifs = 25)
  ggsave(paste(name, reduction,"percentKmer TRA all samples.png"), dpi=300, width=10, height=6)
  compare_top_clone(combined.CR.filtered, "UPN685",10,"CART")
  compare_top_clone(combined.CR.filtered, "UPN775",10,"CART")
  compare_top_clone(combined.CR.filtered, "UPN749",10,"CART")
  compare_top_clone(combined.CR.filtered, "UPN787",10,"CART")
  compare_top_clone(combined.CR.filtered, "UPN746",10,"CART")
  integrated <- combineExpression(combined.CR.filtered, integrated, cloneCall = "aa", group.by = "sample",  proportion = FALSE,  cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
  DimPlot(integrated, reduction = "umap",group.by = "cloneSize") + scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))
  ggsave(paste(name, reduction, "integrated all samples UMAP CR plan B.png"), dpi=300, width=8, height=5)
}
clonal_plots(harmony_inte, "CART", "harmony",combined.CR)

CD4_integrated <- readRDS("CD4 harmony integrated CART cells planB.rds")
var_genes <- VariableFeatures(CD4_integrated)
var_genes <- var_genes[!grepl("TRAV|TRBV|TRDV|TRGV|TRGC|TRAJ|TRAC|TRDC|ENSG|LINC|-AS1|orf",var_genes)]
cluster.markers <- FindAllMarkers(object = CD4_integrated, features = var_genes, only.pos = T, min.pct = 0.3)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
cells_use <- WhichCells(CD4_integrated, downsample = 300)
DoHeatmap(CD4_integrated, cells = cells_use, features = top10$gene,group.by="integrated.cluster",disp.max = 4,raster = TRUE,slot = "scale.data") + scale_fill_gradientn(colors = c("#89cff0", "white", "red"))
ggsave(paste("Heatmap of signature genes in each Harmony integrated cluster", "CD4", "resolution", 0.7, "plan B.svg"), dpi=300,width=12, height=14)
