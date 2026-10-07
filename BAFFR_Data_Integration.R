source("/home/jibzhang/pipeline/Singlecell/Rscript/F_sc.R")
source("/home/jibzhang/pipeline/Singlecell/Rscript/CART_singlecell_functions.R")
library(Matrix)
library(Seurat)
library(reshape2)
library(DoubletFinder)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(scRepertoire)
library(VennDiagram)
library(SingleR)
library(ggridges)
library(ggpubr)
library(clustree)
library(genekitr)
setwd("/net/isi-dcnl/ifs/user_data/lkwak/Data/Jibin/project/BAFFR_singlecell")
qc <- read.table("qc_cutoff.tsv",sep="\t",header=T)
plotdata <- function(qc,qc_measure) {
  subset <- qc %>% filter(var==qc_measure) %>% select(-var,-median_,-mad_)
  melted <- melt(subset,id.vars="Sample")
  ggplot(melted, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot() + 
  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.3, seed = 1234),
             aes(color=factor(Sample)), show.legend = T) +
  xlab("") + ylab(qc_measure) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  theme_ipsum() + guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                        axis.title.y = element_text(hjust = 0.5))
ggsave(paste0("/net/isi-dcnl/ifs/user_data/lkwak/Data/Jibin/project/BAFFR_singlecell/id/",qc_measure,"_boxplot.png"),width=6,height = 6)
}
plotdata(qc,"nCount_RNA")
plotdata(qc,"nFeature_RNA")
plotdata(qc,"percent.mt")
plotdata(qc,"percent.ribo")

#add TCR information and merge and integration individual single cell data
workdir="/net/isi-dcnl/ifs/user_data/lkwak/Data/Jibin/project/BAFFR_singlecell/"
setwd(paste0(workdir,"T_cell_integration"))
meta <- read.table(paste0(workdir,"id/intake_file.tsv"),sep="\t",header=T)
samples <- meta$sample_id
CR <- list()
for (i in 1:length(samples)){
  CR[[i]] <- read.csv(paste0(workdir,samples[i],"/outs/per_sample_outs/",samples[i],"/vdj_t/filtered_contig_annotations.csv")) %>% 
    select(-sample)
  CR[[i]]$raw_clonotype_id <- sub("^([^_]+_){4}", paste0(samples[i],"_"), CR[[i]]$raw_clonotype_id)
  CR[[i]]$raw_consensus_id <- sub("^([^_]+_){4}", paste0(samples[i],"_"), CR[[i]]$raw_consensus_id)
  }
plan("multicore", workers = 4) 
combined.CR <- combineTCR(CR, samples = samples)
combined.CR <- addVariable(combined.CR, variable.name = "Group", variables = c(rep(c("Product","Blood"), 5)))
exportClones(combined.CR,  write.file = TRUE,dir = paste0(workdir,"Integration"),file.name = "allsample.CRclones.csv")

clonalQuant(combined.CR, cloneCall="aa", chain = "both",  scale = F)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("Count of unique clonotypes by aa.png", dpi=300, width=6, height=4) 
clonalQuant(combined.CR, cloneCall="aa", group.by = "Group", scale = F)
ggsave("Count of unique clonotypes by aa in each group.png", dpi=300, width=4, height=4)
clonalLength(combined.CR, cloneCall="aa", chain = "both", scale=F, group.by = "Group") 
ggsave("clonal CDR3 aa length distribution in each group.png", dpi=300, width=5, height=4)
clonalHomeostasis(combined.CR, group.by = "Group",cloneCall = "aa")
ggsave("clonal homeostasis based on aa in each group.png", dpi=300, width=5, height=4)
clonalDiversity(combined.CR, cloneCall = "aa", n.boots = 20, order.by = c("UPN685-P","UPN749-P", "UPN775-P","UPN787-P","UPN746-P", "UPN685-B","UPN749-B","UPN775-B","UPN787-B","UPN746-B"))
ggsave("clonal diversity gene+nt in each sample.png", dpi=300, width=8, height=4)
clonalRarefaction(combined.CR,plot.type = 3,cloneCall="aa",group.by = "Group",hill.numbers = 1,n.boots = 2)
ggsave("clonal Rarefaction gene+nt in each group.png", dpi=300, width=5, height=4)
clonalSizeDistribution(combined.CR, cloneCall = "aa", method= "ward.D2")
ggsave("clonal size distribution distance aa.png", dpi=300, width=6, height=4)

make_clonal_compare <- function(cr_list, patient) {
  sampP <- paste0(patient, "-P")
  sampB <- paste0(patient, "-B")
  clones <- intersect(unique(cr_list[[sampP]][["CTaa"]]),unique(cr_list[[sampB]][["CTaa"]]))
  countP <- cr_list[[sampP]] %>% filter(CTaa %in% clones) %>% dplyr::count(CTaa, name = "nP")
  countB <- cr_list[[sampB]] %>% filter(CTaa %in% clones) %>% dplyr::count(CTaa, name = "nB")
  # Merge and rank by total cells (nP + nB)
  top10_tbl <- full_join(countP, countB, by = "CTaa") %>% mutate(nP = coalesce(nP, 0L), nB = coalesce(nB, 0L), nTotal = nP + nB) %>%
    arrange(desc(nTotal)) %>% slice_head(n = 10)
  clonalCompare(cr_list, clones = top10_tbl$CTaa, order.by = c(sampP, sampB), cloneCall = "aa",graph = "alluvial")
  ggsave(paste("Top10 CDR3 aa clonal comparison between",patient,"Product and Blood.pdf"), dpi=300, width=6, height=6)
}
make_clonal_compare(combined.CR,"UPN749")
make_clonal_compare(combined.CR,"UPN685")
make_clonal_compare(combined.CR,"UPN775")
make_clonal_compare(combined.CR,"UPN787")
make_clonal_compare(combined.CR,"UPN746")

seurat_list2 <- lapply(samples, function(s) {
  seurat_obj <- readRDS(paste0(workdir, s, "/Seurat/plotQC/", s, ".SeuratObj.QC.rds"))
  seurat_obj@meta.data <- seurat_obj@meta.data[, !grepl("^DF.classifications|^pANN", colnames(seurat_obj@meta.data))]
  seurat_obj
})
names(seurat_list2) <- samples
merged <- merge(seurat_list2[[1]], y = seurat_list2[2:length(samples)], add.cell.ids = samples, project = "BAFFR")
combined <- combineExpression(combined.CR, merged, cloneCall="strict", group.by = "sample",proportion = T)
idx <- match(combined@meta.data$orig.ident, meta$sample_id)
combined@meta.data$group <- ifelse(is.na(idx), "Unknown", meta$group[idx])

adt_assay <- combined[["ADT"]]
counts_mat <- GetAssayData(combined, assay = "ADT", layer = "counts")
data_mat   <- GetAssayData(combined, assay = "ADT", layer = "data")
# Collapse EGFR + EGFR-antibody by summing (safe since they're mutually exclusive across samples)
collapse_rows <- function(mat, old_name, new_name) {
  combined_row <- mat[old_name, ] + mat[new_name, ]
  # Remove both old rows, add merged row
  mat <- mat[!rownames(mat) %in% c(old_name, new_name), ]
  mat <- rbind(mat, combined_row)
  rownames(mat)[nrow(mat)] <- new_name
  return(mat)
}
counts_mat <- collapse_rows(counts_mat,"EGFR", "EGFR-antibody")
data_mat   <- collapse_rows(data_mat, "EGFR",  "EGFR-antibody")
# Rebuild the assay
combined[["ADT"]] <- CreateAssayObject(counts = counts_mat)
combined[["ADT"]]@data <- data_mat

combined <- readRDS("merged samples with CR info.rds")
#use previous raw analysis cluster for filtering
BAFFR <- readRDS("T cells harmony integrated samples with CR info.rds")
DimPlot(BAFFR, reduction = "umap.harmony",pt.size = 1, ncol = 2, label = T, group.by="RNA_snn_res.0.1", split.by = "group")
ggsave("Dimplot split by group of merged before filtering.png", dpi=300, width=9, height=4)
FeaturePlot(BAFFR, features = c("CD14","FCGR3A", "CD3D","CD3E","CD3G","ALAS2","CD27","IL2RA", "PPBP", "CD8A", "CD8B", "CD4"), reduction = "umap.harmony",ncol=3)
ggsave("feature plot of merged all samples before filtering.png", dpi=300, width=9, height=9)
Idents(BAFFR) <- "RNA_snn_res.0.1"
BAFFR <- subset(BAFFR, idents = c("0","1","2","3","4","5"))
cells <- colnames(BAFFR)
Tcells <- subset(x = combined, cells = cells)

Tcells <- NormalizeData(Tcells, normalization.method = "LogNormalize", scale.factor = 10000)
Tcells <- FindVariableFeatures(Tcells, selection.method = "vst")
VariableFeatures(Tcells) = Trex::quietTCRgenes(VariableFeatures(Tcells),assay = "RNA")
options(future.globals.maxSize = 10 * 1024^3)
Tcells <- ScaleData(Tcells, vars.to.regress = c("nCount_RNA","percent.mt"))
Tcells <- RunPCA(Tcells)
ElbowPlot(object = Tcells,ndims=50)
ggsave("Elbowplot of 50 dims.png", dpi=300, width=6, height=5)
Tcells <- FindNeighbors(Tcells, dims = 1:20, reduction = "pca")
Tcells <- FindClusters(Tcells, resolution = 0.3, cluster.name = "unintegrated_clusters")
Tcells <- RunUMAP(Tcells,reduction = "pca", dims = 1:20,reduction.name = "umap.unintegrated")
DimPlot(Tcells, reduction = "umap.unintegrated", group.by = c("group", "seurat_clusters"))
ggsave("merged all samples UMAP group.png", dpi=300, width=10, height=5)
#plot to see cluster distribution amont group and samples
DimPlot(Tcells, reduction = "umap.unintegrated",pt.size = 1, ncol = 2, label = T, split.by = "group")
ggsave("Dimplot split by group after merge.png", dpi=300, width=9, height=4)
DimPlot(Tcells, reduction = "umap.unintegrated",pt.size = 1, ncol = 2, label = F, split.by = "orig.ident")
ggsave("Dimplot split by sample after merge.png", dpi=300, width=7, height=20)
#plot CR
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)
DimPlot(Tcells, reduction = "umap.unintegrated",group.by = "cloneSize") + scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))
ggsave("merged all samples UMAP CR.png", dpi=300, width=8, height=5)
saveRDS(Tcells, file = "merged samples with CR info.rds")

#check marker gene expression
features = c("CD4","CD8A","CD8B","CD3D","CD3E","CD3G","CD19","CAR-BAFFR-1","CAR-BAFFR-2","EGFR", "WPRE601","adt_EGFR-antibody")
FeaturePlot(Tcells, features = features, reduction = "umap.unintegrated",ncol=3)
ggsave("merged all samples UMAP cell marker genes.png", dpi=300, width=9, height=12)
DotPlot(Tcells, features = features) + RotatedAxis() 
ggsave("dot plot after all sample merge.png", dpi=300,width=7, height=6,bg="white")
RidgePlot(Tcells, c("adt_EGFR-antibody"),group.by = "orig.ident") + theme(axis.title.y = element_text(hjust = 0.5),axis.title.x = element_text(hjust = 0.5))
ggsave("EGFR Ridgeplot after all sample merge.png",dpi=300,width=9,height=5)
Tcells = CellCycleScoring(object = Tcells, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes)
DimPlot(Tcells,group.by='Phase',reduction="umap.unintegrated")
ggsave("CellCycle_FeaturePlot_all_samples.png",dpi=300,width=5,height=5)

#batch correct ADT with ADTnorm
# Tcells <- NormalizeData(Tcells, normalization.method = "CLR",margin = 2,assay = "ADT")
# cell_x_feature <- Tcells@meta.data
# cell_x_adt <- t(data.frame(Tcells[["ADT"]]@data))
# cell_x_adt_egfr <- cell_x_adt[, "EGFR-antibody", drop = FALSE]
# cell_x_feature <- cell_x_feature %>% dplyr::rename(sample = orig.ident) %>%
#   mutate(batch=sub("-.*","",sample)) %>% select(sample,batch)
# cell_x_adt_norm = ADTnorm(
#   cell_x_adt = cell_x_adt_egfr, exclude_zeroes=T,
#   cell_x_feature = cell_x_feature,
#   save_outpath = "/net/isi-dcnl/ifs/user_data/lkwak/Data/Jibin/project/BAFFR_singlecell/T_cell_integration/ADTnorm",
#   study_name = "EGFR_ADTnorm",peak_type = "midpoint",
#   trimodal_marker = "EGFR-antibody",valley_density_adjust=0.5,
#   brewer_palettes = "Dark2",shoulder_valley=T,
#   save_fig = TRUE,shoulder_valley_slope=-2
# )

Idents(Tcells) <- "orig.ident"
df_adt <- FetchData(Tcells, vars = "adt_EGFR-antibody") %>% mutate(Identity = Idents(Tcells))

# Step 3: summarise per Identity
valley_or_shoulder <- df_adt %>%
  group_by(Identity) %>%
  summarise({
    out <- pick_valley_then_shoulder(
      raw_x = `adt_EGFR-antibody`,
      valley_adjust = 1,
      valley_threshold = 1,
      spar = 0.6, z_thresh = 3, minx = 0.64,
      peak_frac = 0.40, prom_frac = 0.2, min_prom_drop = 0.07
    )
    tibble(event_x = out$event_x, event_y = out$event_y,
           type = out$type, percentage = out$percentage,
           main_peak_x = out$main_peak_x, main_peak_y = out$main_peak_y)
  }, .groups = "drop")

write.table(valley_or_shoulder,"ADT cutoff of each samples.txt", sep="\t",row.names = F,quote=F)

dens_df <- df_adt %>% group_by(Identity) %>%
  reframe({
    d <- density(`adt_EGFR-antibody`, n = 512)   # add bw=... if you want smoother
    tibble(x = d$x, y = d$y)
  })

valley_or_shoulder <- valley_or_shoulder %>% mutate(
  Identity = factor(Identity, levels = c("UPN685-P","UPN749-P", "UPN775-P","UPN787-P", "UPN746-P", "UPN685-B","UPN749-B","UPN775-B","UPN787-B","UPN746-B")),
  y_num = as.numeric(Identity),
  y0 = y_num - 0.45,   # segment half-height; tweak as you like
  y1 = y_num + 0.45)

ggplot(dens_df, aes(x = x, y = Identity, height = y)) +
  geom_ridgeline(fill = "lightblue", alpha = 0.8, size = 0.3) +
  geom_segment(
    data = valley_or_shoulder ,
    aes(x = event_x, xend = event_x, y = y0, yend = y1),
    inherit.aes = FALSE,
    color = "red",
    linewidth = 0.7,linetype = "dashed",alpha = 0.5) +
  labs(x = "ADT EGFR (antibody)", y = NULL)
ggsave("EGFR ridgeplot with cutoff line at the valley.png", dpi=300,width=12, height=7)

#combined the identified cutoff with the adt data
cuts <- valley_or_shoulder %>% select(Identity, event_x)
cell_tab <- df_adt %>% mutate(Cell = rownames(.)) %>% left_join(cuts, by = "Identity") %>% mutate(EGFR_Status = `adt_EGFR-antibody` > event_x)

#get other counts
df <- FetchData(Tcells, vars = c("orig.ident","CAR-BAFFR-1","CAR-BAFFR-2","EGFR","CD4","CD8A","CD8B","WPRE601"))
all_info <- merge(cell_tab,df,by.x="Cell",by.y="row.names" )
all_info <- all_info %>% mutate(
  CD4_pos= CD4>0,CD8_pos=(CD8A > 0 | CD8B > 0),
  CARrna_pos = (`CAR-BAFFR-1` > 0 | `CAR-BAFFR-2` > 0|WPRE601>0),
  CARrna_pos = replace_na(CARrna_pos, FALSE),
  EGFRrna_pos = EGFR>0) %>% mutate(CAR_pos=(CARrna_pos=="TRUE"|EGFRrna_pos=="TRUE"),
             Cell_type = case_when(
               CD4_pos=="TRUE" & CD8_pos=="FALSE" ~ "CD4",
               CD8_pos=="TRUE" & CD4_pos=="FALSE" ~ "CD8"),
             EGFR_pos=(EGFR_Status=="TRUE"| EGFRrna_pos=="TRUE")) %>% 
             mutate(all_pos=(EGFR_Status=="TRUE" & CARrna_pos=="TRUE" & EGFR_pos=="TRUE"),
                    CD4_CAR=(CAR_pos=="TRUE" & CD4_pos=="TRUE" ),
                    CD8_CAR=(CAR_pos=="TRUE" & CD8_pos=="TRUE" ))
signature_counts <- all_info  %>%
            group_by(orig.ident) %>% 
            summarise(n_total = n(), n_CD4_pos = sum(CD4_pos),n_CD8_pos = sum(CD8_pos),
            n_CARrna_pos = sum(CARrna_pos),
            n_EGFRrna_pos = sum(EGFRrna_pos),
            n_EGFRadt_pos = sum(EGFR_Status),
            n_EGFR_pos = sum(EGFR_pos),
            n_CAR_pos = sum(CAR_pos),
            n_CD4_CAR = sum(CD4_CAR),
            n_CD8_CAR = sum(CD8_CAR),
            n_all_pos = sum(all_pos),
            frac_CD4 = n_CD4_pos / n_total,
            frac_CD8 = n_CD8_pos / n_total,
            frac_CARrna = n_CARrna_pos / n_total,
            frac_EGFRrna = n_EGFRrna_pos / n_total,
            frac_EGFRadt = n_EGFRadt_pos / n_total,
            frac_EGFR = n_EGFR_pos / n_total,
            frac_CAR = n_CAR_pos / n_total,
            frac_CD4_CAR = n_CD4_CAR / n_total,
            frac_CD8_CAR = n_CD8_CAR / n_total,
            frac_all = n_all_pos / n_total)
write.table(signature_counts,"CART signature cutoff of each samples plan B.txt",sep="\t",row.names = F,quote=F)
#get counts overlap venndiagram
CARrnapos <- all_info %>% select(Identity,Cell,CARrna_pos) %>% filter(CARrna_pos=="TRUE")
EGFRabpos <- all_info %>% select(Identity,Cell,EGFR_Status)%>% filter(EGFR_Status=="TRUE")
EGFRrnapos <- all_info %>% select(Identity,Cell,EGFRrna_pos)%>% filter(EGFRrna_pos=="TRUE")
find_overlap <- function(sample) {
  idCARrnapos <- CARrnapos %>% filter(Identity == sample)
  idEGFRabpos <- EGFRabpos %>% filter(Identity == sample)
  idEGFRrnapos <- EGFRrnapos %>% filter(Identity == sample)
  venn.diagram(
    x=list(idEGFRabpos$Cell,idEGFRrnapos$Cell,idCARrnapos$Cell),
    category.names = c("EGFR_AB_pos","EGFR_pos","CARrna_pos"),
    fill = c("skyblue", "red", "lightgreen"),
    print.mode = c("raw","percent"),
    filename = paste0(sample,"_CART_signature_overlap.png"), 
    height = 2000 ,
    width = 2200 ,
    resolution = 300,
    alpha = 0.5,
    cat.cex = 1.2,
    cex = 2
  ) 
}
samples <- c("UPN787-B","UPN787-P","UPN775-B","UPN775-P","UPN749-B","UPN749-P","UPN685-B","UPN685-P","UPN746-B","UPN746-P")
for (i in 1:length(samples)) {find_overlap(samples[i])}
#match cell order and add to metadata
all_info <- all_info[match(rownames(Tcells@meta.data), all_info$Cell), ]
identical(rownames(Tcells@meta.data),all_info$Cell)
Tcells@meta.data$CAR_positive <- all_info$CAR_pos
Tcells@meta.data$Cell_type <- all_info$Cell_type
Tcells$TCR_positive <- !is.na(Tcells$cloneSize)
Tcells$patient <- sub("-.*", "", Tcells$orig.ident)
Tcells$disease <- ifelse(Idents(Tcells) %in% c("UPN749-P","UPN749-B","UPN746-B","UPN746-P"), "PD", "CR")
Tcells$subgroup <- paste(Tcells$group, Tcells$disease, sep = "_")

options(future.globals.maxSize = 100 * 1024^3)
integrated <- IntegrateLayers(object = Tcells, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony",verbose = FALSE)
integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:20)
integrated <- FindClusters(integrated, resolution = 0.2, cluster.name = "integrated.cluster")
integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")
DimPlot(integrated,group.by='integrated.cluster',reduction="umap.harmony")
ggsave("UMAP of harmony integrated T cells.png", dpi=300,width=7, height=6)
integrated <- JoinLayers(integrated)
var_genes <- Trex::quietTCRgenes(VariableFeatures(integrated,assay = "RNA"))
var_genes <- var_genes[!grepl("TRAV|TRBV|TRDV|TRGV|TRGC|TRAJ|TRAC|TRDC|ENSG|LINC|-AS1",var_genes)]
cluster.markers <- FindAllMarkers(object = integrated, features = var_genes, only.pos = T, min.pct = 0.2, logfc.threshold = 0.1)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
cells_use <- WhichCells(integrated, downsample = 300)
DoHeatmap(integrated, cells = cells_use, features = top10$gene,group.by="integrated.cluster",disp.max = 4,raster = TRUE,slot = "scale.data") + scale_fill_gradientn(colors = c("#89cff0", "white", "red"))
ggsave("Heatmap of top 10 variable genes in each Harmony integrated cluster in Tcells resolution 0.2.png", dpi=300,width=12, height=18)
saveRDS(integrated,file = "Harmony integrated Tcells with CR info.rds")
integrated <- projectil_analysis(integrated,ref_CD8,"CD8") 
integrated <- projectil_analysis(integrated,ref_CD4,"CD4") 
DimPlot(integrated, reduction = "umap.harmony", group.by = c("orig.ident", "integrated.cluster","subgroup","disease","Cell_type","Phase"))
ggsave(paste("Tcell","CART","integration cluster umap.png"), dpi=300, width=15, height=10)
features = c("TCF7","LEF1","FOXO1","IL7R","CCR7","SELL","BCL2","GZMB","PRF1","IFNG","TOX","NR4A1","NR4A3","PDCD1","HAVCR2","MKI67","TOP2A",
             "LAG3","TIGIT","GNLY","CTSW","SH2D1A","CXCR4","SLAMF6","TBX21","GATA3","ICOS","CD44")
DotPlot(integrated, features = features) + RotatedAxis() + theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste("dot plot after all sample","Tcell", "harmony", "integration.png"), dpi=300,width=9, height=7,bg="white")
#cell signatures for index score
integrated <- AddModuleScore(object=integrated,seed.use = 123, 
                             features = list(prolif_sig,Treg_sig,Tex_sig,CTL_sig,Tem_sig,Temra_sig,Tpex_sig,Tcm_sig,Tnaive_sig,gdT_sig,glyco_sig,OXPHOS_sig,Cytokine1_sig,Cytokine2_sig),
                             name=c("Prolif_score","Treg_score","Tex_score","CTL_score","Tem_score","Temra_score","Tpex_score","Tcm_score","Tn_score","gdT_score","Glyco_score","OXPHOS_score","Cytokine1_socre","Cytokine2_score"))
integrated@meta.data <- integrated@meta.data %>% dplyr::rename(Prolif_score = Prolif_score1,Treg_score = Treg_score2,Tex_score = Tex_score3, CTL_score = CTL_score4,Tem_score = Tem_score5, 
                                                               Temra_score = Temra_score6, Tpex_score = Tpex_score7, Tcm_score = Tcm_score8, Tn_score = Tn_score9, gdT_score = gdT_score10, 
                                                               Glyco_score = Glyco_score11, OXPHOS_score = OXPHOS_score12, Cytokine1_socre=Cytokine1_socre13,Cytokine2_score = Cytokine2_score14)
scores=c("Prolif_score","Treg_score","Tex_score","CTL_score","Tem_score","Temra_score","Tpex_score","Tcm_score","Tn_score","gdT_score","Glyco_score","OXPHOS_score","Cytokine1_socre","Cytokine2_score")
DotPlot(integrated, features = scores,group = "integrated.cluster") + RotatedAxis() 
ggsave("dot plot of index scores for each cluster.png", dpi=300,width=8, height=7,bg="white")
new.cluster.ids <- c("0" = "CD4 Tn", "1" = "Tprol-S", "2" = "Tem-like" ,"3" = "Tprol-M","4" = "Tpex","5" = "Temra","6" = "Tnk", "7" ="Treg")
annotation <- new.cluster.ids[as.character(integrated$integrated.cluster)]
names(annotation) <- colnames(integrated)
integrated <- AddMetaData(object = integrated,metadata = annotation,col.name = 'annotation')
#make dot plot of marker genes
genes <- c("CD4","IL2RA","GATA3","ICOS","MX1","IFIT3","SELL","JUN","BACH2","GZMB","CD8A","CD8B","LEF1","PRF1","FOS","FOSB","CXCR3","BATF3","CCR7","IL7R","TCF7","FASLG",
           "ITGA1", "MCM2","MCM4","MCM7","PCNA","MYB","MCM10","TYMS","STMN1","MKI67","TOP2A","CDK1","HMGB2","HMGB3","HMGN2","BIRC5","EZH2","CCNB2","CCNB1",
           "PLK1","TNF","TOX2", "IFI6","GZMA","CXCR4","CD27","RUNX3","CCL5","KLF2","CD44","GZMK","NKG7","HAVCR2","LAG3","SLAMF6","ZEB2","CCL4","FAS","GNLY","GZMH","TIGIT",
           "FGFBP2","CTLA4","IFNG","KLRG1","CX3CR1","CD83","ADGRG1","TBX21","KLRD1","FOXO1","BCL11B","TOX","RORA","EOMES","FOXP3")
integrated <- readRDS("Tcells harmony integrated CART cells planB.rds")
dp <- DotPlot(integrated, features = genes, group.by = "annotation")
# Step 2: find which cluster each gene is highest in
gene_order <- dp$data %>% group_by(features.plot) %>% slice_max(avg.exp.scaled, n = 1, with_ties = FALSE) %>%
  ungroup() %>% mutate(id = factor(id, levels = c("CD4 Tn","Tem-like","Temra","Tnk","Tpex","Tprol-M","Tprol-S","Treg"))) %>%
  arrange(id, desc(avg.exp.scaled)) %>% pull(features.plot) %>% as.character()
# Step 3: replot with ordered genes
DotPlot(integrated, features = gene_order, group.by = "annotation") + 
  coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
  axis.text = element_text(face = "bold", color = "black",size = 14),
  axis.title = element_text(face = "bold", color = "black", size = 14),
  legend.text = element_text(face = "bold", color = "black", size = 14),
  legend.title = element_text(face = "bold", color = "black", size = 14)) +
  scale_color_gradient2(low = "blue", mid = "lightgrey", high = "red") 
ggsave("Tcell signature genes dotplot for each cluster.svg",dpi = 300, width = 5, height = 15, bg = "white")

Idents(integrated) <- "annotation"
cluster.markers <- FindAllMarkers(object = integrated,only.pos = F, min.pct = 0, logfc.threshold = 0,return.thresh = 1)
entrez=read.table("/net/nfs-irwrsrchnas01/labs/sforman/Seq/CD19CART_singlecell/T_cell_integration/human genes Ensembl115.txt",header = T,sep="\t",stringsAsFactors = F) %>% dplyr::select(Gene.name,NCBI.gene..formerly.Entrezgene..ID) %>% filter(NCBI.gene..formerly.Entrezgene..ID!="NA") %>% distinct()
gsea_analysis <- function(markerlist, group, category, gs, qcutoff, width, height) {
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
  ggsave(paste(group,"dot plot of GSEA",category,"pathways.svg"),dpi=300,width=width, height = height)
}
gs_hallmark <- geneset::getMsigdb(org = "human",category = "H")
gs_kegg <- geneset::getMsigdb(org = "human",category = "C2-CP-KEGG")
gsea_analysis(cluster.markers ,"product","Hallmark",gs_hallmark,0.05,7,9)
gsea_analysis(cluster.markers ,"product","KEGG",gs_kegg,0.05,7,8)

#integration function
data_integration <- function(integrated, name, method, reduction) {
  # ref_ <- readRDS("/net/isi-dcnl/ifs/user_data/lkwak/Data/Jibin/reference/SingleR_Ref/CD8T_human_ref_v1.rds")
  # ref <- as.SingleCellExperiment(ref_, assay = "RNA")
  # pred <- SingleR(test = as.SingleCellExperiment(integrated),ref = ref,labels = ref$PT.annot )
  # integrated$labels <- pred$labels
  # integrated$pruned.labels<- pred$pruned.labels
  # DimPlot(integrated,group.by ="labels",reduction = "umap")
  # ggsave("ProjectTILs reference singleR annotation plan B.png",dpi=300,width=6,height=4)
  DimPlot(integrated, reduction = "umap.harmony",pt.size = 1, ncol = 2, label = F, group.by = "annotation", split.by = "orig.ident")
  ggsave(paste(name,reduction,"integration sample split umap.png"), dpi=300, width=7, height=14)
  Idents(integrated) <- "annotation"
  DimPlot(integrated, reduction = "umap.harmony",pt.size = 1, ncol = 2, label = F, split.by = "Cell_type")
  ggsave(paste(name,reduction,"integration cell type split umap plan B.png"), dpi=300, width=8, height=4)
  DimPlot(integrated, reduction = "umap.harmony",pt.size = 1, ncol = 2, label = F, split.by = "Phase")
  ggsave(paste(name,reduction,"integration phase split umap plan B.png"), dpi=300, width=7, height=8)
  #get the same colors in umap for barplot
  p <- DimPlot(integrated, reduction = "umap.harmony", group.by = "annotation")
  umap_colors <- unique(ggplot_build(p)$data[[1]]$colour)
  cluster_ids <- unique(p[[1]][["data"]][["annotation"]])
  cluster_colors <- setNames(umap_colors, cluster_ids)
  cluster_colors <- cluster_colors[sort(names(cluster_colors))]
  pt <- as.data.frame(table(integrated$annotation, integrated$subgroup,integrated$orig.ident))
  pt$Var1 <- as.character(pt$Var1)
  plot_group_portion(pt,"Product",name, reduction,cluster_colors)
  plot_group_portion(pt,"Blood",name, reduction,cluster_colors)
  plot_sample_portion(pt,"685",name, reduction,cluster_colors)
  plot_sample_portion(pt,"749",name, reduction,cluster_colors)
  plot_sample_portion(pt,"775",name, reduction,cluster_colors)
  plot_sample_portion(pt,"787",name, reduction,cluster_colors)
  plot_sample_portion(pt,"746",name, reduction,cluster_colors)
  colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)
  DimPlot(integrated, reduction = "umap.harmony",group.by = "cloneSize") + scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))
  ggsave(paste(name, reduction, "integrated all samples UMAP CR plan B.png"), dpi=300, width=8, height=5)
  VlnPlot(integrated, c("clonalFrequency"),group.by = "orig.ident")
  ggsave(paste("clonal frequency in each sample after",name, reduction,"integration plan B.png"),dpi=300,width=7,height=7)
  VlnPlot(integrated, c("percent.mt"),group.by = "orig.ident")
  ggsave(paste("mitochondrial gene expression in each sample after",name, reduction, "integration.png"),dpi=300,width=7,height=7)
  VlnPlot(integrated, c("percent.mt"),group.by = "integrated.cluster")
  ggsave(paste("mitochondrial gene expression in each cluster after",name, reduction, "integration.png"),dpi=300,width=7,height=7)
  features = c("TCF7","LEF1","FOXO1","IL7R","CCR7","SELL","BCL2","GZMB","PRF1","IFNG","TOX","NR4A1","NR4A3","PDCD1","HAVCR2","MKI67","TOP2A",
             "LAG3","TIGIT","GNLY","CTSW","SH2D1A","CXCR4","SLAMF6","TBX21","GATA3","ICOS","CD44")
  FeaturePlot(integrated, features = features, reduction = "umap.harmony",ncol=4)
  ggsave(paste(name, reduction, "integrated all samples UMAP cell marker genes plan B.png"), dpi=300, width=12, height=18)
  integrated = CellCycleScoring(object = integrated, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes)
  DimPlot(integrated,group.by='Phase',reduction="umap.harmony")
  ggsave(paste("CellCycle_FeaturePlot_after",name, reduction, "integration plan B.png"),dpi=300,width=5,height=5)
  clonalOccupy(integrated, x.axis = "ident", proportion = TRUE, label = FALSE)
  ggsave(paste("Clonal frequency of each",name, reduction, "integrated cluster.png"),dpi=300,width=9,height=5)
  VlnPlot(integrated, c("CAR-BAFFR-1","CAR-BAFFR-2", "WPRE601","EGFR","adt_EGFR-antibody"),group.by = "subgroup",ncol=2)
  ggsave(paste("CAR features in each group after",name, reduction, "integration.png"),dpi=300,width=7,height=7)
  VlnPlot(integrated, c("CAR-BAFFR-1","CAR-BAFFR-2", "WPRE601","EGFR","adt_EGFR-antibody"),group.by = "integrated.cluster",ncol=2)
  ggsave(paste("CAR features in each cluster after",name, reduction, "integration.png"),dpi=300,width=7,height=7)
  features <- c("CD4","CD8A","CD8B","CD3D","CD3E","CD3G","CD19","CAR-BAFFR-1","CAR-BAFFR-2","WPRE601","EGFR", "adt_EGFR-antibody")
  FeaturePlot(integrated, features = features, reduction = "umap.harmony",ncol=3)
  ggsave(paste("CAR features after",name, reduction, "integration.png"),dpi=300,width=8,height=9)
  FeaturePlot(integrated, features = features, reduction = "umap.harmony",split.by = "group")
  ggsave(paste("split CAR features after",name, reduction, "integration.png"),dpi=300,width=6,height=18)
  CARpos_barplot(integrated$Cell_type,integrated$annotation,name, reduction, "Celltype","cluster")
  CARpos_barplot(integrated$Phase,integrated$annotation,name, reduction, "Cellcycle","cluster")
  CARpos_barplot(integrated$cloneSize,integrated$annotation,name, reduction, "Clonesize","cluster")
  CARpos_barplot(integrated$Cell_type,integrated$subgroup,name, reduction, "Celltype","group")
  CARpos_barplot(integrated$Cell_type,integrated$orig.ident,name, reduction, "Celltype","sample")
  CARpos_barplot(integrated$Phase,integrated$orig.ident,name, reduction, "Cellcycle","sample")
  CARpos_barplot(integrated$Phase,integrated$subgroup,name, reduction, "Cellcycle","group")
  CARpos_barplot(integrated$cloneSize,integrated$orig.ident,name, reduction, "Clonesize","sample")
  CARpos_barplot(integrated$cloneSize,integrated$subgroup,name, reduction, "Clonesize","group")
  CARpos_barplot(integrated$TCR_positive,integrated$orig.ident,name, reduction, "TCR_positive","sample")
  CARpos_barplot(integrated$TCR_positive,integrated$subgroup,name, reduction, "TCR_positive","group")
  CARpos_barplot(integrated$TCR_positive,integrated$annotation,name, reduction, "TCR_positive","cluster")
  CARpos_barplot(harmony_inte$CAR_positive,integrated$annotation,"Tcell", "umap.harmony", "CAR_positive","cluster")
  CARpos_barplot(harmony_inte$CAR_positive,integrated$Cell_type,"Tcell", "umap.harmony", "CAR_positive","celltype")
  CARpos_barplot(harmony_inte$CAR_positive,integrated$subgroup,"Tcell", "umap.harmony", "CAR_positive","subgroup")
  CARpos_barplot(harmony_inte$CAR_positive,integrated$Phase,"Tcell", "umap.harmony", "CAR_positive","Cellcycle")
  
  #plot signature scores
  Idents(integrated) <- "orig.ident"
  product <- subset(integrated, idents = c("UPN787-P","UPN775-P","UPN685-P","UPN749-P","UPN746-P"))
  blood <- subset(integrated, idents = c("UPN787-B","UPN775-B","UPN685-B","UPN749-B","UPN746-B"))
  product$orig.ident <- factor(product$orig.ident,levels=c("UPN787-P","UPN775-P","UPN685-P","UPN749-P","UPN746-P"))
  blood$orig.ident <- factor(blood$orig.ident,levels=c("UPN787-B","UPN775-B","UPN685-B","UPN749-B","UPN746-B"))
  DotPlot(product, features = scores,group = "orig.ident") + RotatedAxis() 
  ggsave("product dot plot of index scores.png", dpi=300,width=8, height=4,bg="white")
  VlnPlot(product, scores,group.by = "disease",pt.size = 0)
  ggsave("product violine plot of index scores per group.png", dpi=300,width=8, height=8,bg="white")
  VlnPlot(product, scores,group.by = "orig.ident",pt.size = 0)
  ggsave("product violine plot of index scores per sample.png", dpi=300,width=10, height=10,bg="white")
  DotPlot(blood, features = scores,group = "orig.ident") + RotatedAxis() 
  ggsave("blood dot plot of index scores.png", dpi=300,width=8, height=4,bg="white")
  VlnPlot(blood, scores,group.by = "disease",pt.size = 0)
  ggsave("blood violine plot of index scores per group.png", dpi=300,width=8, height=8,bg="white")
  VlnPlot(blood, scores,group.by = "orig.ident",pt.size = 0)
  ggsave("blood violine plot of index scores per sample.png", dpi=300,width=10, height=10,bg="white")
  
  Idents(product) <- "annotation"
  Idents(blood) <- "annotation"
  DimPlot(product, reduction = "umap.harmony",pt.size = 1, ncol = 2, label = F, split.by = "subgroup")
  ggsave(paste(name,reduction,"integration product subgroup split umap plan B.png"), dpi=300, width=8, height=4)
  DimPlot(blood, reduction = "umap.harmony",pt.size = 1, ncol = 2, label = F, split.by = "subgroup")
  ggsave(paste(name,reduction,"integration blood subgroup split umap plan B.png"), dpi=300, width=8, height=4)
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
  
clonalOverlay(product, reduction = "umap.harmony",  cutpoint = 1, bins = 10, facet.by = "orig.ident")
ggsave(paste("clonalfrequency overlay in product samples in", name, reduction, "integration plan B.png"),dpi=300,width=8,height=8)
clonalOverlay(blood, reduction = "umap.harmony",  cutpoint = 1, bins = 10, facet.by = "orig.ident")
ggsave(paste("clonalfrequency overlay in blood samples in", name, reduction, "integration plan B.png"),dpi=300,width=8,height=8)

CARpos_barplot(product$Phase,product$subgroup,name, reduction, "Cellcycle","product group")
CARpos_barplot(blood$Phase,blood$subgroup,name, reduction, "Cellcycle","blood group")

StartracDiversity(integrated, type = "group", group.by = "patient")
ggsave(paste0("StartracDiversity in each", reduction, "integrated cluster.png"), dpi=300,width=12, height=9)
saveRDS(integrated, file = paste(name, reduction, "integrated cells planB.rds"))
return(integrated)
}
Tcell_inte <- data_integration(integrated,"Tcells", HarmonyIntegration,"harmony")

cells_keep <- colnames(integrated) 
combined.CR.filtered <- lapply(combined.CR, function(df) {df[df$barcode %in% cells_keep, ]})
compare_top_clone(combined.CR.filtered, "UPN685",10,"Tcell")
compare_top_clone(combined.CR.filtered, "UPN775",10,"Tcell")
compare_top_clone(combined.CR.filtered, "UPN749",10,"Tcell")
compare_top_clone(combined.CR.filtered, "UPN787",10,"Tcell")
compare_top_clone(combined.CR.filtered, "UPN746",10,"Tcell")
clonalProportion(combined.CR.filtered, cloneCall = "aa") 
ggsave(paste("clonal proportion Product and Blood Tcell.png"), dpi=300, width=10, height=6)
clonalHomeostasis(combined.CR.filtered, cloneCall = "aa")
ggsave(paste("clonal homeostasis Product and Blood Tcell.png"), dpi=300, width=10, height=6)

compare_top_clone <- function(cr_list,patient, top_number, name) {
  sampP <- paste0(patient, "-P")
  sampB <- paste0(patient, "-B")
  clonalCompare(cr_list, top.clones = top_number,samples = c(sampP, sampB), order.by = c(sampP, sampB), cloneCall="aa", graph = "alluvial",palette = "viridis")
  ggsave(paste(name,"Top10 CDR3 aa clonal comparison between",patient,"Product and Blood plan B.png"), dpi=300, width=6, height=6)
  clonalScatter(cr_list, cloneCall ="gene", x.axis = sampP, y.axis = sampB,dot.size = "total",graph = "proportion")
  ggsave(paste(name,"aa clonal scatter plot for",patient,"Product and Blood plan B.png"), dpi=300, width=6, height=6)
  vizGenes(combined.CR.filtered[c(sampP, sampB)],x.axis = "TRBV",y.axis = NULL,plot = "barplot") 
  ggsave(paste(name,"TRBV vizGene plot for",patient,"Product and Blood plan B.png"), dpi=300, width=10, height=6)
}

integrated <- readRDS("Tcells harmony integrated CART cells planB.rds")
product <- subset(integrated, idents = c("UPN787-P","UPN775-P","UPN685-P","UPN749-P","UPN746-P"))
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
make_dimplot_with_pecent(product,"product","CAR_positive","annotation", width=6,height=3,ncols=2)

#select samples in this project
Idents(integrated) <- "orig.ident"
#select CART cells for downstream
CART <- subset(integrated, subset = CAR_positive=="TRUE")
cells_keep <- colnames(CART) 
combined.CR.CART <- lapply(combined.CR, function(df) {df[df$barcode %in% cells_keep, ]})
cdr3_stats <- map_df(names(combined.CR.CART),
                     ~ combined.CR.filtered[[.x]] %>% mutate(TCR=paste0(TCR1,"_",TCR2), cdr3_aa=paste0(cdr3_aa1,"_",cdr3_aa2),cdr3_nt=paste0(cdr3_nt1,"_",cdr3_nt2)) %>%
                       summarise(sample = .x,cdr3a_NAs     = sum(is.na(cdr3_aa1)),cdr3a_nonNA  = sum(!is.na(cdr3_aa1)),
                                 cdr3b_NA     = sum(is.na(cdr3_aa2)),cdr3b_nonNA  = sum(!is.na(cdr3_aa2)),
                                 cdr3_aa1_n = n_distinct(cdr3_aa1, na.rm = TRUE),cdr3_aa2_n = n_distinct(cdr3_aa2, na.rm = TRUE),
                                 cdr3_nt1_n = n_distinct(cdr3_nt1, na.rm = TRUE),cdr3_nt2_n = n_distinct(cdr3_nt2, na.rm = TRUE),
                                 TCR1_n = n_distinct(TCR1, na.rm = TRUE),TCR2_n = n_distinct(TCR2, na.rm = TRUE),
                                 cdr3_aa_n = n_distinct(cdr3_aa, na.rm = TRUE),cdr3_nt_n = n_distinct(cdr3_nt, na.rm = TRUE),TCR_n = n_distinct(TCR, na.rm = TRUE)))
write.table(cdr3_stats,"number of clones in each sample CART.csv",sep=",",row.names = F,quote=F)
CART_inte <- data_integration(CART,"CART", HarmonyIntegration,"harmony")

compare_top_clone(combined.CR.CART, "UPN685",10,"CART")
compare_top_clone(combined.CR.CART, "UPN775",10,"CART")
compare_top_clone(combined.CR.CART, "UPN749",10,"CART")
compare_top_clone(combined.CR.CART, "UPN787",10,"CART")
compare_top_clone(combined.CR.CART, "UPN746",10,"CART")

StartracDiversity(harmony_inte, type = "group", group.by = "patient")
ggsave("StartracDiversity in each Harmony integrated cluster.png", dpi=300,width=12, height=9)
clonalBias(harmony_inte, cloneCall = "aa", split.by = "patient", group.by = "integrated.cluster",n.boots = 10, min.expand =5)
ggsave("clonalBias in each Harmony integrated cluster.png", dpi=300,width=12, height=9)
clonalNetwork(harmony_inte, reduction = "umap", group.by = "integrated.cluster",filter.clones = NULL,filter.identity = NULL,cloneCall = "aa")

top5_clones_per_sample <- harmony_inte@meta.data %>% 
  filter(!is.na(CTaa)) %>%                  # keep only cells with a TCR sequence
  dplyr::count(orig.ident, CTaa, name = "n") %>%   # count how many cells per clone per sample
  group_by(orig.ident) %>%                 
  slice_max(order_by = n, n = 5, with_ties = FALSE) %>%  # top 5 per sample
  ungroup()
seq_to_highlight <- unique(top5_clones_per_sample$CTaa)
harmony_inte <- highlightClones(harmony_inte,cloneCall = "aa",sequence  = seq_to_highlight)

samples <- unique(harmony_inte$patient)
for (s in samples) {
  sub_harm <- subset(harmony_inte,subset= patient==s)
  DimPlot(sub_harm,reduction = "umap",group.by  = "highlight",split.by  = "orig.ident")
  ggsave(filename = paste0("UMAP_top5_clones_", s, ".png"), width = 10,height= 4,dpi= 300)
  #seq_to_patient <- top5_clones_per_sample %>% filter(grepl(s,orig.ident))
  #alluvialClones(sub_harm, cloneCall = "aa", y.axes = c("orig.ident","integrated.cluster","Cell_type"), color = unique(seq_to_patient$CTaa))
  #ggsave(filename = paste0("alluvialClones_top5_clones_", s, ".png"), width = 6,height= 4,dpi= 300)
}

#separate CD4 and CD8 cells
CD4 <- subset(integrated,subset=CD4>0)
CD4 <- data_integration(CD4, "CD4", HarmonyIntegration,"harmony")
CD8 <- subset(integrated, subset = CD8A > 0 | CD8B > 0)
CD8 <- data_integration(CD8, "CD8", HarmonyIntegration,"harmony")

# Match sample order
DE_analyse <- function(seu, substring, cells, group ) {
  if (cells %in% c("CD4","CD8")) {
    sub <- subset(seu, subset=Cell_type==cells)
  } else {
    sub <- subset(seu, subset=annotation == cells)
  }
  pb <- make_pseudobulk(sub, group_var = "orig.ident") 
  pb <- pb[-c(1:3),grepl(substring,colnames(pb))]
  meta <- sub@meta.data %>% select(orig.ident, disease) %>% distinct()
  meta <- meta[match(colnames(pb), meta$orig.ident), ] 
  rownames(meta) <- meta$orig.ident
  dds <- DESeqDataSetFromMatrix(countData = pb,colData = meta, design = ~ disease )
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("disease","CR","PD"))
  res <- as.data.frame(res[order(res$padj),])
  write.table(res,paste("DEseq2 DEGs with pseudo bulk of",group, cells, "cells.csv"),sep=",",quote=F)
  vsd <- vst(dds, blind = F)
  write.table(assay(vsd),paste("vsd normalized counts of",group, cells,"cells.csv"), sep=",")
  res$diffexpressed <- "NO"
  res$diffexpressed[res$log2FoldChange > 1 & res$padj < 0.05] <- "UP"
  res$diffexpressed[res$log2FoldChange < -1 & res$padj < 0.05] <- "DOWN"
  res$delabel <- NA
  res$delabel[res$diffexpressed != "NO"] <- rownames(res[res$diffexpressed != "NO",])
  ggplot(data=res, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
    geom_point(size=0.01) + 
    theme(legend.position = "none") +
    ggtitle(paste("CR vs PD DEGs", group, cells, "T cells")) +
    geom_text_repel(size=2,max.overlaps = 5,segment.size = 0.3, seed=40, segment.alpha = 0.5, segment.length = 0.2, segment.color = "grey50", force = 1, min.segment.length = 0.1,box.padding = 0.3,point.padding = 0.1) +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-1, 1), col="red", size=0.3, alpha=0.5) +
    geom_hline(yintercept=-log10(0.05), col="red", size=0.3, alpha=0.5)
  ggsave(paste("volcano plot of DEGs with pseudo bulk of", group, cells, "T cells.pdf"), width=6, height=5, dpi=300)
}

DE_analyse(integrated,"-B","CD4","blood")
DE_analyse(integrated,"-B","CD8","blood")
DE_analyse(integrated,"-P","CD4","product")
DE_analyse(integrated,"-P","CD8","product")
anno <- unique(integrated$annotation)
for (i in c(1:length(unique(integrated$integrated.cluster)))) {
  DE_analyse(integrated,"-P",anno[i],"product")
  DE_analyse(integrated,"-B",anno[i],"blood")
}

harmony_inte <- readRDS("harmony integrated samples.rds")
harmony_inte$integrated.cluster <- factor(
  as.numeric(as.character(harmony_inte$integrated.cluster)),
  levels = sort(unique(as.numeric(as.character(harmony_inte$integrated.cluster))))
)
