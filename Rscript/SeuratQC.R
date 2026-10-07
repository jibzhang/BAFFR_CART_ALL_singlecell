source("/home/jibzhang/pipeline/Singlecell/Rscript/F_sc.R")

#get parameters
args <- commandArgs(trailingOnly = TRUE)
print(args)

rds_in=args[1]
doublet_remove = args[2]
rds_out=args[3]
transform=args[4]

#read rds
obj_raw=readRDS(rds_in)

#if(!dir.exists(dir_raw)){dir.create(dir_raw)}
#if(!dir.exists(dir_qc)){dir.create(dir_qc)}

# df_barcode=get_barcode_df(obj_raw)
# set.seed(10)
# df_barcode_r=df_barcode %>% sample_n(1000)
# obj_raw=subset(obj_raw,cells=df_barcode_r$barcode)
dir_raw=dirname(rds_in)
#remove doublets----------------------------------------------------------------
df_singlelet=get_metadataVar(obj_raw,names(obj_raw@meta.data)[grepl("DF.classifications",names(obj_raw@meta.data))])
df_singlelet=df_singlelet[df_singlelet[,2]=="Singlet",]
write_tsv(df_singlelet,paste0(dir_raw,"/df_singlelet.tsv"))

obj_afterRemovingDoublet=subset(obj_raw,cells=df_singlelet$barcode)
obj_afterRemovingDoublet=add_qc_measurements(obj_afterRemovingDoublet)

if(transform=="SCT"){
  obj_afterRemovingDoublet=SCTransform(obj_afterRemovingDoublet, vars.to.regress = c("percent.mt","nCount_RNA"))
}

if(transform=="basic"){
  obj_afterRemovingDoublet=NormalizeData(obj_afterRemovingDoublet)
  obj_afterRemovingDoublet=FindVariableFeatures(obj_afterRemovingDoublet, selection.method = "vst", nfeatures = 3000)
  VariableFeatures(obj_afterRemovingDoublet) <- grep("TRAV|TRBV|TRDV|TRGV|TRGC|TRAJ|TRAC|TRDC|ENSG|LINC|-AS1|^MT-|^RP[SL]",VariableFeatures(obj_afterRemovingDoublet),value = TRUE, invert = TRUE)
  obj_afterRemovingDoublet=ScaleData(obj_afterRemovingDoublet,vars.to.regress = c("percent.mt","nCount_RNA"))
}

obj_afterRemovingDoublet <- RunPCA(obj_afterRemovingDoublet, verbose = FALSE)
obj_afterRemovingDoublet <- RunUMAP(obj_afterRemovingDoublet, dims = 1:30, verbose = FALSE)
obj_afterRemovingDoublet <- FindNeighbors(obj_afterRemovingDoublet, dims = 1:30, verbose = FALSE)
obj_afterRemovingDoublet <- FindClusters(obj_afterRemovingDoublet, resolution = 0.3,verbose = FALSE)

# obj_afterRemovingDoublet=get_doublet_DoubletFinder(obj_afterRemovingDoublet)
obj_afterRemovingDoublet = CellCycleScoring(object = obj_afterRemovingDoublet, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes)

FeaturePlot(obj_afterRemovingDoublet,features=c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo"),ncol = 2)
ggsave(paste0(dir_raw,"/QC_basic_afterRemovingDoublet_FeaturePlot.png"),width=12,height = 12)
ggsave(paste0(dir_raw,"/QC_basic_afterRemovingDoublet_FeaturePlot.pdf"),width=12,height = 12)

VlnPlot(obj_afterRemovingDoublet,features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo"),ncol = 2)
ggsave(paste0(dir_raw,"/QC_basic_afterRemovingDoublet_VlnPlot.png"),width=12,height = 12)
ggsave(paste0(dir_raw,"/QC_basic_afterRemovingDoublet_VlnPlot.pdf"),width=12,height = 12)

FeatureScatter(obj_afterRemovingDoublet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  FeatureScatter(obj_afterRemovingDoublet, feature1 = "nFeature_RNA", feature2 = "percent.mt")+
  FeatureScatter(obj_afterRemovingDoublet, feature1 = "nCount_RNA", feature2 = "percent.mt")
ggsave(paste0(dir_raw,"/QC_basic_afterRemovingDoublet_FeatureScatter.png"),width=15,height = 5)
ggsave(paste0(dir_raw,"/QC_basic_afterRemovingDoublet_FeatureScatter.pdf"),width=15,height = 5)

FeaturePlot(obj_afterRemovingDoublet,features="S.Score")+
  FeaturePlot(obj_afterRemovingDoublet,features="G2M.Score")+
  DimPlot(obj_afterRemovingDoublet,group.by='Phase')+plot_layout(ncol=2)
ggsave(paste0(dir_raw,"/QC_CellCycle_afterRemovingDoublet_FeaturePlot.png"),width=12,height = 12)
ggsave(paste0(dir_raw,"/QC_CellCycle_afterRemovingDoublet_FeaturePlot.pdf"),width=12,height = 12)

saveRDS(obj_afterRemovingDoublet,file=paste0(dir_raw,"/",doublet_remove))
#make out directory
dir_qc <- paste0(dir_raw,"/plotQC")
dir.create(dir_qc)
#QC ----------------------------------------------------------------
df_cutoff=get_cutoff(obj_afterRemovingDoublet)
write_tsv(df_cutoff,paste0(dir_qc,"df_cutoff.tsv"))

obj_QC1=subset(obj_afterRemovingDoublet, subset = percent.mt < 10)
cat("####################################################################################################")
cat(paste0("number of cells removed because of dead cells: ",ncol(obj_afterRemovingDoublet)-ncol(obj_QC1)))
cat("####################################################################################################")

obj_QC=subset(obj_QC1, subset = nFeature_RNA > 200 & nCount_RNA < df_cutoff$medianAdd3MAD[df_cutoff$var=="nCount_RNA"] )
cat("####################################################################################################")
cat(paste0("number of cells removed because of outlier cells of nFeature: ",ncol(obj_QC1)-ncol(obj_QC)))
cat("####################################################################################################")

obj_QC=add_qc_measurements(obj_QC)
options(future.globals.maxSize = 100 * 1024^3)
if(transform=="SCT"){
  obj_QC=SCTransform(obj_QC, vars.to.regress = c("percent.mt","nCount_RNA"))
}

if(transform=="basic"){
  obj_QC=NormalizeData(obj_QC,normalization.method = "LogNormalize", scale.factor = 1e4)
  obj_QC=FindVariableFeatures(obj_QC, selection.method = "vst", nfeatures = 3000)
  VariableFeatures(obj_QC) <- grep("TRAV|TRBV|TRDV|TRGV|TRGC|TRAJ|TRAC|TRDC|ENSG|LINC|-AS1|^MT-|^RP[SL]",VariableFeatures(obj_QC),value = TRUE, invert = TRUE)
  obj_QC=ScaleData(obj_QC,vars.to.regress = c("percent.mt","nCount_RNA"))
}

obj_QC <- RunPCA(obj_QC, verbose = FALSE)
obj_QC <- RunUMAP(obj_QC, dims = 1:30, verbose = FALSE)
obj_QC <- FindNeighbors(obj_QC, dims = 1:30, verbose = FALSE)
obj_QC <- FindClusters(obj_QC, resolution = 0.3,verbose = FALSE)

obj_QC = CellCycleScoring(object = obj_QC, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes)

embedding <- Embeddings(obj_QC, reduction = "umap")
metadata <- obj_QC@meta.data
cluster<-metadata %>% select(seurat_clusters) 
cluster$CellID <- rownames(cluster)

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
nclust<-length(unique(cluster$seurat_clusters))
color_list <- ggplotColours(n=nclust)
colour<- data.frame(cluster=0:(nclust-1),color_list)
clustercolor<-merge(cluster,colour,by.x = "seurat_clusters",by.y = 'cluster')
clustercolor<-clustercolor[c(2,1,3)]
colnames(clustercolor) <- c("CellID",'cluster','colour')

write.csv(clustercolor, file = paste0(dir_qc,"/clusters.csv"), row.names = F)
write.csv(embedding, file = paste0(dir_qc,"/cell_embeddings.csv"))
write.csv(Cells(obj_QC), file = paste0(dir_qc,"/cellID_obs.csv"), row.names = FALSE)
write.csv(as.matrix(t(GetAssayData(object = obj_QC, layer = "counts"))), paste0(dir_qc,'/raw_counts_after_qc.csv'), sep = ',', row.names = T, col.names = T, quote = F)

FeaturePlot(obj_QC,features=c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo"),ncol = 2)
ggsave(paste0(dir_qc,"/QC_basic_afterRemovingDoublet_FeaturePlot.png"),width=12,height = 12)
ggsave(paste0(dir_qc,"/QC_basic_afterRemovingDoublet_FeaturePlot.pdf"),width=12,height = 12)

VlnPlot(obj_QC,features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo"),ncol = 2)
ggsave(paste0(dir_qc,"/QC_basic_afterRemovingDoublet_VlnPlot.png"),width=12,height = 12)
ggsave(paste0(dir_qc,"/QC_basic_afterRemovingDoublet_VlnPlot.pdf"),width=12,height = 12)

FeatureScatter(obj_QC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  FeatureScatter(obj_QC, feature1 = "nFeature_RNA", feature2 = "percent.mt")+
  FeatureScatter(obj_QC, feature1 = "nCount_RNA", feature2 = "percent.mt")
ggsave(paste0(dir_qc,"/QC_basic_afterRemovingDoublet_FeatureScatter.png"),width=15,height = 5)
ggsave(paste0(dir_qc,"/QC_basic_afterRemovingDoublet_FeatureScatter.pdf"),width=15,height = 5)

FeaturePlot(obj_QC,features="S.Score")+
  FeaturePlot(obj_QC,features="G2M.Score")+
  DimPlot(obj_QC,group.by='Phase')+plot_layout(ncol=2)
ggsave(paste0(dir_qc,"/QC_CellCycle_afterRemovingDoublet_FeaturePlot.png"),width=12,height = 12)
ggsave(paste0(dir_qc,"/QC_CellCycle_afterRemovingDoublet_FeaturePlot.pdf"),width=12,height = 12)

DimPlot(obj_QC,label = TRUE, pt.size = 1,group.by = "seurat_clusters")
ggsave(paste0(dir_qc,"/DimPlot_SeuratCluster.png"),width=11,height = 10)

saveRDS(obj_QC,file=paste0(dir_qc,"/",rds_out))