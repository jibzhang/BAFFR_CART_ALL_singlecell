source("/home/jibzhang/pipeline/Singlecell/Rscript/F_sc.R")

library(Matrix)
library(Seurat)
library(DoubletFinder)

#get parameters
args <- commandArgs(trailingOnly = TRUE)
print(args)

dir_matrix=args[1]
id=args[2]
outdir=args[3]
adt_feature=args[4]
car_feature=args[5]
transform=args[6]

#make out directory
dir.create(outdir)
#setup seurat object and add in the HTO data
data = Read10X(data.dir = dir_matrix,gene.column=2)
rownames(data[["Antibody Capture"]]) <- gsub(pattern = "Anti_EGFR_barcoded_ab_totalseq_C0132_anti_human_", replacement = "",rownames(data[["Antibody Capture"]]))
rownames(data[["Antibody Capture"]]) <- gsub(pattern = "_TotalSeqC", replacement = "",rownames(data[["Antibody Capture"]]))
seurat_object <- CreateSeuratObject(counts = data[["Gene Expression"]], project = id, min.cells = 1, min.features = 200)
seurat_object = add_qc_measurements(seurat_object)
options(future.globals.maxSize = 100 * 1024^3)
if(transform=="SCT"){
  seurat_object=SCTransform(seurat_object,vars.to.regress = c("percent.mt","nCount_RNA"))
}

if(transform=="basic"){
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 3000)
  VariableFeatures(seurat_object) <- grep("TRAV|TRBV|TRDV|TRGV|TRGC|TRAJ|TRAC|TRDC|ENSG|LINC|-AS1|^MT-|^RP[SL]",VariableFeatures(seurat_object),value = TRUE, invert = TRUE)
  seurat_object <- ScaleData(seurat_object,vars.to.regress = c("percent.mt","nCount_RNA"))
}


seurat_Tcells <- subset(x = seurat_object, subset = CD3D>0 & CD3E>0 & CD3G>0)
cat("####################################################################################################")
cat(paste0("number of non Tcells removed: ",ncol(seurat_object)-ncol(seurat_Tcells)))
cat("####################################################################################################")
seurat_object <- seurat_Tcells 

VlnPlot(seurat_object,features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo"),ncol = 2)
ggsave(paste0(outdir,"/",id,"_QC_basic_raw_VlnPlot.png"),width=12,height = 12)
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  FeatureScatter(seurat_object, feature1 = "nFeature_RNA", feature2 = "percent.mt")+
  FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
ggsave(paste0(outdir,"/",id,"_QC_basic_raw_FeatureScatter.png"),width=15,height = 5)

seurat_object[["ADT"]] <- CreateAssayObject(data[["Antibody Capture"]][, colnames(seurat_object), drop = FALSE])
seurat_object <- NormalizeData(seurat_object, assay = "ADT", normalization.method = "CLR")
#FeatureScatter(seurat_object, feature1 = "adt_EGFR-antibody", feature2 = "EGFR", pt.size = 1) + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18))
FeatureScatter(seurat_object, feature1 = adt_feature, feature2 = "EGFR", pt.size = 1) + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18))
ggsave(paste0(outdir,"/",id,"_EGFR_scatterplot.png"),width=6,height = 5)
FeatureScatter(seurat_object, feature1 = adt_feature, feature2 = car_feature, pt.size = 1) + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18))
ggsave(paste0(outdir,"/",id,"_CARadt_scatterplot.png"),width=6,height = 5)
FeatureScatter(seurat_object, feature1 = car_feature, feature2 = "EGFR", pt.size = 1) + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18))
ggsave(paste0(outdir,"/",id,"_CAREGFR_scatterplot.png"),width=6,height = 5)

seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
ElbowPlot(seurat_object, ndims=50,reduction="pca")
ggsave(paste0(outdir,"/",id,"_elbowplot.png"),width=9,height = 5)
FeaturePlot(seurat_object,features=c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo"),ncol = 2)
ggsave(paste0(outdir,"/",id,"_QC_basic_raw_FeaturePlot.png"),width=12,height = 12)

seurat_object <- FindNeighbors(seurat_object, dims = 1:30)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
seurat_object <- RunUMAP(seurat_object, dims = 1:30, verbose = FALSE)
seurat_object = CellCycleScoring(object = seurat_object, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes)
#run double let
if(transform=="SCT"){
  seurat_object=get_doublet_DoubletFinder(seurat_object,sct=T) 
} else {
  seurat_object=get_doublet_DoubletFinder(seurat_object,sct=F) 
}
q
DimPlot(seurat_object,group.by=names(seurat_object@meta.data)[grepl("DF.classifications",names(seurat_object@meta.data))])+
  FeaturePlot(seurat_object,features = names(seurat_object@meta.data)[grepl("pANN",names(seurat_object@meta.data))])+plot_layout(ncol=2)
ggsave(paste0(outdir,"/",id,"_QC_DoubletFinder_raw.png"),width=12,height = 7)

#draw plot
FeaturePlot(seurat_object,features="S.Score")+ FeaturePlot(seurat_object,features="G2M.Score")+DimPlot(seurat_object,group.by='Phase')+plot_layout(ncol=2)
ggsave(paste0(outdir,"/",id,"_QC_CellCycle_raw.png"),width=12,height = 12)

p1 <- FeaturePlot(seurat_object, adt_feature, cols = c("lightgrey", "darkgreen")) + ggtitle("EGFR protein")
p2 <- FeaturePlot(seurat_object, "rna_EGFR") + ggtitle("EGFR RNA") + ggtitle("EGFR RNA")
p3 <- FeaturePlot(seurat_object, car_feature, cols = c("lightgrey", "red")) + ggtitle("CAR RNA")
p1 | p2 | p3
ggsave(paste0(outdir,"/",id,"_CAR_featureplot.png"),width=15,height =6)

VlnPlot(seurat_object, adt_feature)
ggsave(paste0(outdir,"/",id,"_EGFR_violinplot.png"),width=9,height = 4)

saveRDS(seurat_object,file=paste0(outdir,"/",id,".SeuratObj.raw.rds"))
