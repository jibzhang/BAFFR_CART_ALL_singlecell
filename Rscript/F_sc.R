library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(SummarizedExperiment)

# BiocManager::install("affy")

# install.packages("remotes")
# BiocManager::install("SingleR")

# read and write tsv ------------------------
write_tsv=function(indata,outfile){
  write.table(indata,outfile,col.names = T,row.names = F,quote = F,na = "",sep="\t")
}

read_tsv=function(infile){
  read.table(infile,header = T,sep = "\t",stringsAsFactors = F)
}

#get ENSG id -----------------------------------------------------------

hkgenes_human=readLines("/home/jibzhang/pipeline/Singlecell/HouseKeepingGenes.txt")
# 
# #file_all_geneid="/ref_genomes/genome_anno/human/gtf/V102/Homo_sapiens.GRCh38.V102.summary.txt"
# #geneid_all=read.table(file_all_geneid,header = T,stringsAsFactors = F,sep="\t")
# 
# file_geneid="/ref_genomes/genome/human/GRCh38/cellranger/expression/GRCh38/genes/genes.gtf.summary.txt"
# df_geneid=read.table(file_geneid,header = T,stringsAsFactors = F)
# 
# get_ENSGID=function(genes){x1=df_geneid[df_geneid$gene_name %in% genes,];x1$gene_id}
# 
# s_ENSG=get_ENSGID(Seurat::cc.genes$s.genes)
# g2m_ENSG=get_ENSGID(Seurat::cc.genes$g2m.genes)
# MT_ENSG=get_ENSGID(df_geneid$gene_name[grepl("^MT-",df_geneid$gene_name)])
# Ribosomal_ENSG=get_ENSGID(df_geneid$gene_name[grepl("^RP[SL]",df_geneid$gene_name)])
# housekeeping_ENSG=get_ENSGID(hkgenes_human)


#create_seurat_object-------------------------------------------
create_seurat_object=function(data_dir,
                              project_name="Seurat_obj",
                              min.cells=3,
                              min.features=200,
                              gene.column=2){
  print("reading data")
  data = Read10X(data.dir = data_dir,gene.column=gene.column)
  print("create object")
  seurat_object=CreateSeuratObject(counts = data, project = project_name, 
                                   min.cells = min.cells, min.features = min.features)
  return(seurat_object)
}


#Basic -------------------------------------------------------------------------
get_barcode_df=function(object_in){
  data.frame(barcode=row.names(object_in@meta.data),stringsAsFactors = F)
}

get_metadataVar=function(inObj,var){
  df_out=data.frame(barcode=row.names(inObj@meta.data),
                    var=unlist(inObj@meta.data[var]),
                    stringsAsFactors = F)
  names(df_out)[2]=var
  df_out
}

get_embedding_df=function(inobj,reduction_method){
  df_embed = as.data.frame(Embeddings(inobj[[reduction_method]]))
  df_embed$barcode=row.names(df_embed)
  df_embed
}

draw_feature_plot_fromDF=function(objIn,df,varname,label,cols_in,dir_out){
  
  objIn[[label]]=unlist(df[varname])
  
  DimPlot(objIn,group.by = label,label=T,repel = T,
          cols = cols_in[names(cols_in) %in% unlist(objIn@meta.data[label])])
  ggsave(paste0(dir_out,label,".png"),width=11,height = 9)
  ggsave(paste0(dir_out,label,".pdf"),width=11,height = 9)
  
  objIn
}

#For QC ------------------------
add_qc_measurements=function(inObj,housekeepinngGene=hkgenes_human){
  gene_list=rownames(x = inObj)
  MTgene=gene_list[grepl("^MT-",toupper(gene_list))]
  ribosomalGene=gene_list[grepl("^RP[SL]",toupper(gene_list))]
  housekeepinngGene1=gene_list[gene_list %in% housekeepinngGene]
  #add qc measurement
  inObj[['percent.ribo']] = PercentageFeatureSet(inObj,  features = ribosomalGene)
  inObj[['percent.mt']] = PercentageFeatureSet(inObj, features = MTgene)
  inObj[['percent.houseKeeping']] = PercentageFeatureSet(inObj, features = housekeepinngGene1)
  inObj
}

get_cutoff=function(inObj,
                    var_list=c("nCount_RNA","nFeature_RNA","percent.ribo","percent.mt")){
  bind_rows(lapply(var_list, function(var_in){
    x=unlist(inObj@meta.data[var_in])
    df_cutoff=data.frame(
      var=var_in,
      percentile_1=quantile(x,0.01),
      percentile_5=quantile(x,0.05),
      percentile_10=quantile(x,0.10),
      percentile_90=quantile(x,0.90),
      percentile_95=quantile(x,0.95),
      percentile_99=quantile(x,0.99),
      
      median_=median(x),
      mad_=mad(x),
      
      medianAdd3MAD=median(x)+3*mad(x),
      medianMinus3MAD=median(x)-3*mad(x),
      stringsAsFactors = F
    )
    df_cutoff
  }))
}

doubletFinder_new <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE, annotations = NULL) {
  require(Seurat); require(fields); require(KernSmooth)
  
  ## Generate new list of doublet classificatons from existing pANN vector to save time
  if (reuse.pANN != FALSE ) {
    pANN.old <- seu@meta.data[ , reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    return(seu)
  }
  
  if (reuse.pANN == FALSE) {
    ## Make merged real-artifical data
    real.cells <- rownames(seu@meta.data)
    # data <- seu@assays$RNA$counts[, real.cells]
    
    data=Seurat::GetAssayData(seu, layer="counts")[,real.cells]
    
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating",n_doublets,"artificial doublets...",sep=" "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    # Keep track of the types of the simulated doublets
    if(!is.null(annotations)){
      stopifnot(typeof(annotations)=="character")
      stopifnot(length(annotations)==length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    ## Store important pre-processing information
    orig.commands <- seu@commands
    
    ## Pre-process Seurat object
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      
      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets,
                                     normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                                     margin = orig.commands$NormalizeData.RNA@params$margin)
      
      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
                                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
                                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
                                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
                                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
                                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
                                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
                                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
                                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
                                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
                                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
      
      print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets,
                                 features = orig.commands$ScaleData.RNA$features,
                                 model.use = orig.commands$ScaleData.RNA$model.use,
                                 do.scale = orig.commands$ScaleData.RNA$do.scale,
                                 do.center = orig.commands$ScaleData.RNA$do.center,
                                 scale.max = orig.commands$ScaleData.RNA$scale.max,
                                 block.size = orig.commands$ScaleData.RNA$block.size,
                                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
      
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets,
                              features = orig.commands$ScaleData.RNA$features,
                              npcs = length(PCs),
                              rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                              verbose=FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc() # Free up memory
    }
    
    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      
      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)
      
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc()
    }
    
    ## Compute PC distance matrix
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    
    ## Compute pANN
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    if(!is.null(annotations)){
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = length(levels(doublet_types1))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
      if(!is.null(annotations)){
        for(ct in unique(annotations)){
          neighbors_that_are_doublets = neighbors[neighbors>n_real.cells]
          if(length(neighbors_that_are_doublets) > 0){
            neighbor_types[i,] <-
              table( doublet_types1[neighbors_that_are_doublets - n_real.cells] ) +
              table( doublet_types2[neighbors_that_are_doublets - n_real.cells] )
            neighbor_types[i,] <- neighbor_types[i,] / sum( neighbor_types[i,] )
          } else {
            neighbor_types[i,] <- NA
          }
        }
      }
    }
    print("Classifying doublets..")
    classifications <- rep("Singlet",n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN",pN,pK,nExp,sep="_")] <- pANN[rownames(seu@meta.data), 1]
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    if(!is.null(annotations)){
      colnames(neighbor_types) = levels(doublet_types1)
      for(ct in levels(doublet_types1)){
        seu@meta.data[, paste("DF.doublet.contributors",pN,pK,nExp,ct,sep="_")] <- neighbor_types[,ct]
      }
    }
    return(seu)
  }
}

get_doublet_DoubletFinder=function(inObj, sct){
  #get singlet doublet 
  annotations=inObj@meta.data$seurat_clusters
  homotypic.prop <- DoubletFinder::modelHomotypic(annotations)  
  nExp_poi <- round(0.075*nrow(inObj@meta.data))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  #get optimal pk value
  sweep.res <- paramSweep(inObj, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  optimal_pk <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  optimal_pk_value <- as.numeric(as.character(optimal_pk))
  #find doublets
  inObj <- doubletFinder_new(seu=inObj, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj,  sct = T)
  inObj
}

#SC analysis -----------------------------------
SCTransform2cluster=function(objIn,resolution=1,n_pc=50,n_variableGenes=2000,n.neighbors=30){
  print("x1")
  objIn=SCTransform(objIn,variable.features.n = n_variableGenes)
  print("x2")
  objIn = RunPCA(objIn,npcs =n_pc ,verbose = FALSE)
  print("x3")
  objIn = RunUMAP(objIn, dims = 1:n_pc, verbose = FALSE,n.neighbors=n.neighbors)
  print("x4")
  objIn = FindNeighbors(objIn, dims = 1:n_pc, verbose = FALSE)
  print("x5")
  objIn = FindClusters(objIn, verbose = FALSE,resolution = resolution)
  print("x6")
  objIn
}

draw_DimPlot=function(inObj,df_feature,var_key,label,file_out_prefix){
  if (!all(row.names(inObj@meta.data)==unlist(df_feature[1]))){stop("barcode not equal")}
  
  inObj[[label]]=unlist(df_feature[var_key])
  DimPlot(inObj,group.by = label,label=T,repel = T,cols=DiscretePalette(length(unique(unlist(df_feature[var_key])))))
  
  ggsave(paste0(file_out_prefix,".png"),width=10,height = 10)
  ggsave(paste0(file_out_prefix,".pdf"),width=10,height = 10)
  inObj
}


#SingleR Annotation -------------------------

# inObj=obj_QC
# ref_singleR=ref_BALL
# var_group="diag"

create_singleR_ref=function(exp,pheno,var_group){
  pheno1=pheno[var_group]
  
  if(!all(colnames(exp)==row.names(pheno1))){stop("Exp colname not equal to row names of Pheno")}
  
  objOut = SummarizedExperiment(list(logcounts=exp), colData = pheno1)
  objOut
}

# inObj=obj_QC
# ref_singleR=ref_BALL
# var_group="diag"

run_singleR=function(inObj,ref_singleR,var_group,top_variable_features_n=1000,method="classic"){
  #using top variable features
  top_variable_features=head(VariableFeatures(inObj), top_variable_features_n)
  for_singleR=GetAssayData(inObj)
  for_singleR=for_singleR[top_variable_features,]
  
  #using all the features
  #for_singleR=GetAssayData(inObj[["SCT"]], slot = "data")
  
  print(dim(for_singleR))
  
  value_group=unlist(ref_singleR@colData[var_group])
  
  #cell level annotation
  print("cell level annotation")
  out_singleR = SingleR::SingleR(test = for_singleR, ref = ref_singleR,labels = value_group)
  
  #cluster level annotation
  print("cluster level annotation")
  out_singleR_cluster = SingleR::SingleR(test = for_singleR, ref = ref_singleR,labels = value_group,clusters=inObj$seurat_clusters)
  
  #get output 
  df_singleR_out=data.frame(
    barcode=row.names(inObj@meta.data),
    seurat_clusters=inObj$seurat_clusters,
    labels=out_singleR$labels,
    pruned.labels=out_singleR$pruned.labels,
    stringsAsFactors = F) %>%
    left_join(data.frame(
      seurat_clusters=sort(unique(inObj$seurat_clusters)),
      labels_cluster=out_singleR_cluster$labels,
      pruned.labels_cluster=out_singleR_cluster$pruned.labels,
      stringsAsFactors = F
    )
    )
  return(list(out_singleR=out_singleR,out_singleR_cluster=out_singleR_cluster,df_singleR_out=df_singleR_out))
}

run_singleR_SCT=function(inObj,ref_singleR,var_group,top_variable_features_n=1000,method="classic"){
  #using top variable features
  top_variable_features=head(VariableFeatures(inObj), top_variable_features_n)
  for_singleR=GetAssayData(inObj)
  for_singleR=for_singleR[top_variable_features,]
  
  #using all the features
  #for_singleR=GetAssayData(inObj[["SCT"]], slot = "data")

  print(dim(for_singleR))
  
  value_group=unlist(ref_singleR@colData[var_group])
  
  #cell level annotation
  print("cell level annotation")
  out_singleR = SingleR::SingleR(test = for_singleR, ref = ref_singleR,labels = value_group)
  
  #cluster level annotation
  print("cluster level annotation")
  out_singleR_cluster = SingleR::SingleR(test = for_singleR, ref = ref_singleR,labels = value_group,clusters=inObj$seurat_clusters)

  #get output 
  df_singleR_out=data.frame(
    barcode=row.names(inObj@meta.data),
    seurat_clusters=inObj$seurat_clusters,
    labels=out_singleR$labels,
    pruned.labels=out_singleR$pruned.labels,
    stringsAsFactors = F) %>%
    left_join(data.frame(
      seurat_clusters=sort(unique(inObj$seurat_clusters)),
      labels_cluster=out_singleR_cluster$labels,
      pruned.labels_cluster=out_singleR_cluster$pruned.labels,
      stringsAsFactors = F
    )
    )
  return(list(out_singleR=out_singleR,out_singleR_cluster=out_singleR_cluster,df_singleR_out=df_singleR_out))
}

# inObj=obj_QC
# ref_singleR1=ref_HPCA
# ref_singleR2=ref_BlineageRef_30types
# var_group1="label.final"
# var_group2="celltype"
make_clonal_compare <- function(cr_list, patient) {
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
  ggsave(paste("Top10 CDR3 aa clonal comparison between",patient,"Product and Blood.png"), dpi=300, width=6, height=6)
}

keep_top_pct_per_sample <- function(object, feature, sample_col = "sample",
                                    assay = "ADT", top_pct = 0.73) {
  dat <- FetchData(object, vars = c(sample_col, feature), assay = assay)
  colnames(dat) <- c(sample_col, "marker")
  qprob <- 1 - top_pct
  # per-sample logical keep mask
  keep <- ave(dat$marker, dat[[sample_col]],
              FUN = function(v) v >= quantile(v, probs = qprob, na.rm = TRUE)) & is.finite(dat$marker)
  subset(object, cells = rownames(dat)[keep])
}
# # usage:
# obj_top73 <- keep_top_pct_per_sample(obj, feature = "CD3", sample_col = "sample", assay = "ADT", top_pct = 0.73)

first_valley_cut <- function(x, adjust = 1) {
  d <- density(x, adjust = adjust, na.rm = TRUE)
  y <- d$y
  mins <- which(diff(sign(diff(y))) > 0) + 1L
  imode <- which.max(y)
  mins_right <- mins[mins > imode]
  if (length(mins_right) == 0) return(NA_real_)
  d$x[mins_right[1]]
}

find_first_event_after_major_left_peak <- function(x, y,
                                                   spar = 0.65,
                                                   z_thresh = 3,minx=0.6,
                                                   peak_frac = 0.40,
                                                   prom_frac = 0.10,
                                                   min_prom_drop = 0.05) {
  fit <- smooth.spline(x, y, spar = spar)
  ys <- fit$y
  xs <- fit$x
  n  <- length(ys)
  
  is_peak <- function(v) {
    p <- rep(FALSE, length(v))
    if (length(v) >= 3) {
      p[2:(length(v)-1)] <- (v[1:(length(v)-2)] < v[2:(length(v)-1)]) &
        (v[2:(length(v)-1)] > v[3:length(v)])
    }
    p
  }
  peaks_idx <- which(is_peak(ys))
  if (length(peaks_idx) == 0)
    return(data.frame(x=NA, y=NA, type=NA, main_peak_x=NA, main_peak_y=NA))
  
  nearest_min_left  <- function(i, v) {
    j <- as.integer(i)
    while (j > 1L && v[j-1L] <= v[j]) j <- j - 1L
    j
  }
  nearest_min_right <- function(i, v) {
    j <- as.integer(i)
    n <- length(v)
    while (j < n && v[j+1L] <= v[j]) j <- j + 1L
    j
  }
  
  # Either keep vapply with integer FUN.VALUE...
  left_min_idx  <- vapply(peaks_idx,  nearest_min_left,  integer(1), v = ys)
  right_min_idx <- vapply(peaks_idx,  nearest_min_right, integer(1), v = ys)
  
  peak_y    <- ys[peaks_idx]
  globalmax <- max(ys, na.rm = TRUE)
  
  left_min_y    <- ys[left_min_idx]
  right_min_y   <- ys[right_min_idx]
  prom          <- peak_y - pmax(left_min_y, right_min_y)
  
  keep <- (peak_y >= peak_frac * globalmax) & (prom >= prom_frac * globalmax)
  major_candidates <- peaks_idx[keep]
  if (length(major_candidates) == 0) {
    major_candidates <- peaks_idx
  }
  pk_idx <- major_candidates[1]
  pk_x   <- x[pk_idx]; pk_y <- ys[pk_idx]
  
  dy  <- c(NA, diff(ys) / diff(x))
  ddy <- c(NA, diff(dy))
  z   <- (ddy - stats::median(ddy, na.rm=TRUE)) / stats::mad(ddy, na.rm=TRUE)
  
  search_idx <- seq.int(min(pk_idx + 1L, n-1L), n-1L)
  dropped_enough <- ys[search_idx] <= pk_y * (1 - min_prom_drop)
  right_enough <- xs[search_idx] > minx
  
  cand <- search_idx[ which(abs(z[search_idx]) >= z_thresh & dropped_enough & right_enough) ]
  if (length(cand) == 0)
    return(data.frame(x=NA, y=NA, type=NA, main_peak_x=pk_x, main_peak_y=pk_y))
  
  i <- cand[1]
  type <- ifelse(z[i] > 0, "valley", "shoulder")
  data.frame(x = x[i], y = ys[i], type = type, main_peak_x = pk_x, main_peak_y = pk_y)
}

pick_valley_then_shoulder <- function(raw_x,
                                      d_bw = "nrd0", d_n = 512,
                                      valley_adjust = 1,
                                      valley_threshold = 1,
                                      spar = 0.7, z_thresh = 3, minx = 0.6,
                                      peak_frac = 0.40, prom_frac = 0.10,
                                      min_prom_drop = 0.07) {
  # Step 1: try valley
  val <- tryCatch(first_valley_cut(raw_x, adjust = valley_adjust),
                  error = function(e) NA_real_)
  
  d <- density(raw_x, adjust = valley_adjust, n = d_n, bw = d_bw, na.rm = TRUE)
  
  if (!is.na(val) && val < valley_threshold) {
    y <- d$y
    y_at_val <- approx(d$x, y, xout = val)$y
    imode <- which.max(y)
    dx <- diff(d$x)[1]
    area_above <- sum(d$y[d$x > val]) * dx
    return(data.frame(
      event_x = val,
      event_y = y_at_val,
      type = "valley",
      main_peak_x = d$x[imode],
      main_peak_y = y[imode],
      percentage = round(100 * area_above, 2)
    ))
  }
  
  # Step 2: fallback to shoulder
  res <- find_first_event_after_major_left_peak(
    d$x, d$y,
    spar = spar, z_thresh = z_thresh,
    peak_frac = peak_frac, prom_frac = prom_frac,
    minx = minx, min_prom_drop = min_prom_drop
  )
  dx <- diff(d$x)[1]
  area_above <- sum(d$y[d$x > res$x]) * dx
  res$percentage <- round(100 * area_above, 2)
  
  data.frame(event_x = res$x, event_y = res$y, type = res$type,
             main_peak_x = res$main_peak_x, main_peak_y = res$main_peak_y,
             percentage = res$percentage)
}

run_singleR_doubleRef=function(inObj,ref_singleR1,ref_singleR2,var_group1,var_group2){
  #using all the features
  for_singleR=GetAssayData(inObj[["SCT"]], slot = "data")
  
  print(dim(for_singleR))
  
  value_group1=unlist(ref_singleR1@colData[var_group1])
  value_group2=unlist(ref_singleR2@colData[var_group2])
  
  #cell level annotation
  print("cell level annotation")
  out_singleR = SingleR::SingleR(test = for_singleR, 
                                 ref = list(ref1=ref_singleR1,ref2=ref_singleR2),
                                 labels = list(value_group1,value_group2))
  
  #cluster level annotation
  print("cluster level annotation")
  out_singleR_cluster = SingleR::SingleR(test = for_singleR, 
                                         ref = list(ref1=ref_singleR1,ref2=ref_singleR2),
                                         labels = list(value_group1,value_group2),
                                         clusters=inObj$seurat_clusters)
  
  #get output 
  df_singleR_out=data.frame(
    barcode=row.names(inObj@meta.data),
    seurat_clusters=inObj$seurat_clusters,
    labels=out_singleR$labels,
    pruned.labels=out_singleR$pruned.labels,
    stringsAsFactors = F) %>%
    left_join(data.frame(
      seurat_clusters=sort(unique(inObj$seurat_clusters)),
      labels_cluster=out_singleR_cluster$labels,
      pruned.labels_cluster=out_singleR_cluster$pruned.labels,
      stringsAsFactors = F
    )
    )
  return(list(out_singleR=out_singleR,out_singleR_cluster=out_singleR_cluster,df_singleR_out=df_singleR_out))
}

draw_annotation_plot=function(objIn,singleR_obj,sample,label,cols_in,dir_out){
  
  objIn[[label]]=singleR_obj$labels
  objIn[[sample]]=singleR_obj$samples
  DimPlot(objIn,group.by = label,label=T,repel = T, split.by = sample,
          cols = cols_in[names(cols_in) %in% unlist(objIn@meta.data[label])])
  ggsave(paste0(dir_out,label,"split_dimplot.png"),width=11,height = 10)
  ggsave(paste0(dir_out,label,"split_dimplot.pdf"),width=11,height = 10)
  
  objIn
}

# in_singleR_obj=singleR_anno_$out_singleR
# cutoff_precentile=0.7
get_stringent_label=function(in_singleR_obj,cutoff_precentile){
  df_score=as.data.frame(in_singleR_obj$scores)
  
  df_score$barcode=in_singleR_obj@rownames
  
  df_score$pruned.labels=in_singleR_obj$pruned.labels
  
  df_id=df_score["barcode"]
  
  df_score1=df_score %>% reshape2::melt(id.vars=c("barcode","pruned.labels"))  %>%
    group_by(barcode) %>% 
    mutate(median = median(value, na.rm = TRUE),
           delta=value-median) %>%
    filter(delta>0) %>%
    group_by(variable) %>%
    mutate(
      delta_cutoff=quantile(delta,cutoff_precentile)
    ) %>% group_by(barcode) %>%
    filter(value==max(value)) %>%
    mutate(
      label_stringent=ifelse(delta>delta_cutoff,as.character(pruned.labels),NA)
    )
  
  df_id %>% left_join(df_score1)
}

draw_barplot=function(x,file_out_prefix,width=5,height=10){
  x[is.na(x)]="NA";  n=length(x)
  
  df_freq=data.frame(table(x)) %>% arrange(desc(Freq)) %>%
    mutate(Group=as.character(x),
           ratio=Freq/n,
           Percentage=paste0(Group,": ",Freq," (",sprintf("%.2f",ratio*100),"% )"))
  
  ggplot(df_freq,aes(x=ratio,y=Percentage)) + geom_bar(stat="identity") 
  ggsave(paste0(file_out_prefix,".png"),width=width,height = height)
  ggsave(paste0(file_out_prefix,".pdf"),width=width,height = height)
}

#KNN prediction -------------------------

# library(caret)
# indata=df_pca %>% sample_n(500)
# 
# var_x=names(indata)[grepl("PC",names(indata))]
# var_y="pruned.labels"
# KNN_K=5

knn_pred_one=function(indata,var_y,var_x,KNN_K){
  formula_in=formula(paste0(var_y,"~",paste0(var_x,collapse = "+")))
  
  knnFit = train(formula_in, data = indata, method = "knn",
                 preProcess = c("center", "scale"),
                 tuneGrid = expand.grid(k = c(KNN_K)))
  knn_pred=predict(knnFit) %>% as.vector()
  knn_pred
}


#BALL color -------------------------
BALLCol_singleR=c()
{
  BALLCol_singleR["BCL2MYC"]="seagreen2"
  BALLCol_singleR["DUX4"]='grey40'
  BALLCol_singleR["ETV6RUNX1"]="gold2"
  BALLCol_singleR["HLF"]= "skyblue"
  BALLCol_singleR["Hyperdiploid"]="#3E9F32"
  BALLCol_singleR["iAMP21"]="lightslateblue"
  BALLCol_singleR["IKZF1N159Y"]="#CCCC33"
  BALLCol_singleR["KMT2A"]="#1F78B5"
  BALLCol_singleR["Lowhypodiploid"]="#1E90FF"
  BALLCol_singleR["MEF2D"]="#66C2A6"
  BALLCol_singleR["NUTM1"]='black'
  BALLCol_singleR["PAX5ETV6"]="#808000"
  BALLCol_singleR["PAX5P80R"]="orangered"
  BALLCol_singleR["PAX5alt"]="#FFA620"
  BALLCol_singleR["Ph"]="magenta3"
  BALLCol_singleR["TCF3PBX1"]="darkgoldenrod4"
  BALLCol_singleR["Y"]="#E6BEFF"
  BALLCol_singleR["ZEB2CEBPE"]="#D27B1C86"
  BALLCol_singleR["ZNF384"]="#A8DD00"
  
  BALLCol_singleR["Normal"]="red4"
}

subtypeCol=c()
{
  subtypeCol["ETV6-RUNX1"]="gold2"
  subtypeCol["ETV6-RUNX1-like"]="pink"
  subtypeCol["ETV6-RUNX1-sc"]="pink"
  subtypeCol["KMT2A"]="#1F78B5"
  
  subtypeCol["KMT2A_1"]="#1F78B5"
  subtypeCol["KMT2A_2"]="#FF00B6"
  subtypeCol["KMT2A_3"]="#00FF00"
  
  subtypeCol["Ph"]="magenta3"
  subtypeCol["Ph-like"]="red4"
  
  subtypeCol["Ph-major"]="magenta3"
  subtypeCol["Ph-minor"]="red4"
  
  subtypeCol["PAX5alt"]="#FFA620"
  subtypeCol["PAX5-ETV6"]="#808000"
  subtypeCol["PAX5(P80R)"]="orangered"
  subtypeCol["PAX5 P80R"]="orangered"
  
  subtypeCol["PAX5alt-major"]="#FFA620"
  subtypeCol["PAX5alt-minor"]="pink"
  
  subtypeCol["DUX4"]='grey40'
  subtypeCol["TCF3-PBX1"]="darkgoldenrod4"
  subtypeCol["ZNF384"]="#A8DD00"
  subtypeCol["MEF2D"]="#66C2A6"
  subtypeCol["BCL2/MYC"]="seagreen2"
  subtypeCol["NUTM1"]='black'
  subtypeCol["HLF"]= "skyblue"
  subtypeCol["Hyperdiploid"]="#3E9F32"
  subtypeCol["LowHypo"]="#1E90FF"
  subtypeCol["Low hypodiploid"]="#1E90FF"
  subtypeCol["NearHaploid"]='blue3'
  subtypeCol["Near haploid"]='blue3'
  subtypeCol["iAMP21"]="lightslateblue"
  subtypeCol["IKZF1(N159Y)"]="#CCCC33"
  subtypeCol["IKZF1 N159Y"]="#CCCC33"
  subtypeCol["LowHyper"]="cyan"
  subtypeCol["Bother"]='grey75'
  subtypeCol["Low hyperdiploid"]='grey75'
  subtypeCol["CRLF2(non-Ph-like)"]='grey75'
  subtypeCol["KMT2A-like"]='grey75'
  subtypeCol["ZNF384-like"]='grey75'
  subtypeCol["Other"]='grey75'
  subtypeCol["ZEB2/CEBPE"]="#D27B1C86"
  subtypeCol["Y"]="#E6BEFF"
  subtypeCol["CDX2"]="#E6BEFF"
  
  subtypeCol["Unknown"]="#469990"
  subtypeCol["Normal"]="red4"
  subtypeCol["Others"]="grey75"
  subtypeCol["_Prediction"]="red4"
  
  subtypeCol["AML"]="red4"
  subtypeCol["BALL"]="#FFA620"
  subtypeCol["TALL"]='blue3'
  
}

