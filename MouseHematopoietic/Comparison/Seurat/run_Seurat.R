# Apply Seurat 3 to the hematopoietic study.
rm(list=ls())
library(Seurat)
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
set.seed(123)


# Working directory
# setwd("F:/scRNA/code/0601Cpp_BUSseq/BUSseq_implementation_v1/MouseHematopoietic/Comparison/Seurat")
#setwd("/scratch/data01/BUSseq_cpp/BUSseq_implementation_v1/MouseHematopoietic/Comparison/Seurat")
setwd("/your/working/directory/BUSseq_implementation-1.0/MouseHematopoietic/Comparison/Seurat/")

###################
# Load hemat Data #
###################
# Loading hematopoietic count data
countdata <- read.table("../../RawCountData/count_data_hemat_v1.txt") 

dim <- unlist(read.table("../../RawCountData/dim_hemat_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# Load metadata
metadata <- read.table("../../RawCountData/metadata_hemat_v1.txt")
rownames(metadata) <- paste0("Batch_",rep(c(1,2),nb),"_Cell_",c(1:nb[1],1:nb[2]))
gene_list <- unlist(read.table("../../RawCountData/gene_list_hemat_v1.txt",stringsAsFactors = FALSE))

colnames(countdata) <- paste0("Batch_",rep(c(1,2),nb),"_Cell_",c(1:nb[1],1:nb[2]))
rownames(countdata) <- gene_list

##########################################
# Apply Seurat to the Hematopoietic Data #
##########################################
# Referrring to https://satijalab.org/seurat/v3.0/pancreas_integration_label_transfer.html
# Setting the Seurat Object 
hemat <- CreateSeuratObject(countdata, meta.data = metadata)
hemat.list <- SplitObject(hemat, split.by = "Study")

for (i in 1:length(hemat.list)) {
  hemat.list[[i]] <- NormalizeData(hemat.list[[i]], verbose = FALSE)
  hemat.list[[i]] <- FindVariableFeatures(hemat.list[[i]], selection.method = "vst", nfeatures = 3470, 
                                             verbose = FALSE)
}

reference.list <- hemat.list[c("GSE72857", "GSE81682")]
hemat.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, anchor.features = 3470)

hemat.integrated <- IntegrateData(anchorset = hemat.anchors, dims = 1:30)

DefaultAssay(hemat.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
hemat.integrated <- ScaleData(hemat.integrated, verbose = FALSE)
hemat.integrated <- RunPCA(hemat.integrated, npcs = 30, verbose = FALSE)

hemat.integrated <- FindNeighbors(hemat.integrated, reduction = "pca", dims = 1:30, k.param = 20)
hemat.integrated<-FindClusters(hemat.integrated,resolution = 1.5)

w_Seurat <- hemat.integrated$seurat_clusters

Seurat_PCA <- hemat.integrated@reductions$pca@cell.embeddings
Seurat_corrected <- NULL
for(g in 1:G){
  Seurat_corrected <- rbind(Seurat_corrected, hemat.integrated@assays$integrated[g,])
}

#######
# ARI #
#######
ARI_Seurat <- adjustedRandIndex(metadata$CellType,w_Seurat)

#################################
# scatter plots (t-SNE and PCA) #
#################################
if(!dir.exists("Image")){
  dir.create("Image")
}
# coloring for cell types
color.legendF <- c(MEP="orange", GMP="chartreuse4", CMP="magenta", 
                   HSPC="cyan", LTHSC="dodgerblue", MPP="blue", LMPP="light blue", other="grey")
colmatF <- col2rgb(color.legendF) 
colmatA <- colmatF + 100 # A lighter shade.
colmatA[colmatA > 255] <- 255
color.legendA <- setNames(rgb(colmatA[1,], colmatA[2,], colmatA[3,], maxColorValue=255), names(color.legendF))
allcolors <- c(color.legendA[as.character(metadata$CellType[1:nb[1]])],color.legendF[as.character(metadata$CellType[1:nb[2]+nb[1]])])
first.batch <- rep(c(TRUE, FALSE), nb)

plot_by_celltype<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  # The image with legend
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=ifelse(first.batch, 1, 21)[subset], 
       col=ifelse(first.batch, allcolors, "black")[subset],
       bg=allcolors[subset], xlab=xlab, ylab=ylab, main=main) 
  dev.off()
}

batchcolor <- c("lightcoral","lavender")
batch_ind <- rep(1:2,nb)
plot_by_batch<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=ifelse(first.batch, 1, 21)[subset], 
       col=ifelse(first.batch, batchcolor[batch_ind[subset]], "black"),
       bg=batchcolor[batch_ind[subset]], xlab=xlab, ylab=ylab, main=main)#,  xlab="tSNE 1",ylab="tSNE 2")
  dev.off()
}

#################################################
# Draw t-SNE plots by batch and true cell types #
#################################################
set.seed(123)
Seurat_dist <- dist(Seurat_PCA)
all.dists.Seurat <- as.matrix(Seurat_dist)
tsne_Seurat_dist <- Rtsne(all.dists.Seurat, is_distance=TRUE, perplexity = 30)

Seurat_by_celltype<- "Image/tsne_hemat_Seurat_by_celltype.jpeg"
plot_by_celltype(Seurat_by_celltype, tsne_Seurat_dist$Y)

Seurat_by_batch <- "Image/tsne_hemat_Seurat_by_batch.jpeg"
plot_by_batch(Seurat_by_batch, tsne_Seurat_dist$Y)

##################
# Draw PCA plots #
##################
#pca.Seurat <- prcomp(Seurat_PCA, rank=2)
Seu_PCA<- "./Image/pca_hemat_Seurat.jpeg"
plot_by_celltype(Seu_PCA, Seurat_PCA[,1:2], xlab = "PC 1", ylab = "PC 2")


# Store the workspace
# save.image("Seurat_workspace.RData")
save(ARI_Seurat,tsne_Seurat_dist,file = "Seurat_results.RData")