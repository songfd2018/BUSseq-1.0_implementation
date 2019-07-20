# Apply Seurat 3 to the pancreatic study.
# Referring to https://satijalab.org/seurat/v3.0/pancreas_integration_label_transfer.html
rm(list=ls())
library(Seurat)
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
library(ggplot2)
library(cowplot)

set.seed(123)

# coloring
library(scales)
require(WGCNA)
library(RColorBrewer)

###########################
# Load Hematopoietic Data #
###########################
# Working directory
setwd("F:/scRNA/code/0601Cpp_BUSseq/BUSseq_implementation_v1/HumanPancreas/Comparison/Seurat")


# Loading hematopoietic count data
countdata <- read.table("../../RawCountData/count_data_pancreas_v1.txt") 

dim <- unlist(read.table("../../RawCountData/dim_pancreas_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

metadata <- read.table("../../RawCountData/pancreas_metadata.txt")

load("../../RawCountData/pancreas_countdata.RData")
gene_list <- rownames(PancreasCounts[[1]])
cell_list <- c(colnames(PancreasCounts[[1]]),colnames(PancreasCounts[[2]]),colnames(PancreasCounts[[3]]),colnames(PancreasCounts[[4]]))
##########################################
# Apply Seurat to the Hematopoietic Data #
##########################################
# Setting the Seurat Object
colnames(countdata) <- cell_list
rownames(countdata) <- gene_list
rownames(metadata) <- cell_list

pancreas <- CreateSeuratObject(countdata, meta.data = metadata)
pancreas.list <- SplitObject(pancreas, split.by = "Study")


for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2480, 
                                             verbose = FALSE)
}

reference.list <- pancreas.list[c("GSE81076", "GSE85241","GSE86473","E-MTAB-5061")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, anchor.features = 3470)

pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

DefaultAssay(pancreas.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)

# pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
# 
# p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "Study")
# p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "Celltype", label = TRUE, repel = TRUE) + 
#   NoLegend()
# plot_grid(p1, p2)

pancreas.integrated <- FindNeighbors(pancreas.integrated, reduction = "pca", dims = 1:30, k.param = 20)
pancreas.integrated<-FindClusters(pancreas.integrated,resolution = 1.5)

# pancreas.integrated <- RunTSNE(object = pancreas.integrated)
# DimPlot(object = pancreas.integrated, reduction = "tsne")

w_Seurat <- pancreas.integrated$seurat_clusters
reorder <- match(rownames(metadata),names(w_Seurat))
w_Seurat<- w_Seurat[reorder]

adjustedRandIndex(as.numeric(w_Seurat),as.numeric(metadata$CellType))

Seurat_PCA <- pancreas.integrated@reductions$pca@cell.embeddings[reorder,]

Seurat_corrected <- NULL
for(g in 1:G){
  Seurat_corrected <- rbind(Seurat_corrected, pancreas.integrated@assays$integrated[g,])
}

#######
# ARI #
#######
ARI_Seurat <- adjustedRandIndex(metadata$CellType,w_Seurat)

#######################################################################
# Silhouette coefficient by estimated cell type labels of each method #
#######################################################################
# Silhouette coefficients of the top 30 PCs inferred by Seurat
Seurat_dist <- dist(Seurat_PCA)
sil_Seurat <- silhouette(as.integer(w_Seurat), dist = Seurat_dist)
sil_Seurat_true <- silhouette(as.integer(metadata$CellType), dist = Seurat_dist)

#################################
# scatter plots (t-SNE and PCA) #
#################################
if(!dir.exists("Image")){
  dir.create("Image")
}

allcolors<-labels2colors(metadata$CellType)#[allsamples])

allcolors[allcolors=="red"]<-"deeppink"
allcolors[allcolors=="yellow"]<-"orange1"#"darkgoldenrod1"

plot_by_celltype<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  # The image with legend
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=20, 
       col=alpha(allcolors[subset],0.6),
       xlab=xlab, ylab=ylab, main=main) 
  dev.off()
}

colors4<-brewer.pal(4,"Set3")
batch.cols<-colors4[rep(1:B,nb)]
plot_by_batch<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=20, 
       col=alpha(batch.cols[subset],0.6),
       xlab=xlab, ylab=ylab, main=main)#,  xlab="tSNE 1",ylab="tSNE 2")
  dev.off()
}

#################################################
# Draw t-SNE plots by batch and true cell types #
#################################################
set.seed(123)
all.dists.Seurat <- as.matrix(Seurat_dist)
tsne_Seurat_dist <- Rtsne(all.dists.Seurat, is_distance=TRUE, perplexity = 30)

Seurat_by_celltype<- "Image/tsne_pancreas_Seurat_by_celltype.jpeg"
plot_by_celltype(Seurat_by_celltype, tsne_Seurat_dist$Y)

Seurat_by_batch <- "Image/tsne_pancreas_Seurat_by_batch.jpeg"
plot_by_batch(Seurat_by_batch, tsne_Seurat_dist$Y)

##################
# Draw PCA plots #
##################
#pca.Seurat <- prcomp(Seurat_PCA, rank=2)
Seu_PCA<- "Image/pca_pancreas_Seurat.jpeg"
plot_by_celltype(Seu_PCA, Seurat_PCA[,1:2], xlab = "PC 1", ylab = "PC 2")

# legend
pdf(file="Image/legend.pdf", width=10, height=8)
plot(c(-5,5),c(-4,4),type="n", bty="n", axes=FALSE, xlab="", ylab="")
legend(x=-5, y=4, legend=c("GSE81076","GSE85241","GSE86473","E-MTAB-5061"), pch=20, cex=2.5, col=alpha(colors4,0.6), title = "Batch", bty="n")
legend(x=-0, y=4, legend=names(table(metadata$CellType)), pch=20, cex=2.5, col=alpha(c("turquoise","blue","brown","orange1","green","deeppink","black"),0.6), title = "Cell Type", bty="n")
dev.off()

# Store the workspace
save.image("Seurat_workspace.RData")


