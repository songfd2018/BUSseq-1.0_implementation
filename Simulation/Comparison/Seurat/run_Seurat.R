# Apply Seurat 3 to the simulated study.
# Referring to https://satijalab.org/seurat/v3.0/simulation_integration_label_transfer.html
rm(list=ls())
library(Seurat)
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
library(ggplot2)
library(cowplot)

# coloring
library(scales)
require(WGCNA)
library(RColorBrewer)

set.seed(123)
# Working directory
setwd("/your/working/directory/BUSseq_implementation-1.0/Simulation/Comparison/Seurat/")

########################
# Load simulated Data #
########################

# Loading hematopoietic count data
countdata <- read.table("../../RawCountData/count_data_simulation_v4.txt") 

dim <- unlist(read.table("../../RawCountData/dim_simulation_v4.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# Load metadata
metadata <- read.table("../../RawCountData/metadata_simulation_v4.txt")

# Load gene_list
gene_list <- paste0("gene-",1:G)

rownames(countdata) <- gene_list
colnames(countdata) <- rownames(metadata)

##########################################
# Apply Seurat to the Simulated Data #
##########################################
# Setting the Seurat Object
simulation <- CreateSeuratObject(countdata, meta.data = metadata)
simulation.list <- SplitObject(simulation, split.by = "batch")


for (i in 1:length(simulation.list)) {
  simulation.list[[i]] <- NormalizeData(simulation.list[[i]], verbose = FALSE)
  simulation.list[[i]] <- FindVariableFeatures(simulation.list[[i]], selection.method = "vst", nfeatures = G, 
                                             verbose = FALSE)
}

reference.list <- simulation.list[c("Batch_1", "Batch_2","Batch_3","Batch_4")]
simulation.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, anchor.features = G)

simulation.integrated <- IntegrateData(anchorset = simulation.anchors, dims = 1:30)

DefaultAssay(simulation.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
simulation.integrated <- ScaleData(simulation.integrated, verbose = FALSE)
simulation.integrated <- RunPCA(simulation.integrated, npcs = 30, verbose = FALSE)

# simulation.integrated <- RunUMAP(simulation.integrated, reduction = "pca", dims = 1:30)
# 
# p1 <- DimPlot(simulation.integrated, reduction = "umap", group.by = "Study")
# p2 <- DimPlot(simulation.integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + 
#   NoLegend()
# plot_grid(p1, p2)

simulation.integrated <- FindNeighbors(simulation.integrated, reduction = "pca", dims = 1:30, k.param = 20)
simulation.integrated<-FindClusters(simulation.integrated,resolution = 1.5)

# simulation.integrated <- RunTSNE(object = simulation.integrated)
# DimPlot(object = simulation.integrated, reduction = "tsne")

w_Seurat <- simulation.integrated$seurat_clusters
reorder <- match(rownames(metadata),names(w_Seurat))
w_Seurat<- w_Seurat[reorder]

Seurat_PCA <- simulation.integrated@reductions$pca@cell.embeddings[reorder,]

Seurat_corrected <- NULL
for(g in 1:G){
  Seurat_corrected <- rbind(Seurat_corrected, simulation.integrated@assays$integrated[g,])
}

#######
# ARI #
#######
ARI_Seurat <- adjustedRandIndex(metadata$celltype,w_Seurat)

#################################
# scatter plots (t-SNE and PCA) #
#################################
if(!dir.exists("Image")){
  dir.create("Image")
}

#####set cell type colorings
# batch color bar
color_by_batch<-c("#EB4334","#FBBD06","#35AA53","#4586F3")

# cell type color bar
color_by_celltype<-c("#F88A7E", "#FFD87D", "#ABD978","#8097D3","#9C7ACE")


# plot colored by cell type
allcolors<-color_by_celltype[metadata$celltype]#[allsamples])

plot_by_celltype<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  # The image with legend
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(10,10,2,2))
  plot(Y[,1], Y[,2],t="n",xaxt="n",yaxt="n",xlab="", ylab="")
  axis(1,line = 2.5, cex.axis=6,tick = F)#plot the x axis
  axis(2,cex.axis=6,tick = F)#plot the y axis
  mtext(xlab, side=1, line=7.5, cex=6)
  mtext(ylab, side=2, line=5, cex=6)
  points(Y[,1], Y[,2], cex=3,
         pch=20, 
         col=alpha(allcolors[subset],0.6)) 
  dev.off()
}

# plot colored by batch
batch.cols<-color_by_batch[rep(1:B,nb)]
plot_by_batch<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(10,10,2,2))
  plot(Y[,1], Y[,2],t="n",xaxt="n",yaxt="n",xlab="", ylab="")
  axis(1,line = 2.5, cex.axis=6,tick = F)#plot the x axis
  axis(2,cex.axis=6,tick = F)#plot the y axis
  mtext(xlab, side=1, line=7.5, cex=6)
  mtext(ylab, side=2, line=5, cex=6)
  points(Y[,1], Y[,2], cex=3,
         pch=20, 
         col=alpha(batch.cols[subset],0.6)) 
  dev.off()
}

#################################################
# Draw t-SNE plots by batch and true cell types #
#################################################
set.seed(123)
Seurat_dist <- dist(Seurat_PCA)
all.dists.Seurat <- as.matrix(Seurat_dist)
tsne_Seurat_dist <- Rtsne(all.dists.Seurat, is_distance=TRUE, perplexity = 30)

Seurat_by_celltype<- "Image/tsne_simulation_Seurat_by_celltype.jpeg"
plot_by_celltype(Seurat_by_celltype, tsne_Seurat_dist$Y)

Seurat_by_batch <- "Image/tsne_simulation_Seurat_by_batch.jpeg"
plot_by_batch(Seurat_by_batch, tsne_Seurat_dist$Y)

##################
# Draw PCA plots #
##################
pca.Seurat <- prcomp(Seurat_PCA, rank=2)
Seu_PCA<- "Image/pca_simulation_Seurat.jpeg"
plot_by_celltype(Seu_PCA, Seurat_PCA[,1:2], xlab = "PC 1", ylab = "PC 2")

# Store the workspace
save.image("Seurat_workspace.RData")
save(ARI_Seurat,tsne_Seurat_dist,file = "Seurat_results.RData")

