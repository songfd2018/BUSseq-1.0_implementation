# library(devtools)
# install_github('MacoskoLab/liger')
rm(list=ls())
library(liger)
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot

# coloring
library(scales)
require(WGCNA)
library(RColorBrewer)

set.seed(12345)
# Working directory
setwd("/your/working/directory/BUSseq_implementation-1.0/Simulation/Comparison/liger/")

###################
# Load simulation Data #
###################
# Loading simulation count data
countdata <- read.table("../../RawCountData/count_data_simulation_v4.txt") 

# Load dimension
dim <- unlist(read.table("../../RawCountData/dim_simulation_v4.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# Load metadata
metadata <- read.table("../../RawCountData/metadata_simulation_v4.txt")

# Load gene_list
gene_list <- paste0("gene_",1:G)

rownames(countdata) <- gene_list
colnames(countdata) <- rownames(metadata)
simulation.list <- list(Batch1 = countdata[,1:nb[1]], 
                        Batch2 = countdata[,1:nb[2]+nb[1]], 
                        Batch3 = countdata[,1:nb[3]+sum(nb[1:2])], 
                        Batch4 = countdata[,1:nb[4]+sum(nb[1:3])])

################################################################################################################
# Conduct analysis referring to https://github.com/MacoskoLab/liger/blob/master/vignettes/walkthrough_pbmc.pdf #
################################################################################################################
# Create liger object
liger_simulation <- createLiger(simulation.list, remove.missing = FALSE)

# Normalization
liger_simulation <- normalize(liger_simulation)
liger_simulation@var.genes <- gene_list
liger_simulation <- scaleNotCenter(liger_simulation)

# running suggestK on multiple cores can greatly decrease the runtime
# looking for the section of the plot where this metric stops increasing as sharply
k.suggest <- suggestK(liger_simulation, num.cores = 5, gen.new = T, return.data = T, plot.log2 = F)


# Take the lowest objective of three factorizations with different initializations
# Multiple restarts are recommended for initial analyses since iNMF is non-deterministic
liger_simulation <- optimizeALS(liger_simulation, k=30, thresh = 5e-5, nrep = 3)

liger_simulation <- runTSNE(liger_simulation, use.raw = T)
p1 <- plotByDatasetAndCluster(liger_simulation, return.plots = T)
# Plot by dataset
print(p1[[1]])


##########################
# Currently working here #
##########################

# To better integrate the datasets, we perform a quantile alignment step. This process first identifies similarly
# loading cells across datasets by building a similarity graph based on shared factor neighborhoods.
liger_simulation <- quantileAlignSNF(liger_simulation, resolution = 0.4, small.clust.thresh = 20)
w_liger <- liger_simulation@alignment.clusters

#######
# ARI #
#######
ARI_liger <- adjustedRandIndex(metadata$celltype,w_liger)

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
loading_liger <- liger_simulation@H.norm
liger_dist <- dist(loading_liger)
all.dists.liger <- as.matrix(liger_dist)
tsne_liger_dist <- Rtsne(all.dists.liger, is_distance=TRUE, perplexity = 30)

liger_by_celltype<- "Image/tsne_simulation_liger_by_celltype.jpeg"
plot_by_celltype(liger_by_celltype, tsne_liger_dist$Y)

liger_by_batch <- "Image/tsne_simulation_liger_by_batch.jpeg"
plot_by_batch(liger_by_batch, tsne_liger_dist$Y)

# ##################
# # Draw PCA plots #
# ##################
# pca.liger <- prcomp(loading_liger, rank=2)
# liger_PCA<- "Image/pca_simulation_liger.jpeg"
# plot_by_celltype(liger_PCA, pca.liger$x, xlab = "PC 1", ylab = "PC 2")

# Store the workspace
save.image("liger_workspace.RData")
save(ARI_liger,tsne_liger_dist,file = "liger_results.RData")