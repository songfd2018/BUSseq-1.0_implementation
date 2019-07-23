rm(list=ls())
library(mclust)
library(cluster)
library(ggplot2)
library(Rtsne)

library(scales)
require(WGCNA)
library(RColorBrewer)

###################
# Load Workspaces #
###################
# Working directory
setwd("F:/scRNA/code/0601Cpp_BUSseq/BUSseq_implementation_v1/HumanPancreas/")

load("BUSseq/BUSseq_workspace.RData")
comparison_list <- c("liger", "MNN", "Scanorama", "scVI", "Seurat", "ZINBWaVE")
num_comparison <- length(comparison_list)

for(m in 1:num_comparison){
  load(paste0("Comparison/",comparison_list[m],"/",comparison_list[m],"_workspace.RData"))
}

# load the dimentsion information
dim <- unlist(read.table("RawCountData/dim_pancreas_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[1:B + 3]

#############################################################
# Evaluation of the Performance of All the Methods          #
# according to ARI, silhouette coefficients and t-SNE plots #
#############################################################
if(!dir.exists("Results")){
  dir.create("Results")
}

if(!dir.exists("Image")){
  dir.create("Image")
}

#######
# ARI #
#######
ARI_values <- rep(NA, num_comparison + 1)
names(ARI_values) <- c("BUSseq", comparison_list)

# Calculating the ARI
ARI_values[1] <- ARI_BUSseq
ARI_values[2] <- ARI_liger
ARI_values[3] <- ARI_MNN
ARI_values[4] <- ARI_scanorama
ARI_values[5] <- ARI_scVI
ARI_values[6] <- ARI_Seurat
ARI_values[7] <- ARI_ZINBWaVE

write.table(ARI_values, "Results/ARI_values.txt",col.names = F)


###############################################################################
# Silhouette coefficient of tSNE by true cell type labels of each method #
###############################################################################
# Storing the Silhouette coefficients
sil_tSNE_true <- matrix(NA,N,num_comparison+1)
colnames(sil_tSNE_true) <- names(ARI_values)

BUSseq_dist_tSNE <- dist(tsne_BUSseq_dist$Y)
sil_tSNE_true[,1] <- silhouette(as.integer(metadata$CellType), dist = BUSseq_dist_tSNE)[,3]

liger_dist_tSNE <- dist(tsne_liger_dist$Y)
sil_tSNE_true[,2] <- silhouette(as.integer(metadata$CellType), dist = liger_dist_tSNE)[,3]

MNN_dist_tSNE <- dist(tsne_MNN_dist$Y)
sil_tSNE_true[,3] <- silhouette(as.integer(metadata$CellType), dist = MNN_dist_tSNE)[,3]

scanorama_dist_tSNE <- dist(tsne_scanorama_dist$Y)
sil_tSNE_true[,4] <- silhouette(as.integer(metadata$CellType), dist = scanorama_dist_tSNE)[,3]

scVI_dist_tSNE <- dist(tsne_scVI_dist$Y)
sil_tSNE_true[,5] <- silhouette(as.integer(metadata$CellType), dist = scVI_dist_tSNE)[,3]

Seurat_dist_tSNE <- dist(tsne_Seurat_dist$Y)
sil_tSNE_true[,6] <- silhouette(as.integer(metadata$CellType), dist = Seurat_dist_tSNE)[,3]

ZINBWaVE_dist_tSNE <- dist(tsne_ZINBW_dist$Y)
sil_tSNE_true[,7] <- silhouette(as.integer(metadata$CellType), dist = ZINBWaVE_dist_tSNE)[,3]

write.table(sil_tSNE_true,"./Results/Silhouette_coefs.txt",row.names = F)

# Drawing the violin plot
sil_com_tSNE_true<-data.frame(sli_coef=as.vector(sil_tSNE_true),
                              method=factor(rep(colnames(sil_tSNE_true),each=N)))

# png("Image/violin_plot_of_Silhouette_coef_tSNE_true.png",width = 960, height = 720)
# p <- ggplot(sil_com_tSNE_true, aes(x=method, y=sli_coef,fill=method)) +
#   geom_violin() + xlab(NULL) + ylab(NULL)+theme_bw() +
#   theme(text = element_text(size=48),axis.text.x =element_text(size=32,angle = 90), axis.text.y = element_text(size =32),panel.border = element_blank(),panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
# p
# dev.off()

png("Image/boxplot_of_Silhouette_coefs.png",width = 960, height = 720)
p <- ggplot(sil_com_tSNE_true, aes(x=method, y=sli_coef,fill=method)) +
  geom_boxplot() + xlab(NULL) + ylab(NULL)+theme_bw() +
  theme(text = element_text(size=48),axis.text.x =element_text(size=32,angle = 90), axis.text.y = element_text(size =32),panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
p
dev.off()

#################################
# scatter plots (t-SNE and PCA) #
#################################
if(!dir.exists("Image")){
  dir.create("Image")
}

allcolors<-labels2colors(metadata$CellType)#[allsamples])

allcolors[allcolors=="red"]<-"deeppink"
allcolors[allcolors=="yellow"]<-"orange1"#"darkgoldenrod1"

plot_by_celltype<-function(pic_name,Y,subset=NULL,...,xlab = "", ylab = "",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  # The image with legend
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(6,6,4,2))
  plot(Y[,1], Y[,2], cex=3,
       pch=20, 
       col=alpha(allcolors[subset],0.6),
       xaxt="n",yaxt="n",
       xlab=xlab, ylab=ylab, main=main) 
  axis(1,cex.axis=6,line=2.5,tick = F)#plot the x axis
  axis(2,cex.axis=6,tick = F)#plot the y axis
  dev.off()
}

colors4<-brewer.pal(4,"Set3")
batch.cols<-colors4[rep(1:B,nb)]
plot_by_batch<-function(pic_name,Y,subset=NULL,...,xlab = "", ylab = "",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(6,6,4,2))
  plot(Y[,1], Y[,2], cex=3,
       pch=20, 
       col=alpha(batch.cols[subset],0.6),
       xaxt="n",yaxt="n",
       xlab=xlab, ylab=ylab, main=main)#,  xlab="tSNE 1",ylab="tSNE 2")
  axis(1,cex.axis=6,line=2.5,tick = F)#plot the x axis
  axis(2,cex.axis=6,tick = F)#plot the y axis
  dev.off()
}

# Uncorrected #
data_raw <- read.table("RawCountData/count_data_pancreas_v1.txt")

set.seed(123)
all.dists.unc <- as.matrix(dist(log1p(t(data_raw))))
tsne_uncorrected <- Rtsne(all.dists.unc, is_distance=TRUE, perplexity = 30)

unc_by_celltype<- "Image/tsne_pancreas_uncorrected_by_celltype.jpeg"
plot_by_celltype(unc_by_celltype, tsne_uncorrected$Y)

unc_by_batch <- "Image/tsne_pancreas_uncorrected_by_batch.jpeg"
plot_by_batch(unc_by_batch, tsne_uncorrected$Y)

##################
# Draw PCA plots #
##################

# Uncorrected #
pca.unc <- prcomp(log1p(t(data_raw)), rank=2)
unc_PCA<- "Image/pca_pancreas_uncorrected.jpeg"
plot_by_celltype(unc_PCA, pca.unc$x, xlab = "PC 1", ylab = "PC 2")

# legend
jpeg(file="Image/legend.jpeg", width=600, height=480)
plot(c(-5,5),c(-4,4),type="n", bty="n", axes=FALSE, xlab="", ylab="")
legend(x=-5, y=4, legend=c("GSE81076","GSE85241","GSE86473","E-MTAB-5061"), pch=20, cex=2.5, col=alpha(colors4,0.6), title = "Batch", bty="n")
legend(x=-0, y=4, legend=names(table(metadata$CellType)), pch=20, cex=2.5, col=alpha(c("turquoise","blue","brown","orange1","green","deeppink","black"),0.6), title = "Cell Type", bty="n")
dev.off()

# Store the workspace
save.image("comparison_workspace.RData")