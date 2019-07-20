rm(list=ls())
library(mclust)
library(ggplot2)
library(Rtsne)
library(cluster)

###################
# Load Workspaces #
###################
# Working directory
setwd("F:/scRNA/code/0601Cpp_BUSseq/BUSseq_implementation_v1/MouseHematopoietic/")

load("BUSseq/BUSseq_workspace.RData")
comparison_list <- c("liger", "MNN", "Scanorama", "scVI", "Seurat", "ZINBWaVE")
num_comparison <- length(comparison_list)

for(m in 1:num_comparison){
  load(paste0("Comparison/",comparison_list[m],"/",comparison_list[m],"_workspace.RData"))
}

# load the dimentsion information
dim <- unlist(read.table("RawCountData/dim_hemat_v1.txt"))
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

# Loading the true cell type indicators
metadata <- read.table("RawCountData/hemat_metadata.txt")

# Calculating the ARI
ARI_values[1] <- ARI_BUSseq
ARI_values[2] <- ARI_liger
ARI_values[3] <- ARI_MNN
ARI_values[4] <- ARI_scanorama
ARI_values[5] <- ARI_scVI
ARI_values[6] <- ARI_Seurat
ARI_values[7] <- ARI_ZINBWaVE

write.table(ARI_values, "Results/ARI_values.txt",col.names = F)

#######################################################################
# Silhouette coefficient by estimated cell type labels of each method #
#######################################################################
# Storing the Silhouette coefficients
sil_estimated <- matrix(NA,N,num_comparison+1)
colnames(sil_estimated) <- names(ARI_values)

sil_estimated[,1] <- sil_BUSseq[,3]
sil_estimated[,2] <- sil_liger[,3]
sil_estimated[,3] <- sil_MNN[,3]
sil_estimated[,4] <- sil_scanorama[,3]
sil_estimated[,5] <- sil_scVI[,3]
sil_estimated[,6] <- sil_Seurat[,3]
sil_estimated[,7] <- sil_ZINBWaVE[,3]

write.table(sil_estimated,"./Results/Silhouette_coef.txt",row.names = F)

# Drawing the violin plot
sil_com<-data.frame(sli_coef=as.vector(sil_estimated),
                    method=factor(rep(colnames(sil_estimated),each=N)))

png("Image/violin_plot_of_Silhouette_coef_est.png",width = 960, height = 720)
p <- ggplot(sil_com, aes(x=method, y=sli_coef,fill=method)) +
  geom_violin() + xlab(NULL) + ylab(NULL)+theme_bw() +
  theme(text = element_text(size=48),axis.text.x =element_text(size=32,angle = 90), axis.text.y = element_text(size =32),panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
p
dev.off()

##################################################################
# Silhouette coefficient by true cell type labels of each method #
##################################################################
# Storing the Silhouette coefficients
sil_true <- matrix(NA,N,num_comparison+1)
colnames(sil_true) <- names(ARI_values)

sil_true[,1] <- sil_BUSseq_true[,3]
sil_true[,2] <- sil_liger_true[,3]
sil_true[,3] <- sil_MNN_true[,3]
sil_true[,4] <- sil_scanorama_true[,3]
sil_true[,5] <- sil_scVI_true[,3]
sil_true[,6] <- sil_Seurat_true[,3]
sil_true[,7] <- sil_ZINBWaVE_true[,3]

write.table(sil_true,"./Results/Silhouette_coef_true.txt",row.names = F)

# Drawing the violin plot
sil_com_true<-data.frame(sli_coef=as.vector(sil_true),
                    method=factor(rep(colnames(sil_true),each=N)))

png("Image/violin_plot_of_Silhouette_coef_true.png",width = 960, height = 720)
p <- ggplot(sil_com_true, aes(x=method, y=sli_coef,fill=method)) +
  geom_violin() + xlab(NULL) + ylab(NULL)+theme_bw() +
  theme(text = element_text(size=48),axis.text.x =element_text(size=32,angle = 90), axis.text.y = element_text(size =32),panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
p
dev.off()

###############################################################################
# Silhouette coefficient of tSNE by estimated cell type labels of each method #
###############################################################################
# Storing the Silhouette coefficients
sil_tSNE_est <- matrix(NA,N,num_comparison+1)
colnames(sil_tSNE_est) <- names(ARI_values)

BUSseq_dist_tSNE <- dist(tsne_BUSseq_dist$Y)
sil_tSNE_est[,1] <- silhouette(as.integer(w_BUSseq), dist = BUSseq_dist_tSNE)[,3]

liger_dist_tSNE <- dist(tsne_liger_dist$Y)
sil_tSNE_est[,2] <- silhouette(as.integer(w_liger), dist = liger_dist_tSNE)[,3]

MNN_dist_tSNE <- dist(tsne_MNN_dist$Y)
sil_tSNE_est[,3] <- silhouette(as.integer(w_MNN), dist = MNN_dist_tSNE)[,3]

scanorama_dist_tSNE <- dist(tsne_scanorama_dist$Y)
sil_tSNE_est[,4] <- silhouette(as.integer(w_scanorama), dist = scanorama_dist_tSNE)[,3]

scVI_dist_tSNE <- dist(tsne_scVI_dist$Y)
sil_tSNE_est[,5] <- silhouette(as.integer(w_scVI), dist = scVI_dist_tSNE)[,3]

Seurat_dist_tSNE <- dist(tsne_Seurat_dist$Y)
sil_tSNE_est[,6] <- silhouette(as.integer(w_Seurat), dist = Seurat_dist_tSNE)[,3]

ZINBWaVE_dist_tSNE <- dist(tsne_ZINBW_dist$Y)
sil_tSNE_est[,7] <- silhouette(as.integer(w_ZINBWaVE), dist = ZINBWaVE_dist_tSNE)[,3]

write.table(sil_tSNE_est,"./Results/Silhouette_coef_tSNE_est.txt",row.names = F)

# Drawing the violin plot
sil_com_tSNE_est<-data.frame(sli_coef=as.vector(sil_tSNE_est),
                             method=factor(rep(colnames(sil_tSNE_est),each=N)))

png("Image/violin_plot_of_Silhouette_coef_tSNE_est.png",width = 960, height = 720)
p <- ggplot(sil_com_tSNE_est, aes(x=method, y=sli_coef,fill=method)) +
  geom_violin() + xlab(NULL) + ylab(NULL)+theme_bw() +
  theme(text = element_text(size=48),axis.text.x =element_text(size=32,angle = 90), axis.text.y = element_text(size =32),panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
p
dev.off()

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

write.table(sil_tSNE_true,"./Results/Silhouette_coef_tSNE_true.txt",row.names = F)

# Drawing the violin plot
sil_com_tSNE_true<-data.frame(sli_coef=as.vector(sil_tSNE_true),
                             method=factor(rep(colnames(sil_tSNE_true),each=N)))

png("Image/violin_plot_of_Silhouette_coef_tSNE_true.png",width = 960, height = 720)
p <- ggplot(sil_com_tSNE_true, aes(x=method, y=sli_coef,fill=method)) +
  geom_violin() + xlab(NULL) + ylab(NULL)+theme_bw() +
  theme(text = element_text(size=48),axis.text.x =element_text(size=32,angle = 90), axis.text.y = element_text(size =32),panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
p
dev.off()


png("Image/boxplot_of_Silhouette_coef_tSNE_true.png",width = 960, height = 720)
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
# coloring for cell types
color.legendF <- c(MEP="orange", GMP="chartreuse4", CMP="magenta",
                   LTHSC="dodgerblue", MPP="blue", LMPP="light blue", other="grey")
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

# Uncorrected #
data_raw <- NULL
for(b in 1:B){
  data_raw <- cbind(data_raw, HematCounts[[b]])
}

set.seed(123)
all.dists.unc <- as.matrix(dist(log1p(t(data_raw))))
tsne_uncorrected <- Rtsne(all.dists.unc, is_distance=TRUE, perplexity = 30)

unc_by_celltype<- "Image/tsne_hemat_uncorrected_by_celltype.jpeg"
plot_by_celltype(unc_by_celltype, tsne_uncorrected$Y)

unc_by_batch <- "Image/tsne_hemat_uncorrected_by_batch.jpeg"
plot_by_batch(unc_by_batch, tsne_uncorrected$Y)

##################
# Draw PCA plots #
##################

# Uncorrected #
pca.unc <- prcomp(log1p(t(data_raw)), rank=2)
unc_PCA<- "Image/pca_hemat_uncorrected.jpeg"
plot_by_celltype(unc_PCA, pca.unc$x, xlab = "PC 1", ylab = "PC 2")

# legend
jpeg(file="Image/legend.jpeg", width=750, height=600,quality = 100)
plot(c(-5,5),c(-4,4),type="n", bty="n", axes=FALSE, xlab="", ylab="")
legend(x=-5, y=4, legend=c("GSE72857","GSE81682"), pch=c(1,21), cex=2.5, col=c(batchcolor[1],"black"), pt.bg=batchcolor, title = "Batch", bty="n")
legend(x=-1, y=4, legend=names(color.legendF), pch=21, cex=2.5, col="black", pt.bg=color.legendF, title = "GSE81682", bty="n")
legend(x=2, y=4, legend=names(color.legendA)[1:3], pch=1, cex=2.5, col=color.legendA[1:3], pt.bg=color.legendA[1:3], title = "GSE72857", bty="n")
dev.off()

# Store the workspace
save.image("comparison_workspace.RData")