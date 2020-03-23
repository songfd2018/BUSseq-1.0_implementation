# Referring to the source code https://github.com/MarioniLab/MNN2017/blob/master/Haematopoiesis/prepareData.R

# This script prepares data for the hematopoiesis analysis.
# It involves two publicly available datasets.

##########################################
##########################################
rm(list=ls())

setwd("RecoverMNN/")
if(!dir.exists("Image")){
  dir.create("Image")
}

###################
###################
## prepareData.R ##
###################
###################

# Download and read the counts, metadata of Nestorowa et al. 2016
fname <- "../RawData/GSE81682_HTSeq_counts.txt.gz"
if (!file.exists(fname)) { download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81682&format=file&file=GSE81682%5FHTSeq%5Fcounts%2Etxt%2Egz", fname) }
dataF <- read.table(fname, header=TRUE, row.names=1, check.names=FALSE)
dataF <- as.matrix(dataF)
dim(dataF)

fname <- "../RawData/metaF.txt"
if (!file.exists(fname)) { download.file("http://blood.stemcells.cam.ac.uk/data/all_cell_types.txt", fname) }
metaF <- read.table(fname, stringsAsFactors = FALSE, header=TRUE, check.names=FALSE)
metainds <- match(colnames(dataF), rownames(metaF))
missing.meta <- is.na(metainds)
metaF <- metaF[metainds,] # This will contain NA's... which is okay, at this point, to preserve length.

# Defining the cell type based on the metadata.
metatypeF <- rep("other", nrow(metaF))
for (col in rev(colnames(metaF))) { # reverse, so earlier columns end up overwriting later ones.
  chosen <- metaF[,col]==1
  metatypeF[chosen] <- sub("[0-9]?_.*", "", col)
}
metatypeF[metatypeF=="ESLAM"] <- "HSPC"

# Filling in metadata from the cell sorting label, if metadata was missing.
metatypeF[missing.meta] <- sub("_.*", "", colnames(dataF)[missing.meta])
metatypeF[metatypeF=="LT-HSC"] <- "LTHSC"
metatypeF[metatypeF=="Prog"] <- "other"
colnames(dataF)<-metatypeF

# Perform size factor normalization within this data set.
library(scran)
high.abF <- scater::calcAverage(dataF) > 1
clustF <- quickCluster(dataF, method="igraph", subset.row=high.abF)
sizeF <- computeSumFactors(dataF, cluster=clustF, subset.row=high.abF)
dataF2 <- t(t(dataF)/sizeF)

# Cleaning up memory.
gc() 

##########################################
##########################################

# Download and read the counts and meta data of Paul et al. 2015
fname <- "../RawData/GSE72857_umitab.txt.gz"
if (!file.exists(fname)) { download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE72857&format=file&file=GSE72857%5Fumitab%2Etxt%2Egz", fname) }
dataA <- read.table(fname, header=TRUE, row.names=1)
metaA <- read.csv2("../RawData/MAP.csv",sep=",",stringsAsFactors = FALSE, head=TRUE, row.names=1)
dim(dataA)

# Only selecting cells that are in the metadata.
metainds <- match(rownames(metaA), colnames(dataA))
dataA <- dataA[,metainds]
dataA <- as.matrix(dataA)

# Organizing cell type labels.
metatypeA <- character(nrow(metaA))
metatypeA[metaA[,1]<7] <- "ERY"
metatypeA[metaA[,1]>6 & metaA[,1]<12] <- "CMP"
metatypeA[metaA[,1]>11] <- "GMP"
colnames(dataA) <- metatypeA

# Perform size factor normalization within this data set.
high.abA <- scater::calcAverage(dataA) > 1
clustA <- quickCluster(dataA, method="igraph", subset.row=high.abA)
sizeA <- computeSumFactors(dataA, cluster=clustA, subset.row=high.abA)
dataA2 <- t(t(dataA)/sizeA)

# Cleaning up memory.
gc() 

##########################################
##########################################

# Download list of highly variable genes identified by Nestrowa et al. 2016
fname <- "../RawData/coordinates_gene_counts_flow_cytometry.txt.gz"
if (!file.exists(fname)) { download.file("http://blood.stemcells.cam.ac.uk/data/coordinates_gene_counts_flow_cytometry.txt.gz", fname) }
TFs <- read.table(fname, nrows=1, stringsAsFactors=FALSE)
features <- as.character(unlist(TFs))
features <- features[grep("ENSMUS", features)]

# Pull down IDs from BioMaRt.
library(biomaRt)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="www.ensembl.org" )
out <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), values = features, mart = mart,filters = "ensembl_gene_id")

# Selecting features that are HVGs and present in both data sets.
mF <- match(out$ensembl_gene_id, rownames(dataF2))
mA <- pmatch(out$mgi_symbol, rownames(dataA2)) # partial, due to use of concatenated gene symbols.
keep <- !is.na(mF) & !is.na(mA)

dataA3 <- dataA2[mA[keep],]
dataF3 <- dataF2[mF[keep],]
rownames(dataA3) <- rownames(dataF3)

# Rescaling the first dataset to match the coverage of the second.
aveA <- rowMeans(dataA3)
aveF <- rowMeans(dataF3)
dataF3 <- dataF3 * median(aveA/aveF)

# Perform log-transformation and save results to file.
logDataF3 <- log(1 + dataF3)
logDataA3 <- log(1 + dataA3)
save(logDataA3, logDataF3, file="logdataFandA_all.RData")

#######################
#######################
## plotCorrections.R ##
#######################
#######################
# Load the output of "preparedata.R".
load("logdataFandA_all.RData") 
colnames(logDataF3)[colnames(logDataF3)=="other"] <- "Unsorted"
colnames(logDataA3)[colnames(logDataA3)=="ERY"] <- "MEP"
raw.all <- cbind(logDataF3, logDataA3)
first.batch <- rep(c(TRUE, FALSE), c(ncol(logDataF3), ncol(logDataA3)))

# Adding colours.
#base.color <- "grey"
color.legendF <- c(MEP="orange", GMP="chartreuse4", CMP="magenta", 
                   HSPC="cyan", LTHSC="dodgerblue", MPP="blue", LMPP="light blue", Unsorted="grey")
colmatF <- col2rgb(color.legendF) 

colmatA <- colmatF + 100 # A lighter shade.
colmatA[colmatA > 255] <- 255

#colmatF<-colmatF+200
#colmatF[colmatF > 255] <- 255
color.legendA <- setNames(rgb(colmatA[1,], colmatA[2,], colmatA[3,], maxColorValue=255), names(color.legendF))
allcolors <- c(color.legendF[colnames(logDataF3)], color.legendA[colnames(logDataA3)])
batch<-c( rep(1,ncol(logDataF3)),rep(2,ncol(logDataA3)) )

# Only keeping common cell types for PCA.
celltypes <- c(colnames(logDataF3), colnames(logDataA3))
pca.retain <- celltypes %in% c("MEP", "GMP", "CMP") 

# Making a plotting function.
plotFUN <- function(fname, Y, subset=NULL, ..., xlab="tSNE 1",ylab="tSNE 2",main="") {
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  png(fname,width=900,height=700)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=ifelse(first.batch, 21, 1)[subset], 
       col=ifelse(first.batch, "black", allcolors)[subset],
       bg=allcolors[subset], xlab=xlab, ylab=ylab, main=main) 
  dev.off()
}

batchcolor=c("lavender","lightcoral")
plotFUNb <- function(fname, Y, subset=NULL, ...) {
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  png(fname,width=900,height=700)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=ifelse(first.batch, 21, 1)[subset], 
       col=ifelse(first.batch, "black", batchcolor[batch[subset]]),
       bg=batchcolor[batch[subset]], ...)#,  xlab="tSNE 1",ylab="tSNE 2")
  dev.off()
}

######################
# No correction.

X.unc <- raw.all
t.unc <- t(X.unc)

## Generating a t-SNE plot.
require(Rtsne)
set.seed(0)
all.dists.unc <- as.matrix(dist(t.unc))
tsne.unc <- Rtsne(all.dists.unc, is_distance=TRUE, perplexity = 90)
plotFUN("Image/uncFA.png", tsne.unc$Y, main="Uncorrected",  xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("Image/uncFAb.png", tsne.unc$Y, main="Uncorrected",  xlab="tSNE 1",ylab="tSNE 2")

rm(all.dists.unc)
gc()

# Generating a PCA plot.
pca.unc <- prcomp(t.unc[pca.retain,], rank=2)
# pca.unc$x[ (pca.unc$x<(-0.08))]<- (-0.08)
plotFUN("Image/pca_raw.png", pca.unc$x, subset=pca.retain, main="Uncorrected", ylim=c(-0.1, max(pca.unc$x[,2])),  xlab="PC 1",ylab="PC 2")

rm(t.unc)
gc()

###########################################################
###########################################################
## plot t-SNE and pca plot by our preprocessing strategy ##
###########################################################
###########################################################
# Load raw count data 
data_raw <- as.matrix(read.table("../RawCountData/count_data_hemat_v1.txt"))

# Load metadata
metadata <- read.table("../RawCountData/metadata_hemat_v1.txt",stringsAsFactors = FALSE)

# load the dimentsion information
dim <- unlist(read.table("../RawCountData/dim_hemat_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[1:B + 3]

# Load gene lists
gene_list <- unlist(read.table("../RawCountData/gene_list_hemat_v1.txt",stringsAsFactors = FALSE))

if(!dir.exists("Image")){
  dir.create("Image")
}
############
# coloring #
############
raw.all <- data_raw
# raw.all[, 1:nb[2]] <- data_raw[, nb[1] + 1:nb[2]]
# raw.all[, 1:nb[1] + nb[2]] <- data_raw[, 1:nb[1]]

first.batch <- rep(c(FALSE, TRUE), nb)

# Adding colours.
#base.color <- "grey"
color.legendF <- c(MEP="orange", GMP="chartreuse4", CMP="magenta", 
                   HSPC="cyan", LTHSC="dodgerblue", MPP="blue", LMPP="light blue", other="grey")
colmatF <- col2rgb(color.legendF) 

colmatA <- colmatF + 100 # A lighter shade.
colmatA[colmatA > 255] <- 255

#colmatF<-colmatF+200
#colmatF[colmatF > 255] <- 255
color.legendA <- setNames(rgb(colmatA[1,], colmatA[2,], colmatA[3,], maxColorValue=255), names(color.legendF))
allcolors <- c(color.legendA[metadata$CellType[1:nb[1]]], color.legendF[metadata$CellType[nb[1] + 1:nb[2]]])
batch<-rep(1:2, nb)

# Only keeping common cell types for PCA.
pca.retain <- metadata$CellType %in% c("MEP", "GMP", "CMP") 

# Making a plotting function.
plotFUN <- function(fname, Y, subset=NULL, ..., xlab="tSNE 1",ylab="tSNE 2",main="") {
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  png(fname,width=900,height=700)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=ifelse(first.batch, 21, 1)[subset], 
       col=ifelse(first.batch, "black", allcolors)[subset],
       bg=allcolors[subset], xlab=xlab, ylab=ylab, main=main) 
  dev.off()
}

batchcolor=c("lightcoral","lavender")
plotFUNb <- function(fname, Y, subset=NULL, ...) {
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  png(fname,width=900,height=700)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=ifelse(first.batch, 21, 1)[subset], 
       col=ifelse(first.batch, "black", batchcolor[batch[subset]]),
       bg=batchcolor[batch[subset]])#,  xlab="tSNE 1",ylab="tSNE 2")
  dev.off()
}

set.seed(123)
all.dists.unc <- as.matrix(dist(log1p(t(raw.all))))
tsne_uncorrected <- Rtsne(all.dists.unc, is_distance=TRUE, perplexity = 30)
plotFUN("Image/uncorrected_hemat_tsne_our.png", tsne_uncorrected$Y, main="Uncorrected", xlab="tSNE 1", ylab="tSNE 2")
plotFUNb("Image/uncorrected_hemat_tsne_our_batch.png", tsne_uncorrected$Y, main="Uncorrected",  xlab="tSNE 1",ylab="tSNE 2")
# rm(all.dists.unc)
gc()

# Generating a PCA plot.
pca.unc <- prcomp(log1p(t(raw.all[,pca.retain])), rank=2)
plotFUN("Image/uncorrected_hemat_pca_our.png", pca.unc$x, subset=pca.retain, main="Uncorrected", xlab="PC 1",ylab="PC 2")

pca.unc.all <- prcomp(log1p(t(raw.all)), rank=2)
# pca.unc$x[ (pca.unc$x<(-0.08))]<- (-0.08)
plotFUN("Image/uncorrected_hemat_pca_all_our.png", pca.unc.all$x, main="Uncorrected",  xlab="PC 1",ylab="PC 2")

rm(t.unc)
gc()

##############################################
# Draw t-SNE and PCA after our preprocessing #
##############################################
load("F:/scRNA/code/0601Cpp_BUSseq/BUSseq_implementation_v1/MouseHematopoietic/comparison_workspace.RData")
MNN_by_celltype<- "Image/tsne_hemat_uncorrected_by_celltype.jpg"
dat_frame <- data.frame(Var1 = tsne_uncorrected$Y[,1], 
                        Var2 = tsne_uncorrected$Y[,2], 
                        col = celltype_factor)

jpeg(MNN_by_celltype,width = 800, height = 600, quality = 100)
ggplot(dat_frame, aes(x= Var1, y= Var2, colour= col)) +
  geom_point(size=4) + theme_classic() +
  scale_colour_manual(values = alpha(color_by_celltype,0.6)) +
  xlab("tSNE 1") + ylab("tSNE 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 44), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 44),#, angle = 45))
        axis.title=element_text(size=48,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        #legend.position = "none")
        legend.text=element_text(size=28,face="bold"),
        legend.title = element_text(size = 36,face="bold"))
dev.off()

MNN_by_batch <- "Image/tsne_hemat_uncorrected_by_batch.jpg"
dat_frame <- data.frame(Var1 = tsne_uncorrected$Y[,1], 
                        Var2 = tsne_uncorrected$Y[,2], 
                        col = batch_factor)

jpeg(MNN_by_batch,width = 800, height = 600, quality = 100)
ggplot(dat_frame, aes(x= Var1, y= Var2, colour= col)) +
  geom_point(size=4) + theme_classic() +
  scale_colour_manual(values = alpha(color_by_batch,0.6)) +
  xlab("tSNE 1") + ylab("tSNE 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 36), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 36),#, angle = 45))
        axis.title=element_text(size=44,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        #legend.position = "none")
        legend.text=element_text(size=28,face="bold"),
        legend.title = element_text(size = 36,face="bold"))
dev.off()

unc_PCA<- "Image/pca_hemat_uncorrected.jpeg"
dat_frame <- data.frame(Var1 = pca.unc$x[,1], 
                        Var2 = pca.unc$x[,2], 
                        col = celltype_factor)

jpeg(unc_PCA,width = 900, height = 700, quality = 100)
ggplot(dat_frame, aes(x= Var1, y= Var2, colour= col)) +
  geom_point(size=4) + theme_classic() +
  scale_colour_manual(values = alpha(color_by_celltype,0.6)) +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 36), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 36),#, angle = 45))
        axis.title=element_text(size=44,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        #legend.position = "none")
        legend.text=element_text(size=28,face="bold"),
        legend.title = element_text(size = 36,face="bold"))
dev.off()

################################################################################################
################################################################################################
## Plot PCA and t-SNE plots after normalization and use the same 3,470 HVGs as our manuscript ##
################################################################################################
################################################################################################
# Make up for the genes with the changed Gene Symbol
# Selecting features that are HVGs and present in both data sets.
ind_F <- match(gene_list, out$ensembl_gene_id)
interest_symbol <- out[ind_F,]
ind_A <- pmatch(interest_symbol$mgi_symbol,rownames(dataA)) # partial, due to use of concatenated gene symbols.
# Thirteen genes have been renamed in the BioMart database, so they fail to be matched by pmatch.
rename_gene <- which(is.na(ind_A))
print(length(rename_gene))
  
# We add them mutually
# "Dscr3" --> "Vps26c"
# referring to https://asia.ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG00000022898
i <- 1
print(interest_symbol[rename_gene[i],])
ind_A[rename_gene[i]] <- pmatch("Dscr3",rownames(dataA))

# "5630401D24Rik" --> "Eef1aknmt"
# referring to https://asia.ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG00000026694
i <- 2
print(interest_symbol[rename_gene[i],])
ind_A[rename_gene[i]] <- pmatch("Mettl13",rownames(dataA))

# "Ispd" --> "Crppa"
# referring to https://asia.ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG00000043153
i <- 3
print(interest_symbol[rename_gene[i],])
ind_A[rename_gene[i]] <- pmatch("Ispd",rownames(dataA))

# "4933416C03Rik" --> "Taf7l2"
# referring to https://asia.ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG00000074734
i <- 4
print(interest_symbol[rename_gene[i],])
ind_A[rename_gene[i]] <- pmatch("4933416C03Rik",rownames(dataA))

# "Ssfa2" --> "Itprid2"
# referring to https://asia.ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG00000027007
i <- 5
print(interest_symbol[rename_gene[i],])
ind_A[rename_gene[i]] <- pmatch("Ssfa2",rownames(dataA))

# "Tmem206" --> "Pacc1"
# referring to https://asia.ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG00000026627
i <- 6
print(interest_symbol[rename_gene[i],])
ind_A[rename_gene[i]] <- pmatch("Tmem206",rownames(dataA))

# "Fam198a" --> "Gask1a"
# referring to https://asia.ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG00000038233
i <- 7
print(interest_symbol[rename_gene[i],])
ind_A[rename_gene[i]] <- pmatch("Fam198a",rownames(dataA))

# "Narfl" --> "Ciao3"
# referring to https://asia.ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG00000002280
i <- 8
print(interest_symbol[rename_gene[i],])
ind_A[rename_gene[i]] <- pmatch("Narfl",rownames(dataA))

# "Trove2" --> "Ro60"
# referring to https://asia.ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG00000018199
i <- 9
print(interest_symbol[rename_gene[i],])
ind_A[rename_gene[i]] <- pmatch("Trove2",rownames(dataA))

# "Mut" --> "Mmut"
# referring to https://asia.ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG00000023921
i <- 10
print(interest_symbol[rename_gene[i],])
ind_A[rename_gene[i]] <- pmatch("Mut",rownames(dataA))

# "D3Ertd254e" --> "Zfp267"
# referring to https://asia.ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG00000033883
i <- 11
print(interest_symbol[rename_gene[i],])
ind_A[rename_gene[i]] <- pmatch("D3Ertd254e",rownames(dataA))

# "Mum1" --> "Pwwp3a"
# referring to https://asia.ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG00000020156
i <- 12
print(interest_symbol[rename_gene[i],])
ind_A[rename_gene[i]] <- pmatch("Mum1",rownames(dataA))

# "Fbxo18" --> "Fbh1"
# referring to https://asia.ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG00000058594
i <- 13
print(interest_symbol[rename_gene[i],])
ind_A[rename_gene[i]] <- pmatch("Fbxo18",rownames(dataA))

###########################
# Run preprocessing again #
###########################
# Download and read the counts, metadata of Nestorowa et al. 2016
fname <- "../RawData/GSE81682_HTSeq_counts.txt.gz"
# if (!file.exists(fname)) { download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81682&format=file&file=GSE81682%5FHTSeq%5Fcounts%2Etxt%2Egz", fname) }
dataF <- read.table(fname, header=TRUE, row.names=1, check.names=FALSE)
dataF <- as.matrix(dataF)
dim(dataF)

fname <- "../RawData/metaF.txt"
if (!file.exists(fname)) { download.file("http://blood.stemcells.cam.ac.uk/data/all_cell_types.txt", fname) }
metaF <- read.table(fname, stringsAsFactors = FALSE, header=TRUE, check.names=FALSE)
metainds <- match(colnames(dataF), rownames(metaF))
missing.meta <- is.na(metainds)
metaF <- metaF[metainds,] # This will contain NA's... which is okay, at this point, to preserve length.

# Defining the cell type based on the metadata.
metatypeF <- rep("other", nrow(metaF))
for (col in rev(colnames(metaF))) { # reverse, so earlier columns end up overwriting later ones.
  chosen <- metaF[,col]==1
  metatypeF[chosen] <- sub("[0-9]?_.*", "", col)
}
metatypeF[metatypeF=="ESLAM"] <- "HSPC"

# Filling in metadata from the cell sorting label, if metadata was missing.
metatypeF[missing.meta] <- sub("_.*", "", colnames(dataF)[missing.meta])
metatypeF[metatypeF=="LT-HSC"] <- "LTHSC"
metatypeF[metatypeF=="Prog"] <- "other"
colnames(dataF)<-metatypeF

# Perform size factor normalization within this data set.
high.abF <- scater::calcAverage(dataF) > 1
clustF <- quickCluster(dataF, method="igraph", subset.row=high.abF)
sizeF <- computeSumFactors(dataF, cluster=clustF, subset.row=high.abF)
dataF_norm <- t(t(dataF)/sizeF)

# Download and read the counts and meta data of Paul et al. 2015
fname <- "../RawData/GSE72857_umitab.txt.gz"
# if (!file.exists(fname)) { download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE72857&format=file&file=GSE72857%5Fumitab%2Etxt%2Egz", fname) }
dataA <- read.table(fname, header=TRUE, row.names=1)
metaA <- read.csv2("../RawData/MAP.csv",sep=",",stringsAsFactors = FALSE, head=TRUE, row.names=1)
dim(dataA)

# Only selecting cells that are in the metadata.
metainds <- match(rownames(metaA), colnames(dataA))
dataA <- dataA[,metainds]
dataA <- as.matrix(dataA)

# Organizing cell type labels.
metatypeA <- character(nrow(metaA))
metatypeA[metaA[,1]<7] <- "ERY"
metatypeA[metaA[,1]>6 & metaA[,1]<12] <- "CMP"
metatypeA[metaA[,1]>11] <- "GMP"
colnames(dataA) <- metatypeA

# Perform size factor normalization within this data set.
high.abA <- scater::calcAverage(dataA) > 1
clustA <- quickCluster(dataA, method="igraph", subset.row=high.abA)
sizeA <- computeSumFactors(dataA, cluster=clustA, subset.row=high.abA)
dataA_norm <- t(t(dataA)/sizeA)

matched_F <- match(gene_list, rownames(dataF_norm))
dataA_norm2 <- dataA_norm[ind_A,]
dataF_norm2 <- dataF_norm[matched_F,]
rownames(dataA_norm2) <- rownames(dataF_norm2)

# Rescaling the first dataset to match the coverage of the second.
aveA <- rowMeans(dataA_norm2)
aveF <- rowMeans(dataF_norm2)
dataF_norm2 <- dataF_norm2 * median(aveA/aveF)

# Perform log-transformation and save results to file.
logDataF_norm <- log(1 + dataF_norm2)
logDataA_norm <- log(1 + dataA_norm2)

colnames(logDataF_norm)[colnames(logDataF_norm)=="other"] <- "Unsorted"
colnames(logDataA_norm)[colnames(logDataA_norm)=="ERY"] <- "MEP"

############
# Coloring #
############
raw.all <- cbind(logDataF_norm, logDataA_norm)
first.batch <- rep(c(TRUE, FALSE), c(ncol(logDataF_norm), ncol(logDataA_norm)))

# Adding colours.
#base.color <- "grey"
color.legendF <- c(MEP="orange", GMP="chartreuse4", CMP="magenta", 
                   HSPC="cyan", LTHSC="dodgerblue", MPP="blue", LMPP="light blue", Unsorted="grey")
colmatF <- col2rgb(color.legendF) 

colmatA <- colmatF + 100 # A lighter shade.
colmatA[colmatA > 255] <- 255

#colmatF<-colmatF+200
#colmatF[colmatF > 255] <- 255
color.legendA <- setNames(rgb(colmatA[1,], colmatA[2,], colmatA[3,], maxColorValue=255), names(color.legendF))
allcolors <- c(color.legendF[colnames(logDataF_norm)], color.legendA[colnames(logDataA_norm)])
batch<-c( rep(1,ncol(logDataF_norm)),rep(2,ncol(logDataA_norm)) )

# Only keeping common cell types for PCA.
celltypes <- c(colnames(logDataF_norm), colnames(logDataA_norm))
pca.retain <- celltypes %in% c("MEP", "GMP", "CMP") 

# Making a plotting function.
plotFUN <- function(fname, Y, subset=NULL, ..., xlab="tSNE 1",ylab="tSNE 2",main="") {
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  png(fname,width=900,height=700)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=ifelse(first.batch, 21, 1)[subset], 
       col=ifelse(first.batch, "black", allcolors)[subset],
       bg=allcolors[subset], xlab=xlab, ylab=ylab, main=main) 
  dev.off()
}

batchcolor=c("lightcoral","lavender")
plotFUNb <- function(fname, Y, subset=NULL, ...) {
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  png(fname,width=900,height=700)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=ifelse(first.batch, 21, 1)[subset], 
       col=ifelse(first.batch, "black", batchcolor[batch[subset]]),
       bg=batchcolor[batch[subset]])#,  xlab="tSNE 1",ylab="tSNE 2")
  dev.off()
}


##########################
# generating a tSNE plot #
##########################
X.unc <- raw.all
t.unc <- t(X.unc)

## Generating a t-SNE plot.
require(Rtsne)
set.seed(0)
all.dists.unc <- as.matrix(dist(t.unc))
tsne.unc <- Rtsne(all.dists.unc, is_distance=TRUE, perplexity = 90)
plotFUN("Image/uncFA_normalized.png", tsne.unc$Y, main="Uncorrected",  xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("Image/uncFAb_normalized.png", tsne.unc$Y, main="Uncorrected",  xlab="tSNE 1",ylab="tSNE 2")

#rm(all.dists.unc)
gc()

#########################
# Generating a PCA plot #
#########################
pca.unc <- prcomp(t.unc[pca.retain,], rank=2)
# pca.unc$x[ (pca.unc$x<(-0.08))]<- (-0.08)
plotFUN("Image/pca_uncFA_normalized2.png", pca.unc$x, subset=pca.retain, main="Uncorrected", ylim=c(-0.1, max(pca.unc$x[,2])),  xlab="PC 1",ylab="PC 2")

pca.unc.all <- prcomp(t.unc, rank=2)
# pca.unc$x[ (pca.unc$x<(-0.08))]<- (-0.08)
plotFUN("Image/pca_uncFA_all_normalized.png", pca.unc.all$x, main="Uncorrected", ylim=c(-0.1, max(pca.unc$x[,2])),  xlab="PC 1",ylab="PC 2")


#rm(t.unc)
gc()

########################
# Plot by our coloring #
########################
# t-SNE plots
tsne.unc.switch <- tsne.unc$Y
tsne.unc.switch[1:nb[1],] <- tsne.unc$Y[nb[2]+ 1:nb[1],]
tsne.unc.switch[1:nb[2] + nb[1],] <- tsne.unc$Y[1:nb[2],]
MNN_by_celltype<- "Image/tsne_hemat_uncorrected_by_celltype_normalized.jpg"
dat_frame <- data.frame(Var1 = tsne.unc.switch[,1], 
                        Var2 = tsne.unc.switch[,2], 
                        col = celltype_factor)

jpeg(MNN_by_celltype,width = 800, height = 600, quality = 100)
ggplot(dat_frame, aes(x= Var1, y= Var2, colour= col)) +
  geom_point(size=4) + theme_classic() +
  scale_colour_manual(values = alpha(color_by_celltype,0.6)) +
  xlab("tSNE 1") + ylab("tSNE 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 44), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 44),#, angle = 45))
        axis.title=element_text(size=48,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        #legend.position = "none")
        legend.text=element_text(size=28,face="bold"),
        legend.title = element_text(size = 36,face="bold"))
dev.off()

MNN_by_batch <- "Image/tsne_hemat_uncorrected_by_batch_normalized.jpg"
dat_frame <- data.frame(Var1 = tsne.unc.switch[,1], 
                        Var2 = tsne.unc.switch[,2], 
                        col = batch_factor)

jpeg(MNN_by_batch,width = 800, height = 600, quality = 100)
ggplot(dat_frame, aes(x= Var1, y= Var2, colour= col)) +
  geom_point(size=4) + theme_classic() +
  scale_colour_manual(values = alpha(color_by_batch,0.6)) +
  xlab("tSNE 1") + ylab("tSNE 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 36), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 36),#, angle = 45))
        axis.title=element_text(size=44,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        #legend.position = "none")
        legend.text=element_text(size=28,face="bold"),
        legend.title = element_text(size = 36,face="bold"))
dev.off()

# PCA plots
nbpart <- sum(celltype_factor %in% c("CMP","MEP","GMP")) -nb[1]
pca.unc.switchpart <- pca.unc$x
pca.unc.switchpart[1:nb[1],] <- pca.unc$x[1:nb[1] + nbpart,]
pca.unc.switchpart[1:nbpart + nb[1],] <- pca.unc$x[1:nbpart,]

unc_PCA_part<- "Image/pca_part_hemat_uncorrected_normalized.jpg"
dat_frame <- data.frame(Var1 = pca.unc.switchpart[,1], 
                        Var2 = pca.unc.switchpart[,2], 
                        col = celltype_factor[celltype_factor %in% c("CMP","MEP","GMP")])

jpeg(unc_PCA_part,width = 900, height = 700, quality = 100)
ggplot(dat_frame, aes(x= Var1, y= Var2, colour= col)) +
  geom_point(size=4) + theme_classic() +
  scale_colour_manual(values = alpha(color_by_celltype[4:6],0.6)) +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 36), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 36),#, angle = 45))
        axis.title=element_text(size=44,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        #legend.position = "none")
        legend.text=element_text(size=28,face="bold"),
        legend.title = element_text(size = 36,face="bold"))
dev.off()

pca.unc.swtich <- pca.unc.all$x
pca.unc.swtich[1:nb[1],] <- pca.unc.all$x[1:nb[1] + nb[2],]
pca.unc.swtich[1:nb[2] + nb[1],] <- pca.unc.all$x[1:nb[2],]
unc_PCA<- "Image/pca_hemat_uncorrected_normalized.jpg"
dat_frame <- data.frame(Var1 = pca.unc.swtich[,1], 
                        Var2 = pca.unc.swtich[,2], 
                        col = celltype_factor)

jpeg(unc_PCA,width = 900, height = 700, quality = 100)
ggplot(dat_frame, aes(x= Var1, y= Var2, colour= col)) +
  geom_point(size=4) + theme_classic() +
  scale_colour_manual(values = alpha(color_by_celltype,0.6)) +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 36), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 36),#, angle = 45))
        axis.title=element_text(size=44,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        #legend.position = "none")
        legend.text=element_text(size=28,face="bold"),
        legend.title = element_text(size = 36,face="bold"))
dev.off()


save.image("Reproduce_MNN.RData")