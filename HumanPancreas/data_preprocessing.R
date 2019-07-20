# Referring to the source code https://github.com/MarioniLab/MNN2017/tree/master/Pancreas

# This script prepares data for the pancreas analysis.
# It involves four publicly available datasets.
library(scran)
library(biomaRt)
library(limSolve)
library(scater)
library(SingleCellExperiment)
library(cowplot)

##########################################
##########################################
rm(list=ls())

# Setting the working directory
setwd("F:/scRNA/code/0601Cpp_BUSseq/BUSseq_implementation_v1/HumanPancreas")

if(!dir.exists("Data")){
  dir.create("Data")
}

##############
## GSE81076 ##
##############
# Download file from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81076
gse81076 <- './RawData/GSE81076_D2_3_7_10_17.txt.gz'
if (!file.exists(gse81076)) { download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81076/suppl/GSE81076_D2_3_7_10_17.txt.gz", gse81076) }
gse81076.df <- read.table(gse81076, sep='\t', h=T, stringsAsFactors=F)
gse81076.ndim <- dim(gse81076.df)[2]

# construct the meta data from the cell names
donor.names <- unlist(regmatches(colnames(gse81076.df)[2:gse81076.ndim],
                                 gregexpr(pattern="D[0-9]{1,2}", 
                                          colnames(gse81076.df)[2:gse81076.ndim])))

plate.id <- unlist(lapply(strsplit(unlist(lapply(strsplit(head(colnames(gse81076.df)[2:gse81076.ndim]),
                                                          split="D[0-9]{1,2}", perl=T), 
                                                 FUN=function(x) paste0(x[2]))),
                                   fixed=T, split="_"),
                          FUN=function(c) paste0(c[1])))

protocol.id <- rep('CELseq', gse81076.ndim-1)
study.id <- rep('GSE81076', gse81076.ndim-1)
gse81076.meta <- data.frame(list('Donor' = donor.names,
                                 'Plate' = plate.id,
                                 'Protocol' = protocol.id,
                                 'Study' = study.id,
                                 'Sample' = colnames(gse81076.df)[2:gse81076.ndim]))
rownames(gse81076.meta) <- colnames(gse81076.df)[2:gse81076.ndim]
colnames(gse81076.df) <- gsub(colnames(gse81076.df), pattern='X', replacement='gene')

# remove superfluous suffixes from gene IDs
gse81076.df$gene <- gsub(gse81076.df$gene,
                         pattern="__chr[0-9]+", replacement="")
# remove the duplicated gene names
gse81076.df <- gse81076.df[!duplicated(gse81076.df$gene), ]

rownames(gse81076.df) <- gse81076.df$gene

# remove the gene ID column for downstream normalization
gse81076.df <- gse81076.df[, -1]

# remove cells and genes with all 0's
gene_sparsity <- (apply(gse81076.df == 0, MARGIN = 1, sum)/dim(gse81076.df)[2])
keep_genes <- gene_sparsity < 0.9
gse81076.nz <- gse81076.df[keep_genes, ]

cell_sparsity <- apply(gse81076.nz == 0, MARGIN = 2, sum)/dim(gse81076.nz)[1]
keep_cells <- cell_sparsity < 0.8
gse81076.nz <- gse81076.nz[, keep_cells]
gse81076.nz <- apply(gse81076.nz, 2, as.integer)

# use the spike in genes to estimate size factors for normalization
# all values show be non-negative integers
spikes <- grepl(rownames(gse81076.df[keep_genes, ]),
                pattern='ERCC')

sce <- SingleCellExperiment(list(counts = as.matrix(gse81076.nz)))
sce <- calculateQCMetrics(sce, feature_controls=list(Spikes=spikes))

clusters <- quickCluster(sce, get.spikes=TRUE, min.size=120)
sce <- computeSumFactors(sce, sizes=c(10, 20, 40, 60), positive=T,
                         assay.type='counts', clusters=clusters)
summary(sizeFactors(sce))
sce <- normalize(sce)

gse81076.norm <- data.frame(exprs(sce))

gse81076.norm$gene_id <- rownames(gse81076.df[keep_genes, ])
gse81076.norm$gene_id <- gsub(gse81076.norm$gene_id,
                              pattern="__chr[0-9X]+", replacement="")

write.table(gse81076.nz, sep='\t',
            file='./Data/GSE81076_readcount.tsv',
            quote=FALSE, row.names=FALSE, col.names=TRUE)

write.table(gse81076.norm, sep='\t',
            file='./Data/GSE81076_SFnorm.tsv',
            quote=FALSE, row.names=FALSE, col.names=TRUE)

write.table(gse81076.meta, sep="\t",
            file="./Data/GSE81076_metadata.tsv",
            quote=FALSE, row.names=F, col.names=TRUE)

##############
## GSE85241 ##
##############
# clear environment and invoke garbage collector
rm(list=ls())
gc()

# Download file from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85241
gse85241 <- './RawData/GSE85241_cellsystems_dataset_4donors_updated.csv'
if (!file.exists(gse85241)) { download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85241/suppl/GSE85241_cellsystems_dataset_4donors_updated.csv.gz", gse85241) }
gse85241.df <- read.table(gse85241, sep='\t', h=T, stringsAsFactors=F)
gse85241.df$gene_id <- rownames(gse85241.df)

# gene IDs are located in column X for these data
colnames(gse85241.df) <- gsub(colnames(gse85241.df), pattern="X",
                              replacement="gene_id")

donor.id <-unlist(lapply(strsplit(colnames(gse85241.df)[1:(dim(gse85241.df)[2]-1)],
                                  fixed=T, split="."), 
                         FUN=function(x) paste0(x[1])))

plate.id <- unlist(lapply(strsplit(colnames(gse85241.df)[1:(dim(gse85241.df)[2]-1)],
                                   fixed=T, split="."),
                          FUN=function(x) paste0(x[2])))

protocol.id <- rep('CELseq2', dim(gse85241.df)[2]-1)
study.id <- rep('GSE85241', dim(gse85241.df)[2]-1)

gse85241.meta <- data.frame(list('Donor' = donor.id,
                                 'Plate'= plate.id,
                                 'Protocol' = protocol.id,
                                 'Study' = study.id,
                                 'Sample' = colnames(gse85241.df)[1:(dim(gse85241.df)[2]-1)]))

rownames(gse85241.meta) <- colnames(gse85241.df)[1:(dim(gse85241.df)[2]-1)]

write.table(gse85241.meta,
            file="./Data/GSE85241_metadata.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

# set gene IDs as rownames, remove gene ID column
# remove duplicated gene IDs

gse85241.df$gene_id <- gsub(gse85241.df$gene_id,
                            pattern="__chr[0-9X]+", replacement="")
gse85241.df <- gse85241.df[!duplicated(gse85241.df$gene_id), ]
rownames(gse85241.df) <- gse85241.df$gene_id
gse85241.df <- gse85241.df[, 1:(dim(gse85241.df)[2]-1)]

# remove cells and genes with all 0's
gene_sparsity <- (apply(gse85241.df == 0, MARGIN = 1, sum)/dim(gse85241.df)[2])
keep_genes <- gene_sparsity < 0.9
gse85241.nz <- gse85241.df[keep_genes, ]

cell_sparsity <- apply(gse85241.nz == 0, MARGIN = 2, sum)/dim(gse85241.nz)[1]
keep_cells <- cell_sparsity < 0.8
dim(gse85241.nz[, keep_cells])
gse85241.nz <- gse85241.nz[, keep_cells]
gse85241.nz <- apply(gse85241.nz, 2, as.integer)

spikes <- grepl(rownames(gse85241.df[keep_genes, ]),
                pattern='ERCC')
sce <- SingleCellExperiment(list(counts = as.matrix(gse85241.nz)))
sce <- calculateQCMetrics(sce, feature_controls=list(Spikes=spikes))

clusters <- quickCluster(sce, get.spikes=TRUE, min.size=120)
sce <- computeSumFactors(sce, sizes=c(10, 20, 40, 60), positive=T,
                         assay.type='counts', clusters=clusters)
summary(sizeFactors(sce))

sce <- normalize(sce)
gse85241.norm <- data.frame(exprs(sce))
gse85241.norm$gene_id <- rownames(gse85241.df[keep_genes, ])

write.table(gse85241.nz, sep='\t',
            file='./Data/GSE85241_readcount.tsv',
            quote=FALSE, row.names=FALSE, col.names=TRUE)

write.table(gse85241.norm, sep='\t',
            file='./Data/GSE85241_SFnorm.tsv',
            quote=F, row.names=F, col.names=T)

##############
## GSE86473 ##
##############
# clear environment and invoke garbage collector
rm(list=ls())
gc()

# the raw/processed counts table was not available for download from GEO for this data set
# these data were generated by mapping the original fastq's back to mm10, then using featureCounts to quantify
# against mm10 ensembl build 86
# data are contained in ./RawData for each cell type
alpha <- './RawData/alpha-feature_counts.tsv.gz'
beta <- './RawData/beta-feature_counts.tsv.gz'
delta <- './RawData/delta-feature_counts.tsv.gz'
pp <- './RawData/PP-feature_counts.tsv.gz'
# qc_out <- './RawData/pancreas-smarseq2-qcout.tsv'

alpha.df <- read.table(alpha, sep='\t', h=T, stringsAsFactors=F)
beta.df <- read.table(beta, sep='\t', h=T, stringsAsFactors=F)
delta.df <- read.table(delta, sep='\t', h=T, stringsAsFactors=F)
pp.df <- read.table(pp, sep='\t', h=T, stringsAsFactors=F)

# qc.df <- read.table(qc_out, sep='\t', h=F, stringsAsFactors=F)
data.list <- list("alpha"=alpha.df, "beta"=beta.df, "delta"=delta.df,
                  "gamma"=pp.df)

gse86473.df <- Reduce(x=data.list,
                      f=function(x, y) merge(x, y, by='gene_id'))

# remove .dedup suffix
samp.names <- unlist(lapply(strsplit(colnames(gse86473.df), fixed=T,
                                     split="."), 
                            FUN=function(x) paste0(x[1])))
colnames(gse86473.df) <- tolower(samp.names)

gse86473.meta <- read.table('./RawData/GSE86473_experimental_design.tsv',
                            sep='\t', h=T, stringsAsFactors=F)
gse86473.meta$Sample <- tolower(gse86473.meta$Sample)
gse86473.meta$Study <- "GSE86473"

# capitalize first letter of cell labels

gse86473.meta$CellType <- paste(toupper(substr(gse86473.meta$CellType, 1, 1)),
                                substr(gse86473.meta$CellType, 2, nchar(gse86473.meta$CellType)), sep="")

write.table(gse86473.meta,
            "./Data/GSE86473_metadata.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

rownames(gse86473.df) <- gse86473.df$gene_id

# need to map ensemlb gene ids to hgnc symbols to match up with other Pancreas data sets
ensembl <- useEnsembl(biomart='ensembl', dataset='hsapiens_gene_ensembl', GRCh=37)
gene_symbol <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                     filters='ensembl_gene_id', mart=ensembl,
                     values=gse86473.df$gene_id)
dat <- merge(gse86473.df, gene_symbol,
             by.x='gene_id', by.y='ensembl_gene_id')
gene_symbols <- dat$external_gene_name
gse86473.df <- as.data.frame(append(dat[, -1],
                                    list(gene_id=gene_symbols),
                                    after=0))

# unset hgnc_symbol as a factor and remove duplicated IDs
gse86473.df$gene_id <- as.character(gse86473.df$gene_id)
gse86473.df <- gse86473.df[!duplicated(gse86473.df$gene_id), ]

# set the hgnc symbols as rownames and size factor normalize counts table
rownames(gse86473.df) <- gse86473.df$gene_id
gse86473.df <- gse86473.df[, 2:(dim(gse86473.df)[2]-1)]

# remove cells and genes with all 0's
gene_sparsity <- (apply(gse86473.df == 0, MARGIN = 1, sum)/dim(gse86473.df)[2])
keep_genes <- gene_sparsity < 0.9
dim(gse86473.df[keep_genes, ])
gse86473.nz <- gse86473.df[keep_genes, ]

cell_sparsity <- apply(gse86473.nz == 0, MARGIN = 2, sum)/dim(gse86473.nz)[1]
keep_cells <- cell_sparsity < 0.8
dim(gse86473.nz[, keep_cells])
gse86473.nz <- gse86473.nz[, keep_cells]
gse86473.nz <- apply(gse86473.nz, 2, as.integer)

sce <- SingleCellExperiment(list(counts = as.matrix(gse86473.nz)))
sce <- calculateQCMetrics(sce)
clusters <- quickCluster(sce, min.size=120)
sce <- computeSumFactors(sce, sizes=c(10, 20, 40, 60), positive=T,
                         assay.type='counts', clusters=clusters)
summary(sizeFactors((sce)))

sce <- normalize(sce)
gse86473.norm <- data.frame(exprs(sce))
gse86473.norm$gene_id <- rownames(gse86473.df[keep_genes, ])

write.table(gse86473.nz, sep='\t',
            file='./Data/GSE86473_readcount.tsv',
            quote=FALSE, row.names=FALSE, col.names=TRUE)

write.table(gse86473.norm, sep='\t',
            file='./Data/GSE86473_SFnorm.tsv',
            quote=F, row.names=F, col.names=T)

#################
## E-MTAB-5061 ##
#################
# clear environment and invoke garbage collector
rm(list=ls())
gc()

# the download file contains columns of RPKM & counts
# need to pull out just the integer gene count columns
# Please refer to https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5061/ 

# need to unzip
emtab5061 <- "./RawData/pancreas_refseq_rpkms_counts_3514sc.txt"
if(!file.exists(emtab5061)) {
  emtab_combined <- "./RawData/EMTAB5061_rpkm_counts.txt.zip"
  if(!file.exists(emtab_combined)){ download.file("https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5061/files/E-MTAB-5061.processed.1.zip",
                                                  emtab_combined)}
  unzip(emtab_combined, exdir = "./RawData")
}

emtab.df <- read.table(emtab5061,
                       h=FALSE, sep="\t", stringsAsFactors=F)

col.names <- unlist(read.table(emtab5061,
                               h=FALSE, sep="\t", stringsAsFactors=F, comment.char="", nrows = 1))

# first 2 columns are gene symbol and NCBI ID
# because of the way that this table is constructed, there are two cells for each column, but they
# are not contiguous.  Therefore, the sample names need to be read in separately to the counts/rpkm table

emtab5061.df <- emtab.df[, c(1, 3517:dim(emtab.df)[2])]
colnames(emtab5061.df) <- gsub(col.names, pattern="#samples", replacement="gene_id")

emtab.file <- "./RawData/E-MTAB-5061.sdrf.txt"

emtab.sdrf <- read.table(emtab.file,
                         h=TRUE, sep="\t", stringsAsFactors=FALSE)

# construct the appropriate meta data columns, i.e. donor, plate, protocol, study
emtab.meta <- emtab.sdrf[, c("Assay.Name", "Characteristics.cell.type.", "Characteristics.individual.")]
colnames(emtab.meta) <- c("Sample", "CellType", "Donor")
emtab.meta$Study <- "E-MTAB-5061"
emtab.meta$Protocol <- "SmartSeq2"

# remove the marked low quality cells
remove.cells <- unique(emtab.sdrf$Assay.Name[emtab.sdrf$Characteristics.single.cell.well.quality. == "low quality cell"])
emtab5061.df <- emtab5061.df[, !colnames(emtab5061.df) %in% remove.cells]
emtab.meta <- emtab.meta[!emtab.meta$Sample %in% remove.cells, ]
rownames(emtab.meta) <- emtab.meta$Sample

emtab.meta$CellType <- gsub(emtab.meta$CellType,
                            pattern=" cell", replacement="")

emtab.meta$CellType <- paste(toupper(substr(emtab.meta$CellType, 1, 1)),
                             substr(emtab.meta$CellType, 2, nchar(emtab.meta$CellType)), sep="")

write.table(emtab.meta,
            file="./Data/E-MTAB-5061_metadata.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

# remove duplicated gene IDs and set to rownames
emtab5061.df <- emtab5061.df[!duplicated(emtab5061.df$gene_id), ]
rownames(emtab5061.df) <- emtab5061.df$gene_id

emtab5061.df <- emtab5061.df[, -1]

# remove cells and genes with all 0's
gene_sparsity <- (apply(emtab5061.df == 0, MARGIN = 1, sum)/dim(emtab5061.df)[2])
keep_genes <- gene_sparsity < 0.9
dim(emtab5061.df[keep_genes, ])
emtab5061.nz <- emtab5061.df[keep_genes, ]

cell_sparsity <- apply(emtab5061.nz == 0, MARGIN = 2, sum)/dim(emtab5061.nz)[1]
keep_cells <- cell_sparsity < 0.8
dim(emtab5061.nz[, keep_cells])
emtab5061.nz <- emtab5061.nz[, keep_cells]
emtab5061.nz <- apply(emtab5061.nz, 2, as.integer)

spikes <- grepl(x=rownames(emtab5061.df[keep_genes, ]), pattern="ERCC")
sce <- SingleCellExperiment(list(counts = as.matrix(emtab5061.nz)))
sce <- calculateQCMetrics(sce, feature_controls=list(Spikes=spikes))

clusters <- quickCluster(sce, min.size=120)
sce <- computeSumFactors(sce, sizes=c(10, 20, 40, 60), positive=T,
                         assay.type='counts', clusters=clusters)
summary(sizeFactors((sce)))
sce <- normalize(sce)
emtab.norm <- data.frame(exprs(sce))
emtab.norm$gene_id <- rownames(emtab5061.df[keep_genes, ])

write.table(emtab5061.nz, sep='\t',
            file='./Data/E-MTAB-5061_readcount.tsv',
            quote=FALSE, row.names=FALSE, col.names=TRUE)

write.table(emtab.norm,
            file="./Data/E-MTAB-5061_SFnorm.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

# clear the R environment in case this script is directly sourced
rm(list=ls())
gc()


############################
#Find Highly Variable Genes#
############################
find_hvg <- function(dataframe, plot=FALSE, p.threshold=1e-2, return.ranks=FALSE,
                     return.p=FALSE){
  # define a set of highly variable gene for the GFP+ and GFP- separately
  require(MASS)
  require(limSolve)
  require(statmod)
  # assume gene names are in the rows
  # even if they aren't this will still get
  # the input row ordering
  gene.names <- rownames(dataframe)
  means <- rowMeans(dataframe, na.rm = T)
  vars <- apply(dataframe, 1, var, na.rm=T)
  cv2 <- vars/(means^2)
  
  minMeanForFit <- unname(quantile(means[which(cv2 > 0.2)], 0.8))
  
  # select genes with mean value greater than min value for fitting
  # remove values with 1/means == infinite
  recip.means <- 1/means
  recip.means[is.infinite(recip.means)] <- 0
  
  useForFit <- recip.means <= 0.1
  
  # fit with a gamma-distributed GLM
  fit <- glmgam.fit(cbind(a0 = 1, a1tilde=recip.means[!useForFit]), cv2[!useForFit])
  
  # calculate % variance explained by the model fit
  resid.var <- var(fitted.values(fit) - cv2[!useForFit])
  total.var <- var(cv2[!useForFit])
  
  # get fitted values and mean-dispersion dependence line
  a0 <- unname(fit$coefficients["a0"])
  a1 <- unname(fit$coefficients["a1tilde"])
  
  xg <- seq(0, max(means[means != Inf]), length.out=100000)
  vfit <- (a1/xg) + a0
  
  # add confidence intervals
  d.f <- ncol(dataframe) - 1
  
  # rank genes by the significance of their deviation from the fit
  # to call HVGs
  a.fit <- (a1/means) + a0
  varFitRatio <- vars/(a.fit * means^2)
  varOrder <- order(varFitRatio, decreasing=T)
  
  oed <- dataframe[varOrder, ]
  
  if(plot == TRUE){
    smoothScatter(means, cv2, xlab="Mean expression", ylab=expression("CV"^2))
    lines(xg, vfit, col="black", lwd=3 )
    lines(xg, vfit * qchisq(0.975, d.f)/d.f, lty=2, col="black")
    lines(xg, vfit * qchisq(0.025, d.f)/d.f,lty=2,col="black")
    # display the 100 most highly variable genes
    points(means[varOrder[1:100]], cv2[varOrder[1:100]], col='red')
  }
  
  pvals <- pchisq(varFitRatio * d.f, d.f, lower.tail = F)
  pvals[is.na(pvals)] <- 1.0
  adj.pvals <- p.adjust(pvals, method='fdr')
  HVG <- adj.pvals <= p.threshold
  
  if(return.ranks){
    # order p-values, then subset past a threshold
    rank.p <- adj.pvals[order(adj.pvals, decreasing=FALSE)]
    order.names <- gene.names[order(adj.pvals, decreasing=FALSE)]
    thr.p <- rank.p <= p.threshold
    HVG <- order.names[thr.p]
  }
  
  if(return.p){
    HVG <- cbind(HVG, adj.pvals)
  }  
  
  return(HVG)
}

gse81076.norm <- read.table("./Data/GSE81076_SFnorm.tsv",
                            h=TRUE, sep="\t", stringsAsFactors=FALSE)

# remove gene id column for hvg discovery, last column
rownames(gse81076.norm) <- gse81076.norm$gene_id

# output is a data frame of the same order as the input matrix
# need to output the p-value from HVG gene selection for
# combination by meta-analysis later
gse81076.HVG <- find_hvg(dataframe=gse81076.norm[, 1:(dim(gse81076.norm)[2]-1)], 
                         p.threshold=1e-2, plot=FALSE, return.ranks=FALSE,
                         return.p=TRUE)

# select the highly variable genes from the input dataframe column 'gene id'
gse81076.hvg_df <- as.data.frame(gse81076.HVG)
gse81076.hvg_df$gene_id <- rownames(gse81076.hvg_df)
colnames(gse81076.hvg_df) <- c("HVG", "pval", "gene_id")

write.table(gse81076.hvg_df,
            file="./Data/GSE81076-HVG.tsv",
            row.names=FALSE, quote=FALSE, sep="\t")

##############
## GSE85241 ##
##############
gse85241.norm <- read.table("./Data/GSE85241_SFnorm.tsv",
                            h=TRUE, sep="\t", stringsAsFactors=FALSE)

# remove gene id column for hvg discovery, last column
rownames(gse85241.norm) <- gse85241.norm$gene_id


# output is a boolean vector of the same order as the input matrix
gse85241.HVG <- find_hvg(dataframe=gse85241.norm[, 1:(dim(gse85241.norm)[2]-1)], 
                         p.threshold=1e-2, plot=FALSE, return.ranks=FALSE,
                         return.p=TRUE)

# select the highly variable genes from the input dataframe column 'gene id'
gse85241.hvg_df <- as.data.frame(gse85241.HVG)
gse85241.hvg_df$gene_id <- rownames(gse85241.hvg_df)
colnames(gse85241.hvg_df) <- c("HVG", "pval", "gene_id")

write.table(gse85241.hvg_df,
            file="./Data/GSE85241-HVG.tsv",
            row.names=FALSE, quote=FALSE, sep="\t")

##############
## GSE86473 ##
##############
gse86473.norm <- read.table("./Data/GSE86473_SFnorm.tsv",
                            h=TRUE, sep="\t", stringsAsFactors=FALSE)

# remove gene id column for hvg discovery, last column
rownames(gse86473.norm) <- gse86473.norm$gene_id


# output is a boolean vector of the same order as the input matrix
gse86473.HVG <- find_hvg(dataframe=gse86473.norm[, 1:(dim(gse86473.norm)[2]-1)], 
                         p.threshold=1e-2, plot=FALSE, return.ranks=FALSE,
                         return.p=TRUE)

# select the highly variable genes from the input dataframe column 'gene id'
gse86473.hvg_df <- as.data.frame(gse86473.HVG)
gse86473.hvg_df$gene_id <- rownames(gse86473.hvg_df)
colnames(gse86473.hvg_df) <- c("HVG", "pval", "gene_id")

write.table(gse86473.hvg_df,
            file="./Data/GSE86473-HVG.tsv",
            row.names=FALSE, quote=FALSE, sep="\t")

#################
## E-MTAB-5061 ##
#################
emtab5061.norm <- read.table("./Data/E-MTAB-5061_SFnorm.tsv",
                             h=TRUE, sep="\t", stringsAsFactors=FALSE)

# remove gene id column for hvg discovery, last column
rownames(emtab5061.norm) <- emtab5061.norm$gene_id

# output is a boolean vector of the same order as the input matrix
emtab5061.HVG <- find_hvg(dataframe=emtab5061.norm[, 1:(dim(emtab5061.norm)[2]-1)], 
                          p.threshold=1e-2, plot=FALSE, return.ranks=FALSE,
                          return.p=TRUE)

# select the highly variable genes from the input dataframe column 'gene id'
emtab5061.hvg_df <- as.data.frame(emtab5061.HVG)
emtab5061.hvg_df$gene_id <- rownames(emtab5061.hvg_df)
colnames(emtab5061.hvg_df) <- c("HVG", "pval", "gene_id")

write.table(emtab5061.hvg_df,
            file="./Data/E-MTAB-5061-HVG.tsv",
            row.names=FALSE, quote=FALSE, sep="\t")

# clear the environment
rm(list=ls())
gc()

########################################################################
#assigns cell-type labels to cells based on the original implementation#
########################################################################
library(Rtsne)
library(cluster)
library(ggplot2)
############
# GSE81076 #
############
# read in the normalized expression data
gse81076.norm <- read.table("./Data/GSE81076_SFnorm.tsv",
                            h=TRUE, sep="\t", stringsAsFactors=FALSE)

# remove gene id column for hvg discovery, last column
rownames(gse81076.norm) <- gse81076.norm$gene_id

gse81076.hvg_df <- read.table("./Data/GSE81076-HVG.tsv", 
                              sep="\t", h=TRUE, stringsAsFactors=FALSE)

gse81076.hvg <- gse81076.norm[gse81076.norm$gene_id %in% gse81076.hvg_df$gene_id,
                              1:(dim(gse81076.norm)[2]-1)]

gse81076.meta <- read.table("./Data/GSE81076_metadata.tsv",
                            sep="\t", h=TRUE, stringsAsFactors=FALSE)

# get low dimensional embedding of genes by tSNE
gse81076.tsne <- Rtsne(t(gse81076.hvg), perplexity=30)
gse81076.map <- data.frame(gse81076.tsne$Y)
colnames(gse81076.map) <- c("Dim1", "Dim2")
gse81076.map$Sample <- colnames(gse81076.hvg)

# merge metadata and clustering info
gse81076.uber <- merge(gse81076.map, gse81076.meta,
                       by='Sample')

marker_genes <- rownames(gse81076.norm)[grepl(rownames(gse81076.norm), pattern="(^GCG)|(^KRT19)|(^INS)|
                                              |(^SST)|(^PPY)|(^PRSS1)|(^COL1A1)|(^GHRL)|(^ESAM)") ]

marker_exprs <- gse81076.norm[marker_genes, 1:(dim(gse81076.norm)[2]-1)]

set.seed(42)
kmed <- pam(t(gse81076.hvg), 9)
panc.kmed <- data.frame(cbind(kmed$clustering))
colnames(panc.kmed) <- "Kmediods"
panc.kmed$Sample <- rownames(panc.kmed)
panc.kmed$Kmediods <- as.factor(panc.kmed$Kmediods)
panc.meta <- merge(gse81076.uber, panc.kmed, by='Sample')

marker_df <- data.frame(t(marker_exprs))
marker_df$Sample <- rownames(marker_df)

# can we cluster the cells on just the marker genes?
# hierarchical clustering doesn't work so well...

mark.dim <- dim(marker_df)
marker.kmed <- pam(marker_df[, 1:mark.dim[2]-1], 9)
marker.cluster <- data.frame(cbind(marker.kmed$clustering))
colnames(marker.cluster) <- c("markClust")
marker.cluster$Sample <- rownames(marker.cluster)
marker.cluster$markClust <- as.factor(marker.cluster$markClust)

marker_merge <- merge(panc.meta, marker_df, by='Sample')
marker_uber <- merge(marker_merge, marker.cluster, by='Sample')

# uncomment the following plots to see the overlay of each marker gene
# on the tSNE.  This illustrates the cell type label inference.
# mediods 8 = alpha cells
ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=GCG,
                         shape=Kmediods)) +
  geom_point(size=2) + theme_classic() +
  scale_colour_continuous(low='white', high='red') +
  scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))

# mediod 7 = beta cells, INS highest
ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=INS,
                         shape=Kmediods)) +
  #geom_point(data = marker_merge[which(marker_merge$INS>5),], size=2) + 
  geom_point(size=2) + 
  theme_classic() +
  scale_colour_continuous(low='white', high='red') +
  scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))

# mediod 7 = delta cells, SST highest
ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=SST,
                         shape=Kmediods)) +
  geom_point(size=2) + theme_classic() +
  scale_colour_continuous(low='white', high='red') +
  scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))

# mediod 8 = Gamma cells, PPY highest
ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=PPY,
                         shape=Kmediods)) +
  geom_point(size=2) + theme_classic() +
  scale_colour_continuous(low='white', high='red') +
  scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))

# mediod 1, 2 & 6 = acinar cells, PRSS1 highest
ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=PRSS1,
                         shape=Kmediods)) +
  geom_point(size=2) + theme_classic() +
  scale_colour_continuous(low='white', high='red') +
  scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))

# mediod 3, 4 = ductal cells, KRT19 highest
ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=KRT19,
                         shape=Kmediods)) +
  geom_point(size=2) + theme_classic() +
  scale_colour_continuous(low='white', high='red') +
  scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))

# small pop of mediod 9 = mesenchyme cells, COL1A1 highest, same as PP cells
ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=COL1A1,
                         shape=Kmediods)) +
  geom_point(size=2) + theme_classic() +
  scale_colour_continuous(low='white', high='red') +
  scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))


marker_uber$CellType <- ""
marker_uber$CellType[marker_uber$Kmediods %in% c(3, 4, 5)] <- "Ductal"
marker_uber$CellType[marker_uber$Kmediods %in% c(8) | (marker_uber$Kmediods %in% c(5) & marker_uber$GCG >= 7)] <- "Alpha"
marker_uber$CellType[marker_uber$Kmediods %in% c(8) & marker_uber$PPY >= 5] <- "Gamma"
marker_uber$CellType[marker_uber$Kmediods %in% c(7) | (marker_uber$Kmediods %in% c(5) & marker_uber$INS >= 7)] <- "Beta"
marker_uber$CellType[marker_uber$Kmediods %in% c(5, 7) & marker_uber$SST >= 8] <- "Delta"
marker_uber$CellType[marker_uber$Kmediods %in% c(1, 2, 6)] <- "Acinar"
marker_uber$CellType[marker_uber$Kmediods %in% c(9) & marker_uber$COL1A1 >= 2] <- "Mesenchyme"

table(marker_uber$CellType, marker_uber$Kmediods)

# Draw violin plot to show the marker gene experssion levels in each identified cell types  
marker_gene_list <- c("GCG","INS","SST","PPY","PRSS1","KRT19","COL1A1")

# GCG
jpeg(filename = paste0("Image/Label_Celltype/GSE81076_",marker_gene_list[1],"_for_Alpha.jpeg"),width = 1600, height = 600,quality = 100)
p1 <- ggplot(marker_uber[which(marker_uber$CellType!=""),], aes(x = Kmediods, GCG)) + 
  geom_violin(aes(fill = Kmediods))
p2 <- ggplot(marker_uber[which(marker_uber$CellType!=""),], aes(x = CellType, GCG)) + 
  geom_violin(aes(fill = CellType))
plot_grid(p1,p2,nrow = 1)
dev.off()
# INS
jpeg(filename = paste0("Image/Label_Celltype/GSE81076_",marker_gene_list[2],"_for_Beta.jpeg"),width = 1600, height = 600,quality = 100)
p3 <- ggplot(marker_uber[which(marker_uber$CellType!=""),], aes(x = Kmediods, INS)) + 
  geom_violin(aes(fill = Kmediods))
p4 <- ggplot(marker_uber[which(marker_uber$CellType!=""),], aes(x = CellType, INS)) + 
  geom_violin(aes(fill = CellType))
plot_grid(p3,p4,nrow = 1)
dev.off()
# SST
jpeg(filename = paste0("Image/Label_Celltype/GSE81076_",marker_gene_list[3],"_for_Delta.jpeg"),width = 1600, height = 600,quality = 100)
p5 <- ggplot(marker_uber[which(marker_uber$CellType!=""),], aes(x = Kmediods, SST)) + 
  geom_violin(aes(fill = Kmediods))
p6 <- ggplot(marker_uber[which(marker_uber$CellType!=""),], aes(x = CellType, SST)) + 
  geom_violin(aes(fill = CellType))
plot_grid(p5,p6,nrow = 1)
dev.off()
# PPY
jpeg(filename = paste0("Image/Label_Celltype/GSE81076_",marker_gene_list[4],"_for_Gamma.jpeg"),width = 1600, height = 600,quality = 100)
p7 <- ggplot(marker_uber[which(marker_uber$CellType!=""),], aes(x = Kmediods, PPY)) + 
  geom_violin(aes(fill = Kmediods))
p8 <- ggplot(marker_uber[which(marker_uber$CellType!=""),], aes(x = CellType, PPY)) + 
  geom_violin(aes(fill = CellType))
plot_grid(p7,p8,nrow = 1)
dev.off()
# PRSS1
jpeg(filename = paste0("Image/Label_Celltype/GSE81076_",marker_gene_list[5],"_for_Acinar.jpeg"),width = 1600, height = 600,quality = 100)
p9 <- ggplot(marker_uber[which(marker_uber$CellType!=""),], aes(x = Kmediods, PRSS1)) + 
  geom_violin(aes(fill = Kmediods))
p10 <- ggplot(marker_uber[which(marker_uber$CellType!=""),], aes(x = CellType, PRSS1)) + 
  geom_violin(aes(fill = CellType))
plot_grid(p9,p10,nrow = 1)
dev.off()
# KRT19
jpeg(filename = paste0("Image/Label_Celltype/GSE81076_",marker_gene_list[6],"_for_Ductal.jpeg"),width = 1600, height = 600,quality = 100)
p11 <- ggplot(marker_uber[which(marker_uber$CellType!=""),], aes(x = Kmediods, KRT19)) + 
  geom_violin(aes(fill = Kmediods))
p12 <- ggplot(marker_uber[which(marker_uber$CellType!=""),], aes(x = CellType, KRT19)) + 
  geom_violin(aes(fill = CellType))
plot_grid(p11,p12,nrow = 1)
dev.off()
# COL1A1
jpeg(filename = paste0("Image/Label_Celltype/GSE81076_",marker_gene_list[7],"_for_Mesenchyme.jpeg"),width = 1600, height = 600,quality = 100)
p13 <- ggplot(marker_uber[which(marker_uber$CellType!=""),], aes(x = Kmediods, COL1A1)) + 
  geom_violin(aes(fill = Kmediods))
p14 <- ggplot(marker_uber[which(marker_uber$CellType!=""),], aes(x = CellType, COL1A1)) + 
  geom_violin(aes(fill = CellType))
plot_grid(p13,p14,nrow = 1)
dev.off()

jpeg(filename = paste0("Image/Label_Celltype/GSE81076_Labeling.jpeg"),width = 1600, height = 1200,quality = 100)
plot_grid(p2,p4,p6,p8,p10,p12,p14,nrow = 4)
dev.off()

write.table(marker_uber,
            file="./Data/GSE81076_marker_metadata.tsv",
            sep="\t", quote=F, row.names=F, col.names=T)

save.image("Progress/AssignCellType_GSE81076.RData")
############
# GSE85241 #
############
rm(list=ls())
gc()

# read in the normalized expression data
gse85241.norm <- read.table("./Data/GSE85241_SFnorm.tsv",
                            h=TRUE, sep="\t", stringsAsFactors=FALSE)

# remove gene id column for hvg discovery, last column
rownames(gse85241.norm) <- gse85241.norm$gene_id

gse85241.hvg_df <- read.table("./Data/GSE85241-HVG.tsv", 
                              sep="\t", h=TRUE, stringsAsFactors=FALSE)

gse85241.hvg <- gse85241.norm[gse85241.norm$gene_id %in% gse85241.hvg_df$gene_id,
                              1:(dim(gse85241.norm)[2]-1)]

gse85241.meta <- read.table("./Data/GSE85241_metadata.tsv",
                            sep="\t", h=TRUE, stringsAsFactors=FALSE)

# get low dimensional embedding of genes by tSNE
gse85241.tsne <- Rtsne(t(gse85241.hvg), perplexity=30)
gse85241.map <- data.frame(gse85241.tsne$Y)
colnames(gse85241.map) <- c("Dim1", "Dim2")
gse85241.map$Sample <- colnames(gse85241.hvg)

# merge metadata and clustering info
gse85241.uber <- merge(gse85241.map, gse85241.meta,
                       by='Sample')

marker_genes <- rownames(gse85241.norm)[grepl(rownames(gse85241.norm), pattern="(^GCG)|(^KRT19)|(^INS)|
                                              |(^SST)|(^PPY)|(^PRSS1)|(^COL1A1)|(^GHRL)|(^ESAM)") ]

marker_exprs <- gse85241.norm[marker_genes, 1:(dim(gse85241.norm)[2]-1)]

set.seed(42)
kmed <- pam(t(gse85241.hvg), 9)
panc.kmed <- data.frame(cbind(kmed$clustering))
colnames(panc.kmed) <- "Kmediods"
panc.kmed$Sample <- rownames(panc.kmed)
panc.kmed$Kmediods <- as.factor(panc.kmed$Kmediods)
panc.meta <- merge(gse85241.uber, panc.kmed, by='Sample')

marker_df <- data.frame(t(marker_exprs))
marker_df$Sample <- rownames(marker_df)

mark.dim <- dim(marker_df)
marker.kmed <- pam(marker_df[, 1:mark.dim[2]-1], 9)
marker.cluster <- data.frame(cbind(marker.kmed$clustering))
colnames(marker.cluster) <- c("markClust")
marker.cluster$Sample <- rownames(marker.cluster)
marker.cluster$markClust <- as.factor(marker.cluster$markClust)

marker_merge <- merge(panc.meta, marker_df, by='Sample')
marker_uber <- merge(marker_merge, marker.cluster, by='Sample')

# uncomment the following plots to see the overlay of each marker gene
# on the tSNE.  This illustrates the cell type label inference.
# mediods 1 & 3 & 4 = alpha cells
ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=GCG,
                         shape=Kmediods)) +
  #geom_point(data = marker_merge[which(marker_merge$Kmediods==7),], size=2) +
  #geom_point(data = marker_merge[which(marker_merge$GCG>8),], size=2) + 
  geom_point(size=2)+ 
  theme_classic() +
  scale_colour_continuous(low='white', high='red') +
  scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))

# mediod 9 = beta cells, INS highest
ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=INS,
                         shape=Kmediods)) +
  geom_point(size=2) + theme_classic() +
  scale_colour_continuous(low='white', high='red') +
  scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))

# mediod 3 & 4 & 9 = delta cells, SST highest
ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=SST,
                         shape=Kmediods)) +
  #geom_point(data = marker_merge[which(marker_merge$SST>9),], size=2) + 
  geom_point(size=2) + 
  theme_classic() +
  scale_colour_continuous(low='white', high='red') +
  scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))

# mediod 3 = Gamma cells, PPY highest
ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=PPY,
                         shape=Kmediods)) +
  geom_point(size=2) + theme_classic() +
  scale_colour_continuous(low='white', high='red') +
  scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))

# mediod 5 & 8 = acinar cells, PRSS1 highest
ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=PRSS1,
                         shape=Kmediods)) +
  geom_point(size=2) + theme_classic() +
  scale_colour_continuous(low='white', high='red') +
  scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))

# mediod 6 = ductal cells, KRT19 highest
ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=KRT19,
                         shape=Kmediods)) +
  geom_point(size=2) + theme_classic() +
  scale_colour_continuous(low='white', high='red') +
  scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))

# small pop of mediod 2 = mesenchyme cells, COL1A1 highest, same as Gamma cells
ggplot(marker_merge, aes(x=Dim1, y=Dim2, colour=COL1A1,
                         shape=Kmediods)) +
  geom_point(size=2) + theme_classic() +
  scale_colour_continuous(low='white', high='red') +
  scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0))


table(marker_uber$Kmediods, marker_uber$markClust)
marker_uber$CellType <- ""
marker_uber$CellType[marker_uber$Kmediods %in% c(3)] <- "Gamma"
marker_uber$CellType[marker_uber$Kmediods %in%  c(9) | (marker_uber$Kmediods %in% c(3) & marker_merge$INS >= 8) ] <- "Beta"
marker_uber$CellType[marker_uber$Kmediods %in% c(1, 4) | (marker_uber$Kmediods %in% c(3) & marker_merge$GCG >= 8)] <- "Alpha"
marker_uber$CellType[marker_uber$Kmediods %in% c(3, 4, 9) & marker_merge$SST >= 9] <- "Delta"
marker_uber$CellType[marker_uber$Kmediods %in% c(5, 8)] <- "Acinar"
marker_uber$CellType[marker_uber$Kmediods %in% c(6,7) ] <- "Ductal"
marker_uber$CellType[marker_uber$Kmediods %in% c(2)] <- "Mesenchyme"


table(marker_uber$CellType)

# Draw violin plot to show the marker gene experssion levels in each identified cell types  
marker_gene_list <- c("GCG","INS","SST","PPY","PRSS1","KRT19","COL1A1")

# GCG
jpeg(filename = paste0("Image/Label_Celltype/GSE85241_",marker_gene_list[1],"_for_Alpha.jpeg"),width = 1600, height = 600,quality = 100)
p1 <- ggplot(marker_uber, aes(x = Kmediods, GCG)) + 
    geom_violin(aes(fill = Kmediods))
p2 <- ggplot(marker_uber, aes(x = CellType, GCG)) + 
    geom_violin(aes(fill = CellType))
plot_grid(p1,p2,nrow = 1)
dev.off()
# INS
jpeg(filename = paste0("Image/Label_Celltype/GSE85241_",marker_gene_list[2],"_for_Beta.jpeg"),width = 1600, height = 600,quality = 100)
p3 <- ggplot(marker_uber, aes(x = Kmediods, INS)) + 
  geom_violin(aes(fill = Kmediods))
p4 <- ggplot(marker_uber, aes(x = CellType, INS)) + 
  geom_violin(aes(fill = CellType))
plot_grid(p3,p4,nrow = 1)
dev.off()
# SST
jpeg(filename = paste0("Image/Label_Celltype/GSE85241_",marker_gene_list[3],"_for_Delta.jpeg"),width = 1600, height = 600,quality = 100)
p5 <- ggplot(marker_uber, aes(x = Kmediods, SST)) + 
  geom_violin(aes(fill = Kmediods))
p6 <- ggplot(marker_uber, aes(x = CellType, SST)) + 
  geom_violin(aes(fill = CellType))
plot_grid(p5,p6,nrow = 1)
dev.off()
# PPY
jpeg(filename = paste0("Image/Label_Celltype/GSE85241_",marker_gene_list[4],"_for_Gamma.jpeg"),width = 1600, height = 600,quality = 100)
p7 <- ggplot(marker_uber, aes(x = Kmediods, PPY)) + 
  geom_violin(aes(fill = Kmediods))
p8 <- ggplot(marker_uber, aes(x = CellType, PPY)) + 
  geom_violin(aes(fill = CellType))
plot_grid(p7,p8,nrow = 1)
dev.off()
# PRSS1
jpeg(filename = paste0("Image/Label_Celltype/GSE85241_",marker_gene_list[5],"_for_Acinar.jpeg"),width = 1600, height = 600,quality = 100)
p9 <- ggplot(marker_uber, aes(x = Kmediods, PRSS1)) + 
  geom_violin(aes(fill = Kmediods))
p10 <- ggplot(marker_uber, aes(x = CellType, PRSS1)) + 
  geom_violin(aes(fill = CellType))
plot_grid(p9,p10,nrow = 1)
dev.off()
# KRT19
jpeg(filename = paste0("Image/Label_Celltype/GSE85241_",marker_gene_list[6],"_for_Ductal.jpeg"),width = 1600, height = 600,quality = 100)
p11 <- ggplot(marker_uber, aes(x = Kmediods, KRT19)) + 
  geom_violin(aes(fill = Kmediods))
p12 <- ggplot(marker_uber, aes(x = CellType, KRT19)) + 
  geom_violin(aes(fill = CellType))
plot_grid(p11,p12,nrow = 1)
dev.off()
# COL1A1
jpeg(filename = paste0("Image/Label_Celltype/GSE85241_",marker_gene_list[7],"_for_Mesenchyme.jpeg"),width = 1600, height = 600,quality = 100)
p13 <- ggplot(marker_uber, aes(x = Kmediods, COL1A1)) + 
  geom_violin(aes(fill = Kmediods))
p14 <- ggplot(marker_uber, aes(x = CellType, COL1A1)) + 
  geom_violin(aes(fill = CellType))
plot_grid(p13,p14,nrow = 1)
dev.off()

jpeg(filename = paste0("Image/Label_Celltype/GSE85241_Labeling.jpeg"),width = 1600, height = 1200,quality = 100)
plot_grid(p2,p4,p6,p8,p10,p12,p14,nrow = 4)
dev.off()

write.table(marker_uber,
            file="./Data/GSE85241_marker_metadata.tsv",
            sep="\t", quote=F, row.names=F, col.names=T)



save.image("Progress/AssignCellType_GSE85241.RData")

rm(list=ls())
gc()

################################
#PancreasProcessingCorrection.R#
################################
############
# GSE81076 #
############
# do a bit of tidying to make sure expression matrices are all conformable and
# only contain data with cell type assignments

# explicitly read in the study data
datah1 <- read.table("./Data/GSE81076_SFnorm.tsv",
                     h=TRUE, sep="\t", stringsAsFactors=FALSE)

read1 <- read.table("./Data/GSE81076_readcount.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)
# gene IDs are in the last column
genes1 <- datah1$gene_id
rownames(datah1) <- genes1
#datah1 <- datah1[, 1:(dim(datah1)[2]-1)]

hvg1 <- read.table("./Data/GSE81076-HVG.tsv",
                   h=TRUE, sep="\t", stringsAsFactors=F)
HVG1 <- hvg1$gene_id[hvg1$HVG == 1]

meta1 <- read.table("./Data/GSE81076_marker_metadata.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)

# select meta data for which there is expression data
meta1 <- meta1[meta1$Sample %in% colnames(datah1), ]

# standardize the cell labels
## NB there are 112 samples with unassigned cell types, remove these
celltypes1 <- meta1$CellType
no.label1 <- meta1$Sample[meta1$CellType == ""]
datah1 <- datah1[, !colnames(datah1) %in% no.label1]
read1 <- read1[,!colnames(read1) %in% no.label1]
meta1 <- meta1[!meta1$Sample %in% no.label1, ]
celltypes1 <- celltypes1[celltypes1 != ""]
samples1 <- colnames(read1)

# check all dimensions match up
if(dim(datah1)[2] == dim(meta1)[1]) {dim(datah1)[2] == length(celltypes1)}

#generate library id
n1<-ncol(read1)
library_id <- rep(NA,n1)
for(i in 1:n1){
  sample_id<-meta1$Sample[i]
  library_id[i]<-strsplit(sample_id,"_")[[1]][1]
}
meta1<-cbind(meta1,library_id)

############
# GSE85241 #
############
# do a bit of tidying to make sure expression matrices are all conformable and
# only contain data with cell type assignments

# explicitly read in the study data
datah2 <- read.table("./Data/GSE85241_SFnorm.tsv",
                     h=TRUE, sep="\t", stringsAsFactors=FALSE)

read2 <- read.table("./Data/GSE85241_readcount.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)
# gene IDs are in the last column
genes2 <- datah2$gene_id
rownames(datah2) <- genes2
#datah2 <- datah2[, 1:(dim(datah2)[2]-1)]

hvg2 <- read.table("./Data/GSE85241-HVG.tsv",
                   h=TRUE, sep="\t", stringsAsFactors=F)
HVG2 <- hvg2$gene_id[hvg2$HVG == 1]

meta2 <- read.table("./Data/GSE85241_marker_metadata.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)

# select meta data for which there is expression data
meta2 <- meta2[meta2$Sample %in% colnames(datah2), ]

# standardize the cell labels
## NB there are not any samples with unassigned cell types
celltypes2 <- meta2$CellType
no.label2 <- meta2$Sample[meta2$CellType == ""]
datah2 <- datah2[, !colnames(datah2) %in% no.label2]
read2 <- read2[,!colnames(read2) %in% no.label2]
meta2 <- meta2[!meta2$Sample %in% no.label2, ]
celltypes2 <- celltypes2[celltypes2 != ""]
samples2 <- colnames(read2)


# check all dimensions match up
if(dim(datah2)[2] == dim(meta2)[1]) {dim(datah2)[2] == length(celltypes2)}

#generate library id
n2<-ncol(read2)
library_id <- rep(NA,n2)
for(i in 1:n2){
  sample_id<-meta2$Sample[i]
  library_id[i]<-strsplit(sample_id,"_")[[1]][1]
}
meta2<-cbind(meta2,library_id)

############
# GSE86473 #
############
# do a bit of tidying to make sure expression matrices are all conformable and
# only contain data with cell type assignments

# explicitly read in the study data
datah3 <- read.table("./Data/GSE86473_SFnorm.tsv",
                     h=TRUE, sep="\t", stringsAsFactors=FALSE)
read3 <- read.table("./Data/GSE86473_readcount.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)
# gene IDs are in the last column
genes3 <- datah3$gene_id
rownames(datah3) <- genes3
datah3 <- datah3[, 1:(dim(datah3)[2]-1)]

hvg3 <- read.table("./Data/GSE86473-HVG.tsv",
                   h=TRUE, sep="\t", stringsAsFactors=F)
HVG3 <- hvg3$gene_id[hvg3$HVG == 1]

meta3 <- read.table("./Data/GSE86473_metadata.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)

# select meta data for which there is expression data
meta3 <- meta3[meta3$Sample %in% colnames(datah3), ]

# standardize the cell labels
## NB there are not any samples with unassigned cell types
celltypes3 <- meta3$CellType
no.label3 <- meta3$Sample[meta3$CellType == ""]
datah3 <- datah3[, !colnames(datah3) %in% no.label3]
read3 <- read3[,!colnames(read3) %in% no.label3]
meta3 <- meta3[!meta3$Sample %in% no.label3, ]
celltypes3 <- celltypes3[celltypes3 != ""]
samples3 <- colnames(read3)

library_id<-meta3$Plate
meta3<-cbind(meta3,library_id)


###############
# E-MTAB-5061 #
###############
# do a bit of tidying to make sure expression matrices are all conformable and
# only contain data with cell type assignments
# explicitly read in the study data
datah4 <- read.table("./Data/E-MTAB-5061_SFnorm.tsv",
                     h=TRUE, sep="\t", stringsAsFactors=FALSE)

read4 <- read.table("./Data/E-MTAB-5061_readcount.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)
# gene IDs are in the last column
genes4 <- datah4$gene_id
rownames(datah4) <- genes4
datah4 <- datah4[, 1:(dim(datah4)[2]-1)]

hvg4 <- read.table("./Data/E-MTAB-5061-HVG.tsv",
                   h=TRUE, sep="\t", stringsAsFactors=F)
HVG4 <- hvg4$gene_id[hvg4$HVG == 1]

meta4 <- read.table("./Data/E-MTAB-5061_metadata.tsv",
                    h=TRUE, sep="\t", stringsAsFactors=FALSE)

# select meta data for which there is expression data
meta4 <- meta4[meta4$Sample %in% colnames(datah4), ]

# standardize the cell labels
## NB there are not any samples with unassigned cell types
celltypes4 <- meta4$CellType
no.label4 <- meta4$Sample[meta4$CellType == ""]
datah4 <- datah4[, !colnames(datah4) %in% no.label4]
read4 <- read4[,!colnames(read4) %in% no.label4]
meta4 <- meta4[!meta4$Sample %in% no.label4, ]
celltypes4 <- celltypes4[celltypes4 != ""]
samples4 <- colnames(read4)

# check all dimensions match up
if(dim(datah4)[2] == dim(meta4)[1]) {dim(datah4)[2] == length(celltypes4)}

#generate library id
n4<-ncol(read4)
library_id <- rep(NA,n4)
for(i in 1:n4){
  sample_id<-meta4$Sample[i]
  library_id[i]<-strsplit(sample_id,"_")[[1]][1]
}
meta4<-cbind(meta4,library_id)

# create one big meta data frame
all.meta <- do.call(rbind.data.frame, list("b1"=meta1[, c("Sample", "CellType", "Protocol", "Study")],
                                           "b2"=meta2[, c("Sample", "CellType", "Protocol", "Study")],
                                           "b3"=meta3[, c("Sample", "CellType", "Protocol", "Study")],
                                           "b4"=meta4[, c("Sample", "CellType", "Protocol", "Study")]))

# Merge all four data matrices together based on a common set of gene IDs
common.genes <- intersect(genes1, intersect(genes2, intersect(genes3, genes4)))
rownames(read1)<-genes1
rownames(read2)<-genes2
rownames(read3)<-genes3
rownames(read4)<-genes4
r.read1 <- read1[common.genes, ]
r.read2 <- read2[common.genes, ]
r.read3 <- read3[common.genes, ]
r.read4 <- read4[common.genes, ]
# check the orders of names are the same
all(rownames(r.read1) == rownames(r.read2)&rownames(r.read1) == rownames(r.read3)&
      rownames(r.read1) == rownames(r.read4)&rownames(r.read2) == rownames(r.read3)&
      rownames(r.read2) == rownames(r.read4)&rownames(r.read3) == rownames(r.read4))

# merge uncorrected data based on gene IDs
r.read1$gene_id <- rownames(r.read1)
r.read2$gene_id <- rownames(r.read2)
r.read3$gene_id <- rownames(r.read3)
r.read4$gene_id <- rownames(r.read4)

raw.all <- Reduce(x=list("b1"=r.read1, "b2"=r.read2,"b3"=r.read3,"b4"=r.read4),
                  f=function(x, y) merge(x, y, by='gene_id'))

rownames(raw.all) <- raw.all$gene_id
raw.all <- raw.all[, 2:dim(raw.all)[2]]

# geometric mean on v.small fractional numbers might become quite
# unstable.  Could just use take geometric mean on a log scale
# this still uses the intersection of all genes, if it doesn't
# appear in the list of genes at all, set the pval - 1.
all.genes <- unique(c(genes1, genes2, genes3, genes4))

# this is a vector of highly genes that intersect the commonly expressed genes.
all.hvg <- unique(c(HVG1, HVG2, HVG3, HVG4))
hvg_genes <- intersect(common.genes, all.hvg)
raw.read.hvg <- as.matrix(raw.all[rownames(raw.all) %in% hvg_genes, ])

# assign small weird cell types from GSE85241 and E=MTAB-5061 to 'other'
all.meta$CellType[grepl(all.meta$CellType, pattern="PP")] <- "Gamma"
all.meta$CellType[grepl(all.meta$CellType, pattern="Mesenchyme")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Co-ex")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Endo")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Epsi")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Mast")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="MHC")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Uncl")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Not")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="PSC")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="None")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="Stellate")] <- "other"
all.meta$CellType[grepl(all.meta$CellType, pattern="None")] <- "other"

nb <- c(ncol(r.read1)-1, ncol(r.read2)-1,
        ncol(r.read3)-1, ncol(r.read4)-1)

PancreasCounts <- list(GSE81076 = raw.read.hvg[,1:nb[1]],
                       GSE85241 = raw.read.hvg[,1:nb[2] + nb[1]],
                       GSE86473 = raw.read.hvg[,1:nb[3] + sum(nb[1:2])],
                       EMTAB5061 = raw.read.hvg[,1:nb[4] + sum(nb[1:3])])

#store count data matrix
if(!dir.exists("RawCountData")){
  dir.create("RawCountData")
}

save(all.meta, file = "RawCountData/pancreas_metadata.RData")
save(PancreasCounts, file="RawCountData/pancreas_countdata.RData")

pancreas_count <- cbind(PancreasCounts[[1]],PancreasCounts[[2]],PancreasCounts[[3]],PancreasCounts[[4]])
write.table(pancreas_count,file="RawCountData/count_data_pancreas.txt",row.names = FALSE, col.names = FALSE)

# txt version raw count data
mathed_index <- match(colnames(pancreas_count),all.meta$Sample)
metadata <- data.frame(Protocol = all.meta$Protocol[mathed_index], Study = all.meta$Study[mathed_index], CellType = all.meta$CellType[mathed_index])
rownames(metadata) <- all.meta$Sample[mathed_index]

write.table(metadata,file="RawCountData/pancreas_metadata.txt")
#load by metadata <- read.table("RawCountData/pancreas_metadata.txt", header = TRUE)



N <- ncol(pancreas_count)
G <- nrow(pancreas_count)
B <- length(nb)
write.table(c(N,G,B,nb),file="RawCountData/dim_pancreas.txt",row.names = FALSE, col.names = FALSE)

###########
# END