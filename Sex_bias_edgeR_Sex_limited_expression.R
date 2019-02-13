##############################################################################
### install
 
##source("http://bioconductor.org/biocLite.R")
##    biocLite("edgeR")


##############################################################################
### library

library("edgeR")
library("VennDiagram")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(cowplot)
library(stringr)
library(gtable)
library(pheatmap)
library(RColorBrewer)
require(vegan)
library(pvclust)
library(raster)
library("SuperExactTest")

print (sessionInfo())

# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS High Sierra 10.13.3

# Matrix products: default
# BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
 # [1] SuperExactTest_0.99.4 raster_2.5-8          sp_1.2-5              pvclust_2.0-0         vegan_2.4-4           permute_0.9-4         RColorBrewer_1.1-2    pheatmap_1.0.8        gtable_0.2.0          stringr_1.2.0         cowplot_0.8.0        
# [12] lattice_0.20-35       ggplot2_2.2.1         gridExtra_2.3         VennDiagram_1.6.17    futile.logger_1.4.3   edgeR_3.18.1          limma_3.32.10        

# loaded via a namespace (and not attached):
 # [1] Rcpp_0.12.13         cluster_2.0.6        magrittr_1.5         MASS_7.3-47          munsell_0.4.3        colorspace_1.3-2     rlang_0.1.6          plyr_1.8.4           tools_3.4.1          parallel_3.4.1       nlme_3.1-131        
# [12] mgcv_1.8-22          lambda.r_1.2         lazyeval_0.2.1       tibble_1.3.4         Matrix_1.2-11        futile.options_1.0.0 stringi_1.1.5        compiler_3.4.1       scales_0.5.0         locfit_1.5-9.1      
# > 

##############################################################################
# Data

rawdata_Tbi_RBBH <- read.csv("Data/Read_counts/Tbi_Tte_RBBH_orth_counts_Arthropoda+Mixed+NOBLASTHIT_K2edge.counts", check.names=FALSE, stringsAsFactors=FALSE)
rawdata_Tce_RBBH <- read.csv("Data/Read_counts/Tce_Tms_RBBH_orth_counts_Arthropoda+Mixed+NOBLASTHIT_K2edge.counts", check.names=FALSE, stringsAsFactors=FALSE)
rawdata_Tcm_RBBH <- read.csv("Data/Read_counts/Tcm_Tsi_RBBH_orth_counts_Arthropoda+Mixed+NOBLASTHIT_K2edge.counts", check.names=FALSE, stringsAsFactors=FALSE)
rawdata_Tpa_RBBH <- read.csv("Data/Read_counts/Tpa_Tge_RBBH_orth_counts_Arthropoda+Mixed+NOBLASTHIT_K2edge.counts", check.names=FALSE, stringsAsFactors=FALSE)
rawdata_Tps_RBBH <- read.csv("Data/Read_counts/Tps_Tdi_RBBH_orth_counts_Arthropoda+Mixed+NOBLASTHIT_K2edge.counts", check.names=FALSE, stringsAsFactors=FALSE)


head(rawdata_Tbi_RBBH)
head(rawdata_Tce_RBBH)
head(rawdata_Tcm_RBBH)
head(rawdata_Tpa_RBBH)
head(rawdata_Tps_RBBH)

colnames(rawdata_Tbi_RBBH)
length(colnames(rawdata_Tbi_RBBH))


dir.create("Output")
dir.create("Output/DE")
setwd("output/DE")

###### DE analysis DOING for each tissue seperatly 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

## WB
y_WB_RBBH_Tbi_UF <- DGEList(counts=rawdata_Tbi_RBBH[,c(
"Tbi_SF_WB_Md_Re1","Tbi_SF_WB_Md_Re2","Tbi_SF_WB_Md_Re3",
"Tbi_SM_WB_Md_Re1","Tbi_SM_WB_Md_Re2","Tbi_SM_WB_Md_Re3",
"Tte_AF_WB_Vi_Re1","Tte_AF_WB_Vi_Re2","Tte_AF_WB_Vi_Re3")], genes=rawdata_Tbi_RBBH[,1:1])

y_WB_RBBH_Tce_UF <- DGEList(counts=rawdata_Tce_RBBH[,c(
"Tce_SF_WB_Md_Re1","Tce_SF_WB_Md_Re2","Tce_SF_WB_Md_Re3",
"Tce_SM_WB_Md_Re1","Tce_SM_WB_Md_Re2","Tce_SM_WB_Md_Re3",
"Tms_AF_WB_Vi_Re1","Tms_AF_WB_Vi_Re2","Tms_AF_WB_Vi_Re3")], genes=rawdata_Tce_RBBH[,1:1])

y_WB_RBBH_Tcm_UF <- DGEList(counts=rawdata_Tcm_RBBH[,c(
"Tcm_SF_WB_Md_Re1","Tcm_SF_WB_Md_Re2","Tcm_SF_WB_Md_Re3",
"Tcm_SM_WB_Md_Re1","Tcm_SM_WB_Md_Re2","Tcm_SM_WB_Md_Re3",
"Tsi_AF_WB_Vi_Re1","Tsi_AF_WB_Vi_Re2","Tsi_AF_WB_Vi_Re3")], genes=rawdata_Tcm_RBBH[,1:1])

y_WB_RBBH_Tpa_UF <- DGEList(counts=rawdata_Tpa_RBBH[,c(
"Tpa_SF_WB_Md_Re1","Tpa_SF_WB_Md_Re2","Tpa_SF_WB_Md_Re3",
"Tpa_SM_WB_Md_Re1","Tpa_SM_WB_Md_Re2","Tpa_SM_WB_Md_Re3",
"Tge_AF_WB_Vi_Re1","Tge_AF_WB_Vi_Re2","Tge_AF_WB_Vi_Re3")], genes=rawdata_Tpa_RBBH[,1:1])

y_WB_RBBH_Tps_UF <- DGEList(counts=rawdata_Tps_RBBH[,c(
"Tps_SF_WB_Md_Re1","Tps_SF_WB_Md_Re2","Tps_SF_WB_Md_Re3",
"Tps_SM_WB_Md_Re1","Tps_SM_WB_Md_Re2","Tps_SM_WB_Md_Re3",
"Tdi_AF_WB_Vi_Re1","Tdi_AF_WB_Vi_Re2","Tdi_AF_WB_Vi_Re3")], genes=rawdata_Tps_RBBH[,1:1])

## RT
y_RT_RBBH_Tbi_UF <- DGEList(counts=rawdata_Tbi_RBBH[,c(
"Tbi_SF_RT_Md_Re1","Tbi_SF_RT_Md_Re2","Tbi_SF_RT_Md_Re3",
"Tbi_SM_RT_Md_Re1","Tbi_SM_RT_Md_Re2","Tbi_SM_RT_Md_Re3",
"Tte_AF_RT_Vi_Re1","Tte_AF_RT_Vi_Re2","Tte_AF_RT_Vi_Re3")], genes=rawdata_Tbi_RBBH[,1:1])

y_RT_RBBH_Tce_UF <- DGEList(counts=rawdata_Tce_RBBH[,c(
"Tce_SF_RT_Md_Re1","Tce_SF_RT_Md_Re2","Tce_SF_RT_Md_Re3",
"Tce_SM_RT_Md_Re1","Tce_SM_RT_Md_Re2","Tce_SM_RT_Md_Re3",
"Tms_AF_RT_Vi_Re1","Tms_AF_RT_Vi_Re2","Tms_AF_RT_Vi_Re3")], genes=rawdata_Tce_RBBH[,1:1])

y_RT_RBBH_Tcm_UF <- DGEList(counts=rawdata_Tcm_RBBH[,c(
"Tcm_SF_RT_Md_Re1","Tcm_SF_RT_Md_Re2","Tcm_SF_RT_Md_Re3",
"Tcm_SM_RT_Md_Re1","Tcm_SM_RT_Md_Re2","Tcm_SM_RT_Md_Re3",
"Tsi_AF_RT_Vi_Re1","Tsi_AF_RT_Vi_Re2","Tsi_AF_RT_Vi_Re3")], genes=rawdata_Tcm_RBBH[,1:1])

y_RT_RBBH_Tpa_UF <- DGEList(counts=rawdata_Tpa_RBBH[,c(
"Tpa_SF_RT_Md_Re1","Tpa_SF_RT_Md_Re2","Tpa_SF_RT_Md_Re3",
"Tpa_SM_RT_Md_Re1","Tpa_SM_RT_Md_Re2","Tpa_SM_RT_Md_Re3",
"Tge_AF_RT_Vi_Re1","Tge_AF_RT_Vi_Re2","Tge_AF_RT_Vi_Re3")], genes=rawdata_Tpa_RBBH[,1:1])

y_RT_RBBH_Tps_UF <- DGEList(counts=rawdata_Tps_RBBH[,c(
"Tps_SF_RT_Md_Re1","Tps_SF_RT_Md_Re2","Tps_SF_RT_Md_Re3",
"Tps_SM_RT_Md_Re1","Tps_SM_RT_Md_Re2","Tps_SM_RT_Md_Re3",
"Tdi_AF_RT_Vi_Re1","Tdi_AF_RT_Vi_Re2","Tdi_AF_RT_Vi_Re3")], genes=rawdata_Tps_RBBH[,1:1])

## LG
y_LG_RBBH_Tbi_UF <- DGEList(counts=rawdata_Tbi_RBBH[,c(
"Tbi_SF_LG_Md_Re1","Tbi_SF_LG_Md_Re2","Tbi_SF_LG_Md_Re3",
"Tbi_SM_LG_Md_Re1","Tbi_SM_LG_Md_Re2","Tbi_SM_LG_Md_Re3",
"Tte_AF_LG_Vi_Re1","Tte_AF_LG_Vi_Re2","Tte_AF_LG_Vi_Re3")], genes=rawdata_Tbi_RBBH[,1:1])

y_LG_RBBH_Tce_UF <- DGEList(counts=rawdata_Tce_RBBH[,c(
"Tce_SF_LG_Md_Re1","Tce_SF_LG_Md_Re2","Tce_SF_LG_Md_Re3",
"Tce_SM_LG_Md_Re1","Tce_SM_LG_Md_Re2","Tce_SM_LG_Md_Re3",
"Tms_AF_LG_Vi_Re1","Tms_AF_LG_Vi_Re2","Tms_AF_LG_Vi_Re3")], genes=rawdata_Tce_RBBH[,1:1])

y_LG_RBBH_Tcm_UF <- DGEList(counts=rawdata_Tcm_RBBH[,c(
"Tcm_SF_LG_Md_Re1","Tcm_SF_LG_Md_Re2","Tcm_SF_LG_Md_Re3",
"Tcm_SM_LG_Md_Re1","Tcm_SM_LG_Md_Re2","Tcm_SM_LG_Md_Re3",
"Tsi_AF_LG_Vi_Re1","Tsi_AF_LG_Vi_Re2","Tsi_AF_LG_Vi_Re3")], genes=rawdata_Tcm_RBBH[,1:1])

y_LG_RBBH_Tpa_UF <- DGEList(counts=rawdata_Tpa_RBBH[,c(
"Tpa_SF_LG_Md_Re1","Tpa_SF_LG_Md_Re2","Tpa_SF_LG_Md_Re3",
"Tpa_SM_LG_Md_Re1","Tpa_SM_LG_Md_Re2","Tpa_SM_LG_Md_Re3",
"Tge_AF_LG_Vi_Re1","Tge_AF_LG_Vi_Re2","Tge_AF_LG_Vi_Re3")], genes=rawdata_Tpa_RBBH[,1:1])

y_LG_RBBH_Tps_UF <- DGEList(counts=rawdata_Tps_RBBH[,c(
"Tps_SF_LG_Md_Re1","Tps_SF_LG_Md_Re2","Tps_SF_LG_Md_Re3",
"Tps_SM_LG_Md_Re1","Tps_SM_LG_Md_Re2","Tps_SM_LG_Md_Re3",
"Tdi_AF_LG_Vi_Re1","Tdi_AF_LG_Vi_Re2","Tdi_AF_LG_Vi_Re3")], genes=rawdata_Tps_RBBH[,1:1])




##### normalise

norm_no_filter <- function(y,cpm_cut,cut_in_Nsams,gene_lens){
	y$samples$lib.size <- colSums(y$counts) 
	y <- calcNormFactors(y)
	cat("\nLib norm factors\n")
	print(y$samples)	
	cat("\nNumber number of genes / samples after NO filtering\n")
	print(dim(y))
	return(y)
}


y_Tbi_WB_F <- norm_no_filter(y_WB_RBBH_Tbi_UF)
y_Tce_WB_F <- norm_no_filter(y_WB_RBBH_Tce_UF)
y_Tcm_WB_F <- norm_no_filter(y_WB_RBBH_Tcm_UF)
y_Tpa_WB_F <- norm_no_filter(y_WB_RBBH_Tpa_UF)
y_Tps_WB_F <- norm_no_filter(y_WB_RBBH_Tps_UF)

y_Tbi_RT_F <- norm_no_filter(y_RT_RBBH_Tbi_UF)
y_Tce_RT_F <- norm_no_filter(y_RT_RBBH_Tce_UF)
y_Tcm_RT_F <- norm_no_filter(y_RT_RBBH_Tcm_UF)
y_Tpa_RT_F <- norm_no_filter(y_RT_RBBH_Tpa_UF)
y_Tps_RT_F <- norm_no_filter(y_RT_RBBH_Tps_UF)

y_Tbi_LG_F <- norm_no_filter(y_LG_RBBH_Tbi_UF)
y_Tce_LG_F <- norm_no_filter(y_LG_RBBH_Tce_UF)
y_Tcm_LG_F <- norm_no_filter(y_LG_RBBH_Tcm_UF)
y_Tpa_LG_F <- norm_no_filter(y_LG_RBBH_Tpa_UF)
y_Tps_LG_F <- norm_no_filter(y_LG_RBBH_Tps_UF)


#### calc FPKMs


Tbi_WB_F_FPKM <- rpkm(y_Tbi_WB_F, gene.length=rawdata_Tbi_RBBH$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tbi_WB_F_FPKM_1  <- as.data.frame(cbind(rawdata_Tbi_RBBH$Gene_name,Tbi_WB_F_FPKM))
Tce_WB_F_FPKM <- rpkm(y_Tce_WB_F, gene.length=rawdata_Tce_RBBH$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tce_WB_F_FPKM_1  <- as.data.frame(cbind(rawdata_Tce_RBBH$Gene_name,Tce_WB_F_FPKM))
Tcm_WB_F_FPKM <- rpkm(y_Tcm_WB_F, gene.length=rawdata_Tcm_RBBH$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tcm_WB_F_FPKM_1  <- as.data.frame(cbind(rawdata_Tcm_RBBH$Gene_name,Tcm_WB_F_FPKM))
Tpa_WB_F_FPKM <- rpkm(y_Tpa_WB_F, gene.length=rawdata_Tpa_RBBH$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tpa_WB_F_FPKM_1  <- as.data.frame(cbind(rawdata_Tpa_RBBH$Gene_name,Tpa_WB_F_FPKM))
Tps_WB_F_FPKM <- rpkm(y_Tps_WB_F, gene.length=rawdata_Tps_RBBH$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tps_WB_F_FPKM_1  <- as.data.frame(cbind(rawdata_Tps_RBBH$Gene_name,Tps_WB_F_FPKM))

Tbi_RT_F_FPKM <- rpkm(y_Tbi_RT_F, gene.length=rawdata_Tbi_RBBH$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tbi_RT_F_FPKM_1  <- as.data.frame(cbind(rawdata_Tbi_RBBH$Gene_name,Tbi_RT_F_FPKM))
Tce_RT_F_FPKM <- rpkm(y_Tce_RT_F, gene.length=rawdata_Tce_RBBH$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tce_RT_F_FPKM_1  <- as.data.frame(cbind(rawdata_Tce_RBBH$Gene_name,Tce_RT_F_FPKM))
Tcm_RT_F_FPKM <- rpkm(y_Tcm_RT_F, gene.length=rawdata_Tcm_RBBH$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tcm_RT_F_FPKM_1  <- as.data.frame(cbind(rawdata_Tcm_RBBH$Gene_name,Tcm_RT_F_FPKM))
Tpa_RT_F_FPKM <- rpkm(y_Tpa_RT_F, gene.length=rawdata_Tpa_RBBH$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tpa_RT_F_FPKM_1  <- as.data.frame(cbind(rawdata_Tpa_RBBH$Gene_name,Tpa_RT_F_FPKM))
Tps_RT_F_FPKM <- rpkm(y_Tps_RT_F, gene.length=rawdata_Tps_RBBH$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tps_RT_F_FPKM_1  <- as.data.frame(cbind(rawdata_Tps_RBBH$Gene_name,Tps_RT_F_FPKM))

Tbi_LG_F_FPKM <- rpkm(y_Tbi_LG_F, gene.length=rawdata_Tbi_RBBH$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tbi_LG_F_FPKM_1  <- as.data.frame(cbind(rawdata_Tbi_RBBH$Gene_name,Tbi_LG_F_FPKM))
Tce_LG_F_FPKM <- rpkm(y_Tce_LG_F, gene.length=rawdata_Tce_RBBH$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tce_LG_F_FPKM_1  <- as.data.frame(cbind(rawdata_Tce_RBBH$Gene_name,Tce_LG_F_FPKM))
Tcm_LG_F_FPKM <- rpkm(y_Tcm_LG_F, gene.length=rawdata_Tcm_RBBH$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tcm_LG_F_FPKM_1  <- as.data.frame(cbind(rawdata_Tcm_RBBH$Gene_name,Tcm_LG_F_FPKM))
Tpa_LG_F_FPKM <- rpkm(y_Tpa_LG_F, gene.length=rawdata_Tpa_RBBH$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tpa_LG_F_FPKM_1  <- as.data.frame(cbind(rawdata_Tpa_RBBH$Gene_name,Tpa_LG_F_FPKM))
Tps_LG_F_FPKM <- rpkm(y_Tps_LG_F, gene.length=rawdata_Tps_RBBH$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tps_LG_F_FPKM_1  <- as.data.frame(cbind(rawdata_Tps_RBBH$Gene_name,Tps_LG_F_FPKM))


Tbi_WB_F_FPKM_1$Tbi_SF_WB_Md_Re1 <- as.numeric(as.character(Tbi_WB_F_FPKM_1$Tbi_SF_WB_Md_Re1))
Tbi_WB_F_FPKM_1$Tbi_SF_WB_Md_Re2 <- as.numeric(as.character(Tbi_WB_F_FPKM_1$Tbi_SF_WB_Md_Re2))
Tbi_WB_F_FPKM_1$Tbi_SF_WB_Md_Re3 <- as.numeric(as.character(Tbi_WB_F_FPKM_1$Tbi_SF_WB_Md_Re3))
Tbi_WB_F_FPKM_1$Tbi_SM_WB_Md_Re1 <- as.numeric(as.character(Tbi_WB_F_FPKM_1$Tbi_SM_WB_Md_Re1))
Tbi_WB_F_FPKM_1$Tbi_SM_WB_Md_Re2 <- as.numeric(as.character(Tbi_WB_F_FPKM_1$Tbi_SM_WB_Md_Re2))
Tbi_WB_F_FPKM_1$Tbi_SM_WB_Md_Re3 <- as.numeric(as.character(Tbi_WB_F_FPKM_1$Tbi_SM_WB_Md_Re3))
Tbi_WB_F_FPKM_1$Tte_AF_WB_Vi_Re1 <- as.numeric(as.character(Tbi_WB_F_FPKM_1$Tte_AF_WB_Vi_Re1))
Tbi_WB_F_FPKM_1$Tte_AF_WB_Vi_Re2 <- as.numeric(as.character(Tbi_WB_F_FPKM_1$Tte_AF_WB_Vi_Re2))
Tbi_WB_F_FPKM_1$Tte_AF_WB_Vi_Re3 <- as.numeric(as.character(Tbi_WB_F_FPKM_1$Tte_AF_WB_Vi_Re3))

Tce_WB_F_FPKM_1$Tce_SF_WB_Md_Re1 <- as.numeric(as.character(Tce_WB_F_FPKM_1$Tce_SF_WB_Md_Re1))
Tce_WB_F_FPKM_1$Tce_SF_WB_Md_Re2 <- as.numeric(as.character(Tce_WB_F_FPKM_1$Tce_SF_WB_Md_Re2))
Tce_WB_F_FPKM_1$Tce_SF_WB_Md_Re3 <- as.numeric(as.character(Tce_WB_F_FPKM_1$Tce_SF_WB_Md_Re3))
Tce_WB_F_FPKM_1$Tce_SM_WB_Md_Re1 <- as.numeric(as.character(Tce_WB_F_FPKM_1$Tce_SM_WB_Md_Re1))
Tce_WB_F_FPKM_1$Tce_SM_WB_Md_Re2 <- as.numeric(as.character(Tce_WB_F_FPKM_1$Tce_SM_WB_Md_Re2))
Tce_WB_F_FPKM_1$Tce_SM_WB_Md_Re3 <- as.numeric(as.character(Tce_WB_F_FPKM_1$Tce_SM_WB_Md_Re3))
Tce_WB_F_FPKM_1$Tms_AF_WB_Vi_Re1 <- as.numeric(as.character(Tce_WB_F_FPKM_1$Tms_AF_WB_Vi_Re1))
Tce_WB_F_FPKM_1$Tms_AF_WB_Vi_Re2 <- as.numeric(as.character(Tce_WB_F_FPKM_1$Tms_AF_WB_Vi_Re2))
Tce_WB_F_FPKM_1$Tms_AF_WB_Vi_Re3 <- as.numeric(as.character(Tce_WB_F_FPKM_1$Tms_AF_WB_Vi_Re3))

Tcm_WB_F_FPKM_1$Tcm_SF_WB_Md_Re1 <- as.numeric(as.character(Tcm_WB_F_FPKM_1$Tcm_SF_WB_Md_Re1))
Tcm_WB_F_FPKM_1$Tcm_SF_WB_Md_Re2 <- as.numeric(as.character(Tcm_WB_F_FPKM_1$Tcm_SF_WB_Md_Re2))
Tcm_WB_F_FPKM_1$Tcm_SF_WB_Md_Re3 <- as.numeric(as.character(Tcm_WB_F_FPKM_1$Tcm_SF_WB_Md_Re3))
Tcm_WB_F_FPKM_1$Tcm_SM_WB_Md_Re1 <- as.numeric(as.character(Tcm_WB_F_FPKM_1$Tcm_SM_WB_Md_Re1))
Tcm_WB_F_FPKM_1$Tcm_SM_WB_Md_Re2 <- as.numeric(as.character(Tcm_WB_F_FPKM_1$Tcm_SM_WB_Md_Re2))
Tcm_WB_F_FPKM_1$Tcm_SM_WB_Md_Re3 <- as.numeric(as.character(Tcm_WB_F_FPKM_1$Tcm_SM_WB_Md_Re3))
Tcm_WB_F_FPKM_1$Tsi_AF_WB_Vi_Re1 <- as.numeric(as.character(Tcm_WB_F_FPKM_1$Tsi_AF_WB_Vi_Re1))
Tcm_WB_F_FPKM_1$Tsi_AF_WB_Vi_Re2 <- as.numeric(as.character(Tcm_WB_F_FPKM_1$Tsi_AF_WB_Vi_Re2))
Tcm_WB_F_FPKM_1$Tsi_AF_WB_Vi_Re3 <- as.numeric(as.character(Tcm_WB_F_FPKM_1$Tsi_AF_WB_Vi_Re3))

Tpa_WB_F_FPKM_1$Tpa_SF_WB_Md_Re1 <- as.numeric(as.character(Tpa_WB_F_FPKM_1$Tpa_SF_WB_Md_Re1))
Tpa_WB_F_FPKM_1$Tpa_SF_WB_Md_Re2 <- as.numeric(as.character(Tpa_WB_F_FPKM_1$Tpa_SF_WB_Md_Re2))
Tpa_WB_F_FPKM_1$Tpa_SF_WB_Md_Re3 <- as.numeric(as.character(Tpa_WB_F_FPKM_1$Tpa_SF_WB_Md_Re3))
Tpa_WB_F_FPKM_1$Tpa_SM_WB_Md_Re1 <- as.numeric(as.character(Tpa_WB_F_FPKM_1$Tpa_SM_WB_Md_Re1))
Tpa_WB_F_FPKM_1$Tpa_SM_WB_Md_Re2 <- as.numeric(as.character(Tpa_WB_F_FPKM_1$Tpa_SM_WB_Md_Re2))
Tpa_WB_F_FPKM_1$Tpa_SM_WB_Md_Re3 <- as.numeric(as.character(Tpa_WB_F_FPKM_1$Tpa_SM_WB_Md_Re3))
Tpa_WB_F_FPKM_1$Tge_AF_WB_Vi_Re1 <- as.numeric(as.character(Tpa_WB_F_FPKM_1$Tge_AF_WB_Vi_Re1))
Tpa_WB_F_FPKM_1$Tge_AF_WB_Vi_Re2 <- as.numeric(as.character(Tpa_WB_F_FPKM_1$Tge_AF_WB_Vi_Re2))
Tpa_WB_F_FPKM_1$Tge_AF_WB_Vi_Re3 <- as.numeric(as.character(Tpa_WB_F_FPKM_1$Tge_AF_WB_Vi_Re3))

Tps_WB_F_FPKM_1$Tps_SF_WB_Md_Re1 <- as.numeric(as.character(Tps_WB_F_FPKM_1$Tps_SF_WB_Md_Re1))
Tps_WB_F_FPKM_1$Tps_SF_WB_Md_Re2 <- as.numeric(as.character(Tps_WB_F_FPKM_1$Tps_SF_WB_Md_Re2))
Tps_WB_F_FPKM_1$Tps_SF_WB_Md_Re3 <- as.numeric(as.character(Tps_WB_F_FPKM_1$Tps_SF_WB_Md_Re3))
Tps_WB_F_FPKM_1$Tps_SM_WB_Md_Re1 <- as.numeric(as.character(Tps_WB_F_FPKM_1$Tps_SM_WB_Md_Re1))
Tps_WB_F_FPKM_1$Tps_SM_WB_Md_Re2 <- as.numeric(as.character(Tps_WB_F_FPKM_1$Tps_SM_WB_Md_Re2))
Tps_WB_F_FPKM_1$Tps_SM_WB_Md_Re3 <- as.numeric(as.character(Tps_WB_F_FPKM_1$Tps_SM_WB_Md_Re3))
Tps_WB_F_FPKM_1$Tdi_AF_WB_Vi_Re1 <- as.numeric(as.character(Tps_WB_F_FPKM_1$Tdi_AF_WB_Vi_Re1))
Tps_WB_F_FPKM_1$Tdi_AF_WB_Vi_Re2 <- as.numeric(as.character(Tps_WB_F_FPKM_1$Tdi_AF_WB_Vi_Re2))
Tps_WB_F_FPKM_1$Tdi_AF_WB_Vi_Re3 <- as.numeric(as.character(Tps_WB_F_FPKM_1$Tdi_AF_WB_Vi_Re3))

Tbi_RT_F_FPKM_1$Tbi_SF_RT_Md_Re1 <- as.numeric(as.character(Tbi_RT_F_FPKM_1$Tbi_SF_RT_Md_Re1))
Tbi_RT_F_FPKM_1$Tbi_SF_RT_Md_Re2 <- as.numeric(as.character(Tbi_RT_F_FPKM_1$Tbi_SF_RT_Md_Re2))
Tbi_RT_F_FPKM_1$Tbi_SF_RT_Md_Re3 <- as.numeric(as.character(Tbi_RT_F_FPKM_1$Tbi_SF_RT_Md_Re3))
Tbi_RT_F_FPKM_1$Tbi_SM_RT_Md_Re1 <- as.numeric(as.character(Tbi_RT_F_FPKM_1$Tbi_SM_RT_Md_Re1))
Tbi_RT_F_FPKM_1$Tbi_SM_RT_Md_Re2 <- as.numeric(as.character(Tbi_RT_F_FPKM_1$Tbi_SM_RT_Md_Re2))
Tbi_RT_F_FPKM_1$Tbi_SM_RT_Md_Re3 <- as.numeric(as.character(Tbi_RT_F_FPKM_1$Tbi_SM_RT_Md_Re3))
Tbi_RT_F_FPKM_1$Tte_AF_RT_Vi_Re1 <- as.numeric(as.character(Tbi_RT_F_FPKM_1$Tte_AF_RT_Vi_Re1))
Tbi_RT_F_FPKM_1$Tte_AF_RT_Vi_Re2 <- as.numeric(as.character(Tbi_RT_F_FPKM_1$Tte_AF_RT_Vi_Re2))
Tbi_RT_F_FPKM_1$Tte_AF_RT_Vi_Re3 <- as.numeric(as.character(Tbi_RT_F_FPKM_1$Tte_AF_RT_Vi_Re3))

Tce_RT_F_FPKM_1$Tce_SF_RT_Md_Re1 <- as.numeric(as.character(Tce_RT_F_FPKM_1$Tce_SF_RT_Md_Re1))
Tce_RT_F_FPKM_1$Tce_SF_RT_Md_Re2 <- as.numeric(as.character(Tce_RT_F_FPKM_1$Tce_SF_RT_Md_Re2))
Tce_RT_F_FPKM_1$Tce_SF_RT_Md_Re3 <- as.numeric(as.character(Tce_RT_F_FPKM_1$Tce_SF_RT_Md_Re3))
Tce_RT_F_FPKM_1$Tce_SM_RT_Md_Re1 <- as.numeric(as.character(Tce_RT_F_FPKM_1$Tce_SM_RT_Md_Re1))
Tce_RT_F_FPKM_1$Tce_SM_RT_Md_Re2 <- as.numeric(as.character(Tce_RT_F_FPKM_1$Tce_SM_RT_Md_Re2))
Tce_RT_F_FPKM_1$Tce_SM_RT_Md_Re3 <- as.numeric(as.character(Tce_RT_F_FPKM_1$Tce_SM_RT_Md_Re3))
Tce_RT_F_FPKM_1$Tms_AF_RT_Vi_Re1 <- as.numeric(as.character(Tce_RT_F_FPKM_1$Tms_AF_RT_Vi_Re1))
Tce_RT_F_FPKM_1$Tms_AF_RT_Vi_Re2 <- as.numeric(as.character(Tce_RT_F_FPKM_1$Tms_AF_RT_Vi_Re2))
Tce_RT_F_FPKM_1$Tms_AF_RT_Vi_Re3 <- as.numeric(as.character(Tce_RT_F_FPKM_1$Tms_AF_RT_Vi_Re3))

Tcm_RT_F_FPKM_1$Tcm_SF_RT_Md_Re1 <- as.numeric(as.character(Tcm_RT_F_FPKM_1$Tcm_SF_RT_Md_Re1))
Tcm_RT_F_FPKM_1$Tcm_SF_RT_Md_Re2 <- as.numeric(as.character(Tcm_RT_F_FPKM_1$Tcm_SF_RT_Md_Re2))
Tcm_RT_F_FPKM_1$Tcm_SF_RT_Md_Re3 <- as.numeric(as.character(Tcm_RT_F_FPKM_1$Tcm_SF_RT_Md_Re3))
Tcm_RT_F_FPKM_1$Tcm_SM_RT_Md_Re1 <- as.numeric(as.character(Tcm_RT_F_FPKM_1$Tcm_SM_RT_Md_Re1))
Tcm_RT_F_FPKM_1$Tcm_SM_RT_Md_Re2 <- as.numeric(as.character(Tcm_RT_F_FPKM_1$Tcm_SM_RT_Md_Re2))
Tcm_RT_F_FPKM_1$Tcm_SM_RT_Md_Re3 <- as.numeric(as.character(Tcm_RT_F_FPKM_1$Tcm_SM_RT_Md_Re3))
Tcm_RT_F_FPKM_1$Tsi_AF_RT_Vi_Re1 <- as.numeric(as.character(Tcm_RT_F_FPKM_1$Tsi_AF_RT_Vi_Re1))
Tcm_RT_F_FPKM_1$Tsi_AF_RT_Vi_Re2 <- as.numeric(as.character(Tcm_RT_F_FPKM_1$Tsi_AF_RT_Vi_Re2))
Tcm_RT_F_FPKM_1$Tsi_AF_RT_Vi_Re3 <- as.numeric(as.character(Tcm_RT_F_FPKM_1$Tsi_AF_RT_Vi_Re3))

Tpa_RT_F_FPKM_1$Tpa_SF_RT_Md_Re1 <- as.numeric(as.character(Tpa_RT_F_FPKM_1$Tpa_SF_RT_Md_Re1))
Tpa_RT_F_FPKM_1$Tpa_SF_RT_Md_Re2 <- as.numeric(as.character(Tpa_RT_F_FPKM_1$Tpa_SF_RT_Md_Re2))
Tpa_RT_F_FPKM_1$Tpa_SF_RT_Md_Re3 <- as.numeric(as.character(Tpa_RT_F_FPKM_1$Tpa_SF_RT_Md_Re3))
Tpa_RT_F_FPKM_1$Tpa_SM_RT_Md_Re1 <- as.numeric(as.character(Tpa_RT_F_FPKM_1$Tpa_SM_RT_Md_Re1))
Tpa_RT_F_FPKM_1$Tpa_SM_RT_Md_Re2 <- as.numeric(as.character(Tpa_RT_F_FPKM_1$Tpa_SM_RT_Md_Re2))
Tpa_RT_F_FPKM_1$Tpa_SM_RT_Md_Re3 <- as.numeric(as.character(Tpa_RT_F_FPKM_1$Tpa_SM_RT_Md_Re3))
Tpa_RT_F_FPKM_1$Tge_AF_RT_Vi_Re1 <- as.numeric(as.character(Tpa_RT_F_FPKM_1$Tge_AF_RT_Vi_Re1))
Tpa_RT_F_FPKM_1$Tge_AF_RT_Vi_Re2 <- as.numeric(as.character(Tpa_RT_F_FPKM_1$Tge_AF_RT_Vi_Re2))
Tpa_RT_F_FPKM_1$Tge_AF_RT_Vi_Re3 <- as.numeric(as.character(Tpa_RT_F_FPKM_1$Tge_AF_RT_Vi_Re3))

Tps_RT_F_FPKM_1$Tps_SF_RT_Md_Re1 <- as.numeric(as.character(Tps_RT_F_FPKM_1$Tps_SF_RT_Md_Re1))
Tps_RT_F_FPKM_1$Tps_SF_RT_Md_Re2 <- as.numeric(as.character(Tps_RT_F_FPKM_1$Tps_SF_RT_Md_Re2))
Tps_RT_F_FPKM_1$Tps_SF_RT_Md_Re3 <- as.numeric(as.character(Tps_RT_F_FPKM_1$Tps_SF_RT_Md_Re3))
Tps_RT_F_FPKM_1$Tps_SM_RT_Md_Re1 <- as.numeric(as.character(Tps_RT_F_FPKM_1$Tps_SM_RT_Md_Re1))
Tps_RT_F_FPKM_1$Tps_SM_RT_Md_Re2 <- as.numeric(as.character(Tps_RT_F_FPKM_1$Tps_SM_RT_Md_Re2))
Tps_RT_F_FPKM_1$Tps_SM_RT_Md_Re3 <- as.numeric(as.character(Tps_RT_F_FPKM_1$Tps_SM_RT_Md_Re3))
Tps_RT_F_FPKM_1$Tdi_AF_RT_Vi_Re1 <- as.numeric(as.character(Tps_RT_F_FPKM_1$Tdi_AF_RT_Vi_Re1))
Tps_RT_F_FPKM_1$Tdi_AF_RT_Vi_Re2 <- as.numeric(as.character(Tps_RT_F_FPKM_1$Tdi_AF_RT_Vi_Re2))
Tps_RT_F_FPKM_1$Tdi_AF_RT_Vi_Re3 <- as.numeric(as.character(Tps_RT_F_FPKM_1$Tdi_AF_RT_Vi_Re3))

Tbi_LG_F_FPKM_1$Tbi_SF_LG_Md_Re1 <- as.numeric(as.character(Tbi_LG_F_FPKM_1$Tbi_SF_LG_Md_Re1))
Tbi_LG_F_FPKM_1$Tbi_SF_LG_Md_Re2 <- as.numeric(as.character(Tbi_LG_F_FPKM_1$Tbi_SF_LG_Md_Re2))
Tbi_LG_F_FPKM_1$Tbi_SF_LG_Md_Re3 <- as.numeric(as.character(Tbi_LG_F_FPKM_1$Tbi_SF_LG_Md_Re3))
Tbi_LG_F_FPKM_1$Tbi_SM_LG_Md_Re1 <- as.numeric(as.character(Tbi_LG_F_FPKM_1$Tbi_SM_LG_Md_Re1))
Tbi_LG_F_FPKM_1$Tbi_SM_LG_Md_Re2 <- as.numeric(as.character(Tbi_LG_F_FPKM_1$Tbi_SM_LG_Md_Re2))
Tbi_LG_F_FPKM_1$Tbi_SM_LG_Md_Re3 <- as.numeric(as.character(Tbi_LG_F_FPKM_1$Tbi_SM_LG_Md_Re3))
Tbi_LG_F_FPKM_1$Tte_AF_LG_Vi_Re1 <- as.numeric(as.character(Tbi_LG_F_FPKM_1$Tte_AF_LG_Vi_Re1))
Tbi_LG_F_FPKM_1$Tte_AF_LG_Vi_Re2 <- as.numeric(as.character(Tbi_LG_F_FPKM_1$Tte_AF_LG_Vi_Re2))
Tbi_LG_F_FPKM_1$Tte_AF_LG_Vi_Re3 <- as.numeric(as.character(Tbi_LG_F_FPKM_1$Tte_AF_LG_Vi_Re3))

Tce_LG_F_FPKM_1$Tce_SF_LG_Md_Re1 <- as.numeric(as.character(Tce_LG_F_FPKM_1$Tce_SF_LG_Md_Re1))
Tce_LG_F_FPKM_1$Tce_SF_LG_Md_Re2 <- as.numeric(as.character(Tce_LG_F_FPKM_1$Tce_SF_LG_Md_Re2))
Tce_LG_F_FPKM_1$Tce_SF_LG_Md_Re3 <- as.numeric(as.character(Tce_LG_F_FPKM_1$Tce_SF_LG_Md_Re3))
Tce_LG_F_FPKM_1$Tce_SM_LG_Md_Re1 <- as.numeric(as.character(Tce_LG_F_FPKM_1$Tce_SM_LG_Md_Re1))
Tce_LG_F_FPKM_1$Tce_SM_LG_Md_Re2 <- as.numeric(as.character(Tce_LG_F_FPKM_1$Tce_SM_LG_Md_Re2))
Tce_LG_F_FPKM_1$Tce_SM_LG_Md_Re3 <- as.numeric(as.character(Tce_LG_F_FPKM_1$Tce_SM_LG_Md_Re3))
Tce_LG_F_FPKM_1$Tms_AF_LG_Vi_Re1 <- as.numeric(as.character(Tce_LG_F_FPKM_1$Tms_AF_LG_Vi_Re1))
Tce_LG_F_FPKM_1$Tms_AF_LG_Vi_Re2 <- as.numeric(as.character(Tce_LG_F_FPKM_1$Tms_AF_LG_Vi_Re2))
Tce_LG_F_FPKM_1$Tms_AF_LG_Vi_Re3 <- as.numeric(as.character(Tce_LG_F_FPKM_1$Tms_AF_LG_Vi_Re3))

Tcm_LG_F_FPKM_1$Tcm_SF_LG_Md_Re1 <- as.numeric(as.character(Tcm_LG_F_FPKM_1$Tcm_SF_LG_Md_Re1))
Tcm_LG_F_FPKM_1$Tcm_SF_LG_Md_Re2 <- as.numeric(as.character(Tcm_LG_F_FPKM_1$Tcm_SF_LG_Md_Re2))
Tcm_LG_F_FPKM_1$Tcm_SF_LG_Md_Re3 <- as.numeric(as.character(Tcm_LG_F_FPKM_1$Tcm_SF_LG_Md_Re3))
Tcm_LG_F_FPKM_1$Tcm_SM_LG_Md_Re1 <- as.numeric(as.character(Tcm_LG_F_FPKM_1$Tcm_SM_LG_Md_Re1))
Tcm_LG_F_FPKM_1$Tcm_SM_LG_Md_Re2 <- as.numeric(as.character(Tcm_LG_F_FPKM_1$Tcm_SM_LG_Md_Re2))
Tcm_LG_F_FPKM_1$Tcm_SM_LG_Md_Re3 <- as.numeric(as.character(Tcm_LG_F_FPKM_1$Tcm_SM_LG_Md_Re3))
Tcm_LG_F_FPKM_1$Tsi_AF_LG_Vi_Re1 <- as.numeric(as.character(Tcm_LG_F_FPKM_1$Tsi_AF_LG_Vi_Re1))
Tcm_LG_F_FPKM_1$Tsi_AF_LG_Vi_Re2 <- as.numeric(as.character(Tcm_LG_F_FPKM_1$Tsi_AF_LG_Vi_Re2))
Tcm_LG_F_FPKM_1$Tsi_AF_LG_Vi_Re3 <- as.numeric(as.character(Tcm_LG_F_FPKM_1$Tsi_AF_LG_Vi_Re3))

Tpa_LG_F_FPKM_1$Tpa_SF_LG_Md_Re1 <- as.numeric(as.character(Tpa_LG_F_FPKM_1$Tpa_SF_LG_Md_Re1))
Tpa_LG_F_FPKM_1$Tpa_SF_LG_Md_Re2 <- as.numeric(as.character(Tpa_LG_F_FPKM_1$Tpa_SF_LG_Md_Re2))
Tpa_LG_F_FPKM_1$Tpa_SF_LG_Md_Re3 <- as.numeric(as.character(Tpa_LG_F_FPKM_1$Tpa_SF_LG_Md_Re3))
Tpa_LG_F_FPKM_1$Tpa_SM_LG_Md_Re1 <- as.numeric(as.character(Tpa_LG_F_FPKM_1$Tpa_SM_LG_Md_Re1))
Tpa_LG_F_FPKM_1$Tpa_SM_LG_Md_Re2 <- as.numeric(as.character(Tpa_LG_F_FPKM_1$Tpa_SM_LG_Md_Re2))
Tpa_LG_F_FPKM_1$Tpa_SM_LG_Md_Re3 <- as.numeric(as.character(Tpa_LG_F_FPKM_1$Tpa_SM_LG_Md_Re3))
Tpa_LG_F_FPKM_1$Tge_AF_LG_Vi_Re1 <- as.numeric(as.character(Tpa_LG_F_FPKM_1$Tge_AF_LG_Vi_Re1))
Tpa_LG_F_FPKM_1$Tge_AF_LG_Vi_Re2 <- as.numeric(as.character(Tpa_LG_F_FPKM_1$Tge_AF_LG_Vi_Re2))
Tpa_LG_F_FPKM_1$Tge_AF_LG_Vi_Re3 <- as.numeric(as.character(Tpa_LG_F_FPKM_1$Tge_AF_LG_Vi_Re3))

Tps_LG_F_FPKM_1$Tps_SF_LG_Md_Re1 <- as.numeric(as.character(Tps_LG_F_FPKM_1$Tps_SF_LG_Md_Re1))
Tps_LG_F_FPKM_1$Tps_SF_LG_Md_Re2 <- as.numeric(as.character(Tps_LG_F_FPKM_1$Tps_SF_LG_Md_Re2))
Tps_LG_F_FPKM_1$Tps_SF_LG_Md_Re3 <- as.numeric(as.character(Tps_LG_F_FPKM_1$Tps_SF_LG_Md_Re3))
Tps_LG_F_FPKM_1$Tps_SM_LG_Md_Re1 <- as.numeric(as.character(Tps_LG_F_FPKM_1$Tps_SM_LG_Md_Re1))
Tps_LG_F_FPKM_1$Tps_SM_LG_Md_Re2 <- as.numeric(as.character(Tps_LG_F_FPKM_1$Tps_SM_LG_Md_Re2))
Tps_LG_F_FPKM_1$Tps_SM_LG_Md_Re3 <- as.numeric(as.character(Tps_LG_F_FPKM_1$Tps_SM_LG_Md_Re3))
Tps_LG_F_FPKM_1$Tdi_AF_LG_Vi_Re1 <- as.numeric(as.character(Tps_LG_F_FPKM_1$Tdi_AF_LG_Vi_Re1))
Tps_LG_F_FPKM_1$Tdi_AF_LG_Vi_Re2 <- as.numeric(as.character(Tps_LG_F_FPKM_1$Tdi_AF_LG_Vi_Re2))
Tps_LG_F_FPKM_1$Tdi_AF_LG_Vi_Re3 <- as.numeric(as.character(Tps_LG_F_FPKM_1$Tdi_AF_LG_Vi_Re3))

head(Tbi_WB_F_FPKM_1)

write.csv(Tbi_WB_F_FPKM_1, file="Tbi_WB_F_FPKM_RBBH_1.csv", row.names=FALSE)
write.csv(Tce_WB_F_FPKM_1, file="Tce_WB_F_FPKM_RBBH_1.csv", row.names=FALSE)
write.csv(Tcm_WB_F_FPKM_1, file="Tcm_WB_F_FPKM_RBBH_1.csv", row.names=FALSE)
write.csv(Tpa_WB_F_FPKM_1, file="Tpa_WB_F_FPKM_RBBH_1.csv", row.names=FALSE)
write.csv(Tps_WB_F_FPKM_1, file="Tps_WB_F_FPKM_RBBH_1.csv", row.names=FALSE)

write.csv(Tbi_RT_F_FPKM_1, file="Tbi_RT_F_FPKM_RBBH_1.csv", row.names=FALSE)
write.csv(Tce_RT_F_FPKM_1, file="Tce_RT_F_FPKM_RBBH_1.csv", row.names=FALSE)
write.csv(Tcm_RT_F_FPKM_1, file="Tcm_RT_F_FPKM_RBBH_1.csv", row.names=FALSE)
write.csv(Tpa_RT_F_FPKM_1, file="Tpa_RT_F_FPKM_RBBH_1.csv", row.names=FALSE)
write.csv(Tps_RT_F_FPKM_1, file="Tps_RT_F_FPKM_RBBH_1.csv", row.names=FALSE)

write.csv(Tbi_LG_F_FPKM_1, file="Tbi_LG_F_FPKM_RBBH_1.csv", row.names=FALSE)
write.csv(Tce_LG_F_FPKM_1, file="Tce_LG_F_FPKM_RBBH_1.csv", row.names=FALSE)
write.csv(Tcm_LG_F_FPKM_1, file="Tcm_LG_F_FPKM_RBBH_1.csv", row.names=FALSE)
write.csv(Tpa_LG_F_FPKM_1, file="Tpa_LG_F_FPKM_RBBH_1.csv", row.names=FALSE)
write.csv(Tps_LG_F_FPKM_1, file="Tps_LG_F_FPKM_RBBH_1.csv", row.names=FALSE)



################################################################################################################
##### class genes as SL or not

### FPKM_exp = FPKM val to consider EXPRESSED, FPKM_not_exp = FPKM to consider NOT EXPRESSED, Nsams to apply the filter to. 
### e.g get_SL_genes(Tbi_WB_F_FPKM_1, 2, 0, 3) is SL lim genes are genes expressed  ate at least 2 FKPM in 3  libs of one sex and 0 (or less) FKPM in the 3 libs of the other sex


get_SL_genes <- function(FPKM_df,FPKM_exp, FPKM_not_exp, Nsams,sp_w,tiss_w){
	
	### add mean expe
	
	FPKM_df$Avg_FPKM_SF = (FPKM_df[,2] + FPKM_df[,3] + FPKM_df[,4])/3
	FPKM_df$Avg_FPKM_SM = (FPKM_df[,5] + FPKM_df[,6] + FPKM_df[,7])/3
	FPKM_df$Avg_FPKM_AF = (FPKM_df[,8] + FPKM_df[,9] + FPKM_df[,10])/3

	
	## get female-limited genes
	keep_FL <- 
	  rowSums(FPKM_df[,2:4] >= FPKM_exp)     >= Nsams & 
      rowSums(FPKM_df[,5:7] <= FPKM_not_exp) >= Nsams 	
	
	FPKM_df_FL <- FPKM_df[keep_FL,]
	
	## get male-limited genes
	
	keep_ML <- 
	  rowSums(FPKM_df[,2:4] <= FPKM_not_exp) >= Nsams & 
      rowSums(FPKM_df[,5:7] >= FPKM_exp)     >= Nsams 	
	
	FPKM_df_ML <- FPKM_df[keep_ML,]	
	

	#### do wilcoxon
	
	
	ML_N_genes = length(FPKM_df_ML[,1])
	
	if(ML_N_genes == 0){
		ML_out_line = c(sp_w,tiss_w,"ML",NA, NA, NA, NA, 0)
		
	} else {
	
		ML_med_SM = median(FPKM_df_ML$Avg_FPKM_SM)
		ML_med_AF = median(FPKM_df_ML$Avg_FPKM_AF)
		ML_med_SF = median(FPKM_df_ML$Avg_FPKM_SF)
		ML_wil <- wilcox.test(FPKM_df_ML$Avg_FPKM_SM, FPKM_df_ML$Avg_FPKM_AF, paired = TRUE)$p.value
	
		ML_out_line = c(sp_w,tiss_w,"ML",ML_med_SM, ML_med_SF, ML_med_AF, ML_wil, ML_N_genes)
	}
	
	
	
	FL_N_genes = length(FPKM_df_FL[,1])
	
	if(FL_N_genes == 0){
		FL_out_line = c(sp_w,tiss_w,"FL",NA, NA, NA, NA, 0)
		
	} else {
	
		FL_med_SM = median(FPKM_df_FL$Avg_FPKM_SM)
		FL_med_AF = median(FPKM_df_FL$Avg_FPKM_AF)
		FL_med_SF = median(FPKM_df_FL$Avg_FPKM_SF)
		FL_wil <- wilcox.test(FPKM_df_FL$Avg_FPKM_SM, FPKM_df_FL$Avg_FPKM_AF, paired = TRUE)$p.value
	
		FL_out_line = c(sp_w,tiss_w,"FL",FL_med_SM, FL_med_SF, FL_med_AF, FL_wil, FL_N_genes)
	}
		
	
	wilcox_out_df <- as.data.frame(rbind(FL_out_line,ML_out_line), row.names = F)
	colnames(wilcox_out_df) <- c("sp","tiss","lim_type", "SM_exp", "SF_exp", "AF_exp", "wil_p", "N_genes")
	
	out_df_list <- list("FL_genes" = FPKM_df_FL, "ML_genes" = FPKM_df_ML, "wilcox_df" = wilcox_out_df)	
	return(out_df_list)
}



FKPM_exp = 2
FKPM_notexp = 0
N_libs = 3

Tbi_WB_FL <- get_SL_genes(Tbi_WB_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tbi", "WB")$FL_genes
Tbi_WB_ML <- get_SL_genes(Tbi_WB_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tbi", "WB")$ML_genes
Tce_WB_FL <- get_SL_genes(Tce_WB_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tce", "WB")$FL_genes
Tce_WB_ML <- get_SL_genes(Tce_WB_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tce", "WB")$ML_genes
Tcm_WB_FL <- get_SL_genes(Tcm_WB_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tcm", "WB")$FL_genes
Tcm_WB_ML <- get_SL_genes(Tcm_WB_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tcm", "WB")$ML_genes
Tpa_WB_FL <- get_SL_genes(Tpa_WB_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tpa", "WB")$FL_genes
Tpa_WB_ML <- get_SL_genes(Tpa_WB_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tpa", "WB")$ML_genes
Tps_WB_FL <- get_SL_genes(Tps_WB_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tps", "WB")$FL_genes
Tps_WB_ML <- get_SL_genes(Tps_WB_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tps", "WB")$ML_genes

Tbi_RT_FL <- get_SL_genes(Tbi_RT_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tbi", "RT")$FL_genes
Tbi_RT_ML <- get_SL_genes(Tbi_RT_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tbi", "RT")$ML_genes
Tce_RT_FL <- get_SL_genes(Tce_RT_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tce", "RT")$FL_genes
Tce_RT_ML <- get_SL_genes(Tce_RT_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tce", "RT")$ML_genes
Tcm_RT_FL <- get_SL_genes(Tcm_RT_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tcm", "RT")$FL_genes
Tcm_RT_ML <- get_SL_genes(Tcm_RT_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tcm", "RT")$ML_genes
Tpa_RT_FL <- get_SL_genes(Tpa_RT_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tpa", "RT")$FL_genes
Tpa_RT_ML <- get_SL_genes(Tpa_RT_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tpa", "RT")$ML_genes
Tps_RT_FL <- get_SL_genes(Tps_RT_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tps", "RT")$FL_genes
Tps_RT_ML <- get_SL_genes(Tps_RT_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tps", "RT")$ML_genes

Tbi_LG_FL <- get_SL_genes(Tbi_LG_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tbi", "LG")$FL_genes
Tbi_LG_ML <- get_SL_genes(Tbi_LG_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tbi", "LG")$ML_genes
Tce_LG_FL <- get_SL_genes(Tce_LG_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tce", "LG")$FL_genes
Tce_LG_ML <- get_SL_genes(Tce_LG_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tce", "LG")$ML_genes
Tcm_LG_FL <- get_SL_genes(Tcm_LG_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tcm", "LG")$FL_genes
Tcm_LG_ML <- get_SL_genes(Tcm_LG_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tcm", "LG")$ML_genes
Tpa_LG_FL <- get_SL_genes(Tpa_LG_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tpa", "LG")$FL_genes
Tpa_LG_ML <- get_SL_genes(Tpa_LG_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tpa", "LG")$ML_genes
Tps_LG_FL <- get_SL_genes(Tps_LG_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tps", "LG")$FL_genes
Tps_LG_ML <- get_SL_genes(Tps_LG_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tps", "LG")$ML_genes

make_long_table_ML <- function(df_ML,sp_u){
	print(length(df_ML[,1]))
	df_ML_t1 = as.data.frame(cbind(
		c(df_ML$Avg_FPKM_SM, df_ML$Avg_FPKM_AF),
		c(as.character(df_ML$V1), as.character(df_ML$V1)),
		c(rep("SM", length(df_ML[,1])),  rep("AF", length(df_ML[,1]))),
		c(rep(sp_u, length(df_ML[,1]) * 2 ))))
	
	colnames(df_ML_t1)	<- c("Avg_FPKM","genename","sex","sp")
	df_ML_t1$sex_ord <- ordered(df_ML_t1$sex, levels = c("SM", "AF"))
	
	df_ML_t1$Avg_FPKM <- as.numeric(as.character(df_ML_t1$Avg_FPKM))
	return(df_ML_t1)

}

make_long_table_ML_forfemales <- function(df_ML,sp_u){
	print(length(df_ML[,1]))
	df_ML_t1 = as.data.frame(cbind(
		c(df_ML$Avg_FPKM_SF, df_ML$Avg_FPKM_AF),
		c(as.character(df_ML$V1), as.character(df_ML$V1)),
		c(rep("SF", length(df_ML[,1])),  rep("AF", length(df_ML[,1]))),
		c(rep(sp_u, length(df_ML[,1]) * 2 ))))
	
	colnames(df_ML_t1)	<- c("Avg_FPKM","genename","sex","sp")
	df_ML_t1$sex_ord <- ordered(df_ML_t1$sex, levels = c("SF", "AF"))
	
	df_ML_t1$Avg_FPKM <- as.numeric(as.character(df_ML_t1$Avg_FPKM))
	return(df_ML_t1)

}


make_long_table_FL <- function(df_FL,sp_u){
	print(length(df_FL[,1]))
	df_FL_t1 = as.data.frame(cbind(
		c(df_FL$Avg_FPKM_SF, df_FL$Avg_FPKM_AF),
		c(as.character(df_FL$V1), as.character(df_FL$V1)),
		c(rep("SF", length(df_FL[,1])),  rep("AF", length(df_FL[,1]))),
		c(rep(sp_u, length(df_FL[,1]) * 2 ))))
	
	colnames(df_FL_t1)	<- c("Avg_FPKM","genename","sex","sp")
	df_FL_t1$sex_ord <- ordered(df_FL_t1$sex, levels = c("SF", "AF"))	
	df_FL_t1$Avg_FPKM <- as.numeric(as.character(df_FL_t1$Avg_FPKM))
	return(df_FL_t1)
}




Tbi_WB_ML_l <- make_long_table_ML(Tbi_WB_ML, "Tbi")
Tbi_WB_FL_l <- make_long_table_FL(Tbi_WB_FL, "Tbi")
Tce_WB_ML_l <- make_long_table_ML(Tce_WB_ML, "Tce")
Tce_WB_FL_l <- make_long_table_FL(Tce_WB_FL, "Tce")
Tcm_WB_ML_l <- make_long_table_ML(Tcm_WB_ML, "Tcm")
Tcm_WB_FL_l <- make_long_table_FL(Tcm_WB_FL, "Tcm")
Tpa_WB_ML_l <- make_long_table_ML(Tpa_WB_ML, "Tpa")
Tpa_WB_FL_l <- make_long_table_FL(Tpa_WB_FL, "Tpa")
Tps_WB_ML_l <- make_long_table_ML(Tps_WB_ML, "Tps")
Tps_WB_FL_l <- make_long_table_FL(Tps_WB_FL, "Tps")

Tbi_RT_ML_l <- make_long_table_ML(Tbi_RT_ML, "Tbi")
Tbi_RT_FL_l <- make_long_table_FL(Tbi_RT_FL, "Tbi")
Tce_RT_ML_l <- make_long_table_ML(Tce_RT_ML, "Tce")
Tce_RT_FL_l <- make_long_table_FL(Tce_RT_FL, "Tce")
Tcm_RT_ML_l <- make_long_table_ML(Tcm_RT_ML, "Tcm")
Tcm_RT_FL_l <- make_long_table_FL(Tcm_RT_FL, "Tcm")
Tpa_RT_ML_l <- make_long_table_ML(Tpa_RT_ML, "Tpa")
Tpa_RT_FL_l <- make_long_table_FL(Tpa_RT_FL, "Tpa")
Tps_RT_ML_l <- make_long_table_ML(Tps_RT_ML, "Tps")
Tps_RT_FL_l <- make_long_table_FL(Tps_RT_FL, "Tps")

Tbi_LG_ML_l <- make_long_table_ML(Tbi_LG_ML, "Tbi")
Tbi_LG_FL_l <- make_long_table_FL(Tbi_LG_FL, "Tbi")
Tce_LG_ML_l <- make_long_table_ML(Tce_LG_ML, "Tce")
Tce_LG_FL_l <- make_long_table_FL(Tce_LG_FL, "Tce")
Tcm_LG_ML_l <- make_long_table_ML(Tcm_LG_ML, "Tcm")
Tcm_LG_FL_l <- make_long_table_FL(Tcm_LG_FL, "Tcm")
Tpa_LG_ML_l <- make_long_table_ML(Tpa_LG_ML, "Tpa")
Tpa_LG_FL_l <- make_long_table_FL(Tpa_LG_FL, "Tpa")
Tps_LG_ML_l <- make_long_table_ML(Tps_LG_ML, "Tps")
Tps_LG_FL_l <- make_long_table_FL(Tps_LG_FL, "Tps")




### male lim - just female vals
Tbi_WB_ML_ff_l <- make_long_table_ML_forfemales(Tbi_WB_ML, "Tbi")
Tce_WB_ML_ff_l <- make_long_table_ML_forfemales(Tce_WB_ML, "Tce")
Tcm_WB_ML_ff_l <- make_long_table_ML_forfemales(Tcm_WB_ML, "Tcm")
Tpa_WB_ML_ff_l <- make_long_table_ML_forfemales(Tpa_WB_ML, "Tpa")
Tps_WB_ML_ff_l <- make_long_table_ML_forfemales(Tps_WB_ML, "Tps")

Tbi_RT_ML_ff_l <- make_long_table_ML_forfemales(Tbi_RT_ML, "Tbi")
Tce_RT_ML_ff_l <- make_long_table_ML_forfemales(Tce_RT_ML, "Tce")
Tcm_RT_ML_ff_l <- make_long_table_ML_forfemales(Tcm_RT_ML, "Tcm")
Tpa_RT_ML_ff_l <- make_long_table_ML_forfemales(Tpa_RT_ML, "Tpa")
Tps_RT_ML_ff_l <- make_long_table_ML_forfemales(Tps_RT_ML, "Tps")

Tbi_LG_ML_ff_l <- make_long_table_ML_forfemales(Tbi_LG_ML, "Tbi")
Tce_LG_ML_ff_l <- make_long_table_ML_forfemales(Tce_LG_ML, "Tce")
Tcm_LG_ML_ff_l <- make_long_table_ML_forfemales(Tcm_LG_ML, "Tcm")
Tpa_LG_ML_ff_l <- make_long_table_ML_forfemales(Tpa_LG_ML, "Tpa")
Tps_LG_ML_ff_l <- make_long_table_ML_forfemales(Tps_LG_ML, "Tps")




#### plot boxplot, line plot and heat maps WITH OUTLIER

plot_box_and_line_2 <- function(ML_df, FL_df, ML_ff_df){

	### some fake plots for when there are no ML / FL genes

	fak_df_FL <- as.data.frame(
	cbind(
	c(7,10,2,6),
	c("g1","g2","g1","g2"),
	c("SF", "SF", "AF","AF"),
	c("SpF","SpF","SpF","SpF")
	))
	colnames(fak_df_FL ) <-  c("Avg_FPKM","genename","sex","sp")
	fak_df_FL$sex_ord <- ordered(fak_df_FL$sex, levels = c("SF", "AF"))	
	fak_df_FL$Avg_FPKM <- as.numeric(as.character(fak_df_FL$Avg_FPKM))

	fak_df_ML <- as.data.frame(
	cbind(
	c(7,10,2,6),
	c("g1","g2","g1","g2"),
	c("SM", "SM", "AF","AF"),
	c("SpF","SpF","SpF","SpF")
	))
	colnames(fak_df_ML ) <-  c("Avg_FPKM","genename","sex","sp")
	fak_df_ML$sex_ord <- ordered(fak_df_ML$sex, levels = c("SM", "AF"))	
	fak_df_ML$Avg_FPKM <- as.numeric(as.character(fak_df_ML$Avg_FPKM))


	FAKE_FL_P1 <- ggplot(fak_df_FL, aes(y = Avg_FPKM , x = sex_ord, color = "black")) +
				scale_color_manual(values=c("white")) +
 				geom_line(aes(group = genename)) +
  				geom_point(colour="white", size=0.5) +
				ylab ("FPKM") +
				xlab ("Sex")  
				
	FAKE_FL_P2 <- FAKE_FL_P1 + theme(legend.position='none', text = element_text(size=14),axis.text.x = element_text(colour="black",size=12), axis.text.y = element_text(colour="black",size=12))	


	FAKE_FL_PA <- ggplot(fak_df_FL, aes(sex_ord,Avg_FPKM)) + 
				theme_classic() +
				geom_boxplot(aes(fill = factor(sex_ord)),position=position_dodge(0.6), width = 0, outlier.size = 0, outlier.shape = NA, colour = "white") +
				#coord_cartesian(ylim=c(0,80)) +
				scale_fill_manual(values=c("firebrick2", "darkorchid2")) +
				ylab ("FPKM") +
				xlab ("Sex") 	

	FAKE_FL_PB <- FAKE_FL_PA + theme(legend.position='none',text = element_text(size=14),axis.text.x = element_text(colour="black",size=12), axis.text.y = element_text(colour="black",size=12))


	FAKE_ML_P1 <- ggplot(fak_df_ML, aes(y = Avg_FPKM , x = sex_ord, color = "black")) +
				scale_color_manual(values=c("white")) +
 				geom_line(aes(group = genename)) +
  				geom_point(colour="white", size=0.5) +
				ylab ("FPKM") +
				xlab ("Sex")  
				
	FAKE_ML_P2 <- FAKE_ML_P1 + theme(legend.position='none', text = element_text(size=14),axis.text.x = element_text(colour="black",size=12), axis.text.y = element_text(colour="black",size=12))	


	FAKE_ML_PA <- ggplot(fak_df_ML, aes(sex_ord,Avg_FPKM)) + 
				theme_classic() +
				geom_boxplot(aes(fill = factor(sex_ord)),position=position_dodge(0.6), width = 0, outlier.size = 0, outlier.shape = NA, colour = "white") +
				#coord_cartesian(ylim=c(0,80)) +
				scale_fill_manual(values=c("firebrick2", "darkorchid2")) +
				ylab ("FPKM") +
				xlab ("Sex") 	

	FAKE_ML_PB <- FAKE_ML_PA + theme(legend.position='none',text = element_text(size=14),axis.text.x = element_text(colour="black",size=12), axis.text.y = element_text(colour="black",size=12))


	### male lim

	if (length(ML_df[,1]) == 0) {
 		print("No male-limited genes")
  
		ML_P2 <- FAKE_ML_P2
		ML_PB <- FAKE_ML_PB 
		ML_ff_P2 <- FAKE_FL_P2
		ML_ff_PB <- FAKE_FL_PB 		
  
	} else {
		print("Some male-limited genes")

		ML_P1 <- ggplot(ML_df, aes(y = Avg_FPKM , x = sex_ord, color = "black")) +
				scale_color_manual(values=c("grey28")) +
 				geom_line(aes(group = genename)) +
  				geom_point(colour="black", size=0.5) +
				ylab ("FPKM") +
				xlab ("Sex")  
		ML_P2 = 	ML_P1 + theme(legend.position='none', text = element_text(size=14),axis.text.x = element_text(colour="black",size=12), axis.text.y = element_text(colour="black",size=12))	
	
		ML_PA <- ggplot(ML_df, aes(sex_ord,Avg_FPKM)) + 
				theme_classic() +
				geom_boxplot(aes(fill = factor(sex_ord)),position=position_dodge(0.6), width = 0.5, outlier.size = 0.7) +
				#coord_cartesian(ylim=c(0,80)) +
				scale_fill_manual(values=c("royalblue2", "darkorchid2")) +
				ylab ("FPKM") +
				xlab ("Sex") 	

		ML_PB <- ML_PA + theme(legend.position='none',text = element_text(size=14),axis.text.x = element_text(colour="black",size=12), axis.text.y = element_text(colour="black",size=12))


		print("Some male-limited genes")

		ML_ff_P1 <- ggplot(ML_ff_df, aes(y = Avg_FPKM , x = sex_ord, color = "black")) +
				scale_color_manual(values=c("grey28")) +
 				geom_line(aes(group = genename)) +
  				geom_point(colour="black", size=0.5) +
				ylab ("FPKM") +
				xlab ("Sex")  
		ML_ff_P2 = 	ML_ff_P1 + theme(legend.position='none', text = element_text(size=14),axis.text.x = element_text(colour="black",size=12), axis.text.y = element_text(colour="black",size=12))	
	
		ML_ff_PA <- ggplot(ML_ff_df, aes(sex_ord,Avg_FPKM)) + 
				theme_classic() +
				geom_boxplot(aes(fill = factor(sex_ord)),position=position_dodge(0.6), width = 0.5, outlier.size = 0.7) +
				#coord_cartesian(ylim=c(0,80)) +
				scale_fill_manual(values=c("royalblue2", "darkorchid2")) +
				ylab ("FPKM") +
				xlab ("Sex") 	

		ML_ff_PB <- ML_ff_PA + theme(legend.position='none',text = element_text(size=14),axis.text.x = element_text(colour="black",size=12), axis.text.y = element_text(colour="black",size=12))

	}

	##### Female lim

	if (length(FL_df[,1]) == 0) {
 		print("No female-limited genes")
  
		FL_P2 <- FAKE_FL_P2
		FL_PB <- FAKE_FL_PB 
  
	} else {
		print("Some female-limited genes")
		
		FL_P1 <- ggplot(FL_df, aes(y = Avg_FPKM , x = sex_ord, color = "black")) +
				scale_color_manual(values=c("grey28")) +
 				geom_line(aes(group = genename)) +
  				geom_point(colour="black", size=0.5) +
				ylab ("FPKM") +
				xlab ("Sex") 	
		FL_P2 = 	FL_P1 + theme(legend.position='none', text = element_text(size=14),axis.text.x = element_text(colour="black",size=12), axis.text.y = element_text(colour="black",size=12))	
	
		FL_PA <- ggplot(FL_df, aes(sex_ord,Avg_FPKM)) + 
				theme_classic() +
				geom_boxplot(aes(fill = factor(sex_ord)),position=position_dodge(0.6), width = 0.5, outlier.size = 0.7) +
				#coord_cartesian(ylim=c(0,80)) +
				scale_fill_manual(values=c("firebrick2", "darkorchid2")) +
				ylab ("FPKM") +
				xlab ("Sex") 	

		FL_PB <- FL_PA + theme(legend.position='none',text = element_text(size=14),axis.text.x = element_text(colour="black",size=12), axis.text.y = element_text(colour="black",size=12))

	}
	
	out_plots <- list("ML_line" = ML_P2, "ML_box" = ML_PB, "FL_line" = FL_P2, "FL_box" = FL_PB, "ML_ff_line" = ML_ff_P2, "ML_ff_box" = ML_ff_PB)
	return(out_plots )
}


Tbi_WB_SL_plots_2 <- plot_box_and_line_2(Tbi_WB_ML_l,Tbi_WB_FL_l,Tbi_WB_ML_ff_l)
Tce_WB_SL_plots_2 <- plot_box_and_line_2(Tce_WB_ML_l,Tce_WB_FL_l,Tce_WB_ML_ff_l)
Tcm_WB_SL_plots_2 <- plot_box_and_line_2(Tcm_WB_ML_l,Tcm_WB_FL_l,Tcm_WB_ML_ff_l)
Tpa_WB_SL_plots_2 <- plot_box_and_line_2(Tpa_WB_ML_l,Tpa_WB_FL_l,Tpa_WB_ML_ff_l)
Tps_WB_SL_plots_2 <- plot_box_and_line_2(Tps_WB_ML_l,Tps_WB_FL_l,Tps_WB_ML_ff_l)

Tbi_RT_SL_plots_2 <- plot_box_and_line_2(Tbi_RT_ML_l,Tbi_RT_FL_l,Tbi_RT_ML_ff_l)
Tce_RT_SL_plots_2 <- plot_box_and_line_2(Tce_RT_ML_l,Tce_RT_FL_l,Tce_RT_ML_ff_l)
Tcm_RT_SL_plots_2 <- plot_box_and_line_2(Tcm_RT_ML_l,Tcm_RT_FL_l,Tcm_RT_ML_ff_l)
Tpa_RT_SL_plots_2 <- plot_box_and_line_2(Tpa_RT_ML_l,Tpa_RT_FL_l,Tpa_RT_ML_ff_l)
Tps_RT_SL_plots_2 <- plot_box_and_line_2(Tps_RT_ML_l,Tps_RT_FL_l,Tps_RT_ML_ff_l)

Tbi_LG_SL_plots_2 <- plot_box_and_line_2(Tbi_LG_ML_l,Tbi_LG_FL_l,Tbi_LG_ML_ff_l)
Tce_LG_SL_plots_2 <- plot_box_and_line_2(Tce_LG_ML_l,Tce_LG_FL_l,Tce_LG_ML_ff_l)
Tcm_LG_SL_plots_2 <- plot_box_and_line_2(Tcm_LG_ML_l,Tcm_LG_FL_l,Tcm_LG_ML_ff_l)
Tpa_LG_SL_plots_2 <- plot_box_and_line_2(Tpa_LG_ML_l,Tpa_LG_FL_l,Tpa_LG_ML_ff_l)
Tps_LG_SL_plots_2 <- plot_box_and_line_2(Tps_LG_ML_l,Tps_LG_FL_l,Tps_LG_ML_ff_l)


All_sp_plots_2_WB <- plot_grid(
Tbi_WB_SL_plots_2$FL_line, Tbi_WB_SL_plots_2$FL_box,Tbi_WB_SL_plots_2$ML_ff_line, Tbi_WB_SL_plots_2$ML_ff_box,Tbi_WB_SL_plots_2$ML_line, Tbi_WB_SL_plots_2$ML_box,
Tce_WB_SL_plots_2$FL_line, Tce_WB_SL_plots_2$FL_box,Tce_WB_SL_plots_2$ML_ff_line, Tce_WB_SL_plots_2$ML_ff_box,Tce_WB_SL_plots_2$ML_line, Tce_WB_SL_plots_2$ML_box,
Tps_WB_SL_plots_2$FL_line, Tps_WB_SL_plots_2$FL_box,Tps_WB_SL_plots_2$ML_ff_line, Tps_WB_SL_plots_2$ML_ff_box,Tps_WB_SL_plots_2$ML_line, Tps_WB_SL_plots_2$ML_box,
Tcm_WB_SL_plots_2$FL_line, Tcm_WB_SL_plots_2$FL_box,Tcm_WB_SL_plots_2$ML_ff_line, Tcm_WB_SL_plots_2$ML_ff_box,Tcm_WB_SL_plots_2$ML_line, Tcm_WB_SL_plots_2$ML_box,
Tpa_WB_SL_plots_2$FL_line, Tpa_WB_SL_plots_2$FL_box,Tpa_WB_SL_plots_2$ML_ff_line, Tpa_WB_SL_plots_2$ML_ff_box,Tpa_WB_SL_plots_2$ML_line, Tpa_WB_SL_plots_2$ML_box,
ncol = 6, nrow = 5, align =  "h", rel_widths = c(1, 0.6))
title_WB <- ggdraw() + draw_label(paste("WB |","FKPM_exp = ",FKPM_exp, "| FKPM_notexp = ",FKPM_notexp, "| N_libs = ",N_libs) , fontface = 'bold')
All_sp_plots_2_WB_t <- plot_grid(title_WB, All_sp_plots_2_WB , ncol = 1, rel_heights = c(0.04, 1)) # rel_heights values control title margins

All_sp_plots_2_RT <- plot_grid(
Tbi_RT_SL_plots_2$FL_line, Tbi_RT_SL_plots_2$FL_box,Tbi_RT_SL_plots_2$ML_ff_line, Tbi_RT_SL_plots_2$ML_ff_box,Tbi_RT_SL_plots_2$ML_line, Tbi_RT_SL_plots_2$ML_box,
Tce_RT_SL_plots_2$FL_line, Tce_RT_SL_plots_2$FL_box,Tce_RT_SL_plots_2$ML_ff_line, Tce_RT_SL_plots_2$ML_ff_box,Tce_RT_SL_plots_2$ML_line, Tce_RT_SL_plots_2$ML_box,
Tps_RT_SL_plots_2$FL_line, Tps_RT_SL_plots_2$FL_box,Tps_RT_SL_plots_2$ML_ff_line, Tps_RT_SL_plots_2$ML_ff_box,Tps_RT_SL_plots_2$ML_line, Tps_RT_SL_plots_2$ML_box,
Tcm_RT_SL_plots_2$FL_line, Tcm_RT_SL_plots_2$FL_box,Tcm_RT_SL_plots_2$ML_ff_line, Tcm_RT_SL_plots_2$ML_ff_box,Tcm_RT_SL_plots_2$ML_line, Tcm_RT_SL_plots_2$ML_box,
Tpa_RT_SL_plots_2$FL_line, Tpa_RT_SL_plots_2$FL_box,Tpa_RT_SL_plots_2$ML_ff_line, Tpa_RT_SL_plots_2$ML_ff_box,Tpa_RT_SL_plots_2$ML_line, Tpa_RT_SL_plots_2$ML_box,
ncol = 6, nrow = 5, align =  "hv", rel_widths = c(1, 0.6))
title_RT <- ggdraw() + draw_label(paste("RT |","FKPM_exp = ",FKPM_exp, "| FKPM_notexp = ",FKPM_notexp, "| N_libs = ",N_libs) , fontface = 'bold')
All_sp_plots_2_RT_t <- plot_grid(title_RT, All_sp_plots_2_RT , ncol = 1, rel_heights = c(0.04, 1)) # rel_heights values control title margins


All_sp_plots_2_LG <- plot_grid(
Tbi_LG_SL_plots_2$FL_line, Tbi_LG_SL_plots_2$FL_box,Tbi_LG_SL_plots_2$ML_ff_line, Tbi_LG_SL_plots_2$ML_ff_box,Tbi_LG_SL_plots_2$ML_line, Tbi_LG_SL_plots_2$ML_box,
Tce_LG_SL_plots_2$FL_line, Tce_LG_SL_plots_2$FL_box,Tce_LG_SL_plots_2$ML_ff_line, Tce_LG_SL_plots_2$ML_ff_box,Tce_LG_SL_plots_2$ML_line, Tce_LG_SL_plots_2$ML_box,
Tps_LG_SL_plots_2$FL_line, Tps_LG_SL_plots_2$FL_box,Tps_LG_SL_plots_2$ML_ff_line, Tps_LG_SL_plots_2$ML_ff_box,Tps_LG_SL_plots_2$ML_line, Tps_LG_SL_plots_2$ML_box,
Tcm_LG_SL_plots_2$FL_line, Tcm_LG_SL_plots_2$FL_box,Tcm_LG_SL_plots_2$ML_ff_line, Tcm_LG_SL_plots_2$ML_ff_box,Tcm_LG_SL_plots_2$ML_line, Tcm_LG_SL_plots_2$ML_box,
Tpa_LG_SL_plots_2$FL_line, Tpa_LG_SL_plots_2$FL_box,Tpa_LG_SL_plots_2$ML_ff_line, Tpa_LG_SL_plots_2$ML_ff_box,Tpa_LG_SL_plots_2$ML_line, Tpa_LG_SL_plots_2$ML_box,
ncol = 6, nrow = 5, align =  "hv", rel_widths = c(1, 0.6))
title_LG <- ggdraw() + draw_label(paste("LG |","FKPM_exp = ",FKPM_exp, "| FKPM_notexp = ",FKPM_notexp, "| N_libs = ",N_libs) , fontface = 'bold')
All_sp_plots_2_LG_t <- plot_grid(title_LG, All_sp_plots_2_LG , ncol = 1, rel_heights = c(0.04, 1)) # rel_heights values control title margins


### output
png(filename = paste("All_sp_plots_WB_2","FKPM_exp",FKPM_exp, "FKPM_notexp",FKPM_notexp, "N_libs",N_libs,".png", sep = "_"), width = 12, height = 16, units = "in", bg = "white", res = 300)
All_sp_plots_2_WB_t
dev.off()
getwd() ## where has my plot gone....

png(filename = paste("All_sp_plots_RT_2","FKPM_exp",FKPM_exp, "FKPM_notexp",FKPM_notexp, "N_libs",N_libs,".png", sep = "_"), width = 12, height = 16, units = "in", bg = "white", res = 300)
All_sp_plots_2_RT_t
dev.off()
getwd() ## where has my plot gone....

png(filename = paste("All_sp_plots_LG_2","FKPM_exp",FKPM_exp, "FKPM_notexp",FKPM_notexp, "N_libs",N_libs,".png", sep = "_"), width = 12, height = 16, units = "in", bg = "white", res = 300)
All_sp_plots_2_LG_t
dev.off()
getwd() ## where has my plot gone....



######################################################################################################################################################################## STATS

SL_wilcox <- as.data.frame(
rbind(
get_SL_genes(Tbi_WB_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tbi", "WB")$wilcox_df,
get_SL_genes(Tce_WB_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tce", "WB")$wilcox_df,
get_SL_genes(Tps_WB_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tps", "WB")$wilcox_df,
get_SL_genes(Tcm_WB_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tcm", "WB")$wilcox_df,
get_SL_genes(Tpa_WB_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tpa", "WB")$wilcox_df, 
get_SL_genes(Tbi_RT_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tbi", "RT")$wilcox_df,
get_SL_genes(Tce_RT_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tce", "RT")$wilcox_df,
get_SL_genes(Tps_RT_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tps", "RT")$wilcox_df,
get_SL_genes(Tcm_RT_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tcm", "RT")$wilcox_df,
get_SL_genes(Tpa_RT_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tpa", "RT")$wilcox_df, 
get_SL_genes(Tbi_LG_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tbi", "LG")$wilcox_df,
get_SL_genes(Tce_LG_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tce", "LG")$wilcox_df,
get_SL_genes(Tps_LG_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tps", "LG")$wilcox_df,
get_SL_genes(Tcm_LG_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tcm", "LG")$wilcox_df,
get_SL_genes(Tpa_LG_F_FPKM_1, FKPM_exp, FKPM_notexp, N_libs, "Tpa", "LG")$wilcox_df
))

SL_wilcox$wil_p <- as.numeric(as.character((SL_wilcox$wil_p)))
SL_wilcox$wil_FDR <- p.adjust(SL_wilcox$wil_p, method = "BH")

write.csv(SL_wilcox, file= paste("SL_wilcox","FKPM_exp",FKPM_exp, "FKPM_notexp",FKPM_notexp, "N_libs",N_libs,".csv", sep = "_"), row.names=FALSE)


########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), "sex_bias_edgeR_Sex_limited_expression.R_sessionInfo.txt")


