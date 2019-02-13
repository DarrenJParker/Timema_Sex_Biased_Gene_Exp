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
### sexuals 

rawdata_Tbi_longest_iso <- read.csv("Data/Read_counts/to_Tbi_longest_iso_Arthropoda+Mixed+NOBLASTHIT_K2edge.counts", check.names=FALSE, stringsAsFactors=FALSE)
rawdata_Tce_longest_iso <- read.csv("Data/Read_counts/to_Tce_longest_iso_Arthropoda+Mixed+NOBLASTHIT_K2edge.counts", check.names=FALSE, stringsAsFactors=FALSE)
rawdata_Tcm_longest_iso <- read.csv("Data/Read_counts/to_Tcm_longest_iso_Arthropoda+Mixed+NOBLASTHIT_K2edge.counts", check.names=FALSE, stringsAsFactors=FALSE)
rawdata_Tpa_longest_iso <- read.csv("Data/Read_counts/to_Tpa_longest_iso_Arthropoda+Mixed+NOBLASTHIT_K2edge.counts", check.names=FALSE, stringsAsFactors=FALSE)
rawdata_Tps_longest_iso <- read.csv("Data/Read_counts/to_Tps_longest_iso_Arthropoda+Mixed+NOBLASTHIT_K2edge.counts", check.names=FALSE, stringsAsFactors=FALSE)

head(rawdata_Tbi_longest_iso)
head(rawdata_Tce_longest_iso)
head(rawdata_Tcm_longest_iso)
head(rawdata_Tpa_longest_iso)
head(rawdata_Tps_longest_iso)

colnames(rawdata_Tbi_longest_iso)
length(colnames(rawdata_Tbi_longest_iso))


dir.create("Output")
dir.create("Output/DE_sexual_ref")

setwd("Output/DE_sexual_ref")


###### get cpm and FPKMs

#### WB

## for SB
y_WB_longest_iso_Tbi_UF <- DGEList(counts=rawdata_Tbi_longest_iso[,c(
"Tbi_SF_WB_Md_Re1_toTbi","Tbi_SF_WB_Md_Re2_toTbi","Tbi_SF_WB_Md_Re3_toTbi",
"Tbi_SM_WB_Md_Re1_toTbi","Tbi_SM_WB_Md_Re2_toTbi","Tbi_SM_WB_Md_Re3_toTbi",
"Tte_AF_WB_Vi_Re1_toTbi","Tte_AF_WB_Vi_Re2_toTbi","Tte_AF_WB_Vi_Re3_toTbi"
)], genes=rawdata_Tbi_longest_iso[,1:1])

y_WB_longest_iso_Tce_UF <- DGEList(counts=rawdata_Tce_longest_iso[,c(
"Tce_SF_WB_Md_Re1_toTce","Tce_SF_WB_Md_Re2_toTce","Tce_SF_WB_Md_Re3_toTce",
"Tce_SM_WB_Md_Re1_toTce","Tce_SM_WB_Md_Re2_toTce","Tce_SM_WB_Md_Re3_toTce",
"Tms_AF_WB_Vi_Re1_toTce","Tms_AF_WB_Vi_Re2_toTce","Tms_AF_WB_Vi_Re3_toTce"
)], genes=rawdata_Tce_longest_iso[,1:1])

y_WB_longest_iso_Tcm_UF <- DGEList(counts=rawdata_Tcm_longest_iso[,c(
"Tcm_SF_WB_Md_Re1_toTcm","Tcm_SF_WB_Md_Re2_toTcm","Tcm_SF_WB_Md_Re3_toTcm",
"Tcm_SM_WB_Md_Re1_toTcm","Tcm_SM_WB_Md_Re2_toTcm","Tcm_SM_WB_Md_Re3_toTcm",
"Tsi_AF_WB_Vi_Re1_toTcm","Tsi_AF_WB_Vi_Re2_toTcm","Tsi_AF_WB_Vi_Re3_toTcm"
)], genes=rawdata_Tcm_longest_iso[,1:1])

y_WB_longest_iso_Tpa_UF <- DGEList(counts=rawdata_Tpa_longest_iso[,c(
"Tpa_SF_WB_Md_Re1_toTpa","Tpa_SF_WB_Md_Re2_toTpa","Tpa_SF_WB_Md_Re3_toTpa",
"Tpa_SM_WB_Md_Re1_toTpa","Tpa_SM_WB_Md_Re2_toTpa","Tpa_SM_WB_Md_Re3_toTpa",
"Tge_AF_WB_Vi_Re1_toTpa","Tge_AF_WB_Vi_Re2_toTpa","Tge_AF_WB_Vi_Re3_toTpa"
)], genes=rawdata_Tpa_longest_iso[,1:1])

y_WB_longest_iso_Tps_UF <- DGEList(counts=rawdata_Tps_longest_iso[,c(
"Tps_SF_WB_Md_Re1_toTps","Tps_SF_WB_Md_Re2_toTps","Tps_SF_WB_Md_Re3_toTps",
"Tps_SM_WB_Md_Re1_toTps","Tps_SM_WB_Md_Re2_toTps","Tps_SM_WB_Md_Re3_toTps",
"Tdi_AF_WB_Vi_Re1_toTps","Tdi_AF_WB_Vi_Re2_toTps","Tdi_AF_WB_Vi_Re3_toTps"
)], genes=rawdata_Tps_longest_iso[,1:1])

#### RT

## for SB
y_RT_longest_iso_Tbi_UF <- DGEList(counts=rawdata_Tbi_longest_iso[,c(
"Tbi_SF_RT_Md_Re1_toTbi","Tbi_SF_RT_Md_Re2_toTbi","Tbi_SF_RT_Md_Re3_toTbi",
"Tbi_SM_RT_Md_Re1_toTbi","Tbi_SM_RT_Md_Re2_toTbi","Tbi_SM_RT_Md_Re3_toTbi",
"Tte_AF_RT_Vi_Re1_toTbi","Tte_AF_RT_Vi_Re2_toTbi","Tte_AF_RT_Vi_Re3_toTbi"
)], genes=rawdata_Tbi_longest_iso[,1:1])

y_RT_longest_iso_Tce_UF <- DGEList(counts=rawdata_Tce_longest_iso[,c(
"Tce_SF_RT_Md_Re1_toTce","Tce_SF_RT_Md_Re2_toTce","Tce_SF_RT_Md_Re3_toTce",
"Tce_SM_RT_Md_Re1_toTce","Tce_SM_RT_Md_Re2_toTce","Tce_SM_RT_Md_Re3_toTce",
"Tms_AF_RT_Vi_Re1_toTce","Tms_AF_RT_Vi_Re2_toTce","Tms_AF_RT_Vi_Re3_toTce"
)], genes=rawdata_Tce_longest_iso[,1:1])

y_RT_longest_iso_Tcm_UF <- DGEList(counts=rawdata_Tcm_longest_iso[,c(
"Tcm_SF_RT_Md_Re1_toTcm","Tcm_SF_RT_Md_Re2_toTcm","Tcm_SF_RT_Md_Re3_toTcm",
"Tcm_SM_RT_Md_Re1_toTcm","Tcm_SM_RT_Md_Re2_toTcm","Tcm_SM_RT_Md_Re3_toTcm",
"Tsi_AF_RT_Vi_Re1_toTcm","Tsi_AF_RT_Vi_Re2_toTcm","Tsi_AF_RT_Vi_Re3_toTcm"
)], genes=rawdata_Tcm_longest_iso[,1:1])

y_RT_longest_iso_Tpa_UF <- DGEList(counts=rawdata_Tpa_longest_iso[,c(
"Tpa_SF_RT_Md_Re1_toTpa","Tpa_SF_RT_Md_Re2_toTpa","Tpa_SF_RT_Md_Re3_toTpa",
"Tpa_SM_RT_Md_Re1_toTpa","Tpa_SM_RT_Md_Re2_toTpa","Tpa_SM_RT_Md_Re3_toTpa",
"Tge_AF_RT_Vi_Re1_toTpa","Tge_AF_RT_Vi_Re2_toTpa","Tge_AF_RT_Vi_Re3_toTpa"
)], genes=rawdata_Tpa_longest_iso[,1:1])

y_RT_longest_iso_Tps_UF <- DGEList(counts=rawdata_Tps_longest_iso[,c(
"Tps_SF_RT_Md_Re1_toTps","Tps_SF_RT_Md_Re2_toTps","Tps_SF_RT_Md_Re3_toTps",
"Tps_SM_RT_Md_Re1_toTps","Tps_SM_RT_Md_Re2_toTps","Tps_SM_RT_Md_Re3_toTps",
"Tdi_AF_RT_Vi_Re1_toTps","Tdi_AF_RT_Vi_Re2_toTps","Tdi_AF_RT_Vi_Re3_toTps"
)], genes=rawdata_Tps_longest_iso[,1:1])


#### LG

## for SB
y_LG_longest_iso_Tbi_UF <- DGEList(counts=rawdata_Tbi_longest_iso[,c(
"Tbi_SF_LG_Md_Re1_toTbi","Tbi_SF_LG_Md_Re2_toTbi","Tbi_SF_LG_Md_Re3_toTbi",
"Tbi_SM_LG_Md_Re1_toTbi","Tbi_SM_LG_Md_Re2_toTbi","Tbi_SM_LG_Md_Re3_toTbi",
"Tte_AF_LG_Vi_Re1_toTbi","Tte_AF_LG_Vi_Re2_toTbi","Tte_AF_LG_Vi_Re3_toTbi"
)], genes=rawdata_Tbi_longest_iso[,1:1])

y_LG_longest_iso_Tce_UF <- DGEList(counts=rawdata_Tce_longest_iso[,c(
"Tce_SF_LG_Md_Re1_toTce","Tce_SF_LG_Md_Re2_toTce","Tce_SF_LG_Md_Re3_toTce",
"Tce_SM_LG_Md_Re1_toTce","Tce_SM_LG_Md_Re2_toTce","Tce_SM_LG_Md_Re3_toTce",
"Tms_AF_LG_Vi_Re1_toTce","Tms_AF_LG_Vi_Re2_toTce","Tms_AF_LG_Vi_Re3_toTce"
)], genes=rawdata_Tce_longest_iso[,1:1])

y_LG_longest_iso_Tcm_UF <- DGEList(counts=rawdata_Tcm_longest_iso[,c(
"Tcm_SF_LG_Md_Re1_toTcm","Tcm_SF_LG_Md_Re2_toTcm","Tcm_SF_LG_Md_Re3_toTcm",
"Tcm_SM_LG_Md_Re1_toTcm","Tcm_SM_LG_Md_Re2_toTcm","Tcm_SM_LG_Md_Re3_toTcm",
"Tsi_AF_LG_Vi_Re1_toTcm","Tsi_AF_LG_Vi_Re2_toTcm","Tsi_AF_LG_Vi_Re3_toTcm"
)], genes=rawdata_Tcm_longest_iso[,1:1])

y_LG_longest_iso_Tpa_UF <- DGEList(counts=rawdata_Tpa_longest_iso[,c(
"Tpa_SF_LG_Md_Re1_toTpa","Tpa_SF_LG_Md_Re2_toTpa","Tpa_SF_LG_Md_Re3_toTpa",
"Tpa_SM_LG_Md_Re1_toTpa","Tpa_SM_LG_Md_Re2_toTpa","Tpa_SM_LG_Md_Re3_toTpa",
"Tge_AF_LG_Vi_Re1_toTpa","Tge_AF_LG_Vi_Re2_toTpa","Tge_AF_LG_Vi_Re3_toTpa"
)], genes=rawdata_Tpa_longest_iso[,1:1])

y_LG_longest_iso_Tps_UF <- DGEList(counts=rawdata_Tps_longest_iso[,c(
"Tps_SF_LG_Md_Re1_toTps","Tps_SF_LG_Md_Re2_toTps","Tps_SF_LG_Md_Re3_toTps",
"Tps_SM_LG_Md_Re1_toTps","Tps_SM_LG_Md_Re2_toTps","Tps_SM_LG_Md_Re3_toTps",
"Tdi_AF_LG_Vi_Re1_toTps","Tdi_AF_LG_Vi_Re2_toTps","Tdi_AF_LG_Vi_Re3_toTps"
)], genes=rawdata_Tps_longest_iso[,1:1])





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


y_Tbi_WB_N <- norm_no_filter(y_WB_longest_iso_Tbi_UF)
y_Tce_WB_N <- norm_no_filter(y_WB_longest_iso_Tce_UF)
y_Tcm_WB_N <- norm_no_filter(y_WB_longest_iso_Tcm_UF)
y_Tpa_WB_N <- norm_no_filter(y_WB_longest_iso_Tpa_UF)
y_Tps_WB_N <- norm_no_filter(y_WB_longest_iso_Tps_UF)

y_Tbi_RT_N <- norm_no_filter(y_RT_longest_iso_Tbi_UF)
y_Tce_RT_N <- norm_no_filter(y_RT_longest_iso_Tce_UF)
y_Tcm_RT_N <- norm_no_filter(y_RT_longest_iso_Tcm_UF)
y_Tpa_RT_N <- norm_no_filter(y_RT_longest_iso_Tpa_UF)
y_Tps_RT_N <- norm_no_filter(y_RT_longest_iso_Tps_UF)

y_Tbi_LG_N <- norm_no_filter(y_LG_longest_iso_Tbi_UF)
y_Tce_LG_N <- norm_no_filter(y_LG_longest_iso_Tce_UF)
y_Tcm_LG_N <- norm_no_filter(y_LG_longest_iso_Tcm_UF)
y_Tpa_LG_N <- norm_no_filter(y_LG_longest_iso_Tpa_UF)
y_Tps_LG_N <- norm_no_filter(y_LG_longest_iso_Tps_UF)



#### calc FPKMs


Tbi_WB_N_FPKM <- rpkm(y_Tbi_WB_N, gene.length=rawdata_Tbi_longest_iso$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tbi_WB_N_FPKM_1  <- as.data.frame(cbind(rawdata_Tbi_longest_iso$Gene_name,Tbi_WB_N_FPKM))
Tce_WB_N_FPKM <- rpkm(y_Tce_WB_N, gene.length=rawdata_Tce_longest_iso$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tce_WB_N_FPKM_1  <- as.data.frame(cbind(rawdata_Tce_longest_iso$Gene_name,Tce_WB_N_FPKM))
Tcm_WB_N_FPKM <- rpkm(y_Tcm_WB_N, gene.length=rawdata_Tcm_longest_iso$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tcm_WB_N_FPKM_1  <- as.data.frame(cbind(rawdata_Tcm_longest_iso$Gene_name,Tcm_WB_N_FPKM))
Tpa_WB_N_FPKM <- rpkm(y_Tpa_WB_N, gene.length=rawdata_Tpa_longest_iso$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tpa_WB_N_FPKM_1  <- as.data.frame(cbind(rawdata_Tpa_longest_iso$Gene_name,Tpa_WB_N_FPKM))
Tps_WB_N_FPKM <- rpkm(y_Tps_WB_N, gene.length=rawdata_Tps_longest_iso$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tps_WB_N_FPKM_1  <- as.data.frame(cbind(rawdata_Tps_longest_iso$Gene_name,Tps_WB_N_FPKM))

Tbi_RT_N_FPKM <- rpkm(y_Tbi_RT_N, gene.length=rawdata_Tbi_longest_iso$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tbi_RT_N_FPKM_1  <- as.data.frame(cbind(rawdata_Tbi_longest_iso$Gene_name,Tbi_RT_N_FPKM))
Tce_RT_N_FPKM <- rpkm(y_Tce_RT_N, gene.length=rawdata_Tce_longest_iso$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tce_RT_N_FPKM_1  <- as.data.frame(cbind(rawdata_Tce_longest_iso$Gene_name,Tce_RT_N_FPKM))
Tcm_RT_N_FPKM <- rpkm(y_Tcm_RT_N, gene.length=rawdata_Tcm_longest_iso$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tcm_RT_N_FPKM_1  <- as.data.frame(cbind(rawdata_Tcm_longest_iso$Gene_name,Tcm_RT_N_FPKM))
Tpa_RT_N_FPKM <- rpkm(y_Tpa_RT_N, gene.length=rawdata_Tpa_longest_iso$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tpa_RT_N_FPKM_1  <- as.data.frame(cbind(rawdata_Tpa_longest_iso$Gene_name,Tpa_RT_N_FPKM))
Tps_RT_N_FPKM <- rpkm(y_Tps_RT_N, gene.length=rawdata_Tps_longest_iso$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tps_RT_N_FPKM_1  <- as.data.frame(cbind(rawdata_Tps_longest_iso$Gene_name,Tps_RT_N_FPKM))

Tbi_LG_N_FPKM <- rpkm(y_Tbi_LG_N, gene.length=rawdata_Tbi_longest_iso$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tbi_LG_N_FPKM_1  <- as.data.frame(cbind(rawdata_Tbi_longest_iso$Gene_name,Tbi_LG_N_FPKM))
Tce_LG_N_FPKM <- rpkm(y_Tce_LG_N, gene.length=rawdata_Tce_longest_iso$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tce_LG_N_FPKM_1  <- as.data.frame(cbind(rawdata_Tce_longest_iso$Gene_name,Tce_LG_N_FPKM))
Tcm_LG_N_FPKM <- rpkm(y_Tcm_LG_N, gene.length=rawdata_Tcm_longest_iso$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tcm_LG_N_FPKM_1  <- as.data.frame(cbind(rawdata_Tcm_longest_iso$Gene_name,Tcm_LG_N_FPKM))
Tpa_LG_N_FPKM <- rpkm(y_Tpa_LG_N, gene.length=rawdata_Tpa_longest_iso$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tpa_LG_N_FPKM_1  <- as.data.frame(cbind(rawdata_Tpa_longest_iso$Gene_name,Tpa_LG_N_FPKM))
Tps_LG_N_FPKM <- rpkm(y_Tps_LG_N, gene.length=rawdata_Tps_longest_iso$gene_len, normalized.lib.sizes=TRUE, log=FALSE)
Tps_LG_N_FPKM_1  <- as.data.frame(cbind(rawdata_Tps_longest_iso$Gene_name,Tps_LG_N_FPKM))


###### export as csvs

write.csv(Tbi_WB_N_FPKM_1, file="Tbi_WB_N_FPKM_longest_iso_1.csv", row.names=FALSE)
write.csv(Tce_WB_N_FPKM_1, file="Tce_WB_N_FPKM_longest_iso_1.csv", row.names=FALSE)
write.csv(Tcm_WB_N_FPKM_1, file="Tcm_WB_N_FPKM_longest_iso_1.csv", row.names=FALSE)
write.csv(Tpa_WB_N_FPKM_1, file="Tpa_WB_N_FPKM_longest_iso_1.csv", row.names=FALSE)
write.csv(Tps_WB_N_FPKM_1, file="Tps_WB_N_FPKM_longest_iso_1.csv", row.names=FALSE)

write.csv(Tbi_RT_N_FPKM_1, file="Tbi_RT_N_FPKM_longest_iso_1.csv", row.names=FALSE)
write.csv(Tce_RT_N_FPKM_1, file="Tce_RT_N_FPKM_longest_iso_1.csv", row.names=FALSE)
write.csv(Tcm_RT_N_FPKM_1, file="Tcm_RT_N_FPKM_longest_iso_1.csv", row.names=FALSE)
write.csv(Tpa_RT_N_FPKM_1, file="Tpa_RT_N_FPKM_longest_iso_1.csv", row.names=FALSE)
write.csv(Tps_RT_N_FPKM_1, file="Tps_RT_N_FPKM_longest_iso_1.csv", row.names=FALSE)

write.csv(Tbi_LG_N_FPKM_1, file="Tbi_LG_N_FPKM_longest_iso_1.csv", row.names=FALSE)
write.csv(Tce_LG_N_FPKM_1, file="Tce_LG_N_FPKM_longest_iso_1.csv", row.names=FALSE)
write.csv(Tcm_LG_N_FPKM_1, file="Tcm_LG_N_FPKM_longest_iso_1.csv", row.names=FALSE)
write.csv(Tpa_LG_N_FPKM_1, file="Tpa_LG_N_FPKM_longest_iso_1.csv", row.names=FALSE)
write.csv(Tps_LG_N_FPKM_1, file="Tps_LG_N_FPKM_longest_iso_1.csv", row.names=FALSE)





###### DE analysis DOING for each tissue seperatly 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

#### WB

## for SB
y_WB_SB_longest_iso_Tbi_UF <- DGEList(counts=rawdata_Tbi_longest_iso[,c("Tbi_SF_WB_Md_Re1_toTbi","Tbi_SF_WB_Md_Re2_toTbi","Tbi_SF_WB_Md_Re3_toTbi","Tbi_SM_WB_Md_Re1_toTbi","Tbi_SM_WB_Md_Re2_toTbi","Tbi_SM_WB_Md_Re3_toTbi")], genes=rawdata_Tbi_longest_iso[,1:1])
y_WB_SB_longest_iso_Tce_UF <- DGEList(counts=rawdata_Tce_longest_iso[,c("Tce_SF_WB_Md_Re1_toTce","Tce_SF_WB_Md_Re2_toTce","Tce_SF_WB_Md_Re3_toTce","Tce_SM_WB_Md_Re1_toTce","Tce_SM_WB_Md_Re2_toTce","Tce_SM_WB_Md_Re3_toTce")], genes=rawdata_Tce_longest_iso[,1:1])
y_WB_SB_longest_iso_Tcm_UF <- DGEList(counts=rawdata_Tcm_longest_iso[,c("Tcm_SF_WB_Md_Re1_toTcm","Tcm_SF_WB_Md_Re2_toTcm","Tcm_SF_WB_Md_Re3_toTcm","Tcm_SM_WB_Md_Re1_toTcm","Tcm_SM_WB_Md_Re2_toTcm","Tcm_SM_WB_Md_Re3_toTcm")], genes=rawdata_Tcm_longest_iso[,1:1])
y_WB_SB_longest_iso_Tpa_UF <- DGEList(counts=rawdata_Tpa_longest_iso[,c("Tpa_SF_WB_Md_Re1_toTpa","Tpa_SF_WB_Md_Re2_toTpa","Tpa_SF_WB_Md_Re3_toTpa","Tpa_SM_WB_Md_Re1_toTpa","Tpa_SM_WB_Md_Re2_toTpa","Tpa_SM_WB_Md_Re3_toTpa")], genes=rawdata_Tpa_longest_iso[,1:1])
y_WB_SB_longest_iso_Tps_UF <- DGEList(counts=rawdata_Tps_longest_iso[,c("Tps_SF_WB_Md_Re1_toTps","Tps_SF_WB_Md_Re2_toTps","Tps_SF_WB_Md_Re3_toTps","Tps_SM_WB_Md_Re1_toTps","Tps_SM_WB_Md_Re2_toTps","Tps_SM_WB_Md_Re3_toTps")], genes=rawdata_Tps_longest_iso[,1:1])

## for sex-asex
y_WB_SA_longest_iso_Tbi_UF <- DGEList(counts=rawdata_Tbi_longest_iso[,c("Tbi_SF_WB_Md_Re1_toTbi","Tbi_SF_WB_Md_Re2_toTbi","Tbi_SF_WB_Md_Re3_toTbi","Tte_AF_WB_Vi_Re1_toTbi","Tte_AF_WB_Vi_Re2_toTbi","Tte_AF_WB_Vi_Re3_toTbi")], genes=rawdata_Tbi_longest_iso[,1:1])
y_WB_SA_longest_iso_Tce_UF <- DGEList(counts=rawdata_Tce_longest_iso[,c("Tce_SF_WB_Md_Re1_toTce","Tce_SF_WB_Md_Re2_toTce","Tce_SF_WB_Md_Re3_toTce","Tms_AF_WB_Vi_Re1_toTce","Tms_AF_WB_Vi_Re2_toTce","Tms_AF_WB_Vi_Re3_toTce")], genes=rawdata_Tce_longest_iso[,1:1])
y_WB_SA_longest_iso_Tcm_UF <- DGEList(counts=rawdata_Tcm_longest_iso[,c("Tcm_SF_WB_Md_Re1_toTcm","Tcm_SF_WB_Md_Re2_toTcm","Tcm_SF_WB_Md_Re3_toTcm","Tsi_AF_WB_Vi_Re1_toTcm","Tsi_AF_WB_Vi_Re2_toTcm","Tsi_AF_WB_Vi_Re3_toTcm")], genes=rawdata_Tcm_longest_iso[,1:1])
y_WB_SA_longest_iso_Tpa_UF <- DGEList(counts=rawdata_Tpa_longest_iso[,c("Tpa_SF_WB_Md_Re1_toTpa","Tpa_SF_WB_Md_Re2_toTpa","Tpa_SF_WB_Md_Re3_toTpa","Tge_AF_WB_Vi_Re1_toTpa","Tge_AF_WB_Vi_Re2_toTpa","Tge_AF_WB_Vi_Re3_toTpa")], genes=rawdata_Tpa_longest_iso[,1:1])
y_WB_SA_longest_iso_Tps_UF <- DGEList(counts=rawdata_Tps_longest_iso[,c("Tps_SF_WB_Md_Re1_toTps","Tps_SF_WB_Md_Re2_toTps","Tps_SF_WB_Md_Re3_toTps","Tdi_AF_WB_Vi_Re1_toTps","Tdi_AF_WB_Vi_Re2_toTps","Tdi_AF_WB_Vi_Re3_toTps")], genes=rawdata_Tps_longest_iso[,1:1])

#### RT

## for SB
y_RT_SB_longest_iso_Tbi_UF <- DGEList(counts=rawdata_Tbi_longest_iso[,c("Tbi_SF_RT_Md_Re1_toTbi","Tbi_SF_RT_Md_Re2_toTbi","Tbi_SF_RT_Md_Re3_toTbi","Tbi_SM_RT_Md_Re1_toTbi","Tbi_SM_RT_Md_Re2_toTbi","Tbi_SM_RT_Md_Re3_toTbi")], genes=rawdata_Tbi_longest_iso[,1:1])
y_RT_SB_longest_iso_Tce_UF <- DGEList(counts=rawdata_Tce_longest_iso[,c("Tce_SF_RT_Md_Re1_toTce","Tce_SF_RT_Md_Re2_toTce","Tce_SF_RT_Md_Re3_toTce","Tce_SM_RT_Md_Re1_toTce","Tce_SM_RT_Md_Re2_toTce","Tce_SM_RT_Md_Re3_toTce")], genes=rawdata_Tce_longest_iso[,1:1])
y_RT_SB_longest_iso_Tcm_UF <- DGEList(counts=rawdata_Tcm_longest_iso[,c("Tcm_SF_RT_Md_Re1_toTcm","Tcm_SF_RT_Md_Re2_toTcm","Tcm_SF_RT_Md_Re3_toTcm","Tcm_SM_RT_Md_Re1_toTcm","Tcm_SM_RT_Md_Re2_toTcm","Tcm_SM_RT_Md_Re3_toTcm")], genes=rawdata_Tcm_longest_iso[,1:1])
y_RT_SB_longest_iso_Tpa_UF <- DGEList(counts=rawdata_Tpa_longest_iso[,c("Tpa_SF_RT_Md_Re1_toTpa","Tpa_SF_RT_Md_Re2_toTpa","Tpa_SF_RT_Md_Re3_toTpa","Tpa_SM_RT_Md_Re1_toTpa","Tpa_SM_RT_Md_Re2_toTpa","Tpa_SM_RT_Md_Re3_toTpa")], genes=rawdata_Tpa_longest_iso[,1:1])
y_RT_SB_longest_iso_Tps_UF <- DGEList(counts=rawdata_Tps_longest_iso[,c("Tps_SF_RT_Md_Re1_toTps","Tps_SF_RT_Md_Re2_toTps","Tps_SF_RT_Md_Re3_toTps","Tps_SM_RT_Md_Re1_toTps","Tps_SM_RT_Md_Re2_toTps","Tps_SM_RT_Md_Re3_toTps")], genes=rawdata_Tps_longest_iso[,1:1])

## for sex-asex
y_RT_SA_longest_iso_Tbi_UF <- DGEList(counts=rawdata_Tbi_longest_iso[,c("Tbi_SF_RT_Md_Re1_toTbi","Tbi_SF_RT_Md_Re2_toTbi","Tbi_SF_RT_Md_Re3_toTbi","Tte_AF_RT_Vi_Re1_toTbi","Tte_AF_RT_Vi_Re2_toTbi","Tte_AF_RT_Vi_Re3_toTbi")], genes=rawdata_Tbi_longest_iso[,1:1])
y_RT_SA_longest_iso_Tce_UF <- DGEList(counts=rawdata_Tce_longest_iso[,c("Tce_SF_RT_Md_Re1_toTce","Tce_SF_RT_Md_Re2_toTce","Tce_SF_RT_Md_Re3_toTce","Tms_AF_RT_Vi_Re1_toTce","Tms_AF_RT_Vi_Re2_toTce","Tms_AF_RT_Vi_Re3_toTce")], genes=rawdata_Tce_longest_iso[,1:1])
y_RT_SA_longest_iso_Tcm_UF <- DGEList(counts=rawdata_Tcm_longest_iso[,c("Tcm_SF_RT_Md_Re1_toTcm","Tcm_SF_RT_Md_Re2_toTcm","Tcm_SF_RT_Md_Re3_toTcm","Tsi_AF_RT_Vi_Re1_toTcm","Tsi_AF_RT_Vi_Re2_toTcm","Tsi_AF_RT_Vi_Re3_toTcm")], genes=rawdata_Tcm_longest_iso[,1:1])
y_RT_SA_longest_iso_Tpa_UF <- DGEList(counts=rawdata_Tpa_longest_iso[,c("Tpa_SF_RT_Md_Re1_toTpa","Tpa_SF_RT_Md_Re2_toTpa","Tpa_SF_RT_Md_Re3_toTpa","Tge_AF_RT_Vi_Re1_toTpa","Tge_AF_RT_Vi_Re2_toTpa","Tge_AF_RT_Vi_Re3_toTpa")], genes=rawdata_Tpa_longest_iso[,1:1])
y_RT_SA_longest_iso_Tps_UF <- DGEList(counts=rawdata_Tps_longest_iso[,c("Tps_SF_RT_Md_Re1_toTps","Tps_SF_RT_Md_Re2_toTps","Tps_SF_RT_Md_Re3_toTps","Tdi_AF_RT_Vi_Re1_toTps","Tdi_AF_RT_Vi_Re2_toTps","Tdi_AF_RT_Vi_Re3_toTps")], genes=rawdata_Tps_longest_iso[,1:1])


#### LG

## for SB
y_LG_SB_longest_iso_Tbi_UF <- DGEList(counts=rawdata_Tbi_longest_iso[,c("Tbi_SF_LG_Md_Re1_toTbi","Tbi_SF_LG_Md_Re2_toTbi","Tbi_SF_LG_Md_Re3_toTbi","Tbi_SM_LG_Md_Re1_toTbi","Tbi_SM_LG_Md_Re2_toTbi","Tbi_SM_LG_Md_Re3_toTbi")], genes=rawdata_Tbi_longest_iso[,1:1])
y_LG_SB_longest_iso_Tce_UF <- DGEList(counts=rawdata_Tce_longest_iso[,c("Tce_SF_LG_Md_Re1_toTce","Tce_SF_LG_Md_Re2_toTce","Tce_SF_LG_Md_Re3_toTce","Tce_SM_LG_Md_Re1_toTce","Tce_SM_LG_Md_Re2_toTce","Tce_SM_LG_Md_Re3_toTce")], genes=rawdata_Tce_longest_iso[,1:1])
y_LG_SB_longest_iso_Tcm_UF <- DGEList(counts=rawdata_Tcm_longest_iso[,c("Tcm_SF_LG_Md_Re1_toTcm","Tcm_SF_LG_Md_Re2_toTcm","Tcm_SF_LG_Md_Re3_toTcm","Tcm_SM_LG_Md_Re1_toTcm","Tcm_SM_LG_Md_Re2_toTcm","Tcm_SM_LG_Md_Re3_toTcm")], genes=rawdata_Tcm_longest_iso[,1:1])
y_LG_SB_longest_iso_Tpa_UF <- DGEList(counts=rawdata_Tpa_longest_iso[,c("Tpa_SF_LG_Md_Re1_toTpa","Tpa_SF_LG_Md_Re2_toTpa","Tpa_SF_LG_Md_Re3_toTpa","Tpa_SM_LG_Md_Re1_toTpa","Tpa_SM_LG_Md_Re2_toTpa","Tpa_SM_LG_Md_Re3_toTpa")], genes=rawdata_Tpa_longest_iso[,1:1])
y_LG_SB_longest_iso_Tps_UF <- DGEList(counts=rawdata_Tps_longest_iso[,c("Tps_SF_LG_Md_Re1_toTps","Tps_SF_LG_Md_Re2_toTps","Tps_SF_LG_Md_Re3_toTps","Tps_SM_LG_Md_Re1_toTps","Tps_SM_LG_Md_Re2_toTps","Tps_SM_LG_Md_Re3_toTps")], genes=rawdata_Tps_longest_iso[,1:1])

## for sex-asex
y_LG_SA_longest_iso_Tbi_UF <- DGEList(counts=rawdata_Tbi_longest_iso[,c("Tbi_SF_LG_Md_Re1_toTbi","Tbi_SF_LG_Md_Re2_toTbi","Tbi_SF_LG_Md_Re3_toTbi","Tte_AF_LG_Vi_Re1_toTbi","Tte_AF_LG_Vi_Re2_toTbi","Tte_AF_LG_Vi_Re3_toTbi")], genes=rawdata_Tbi_longest_iso[,1:1])
y_LG_SA_longest_iso_Tce_UF <- DGEList(counts=rawdata_Tce_longest_iso[,c("Tce_SF_LG_Md_Re1_toTce","Tce_SF_LG_Md_Re2_toTce","Tce_SF_LG_Md_Re3_toTce","Tms_AF_LG_Vi_Re1_toTce","Tms_AF_LG_Vi_Re2_toTce","Tms_AF_LG_Vi_Re3_toTce")], genes=rawdata_Tce_longest_iso[,1:1])
y_LG_SA_longest_iso_Tcm_UF <- DGEList(counts=rawdata_Tcm_longest_iso[,c("Tcm_SF_LG_Md_Re1_toTcm","Tcm_SF_LG_Md_Re2_toTcm","Tcm_SF_LG_Md_Re3_toTcm","Tsi_AF_LG_Vi_Re1_toTcm","Tsi_AF_LG_Vi_Re2_toTcm","Tsi_AF_LG_Vi_Re3_toTcm")], genes=rawdata_Tcm_longest_iso[,1:1])
y_LG_SA_longest_iso_Tpa_UF <- DGEList(counts=rawdata_Tpa_longest_iso[,c("Tpa_SF_LG_Md_Re1_toTpa","Tpa_SF_LG_Md_Re2_toTpa","Tpa_SF_LG_Md_Re3_toTpa","Tge_AF_LG_Vi_Re1_toTpa","Tge_AF_LG_Vi_Re2_toTpa","Tge_AF_LG_Vi_Re3_toTpa")], genes=rawdata_Tpa_longest_iso[,1:1])
y_LG_SA_longest_iso_Tps_UF <- DGEList(counts=rawdata_Tps_longest_iso[,c("Tps_SF_LG_Md_Re1_toTps","Tps_SF_LG_Md_Re2_toTps","Tps_SF_LG_Md_Re3_toTps","Tdi_AF_LG_Vi_Re1_toTps","Tdi_AF_LG_Vi_Re2_toTps","Tdi_AF_LG_Vi_Re3_toTps")], genes=rawdata_Tps_longest_iso[,1:1])






#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### filter and normalisation (code)


filt_and_norm <- function(y,cpm_cut,cut_in_Nsams){
	cat("\nNumber number of genes / samples in orig data\n")
	print(dim(y)) ### number of genes / samples
	head(cpm(y)) 
	keep <- rowSums(cpm(y)>cpm_cut) >= cut_in_Nsams 
	y <- y[keep,]

	y$samples$lib.size <- colSums(y$counts) # if filter recalc lib size
	y <- calcNormFactors(y)
	cat("\nLib norm factors\n")
	print(y$samples)	
	cat("\nNumber number of genes / samples after filtering\n")
	print(dim(y))
	return(y)
}

filt_and_norm_maj <- function(y,cpm_cut,cut_in_Nsams){
	
	cat("\nNumber number of genes / samples in orig data\n")
	print(dim(y)) ### number of genes / samples
	head(cpm(y)) 
	keep <- 
	  rowSums(cpm(y[,1:3])> cpm_cut) >= cut_in_Nsams & 
      rowSums(cpm(y[,4:6])> cpm_cut) >= cut_in_Nsams 
      
	y <- y[keep,]

	y$samples$lib.size <- colSums(y$counts) # if filter recalc lib size
	y <- calcNormFactors(y)
	cat("\nLib norm factors\n")
	print(y$samples)	
	cat("\nNumber number of genes / samples after filtering\n")
	print(dim(y))
	
	return(y)
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
### DE

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### filter and normalisation (samples)
## at least 0.5 cpm  in 2/3 reps in all samps (i.e. NO sex- or sp- specific genes)



#### WB

y_WB_SB_longest_iso_Tbi <- filt_and_norm_maj(y_WB_SB_longest_iso_Tbi_UF,0.5, 2)
y_WB_SB_longest_iso_Tce <- filt_and_norm_maj(y_WB_SB_longest_iso_Tce_UF,0.5, 2)
y_WB_SB_longest_iso_Tcm <- filt_and_norm_maj(y_WB_SB_longest_iso_Tcm_UF,0.5, 2)
y_WB_SB_longest_iso_Tpa <- filt_and_norm_maj(y_WB_SB_longest_iso_Tpa_UF,0.5, 2)
y_WB_SB_longest_iso_Tps <- filt_and_norm_maj(y_WB_SB_longest_iso_Tps_UF,0.5, 2)


y_WB_SA_longest_iso_Tbi <- filt_and_norm_maj(y_WB_SA_longest_iso_Tbi_UF,0.5, 2)
y_WB_SA_longest_iso_Tce <- filt_and_norm_maj(y_WB_SA_longest_iso_Tce_UF,0.5, 2)
y_WB_SA_longest_iso_Tcm <- filt_and_norm_maj(y_WB_SA_longest_iso_Tcm_UF,0.5, 2)
y_WB_SA_longest_iso_Tpa <- filt_and_norm_maj(y_WB_SA_longest_iso_Tpa_UF,0.5, 2)
y_WB_SA_longest_iso_Tps <- filt_and_norm_maj(y_WB_SA_longest_iso_Tps_UF,0.5, 2)


#### RT

y_RT_SB_longest_iso_Tbi <- filt_and_norm_maj(y_RT_SB_longest_iso_Tbi_UF,0.5, 2)
y_RT_SB_longest_iso_Tce <- filt_and_norm_maj(y_RT_SB_longest_iso_Tce_UF,0.5, 2)
y_RT_SB_longest_iso_Tcm <- filt_and_norm_maj(y_RT_SB_longest_iso_Tcm_UF,0.5, 2)
y_RT_SB_longest_iso_Tpa <- filt_and_norm_maj(y_RT_SB_longest_iso_Tpa_UF,0.5, 2)
y_RT_SB_longest_iso_Tps <- filt_and_norm_maj(y_RT_SB_longest_iso_Tps_UF,0.5, 2)

y_RT_SA_longest_iso_Tbi <- filt_and_norm_maj(y_RT_SA_longest_iso_Tbi_UF,0.5, 2)
y_RT_SA_longest_iso_Tce <- filt_and_norm_maj(y_RT_SA_longest_iso_Tce_UF,0.5, 2)
y_RT_SA_longest_iso_Tcm <- filt_and_norm_maj(y_RT_SA_longest_iso_Tcm_UF,0.5, 2)
y_RT_SA_longest_iso_Tpa <- filt_and_norm_maj(y_RT_SA_longest_iso_Tpa_UF,0.5, 2)
y_RT_SA_longest_iso_Tps <- filt_and_norm_maj(y_RT_SA_longest_iso_Tps_UF,0.5, 2)

#### LG

y_LG_SB_longest_iso_Tbi <- filt_and_norm_maj(y_LG_SB_longest_iso_Tbi_UF,0.5, 2)
y_LG_SB_longest_iso_Tce <- filt_and_norm_maj(y_LG_SB_longest_iso_Tce_UF,0.5, 2)
y_LG_SB_longest_iso_Tcm <- filt_and_norm_maj(y_LG_SB_longest_iso_Tcm_UF,0.5, 2)
y_LG_SB_longest_iso_Tpa <- filt_and_norm_maj(y_LG_SB_longest_iso_Tpa_UF,0.5, 2)
y_LG_SB_longest_iso_Tps <- filt_and_norm_maj(y_LG_SB_longest_iso_Tps_UF,0.5, 2)

y_LG_SA_longest_iso_Tbi <- filt_and_norm_maj(y_LG_SA_longest_iso_Tbi_UF,0.5, 2)
y_LG_SA_longest_iso_Tce <- filt_and_norm_maj(y_LG_SA_longest_iso_Tce_UF,0.5, 2)
y_LG_SA_longest_iso_Tcm <- filt_and_norm_maj(y_LG_SA_longest_iso_Tcm_UF,0.5, 2)
y_LG_SA_longest_iso_Tpa <- filt_and_norm_maj(y_LG_SA_longest_iso_Tpa_UF,0.5, 2)
y_LG_SA_longest_iso_Tps <- filt_and_norm_maj(y_LG_SA_longest_iso_Tps_UF,0.5, 2)


#################################################################################################
### get DE_gene_table (code)

## returns full table and sig DE genes as a vector
get_DE_genes <- function(fita,FDRa){
	TT1 = topTags(fita, n =3000000000)
	TT2 = TT1$table
	temp_sig <- subset(TT2, TT2$FDR <= FDRa)	
	sig_genes = temp_sig$genes
	sig_logFC = temp_sig$logFC
	
	N_sig_genes <- length(sig_genes)
	cat("Number of sig genes: ", N_sig_genes )
	
	r_list <- list("table" = TT2, "S_gene_list" = sig_genes, "S_logFC_list" = sig_logFC )
	return(r_list)
}


##################################################################################################
####### get SB_genes

####### design matrix

## using 0M for male and 1F for female samples

sex_SB = factor(c(
"1F","1F","1F","0M","0M","0M"
))

data.frame(Sample=colnames(y_WB_SB_longest_iso_Tbi),sex_SB)


### design matrix WB 
design_WB_SB_longest_iso_Tbi <- model.matrix(~sex_SB) 
rownames(design_WB_SB_longest_iso_Tbi) <- colnames(y_WB_SB_longest_iso_Tbi)
design_WB_SB_longest_iso_Tbi 

### design matrix RT 
design_RT_SB_longest_iso_Tbi <- model.matrix(~sex_SB) 
rownames(design_RT_SB_longest_iso_Tbi ) <- colnames(y_RT_SB_longest_iso_Tbi)
design_RT_SB_longest_iso_Tbi 

### design matrix LG 
design_LG_SB_longest_iso_Tbi  <- model.matrix(~sex_SB ) 
rownames(design_LG_SB_longest_iso_Tbi ) <- colnames(y_LG_SB_longest_iso_Tbi )
design_LG_SB_longest_iso_Tbi 

###### Est dispersion (samples)

y_WB_SB_longest_iso_Tbi  <- estimateDisp(y_WB_SB_longest_iso_Tbi, design_WB_SB_longest_iso_Tbi)
y_RT_SB_longest_iso_Tbi  <- estimateDisp(y_RT_SB_longest_iso_Tbi, design_RT_SB_longest_iso_Tbi)
y_LG_SB_longest_iso_Tbi  <- estimateDisp(y_LG_SB_longest_iso_Tbi, design_LG_SB_longest_iso_Tbi)

### design matrix WB 
design_WB_SB_longest_iso_Tce <- model.matrix(~sex_SB) 
rownames(design_WB_SB_longest_iso_Tce) <- colnames(y_WB_SB_longest_iso_Tce)
design_WB_SB_longest_iso_Tce 

### design matrix RT 
design_RT_SB_longest_iso_Tce <- model.matrix(~sex_SB) 
rownames(design_RT_SB_longest_iso_Tce ) <- colnames(y_RT_SB_longest_iso_Tce)
design_RT_SB_longest_iso_Tce 

### design matrix LG 
design_LG_SB_longest_iso_Tce  <- model.matrix(~sex_SB ) 
rownames(design_LG_SB_longest_iso_Tce ) <- colnames(y_LG_SB_longest_iso_Tce )
design_LG_SB_longest_iso_Tce 

###### Est dispersion (samples)

y_WB_SB_longest_iso_Tce  <- estimateDisp(y_WB_SB_longest_iso_Tce, design_WB_SB_longest_iso_Tce)
y_RT_SB_longest_iso_Tce  <- estimateDisp(y_RT_SB_longest_iso_Tce, design_RT_SB_longest_iso_Tce)
y_LG_SB_longest_iso_Tce  <- estimateDisp(y_LG_SB_longest_iso_Tce, design_LG_SB_longest_iso_Tce)


### design matrix WB 
design_WB_SB_longest_iso_Tcm <- model.matrix(~sex_SB) 
rownames(design_WB_SB_longest_iso_Tcm) <- colnames(y_WB_SB_longest_iso_Tcm)
design_WB_SB_longest_iso_Tcm 

### design matrix RT 
design_RT_SB_longest_iso_Tcm <- model.matrix(~sex_SB) 
rownames(design_RT_SB_longest_iso_Tcm ) <- colnames(y_RT_SB_longest_iso_Tcm)
design_RT_SB_longest_iso_Tcm 


### design matrix LG 
design_LG_SB_longest_iso_Tcm  <- model.matrix(~sex_SB ) 
rownames(design_LG_SB_longest_iso_Tcm ) <- colnames(y_LG_SB_longest_iso_Tcm )
design_LG_SB_longest_iso_Tcm 

###### Est dispersion (samples)

y_WB_SB_longest_iso_Tcm  <- estimateDisp(y_WB_SB_longest_iso_Tcm, design_WB_SB_longest_iso_Tcm)
y_RT_SB_longest_iso_Tcm  <- estimateDisp(y_RT_SB_longest_iso_Tcm, design_RT_SB_longest_iso_Tcm)
y_LG_SB_longest_iso_Tcm  <- estimateDisp(y_LG_SB_longest_iso_Tcm, design_LG_SB_longest_iso_Tcm)


### design matrix WB 
design_WB_SB_longest_iso_Tpa <- model.matrix(~sex_SB) 
rownames(design_WB_SB_longest_iso_Tpa) <- colnames(y_WB_SB_longest_iso_Tpa)
design_WB_SB_longest_iso_Tpa 

### design matrix RT 
design_RT_SB_longest_iso_Tpa <- model.matrix(~sex_SB) 
rownames(design_RT_SB_longest_iso_Tpa ) <- colnames(y_RT_SB_longest_iso_Tpa)
design_RT_SB_longest_iso_Tpa 

### design matrix LG 
design_LG_SB_longest_iso_Tpa  <- model.matrix(~sex_SB ) 
rownames(design_LG_SB_longest_iso_Tpa ) <- colnames(y_LG_SB_longest_iso_Tpa )
design_LG_SB_longest_iso_Tpa 

###### Est dispersion (samples)

y_WB_SB_longest_iso_Tpa  <- estimateDisp(y_WB_SB_longest_iso_Tpa, design_WB_SB_longest_iso_Tpa)
y_RT_SB_longest_iso_Tpa  <- estimateDisp(y_RT_SB_longest_iso_Tpa, design_RT_SB_longest_iso_Tpa)
y_LG_SB_longest_iso_Tpa  <- estimateDisp(y_LG_SB_longest_iso_Tpa, design_LG_SB_longest_iso_Tpa)


### design matrix WB 
design_WB_SB_longest_iso_Tps <- model.matrix(~sex_SB) 
rownames(design_WB_SB_longest_iso_Tps) <- colnames(y_WB_SB_longest_iso_Tps)
design_WB_SB_longest_iso_Tps 

### design matrix RT 
design_RT_SB_longest_iso_Tps <- model.matrix(~sex_SB) 
rownames(design_RT_SB_longest_iso_Tps ) <- colnames(y_RT_SB_longest_iso_Tps)
design_RT_SB_longest_iso_Tps 

### design matrix LG 
design_LG_SB_longest_iso_Tps  <- model.matrix(~sex_SB ) 
rownames(design_LG_SB_longest_iso_Tps ) <- colnames(y_LG_SB_longest_iso_Tps )
design_LG_SB_longest_iso_Tps 

###### Est dispersion (samples)

y_WB_SB_longest_iso_Tps  <- estimateDisp(y_WB_SB_longest_iso_Tps, design_WB_SB_longest_iso_Tps)
y_RT_SB_longest_iso_Tps  <- estimateDisp(y_RT_SB_longest_iso_Tps, design_RT_SB_longest_iso_Tps)
y_LG_SB_longest_iso_Tps  <- estimateDisp(y_LG_SB_longest_iso_Tps, design_LG_SB_longest_iso_Tps)






#####################################################################################################################################################################
###### fit glm 

fit_WB_SB_longest_iso_Tbi <- glmFit(y_WB_SB_longest_iso_Tbi, design_WB_SB_longest_iso_Tbi,robust=TRUE)
fit_RT_SB_longest_iso_Tbi <- glmFit(y_RT_SB_longest_iso_Tbi, design_RT_SB_longest_iso_Tbi,robust=TRUE)
fit_LG_SB_longest_iso_Tbi <- glmFit(y_LG_SB_longest_iso_Tbi, design_LG_SB_longest_iso_Tbi,robust=TRUE)

fit_WB_SB_longest_iso_Tce <- glmFit(y_WB_SB_longest_iso_Tce, design_WB_SB_longest_iso_Tce,robust=TRUE)
fit_RT_SB_longest_iso_Tce <- glmFit(y_RT_SB_longest_iso_Tce, design_RT_SB_longest_iso_Tce,robust=TRUE)
fit_LG_SB_longest_iso_Tce <- glmFit(y_LG_SB_longest_iso_Tce, design_LG_SB_longest_iso_Tce,robust=TRUE)

fit_WB_SB_longest_iso_Tcm <- glmFit(y_WB_SB_longest_iso_Tcm, design_WB_SB_longest_iso_Tcm,robust=TRUE)
fit_RT_SB_longest_iso_Tcm <- glmFit(y_RT_SB_longest_iso_Tcm, design_RT_SB_longest_iso_Tcm,robust=TRUE)
fit_LG_SB_longest_iso_Tcm <- glmFit(y_LG_SB_longest_iso_Tcm, design_LG_SB_longest_iso_Tcm,robust=TRUE)

fit_WB_SB_longest_iso_Tpa <- glmFit(y_WB_SB_longest_iso_Tpa, design_WB_SB_longest_iso_Tpa,robust=TRUE)
fit_RT_SB_longest_iso_Tpa <- glmFit(y_RT_SB_longest_iso_Tpa, design_RT_SB_longest_iso_Tpa,robust=TRUE)
fit_LG_SB_longest_iso_Tpa <- glmFit(y_LG_SB_longest_iso_Tpa, design_LG_SB_longest_iso_Tpa,robust=TRUE)

fit_WB_SB_longest_iso_Tps <- glmFit(y_WB_SB_longest_iso_Tps, design_WB_SB_longest_iso_Tps,robust=TRUE)
fit_RT_SB_longest_iso_Tps <- glmFit(y_RT_SB_longest_iso_Tps, design_RT_SB_longest_iso_Tps,robust=TRUE)
fit_LG_SB_longest_iso_Tps <- glmFit(y_LG_SB_longest_iso_Tps, design_LG_SB_longest_iso_Tps,robust=TRUE)

colnames(fit_LG_SB_longest_iso_Tps)
colnames(fit_WB_SB_longest_iso_Tbi)

# [1] "(Intercept)" "sex_SB1F"   


## get SB genes ### NOTE +ve FC = higher experssion in females (female-biased)

fit_WB_SB_longest_iso_Tbi_c <- glmLRT(fit_WB_SB_longest_iso_Tbi,coef=2)
fit_RT_SB_longest_iso_Tbi_c <- glmLRT(fit_RT_SB_longest_iso_Tbi,coef=2)
fit_LG_SB_longest_iso_Tbi_c <- glmLRT(fit_LG_SB_longest_iso_Tbi,coef=2)

TTT_WB_SB_longest_iso_Tbi_sex_bias <- get_DE_genes(fit_WB_SB_longest_iso_Tbi_c , 0.05)
TTT_RT_SB_longest_iso_Tbi_sex_bias <- get_DE_genes(fit_RT_SB_longest_iso_Tbi_c , 0.05)
TTT_LG_SB_longest_iso_Tbi_sex_bias <- get_DE_genes(fit_LG_SB_longest_iso_Tbi_c , 0.05)

fit_WB_SB_longest_iso_Tce_c <- glmLRT(fit_WB_SB_longest_iso_Tce,coef=2)
fit_RT_SB_longest_iso_Tce_c <- glmLRT(fit_RT_SB_longest_iso_Tce,coef=2)
fit_LG_SB_longest_iso_Tce_c <- glmLRT(fit_LG_SB_longest_iso_Tce,coef=2)

TTT_WB_SB_longest_iso_Tce_sex_bias <- get_DE_genes(fit_WB_SB_longest_iso_Tce_c , 0.05)
TTT_RT_SB_longest_iso_Tce_sex_bias <- get_DE_genes(fit_RT_SB_longest_iso_Tce_c , 0.05)
TTT_LG_SB_longest_iso_Tce_sex_bias <- get_DE_genes(fit_LG_SB_longest_iso_Tce_c , 0.05)

fit_WB_SB_longest_iso_Tcm_c <- glmLRT(fit_WB_SB_longest_iso_Tcm,coef=2)
fit_RT_SB_longest_iso_Tcm_c <- glmLRT(fit_RT_SB_longest_iso_Tcm,coef=2)
fit_LG_SB_longest_iso_Tcm_c <- glmLRT(fit_LG_SB_longest_iso_Tcm,coef=2)

TTT_WB_SB_longest_iso_Tcm_sex_bias <- get_DE_genes(fit_WB_SB_longest_iso_Tcm_c , 0.05)
TTT_RT_SB_longest_iso_Tcm_sex_bias <- get_DE_genes(fit_RT_SB_longest_iso_Tcm_c , 0.05)
TTT_LG_SB_longest_iso_Tcm_sex_bias <- get_DE_genes(fit_LG_SB_longest_iso_Tcm_c , 0.05)


fit_WB_SB_longest_iso_Tpa_c <- glmLRT(fit_WB_SB_longest_iso_Tpa,coef=2)
fit_RT_SB_longest_iso_Tpa_c <- glmLRT(fit_RT_SB_longest_iso_Tpa,coef=2)
fit_LG_SB_longest_iso_Tpa_c <- glmLRT(fit_LG_SB_longest_iso_Tpa,coef=2)

TTT_WB_SB_longest_iso_Tpa_sex_bias <- get_DE_genes(fit_WB_SB_longest_iso_Tpa_c , 0.05)
TTT_RT_SB_longest_iso_Tpa_sex_bias <- get_DE_genes(fit_RT_SB_longest_iso_Tpa_c , 0.05)
TTT_LG_SB_longest_iso_Tpa_sex_bias <- get_DE_genes(fit_LG_SB_longest_iso_Tpa_c , 0.05)


fit_WB_SB_longest_iso_Tps_c <- glmLRT(fit_WB_SB_longest_iso_Tps,coef=2)
fit_RT_SB_longest_iso_Tps_c <- glmLRT(fit_RT_SB_longest_iso_Tps,coef=2)
fit_LG_SB_longest_iso_Tps_c <- glmLRT(fit_LG_SB_longest_iso_Tps,coef=2)

TTT_WB_SB_longest_iso_Tps_sex_bias <- get_DE_genes(fit_WB_SB_longest_iso_Tps_c , 0.05)
TTT_RT_SB_longest_iso_Tps_sex_bias <- get_DE_genes(fit_RT_SB_longest_iso_Tps_c , 0.05)
TTT_LG_SB_longest_iso_Tps_sex_bias <- get_DE_genes(fit_LG_SB_longest_iso_Tps_c , 0.05)


####### export full tables

write.table(TTT_WB_SB_longest_iso_Tbi_sex_bias$table, "TTT_lrt_Tbi_sex_bias_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SB_longest_iso_Tbi_sex_bias$table, "TTT_lrt_Tbi_sex_bias_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SB_longest_iso_Tbi_sex_bias$table, "TTT_lrt_Tbi_sex_bias_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)

write.table(TTT_WB_SB_longest_iso_Tce_sex_bias$table, "TTT_lrt_Tce_sex_bias_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SB_longest_iso_Tce_sex_bias$table, "TTT_lrt_Tce_sex_bias_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SB_longest_iso_Tce_sex_bias$table, "TTT_lrt_Tce_sex_bias_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)

write.table(TTT_WB_SB_longest_iso_Tcm_sex_bias$table, "TTT_lrt_Tcm_sex_bias_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SB_longest_iso_Tcm_sex_bias$table, "TTT_lrt_Tcm_sex_bias_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SB_longest_iso_Tcm_sex_bias$table, "TTT_lrt_Tcm_sex_bias_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)

write.table(TTT_WB_SB_longest_iso_Tpa_sex_bias$table, "TTT_lrt_Tpa_sex_bias_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SB_longest_iso_Tpa_sex_bias$table, "TTT_lrt_Tpa_sex_bias_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SB_longest_iso_Tpa_sex_bias$table, "TTT_lrt_Tpa_sex_bias_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)

write.table(TTT_WB_SB_longest_iso_Tps_sex_bias$table, "TTT_lrt_Tps_sex_bias_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SB_longest_iso_Tps_sex_bias$table, "TTT_lrt_Tps_sex_bias_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SB_longest_iso_Tps_sex_bias$table, "TTT_lrt_Tps_sex_bias_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)




##################################################################################################
####### get SA_genes

####### design matrix

## using 0SF for sexual female and 1AF for asexual female samples

RM_SA = factor(c(
"0SF","0SF","0SF","1AF","1AF","1AF"
))

data.frame(Sample=colnames(y_WB_SA_longest_iso_Tbi),RM_SA)


### design matrix WB 
design_WB_SA_longest_iso_Tbi <- model.matrix(~RM_SA) 
rownames(design_WB_SA_longest_iso_Tbi) <- colnames(y_WB_SA_longest_iso_Tbi)
design_WB_SA_longest_iso_Tbi 

### design matrix RT 
design_RT_SA_longest_iso_Tbi <- model.matrix(~RM_SA) 
rownames(design_RT_SA_longest_iso_Tbi ) <- colnames(y_RT_SA_longest_iso_Tbi)
design_RT_SA_longest_iso_Tbi 

### design matrix LG 
design_LG_SA_longest_iso_Tbi  <- model.matrix(~RM_SA ) 
rownames(design_LG_SA_longest_iso_Tbi ) <- colnames(y_LG_SA_longest_iso_Tbi )
design_LG_SA_longest_iso_Tbi 

###### Est dispersion (samples)

y_WB_SA_longest_iso_Tbi  <- estimateDisp(y_WB_SA_longest_iso_Tbi, design_WB_SA_longest_iso_Tbi)
y_RT_SA_longest_iso_Tbi  <- estimateDisp(y_RT_SA_longest_iso_Tbi, design_RT_SA_longest_iso_Tbi)
y_LG_SA_longest_iso_Tbi  <- estimateDisp(y_LG_SA_longest_iso_Tbi, design_LG_SA_longest_iso_Tbi)

### design matrix WB 
design_WB_SA_longest_iso_Tce <- model.matrix(~RM_SA) 
rownames(design_WB_SA_longest_iso_Tce) <- colnames(y_WB_SA_longest_iso_Tce)
design_WB_SA_longest_iso_Tce 

### design matrix RT 
design_RT_SA_longest_iso_Tce <- model.matrix(~RM_SA) 
rownames(design_RT_SA_longest_iso_Tce ) <- colnames(y_RT_SA_longest_iso_Tce)
design_RT_SA_longest_iso_Tce 


### design matrix LG 
design_LG_SA_longest_iso_Tce  <- model.matrix(~RM_SA ) 
rownames(design_LG_SA_longest_iso_Tce ) <- colnames(y_LG_SA_longest_iso_Tce )
design_LG_SA_longest_iso_Tce 

###### Est dispersion (samples)

y_WB_SA_longest_iso_Tce  <- estimateDisp(y_WB_SA_longest_iso_Tce, design_WB_SA_longest_iso_Tce)
y_RT_SA_longest_iso_Tce  <- estimateDisp(y_RT_SA_longest_iso_Tce, design_RT_SA_longest_iso_Tce)
y_LG_SA_longest_iso_Tce  <- estimateDisp(y_LG_SA_longest_iso_Tce, design_LG_SA_longest_iso_Tce)

### design matrix WB 
design_WB_SA_longest_iso_Tcm <- model.matrix(~RM_SA) 
rownames(design_WB_SA_longest_iso_Tcm) <- colnames(y_WB_SA_longest_iso_Tcm)
design_WB_SA_longest_iso_Tcm 

### design matrix RT 
design_RT_SA_longest_iso_Tcm <- model.matrix(~RM_SA) 
rownames(design_RT_SA_longest_iso_Tcm ) <- colnames(y_RT_SA_longest_iso_Tcm)
design_RT_SA_longest_iso_Tcm 

### design matrix LG 
design_LG_SA_longest_iso_Tcm  <- model.matrix(~RM_SA ) 
rownames(design_LG_SA_longest_iso_Tcm ) <- colnames(y_LG_SA_longest_iso_Tcm )
design_LG_SA_longest_iso_Tcm 

###### Est dispersion (samples)

y_WB_SA_longest_iso_Tcm  <- estimateDisp(y_WB_SA_longest_iso_Tcm, design_WB_SA_longest_iso_Tcm)
y_RT_SA_longest_iso_Tcm  <- estimateDisp(y_RT_SA_longest_iso_Tcm, design_RT_SA_longest_iso_Tcm)
y_LG_SA_longest_iso_Tcm  <- estimateDisp(y_LG_SA_longest_iso_Tcm, design_LG_SA_longest_iso_Tcm)


### design matrix WB 
design_WB_SA_longest_iso_Tpa <- model.matrix(~RM_SA) 
rownames(design_WB_SA_longest_iso_Tpa) <- colnames(y_WB_SA_longest_iso_Tpa)
design_WB_SA_longest_iso_Tpa 

### design matrix RT 
design_RT_SA_longest_iso_Tpa <- model.matrix(~RM_SA) 
rownames(design_RT_SA_longest_iso_Tpa ) <- colnames(y_RT_SA_longest_iso_Tpa)
design_RT_SA_longest_iso_Tpa 

### design matrix LG 
design_LG_SA_longest_iso_Tpa  <- model.matrix(~RM_SA ) 
rownames(design_LG_SA_longest_iso_Tpa ) <- colnames(y_LG_SA_longest_iso_Tpa )
design_LG_SA_longest_iso_Tpa 

###### Est dispersion (samples)

y_WB_SA_longest_iso_Tpa  <- estimateDisp(y_WB_SA_longest_iso_Tpa, design_WB_SA_longest_iso_Tpa)
y_RT_SA_longest_iso_Tpa  <- estimateDisp(y_RT_SA_longest_iso_Tpa, design_RT_SA_longest_iso_Tpa)
y_LG_SA_longest_iso_Tpa  <- estimateDisp(y_LG_SA_longest_iso_Tpa, design_LG_SA_longest_iso_Tpa)


### design matrix WB 
design_WB_SA_longest_iso_Tps <- model.matrix(~RM_SA) 
rownames(design_WB_SA_longest_iso_Tps) <- colnames(y_WB_SA_longest_iso_Tps)
design_WB_SA_longest_iso_Tps 

### design matrix RT 
design_RT_SA_longest_iso_Tps <- model.matrix(~RM_SA) 
rownames(design_RT_SA_longest_iso_Tps ) <- colnames(y_RT_SA_longest_iso_Tps)
design_RT_SA_longest_iso_Tps 

### design matrix LG 
design_LG_SA_longest_iso_Tps  <- model.matrix(~RM_SA ) 
rownames(design_LG_SA_longest_iso_Tps ) <- colnames(y_LG_SA_longest_iso_Tps )
design_LG_SA_longest_iso_Tps 

###### Est dispersion (samples)

y_WB_SA_longest_iso_Tps  <- estimateDisp(y_WB_SA_longest_iso_Tps, design_WB_SA_longest_iso_Tps)
y_RT_SA_longest_iso_Tps  <- estimateDisp(y_RT_SA_longest_iso_Tps, design_RT_SA_longest_iso_Tps)
y_LG_SA_longest_iso_Tps  <- estimateDisp(y_LG_SA_longest_iso_Tps, design_LG_SA_longest_iso_Tps)







#####################################################################################################################################################################
###### fit glm 

fit_WB_SA_longest_iso_Tbi <- glmFit(y_WB_SA_longest_iso_Tbi, design_WB_SA_longest_iso_Tbi,robust=TRUE)
fit_RT_SA_longest_iso_Tbi <- glmFit(y_RT_SA_longest_iso_Tbi, design_RT_SA_longest_iso_Tbi,robust=TRUE)
fit_LG_SA_longest_iso_Tbi <- glmFit(y_LG_SA_longest_iso_Tbi, design_LG_SA_longest_iso_Tbi,robust=TRUE)

fit_WB_SA_longest_iso_Tce <- glmFit(y_WB_SA_longest_iso_Tce, design_WB_SA_longest_iso_Tce,robust=TRUE)
fit_RT_SA_longest_iso_Tce <- glmFit(y_RT_SA_longest_iso_Tce, design_RT_SA_longest_iso_Tce,robust=TRUE)
fit_LG_SA_longest_iso_Tce <- glmFit(y_LG_SA_longest_iso_Tce, design_LG_SA_longest_iso_Tce,robust=TRUE)

fit_WB_SA_longest_iso_Tcm <- glmFit(y_WB_SA_longest_iso_Tcm, design_WB_SA_longest_iso_Tcm,robust=TRUE)
fit_RT_SA_longest_iso_Tcm <- glmFit(y_RT_SA_longest_iso_Tcm, design_RT_SA_longest_iso_Tcm,robust=TRUE)
fit_LG_SA_longest_iso_Tcm <- glmFit(y_LG_SA_longest_iso_Tcm, design_LG_SA_longest_iso_Tcm,robust=TRUE)

fit_WB_SA_longest_iso_Tpa <- glmFit(y_WB_SA_longest_iso_Tpa, design_WB_SA_longest_iso_Tpa,robust=TRUE)
fit_RT_SA_longest_iso_Tpa <- glmFit(y_RT_SA_longest_iso_Tpa, design_RT_SA_longest_iso_Tpa,robust=TRUE)
fit_LG_SA_longest_iso_Tpa <- glmFit(y_LG_SA_longest_iso_Tpa, design_LG_SA_longest_iso_Tpa,robust=TRUE)

fit_WB_SA_longest_iso_Tps <- glmFit(y_WB_SA_longest_iso_Tps, design_WB_SA_longest_iso_Tps,robust=TRUE)
fit_RT_SA_longest_iso_Tps <- glmFit(y_RT_SA_longest_iso_Tps, design_RT_SA_longest_iso_Tps,robust=TRUE)
fit_LG_SA_longest_iso_Tps <- glmFit(y_LG_SA_longest_iso_Tps, design_LG_SA_longest_iso_Tps,robust=TRUE)

colnames(fit_LG_SA_longest_iso_Tps)
colnames(fit_WB_SA_longest_iso_Tbi)

# [1] "(Intercept)" "RM_SA1AF"     


## get Sex-asex genes ### NOTE +ve FC = higher experssion in ASEX

fit_WB_SA_longest_iso_Tbi_c <- glmLRT(fit_WB_SA_longest_iso_Tbi,coef=2)
fit_RT_SA_longest_iso_Tbi_c <- glmLRT(fit_RT_SA_longest_iso_Tbi,coef=2)
fit_LG_SA_longest_iso_Tbi_c <- glmLRT(fit_LG_SA_longest_iso_Tbi,coef=2)

TTT_WB_SA_longest_iso_Tbi_sex_asex <- get_DE_genes(fit_WB_SA_longest_iso_Tbi_c , 0.05)
TTT_RT_SA_longest_iso_Tbi_sex_asex <- get_DE_genes(fit_RT_SA_longest_iso_Tbi_c , 0.05)
TTT_LG_SA_longest_iso_Tbi_sex_asex <- get_DE_genes(fit_LG_SA_longest_iso_Tbi_c , 0.05)

fit_WB_SA_longest_iso_Tce_c <- glmLRT(fit_WB_SA_longest_iso_Tce,coef=2)
fit_RT_SA_longest_iso_Tce_c <- glmLRT(fit_RT_SA_longest_iso_Tce,coef=2)
fit_LG_SA_longest_iso_Tce_c <- glmLRT(fit_LG_SA_longest_iso_Tce,coef=2)

TTT_WB_SA_longest_iso_Tce_sex_asex <- get_DE_genes(fit_WB_SA_longest_iso_Tce_c , 0.05)
TTT_RT_SA_longest_iso_Tce_sex_asex <- get_DE_genes(fit_RT_SA_longest_iso_Tce_c , 0.05)
TTT_LG_SA_longest_iso_Tce_sex_asex <- get_DE_genes(fit_LG_SA_longest_iso_Tce_c , 0.05)

fit_WB_SA_longest_iso_Tcm_c <- glmLRT(fit_WB_SA_longest_iso_Tcm,coef=2)
fit_RT_SA_longest_iso_Tcm_c <- glmLRT(fit_RT_SA_longest_iso_Tcm,coef=2)
fit_LG_SA_longest_iso_Tcm_c <- glmLRT(fit_LG_SA_longest_iso_Tcm,coef=2)

TTT_WB_SA_longest_iso_Tcm_sex_asex <- get_DE_genes(fit_WB_SA_longest_iso_Tcm_c , 0.05)
TTT_RT_SA_longest_iso_Tcm_sex_asex <- get_DE_genes(fit_RT_SA_longest_iso_Tcm_c , 0.05)
TTT_LG_SA_longest_iso_Tcm_sex_asex <- get_DE_genes(fit_LG_SA_longest_iso_Tcm_c , 0.05)

fit_WB_SA_longest_iso_Tpa_c <- glmLRT(fit_WB_SA_longest_iso_Tpa,coef=2)
fit_RT_SA_longest_iso_Tpa_c <- glmLRT(fit_RT_SA_longest_iso_Tpa,coef=2)
fit_LG_SA_longest_iso_Tpa_c <- glmLRT(fit_LG_SA_longest_iso_Tpa,coef=2)

TTT_WB_SA_longest_iso_Tpa_sex_asex <- get_DE_genes(fit_WB_SA_longest_iso_Tpa_c , 0.05)
TTT_RT_SA_longest_iso_Tpa_sex_asex <- get_DE_genes(fit_RT_SA_longest_iso_Tpa_c , 0.05)
TTT_LG_SA_longest_iso_Tpa_sex_asex <- get_DE_genes(fit_LG_SA_longest_iso_Tpa_c , 0.05)

fit_WB_SA_longest_iso_Tps_c <- glmLRT(fit_WB_SA_longest_iso_Tps,coef=2)
fit_RT_SA_longest_iso_Tps_c <- glmLRT(fit_RT_SA_longest_iso_Tps,coef=2)
fit_LG_SA_longest_iso_Tps_c <- glmLRT(fit_LG_SA_longest_iso_Tps,coef=2)

TTT_WB_SA_longest_iso_Tps_sex_asex <- get_DE_genes(fit_WB_SA_longest_iso_Tps_c , 0.05)
TTT_RT_SA_longest_iso_Tps_sex_asex <- get_DE_genes(fit_RT_SA_longest_iso_Tps_c , 0.05)
TTT_LG_SA_longest_iso_Tps_sex_asex <- get_DE_genes(fit_LG_SA_longest_iso_Tps_c , 0.05)

####### export full tables

write.table(TTT_WB_SA_longest_iso_Tbi_sex_asex$table, "TTT_lrt_Tbi_sex_asex_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SA_longest_iso_Tbi_sex_asex$table, "TTT_lrt_Tbi_sex_asex_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SA_longest_iso_Tbi_sex_asex$table, "TTT_lrt_Tbi_sex_asex_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)

write.table(TTT_WB_SA_longest_iso_Tce_sex_asex$table, "TTT_lrt_Tce_sex_asex_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SA_longest_iso_Tce_sex_asex$table, "TTT_lrt_Tce_sex_asex_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SA_longest_iso_Tce_sex_asex$table, "TTT_lrt_Tce_sex_asex_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)

write.table(TTT_WB_SA_longest_iso_Tcm_sex_asex$table, "TTT_lrt_Tcm_sex_asex_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SA_longest_iso_Tcm_sex_asex$table, "TTT_lrt_Tcm_sex_asex_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SA_longest_iso_Tcm_sex_asex$table, "TTT_lrt_Tcm_sex_asex_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)

write.table(TTT_WB_SA_longest_iso_Tpa_sex_asex$table, "TTT_lrt_Tpa_sex_asex_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SA_longest_iso_Tpa_sex_asex$table, "TTT_lrt_Tpa_sex_asex_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SA_longest_iso_Tpa_sex_asex$table, "TTT_lrt_Tpa_sex_asex_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)

write.table(TTT_WB_SA_longest_iso_Tps_sex_asex$table, "TTT_lrt_Tps_sex_asex_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SA_longest_iso_Tps_sex_asex$table, "TTT_lrt_Tps_sex_asex_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SA_longest_iso_Tps_sex_asex$table, "TTT_lrt_Tps_sex_asex_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)








#######################################################################################################################
##### summarise the data a bit

## N sex-biased genes

N_sex_bias_genes <- as.data.frame(c(
length(TTT_WB_SB_longest_iso_Tbi_sex_bias$S_gene_list),
length(TTT_WB_SB_longest_iso_Tce_sex_bias$S_gene_list),
length(TTT_WB_SB_longest_iso_Tcm_sex_bias$S_gene_list),
length(TTT_WB_SB_longest_iso_Tpa_sex_bias$S_gene_list),
length(TTT_WB_SB_longest_iso_Tps_sex_bias$S_gene_list),
length(TTT_RT_SB_longest_iso_Tbi_sex_bias$S_gene_list),
length(TTT_RT_SB_longest_iso_Tce_sex_bias$S_gene_list),
length(TTT_RT_SB_longest_iso_Tcm_sex_bias$S_gene_list),
length(TTT_RT_SB_longest_iso_Tpa_sex_bias$S_gene_list),
length(TTT_RT_SB_longest_iso_Tps_sex_bias$S_gene_list),
length(TTT_LG_SB_longest_iso_Tbi_sex_bias$S_gene_list),
length(TTT_LG_SB_longest_iso_Tce_sex_bias$S_gene_list),
length(TTT_LG_SB_longest_iso_Tcm_sex_bias$S_gene_list),
length(TTT_LG_SB_longest_iso_Tpa_sex_bias$S_gene_list),
length(TTT_LG_SB_longest_iso_Tps_sex_bias$S_gene_list)
))


N_sex_bias_genes <- as.data.frame(cbind(N_sex_bias_genes,c(
length(TTT_WB_SB_longest_iso_Tbi_sex_bias$table[,1]),
length(TTT_WB_SB_longest_iso_Tce_sex_bias$table[,1]),
length(TTT_WB_SB_longest_iso_Tcm_sex_bias$table[,1]),
length(TTT_WB_SB_longest_iso_Tpa_sex_bias$table[,1]),
length(TTT_WB_SB_longest_iso_Tps_sex_bias$table[,1]),
length(TTT_RT_SB_longest_iso_Tbi_sex_bias$table[,1]),
length(TTT_RT_SB_longest_iso_Tce_sex_bias$table[,1]),
length(TTT_RT_SB_longest_iso_Tcm_sex_bias$table[,1]),
length(TTT_RT_SB_longest_iso_Tpa_sex_bias$table[,1]),
length(TTT_RT_SB_longest_iso_Tps_sex_bias$table[,1]),
length(TTT_LG_SB_longest_iso_Tbi_sex_bias$table[,1]),
length(TTT_LG_SB_longest_iso_Tce_sex_bias$table[,1]),
length(TTT_LG_SB_longest_iso_Tcm_sex_bias$table[,1]),
length(TTT_LG_SB_longest_iso_Tpa_sex_bias$table[,1]),
length(TTT_LG_SB_longest_iso_Tps_sex_bias$table[,1])
)))

names(N_sex_bias_genes) <- c("N_sex_bias_genes", "N_genes")

N_sex_bias_genes$sp <- c("Tbi","Tce","Tcm","Tpa","Tps","Tbi","Tce","Tcm","Tpa","Tps","Tbi","Tce","Tcm","Tpa","Tps")
N_sex_bias_genes$tiss <- c("WB","WB","WB","WB","WB","RT","RT","RT","RT","RT","LG","LG","LG","LG","LG")
N_sex_bias_genes$group <- paste(N_sex_bias_genes$sp, N_sex_bias_genes$tiss, sep = "_")

N_sex_bias_genes$group_ordered <-  ordered(N_sex_bias_genes$group, levels=c(
"Tbi_WB","Tbi_RT","Tbi_LG",
"Tce_WB","Tce_RT","Tce_LG",
"Tps_WB","Tps_RT","Tps_LG",
"Tcm_WB","Tcm_RT","Tcm_LG",
"Tpa_WB","Tpa_RT","Tpa_LG"
))


N_sex_bias_genes$perc_sex_bias_genes <- N_sex_bias_genes$N_sex_bias_genes / N_sex_bias_genes$N_genes * 100


write.csv(N_sex_bias_genes, "N_sex_bias_genes_longest_iso.csv")


p3 <- ggplot(N_sex_bias_genes, aes(factor(group_ordered), N_sex_bias_genes, fill = tiss )) + 
	geom_bar(stat="identity", position = "dodge") + 
	theme_bw() +
	xlab ("Species pair") + 
	ylab ("Number of sex biased genes, FDR < 0.05")  +
	scale_y_continuous(expand = c(0,0)) + 
	coord_cartesian(ylim=c(-0,10000)) 
	
p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
   scale_fill_manual(values=c("firebrick3", "#56B4E9", "black"))


png(filename = "N_sex_bias_genes.png",
width = 6, height = 7, units = "in", pointsize = 12,
bg = "white", res = 300)
p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
   scale_fill_manual(values=c("firebrick3", "#56B4E9", "black"))
dev.off()
getwd() ## where has my plot gone....?


p3a <- ggplot(N_sex_bias_genes, aes(factor(group_ordered), perc_sex_bias_genes, fill = tiss )) + 
	geom_bar(stat="identity", position = "dodge") + 
	theme_bw() +
	xlab ("Species pair") + 
	ylab ("Percentage of sex biased genes, FDR < 0.05")  +
	scale_y_continuous(expand = c(0,0)) + 
	coord_cartesian(ylim=c(-0,60)) 
	
p3a + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
   scale_fill_manual(values=c("firebrick3", "#56B4E9", "black"))


png(filename = "Perc_sex_bias_genes.png",
width = 6, height = 7, units = "in", pointsize = 12,
bg = "white", res = 300)
p3a + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
   scale_fill_manual(values=c("firebrick3", "#56B4E9", "black"))
dev.off()
getwd() ## where has my plot gone....?



########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), "sex_bias_edgeR_sexual_ref.R_sessionInfo.txt")





