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

setwd("Output/DE")

###### DE analysis DOING for each tissue seperatly 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

#### WB

## for SB
y_WB_SB_RBBH_Tbi_UF <- DGEList(counts=rawdata_Tbi_RBBH[,c("Tbi_SF_WB_Md_Re1","Tbi_SF_WB_Md_Re2","Tbi_SF_WB_Md_Re3","Tbi_SM_WB_Md_Re1","Tbi_SM_WB_Md_Re2","Tbi_SM_WB_Md_Re3")], genes=rawdata_Tbi_RBBH[,1:1])
y_WB_SB_RBBH_Tce_UF <- DGEList(counts=rawdata_Tce_RBBH[,c("Tce_SF_WB_Md_Re1","Tce_SF_WB_Md_Re2","Tce_SF_WB_Md_Re3","Tce_SM_WB_Md_Re1","Tce_SM_WB_Md_Re2","Tce_SM_WB_Md_Re3")], genes=rawdata_Tce_RBBH[,1:1])
y_WB_SB_RBBH_Tcm_UF <- DGEList(counts=rawdata_Tcm_RBBH[,c("Tcm_SF_WB_Md_Re1","Tcm_SF_WB_Md_Re2","Tcm_SF_WB_Md_Re3","Tcm_SM_WB_Md_Re1","Tcm_SM_WB_Md_Re2","Tcm_SM_WB_Md_Re3")], genes=rawdata_Tcm_RBBH[,1:1])
y_WB_SB_RBBH_Tpa_UF <- DGEList(counts=rawdata_Tpa_RBBH[,c("Tpa_SF_WB_Md_Re1","Tpa_SF_WB_Md_Re2","Tpa_SF_WB_Md_Re3","Tpa_SM_WB_Md_Re1","Tpa_SM_WB_Md_Re2","Tpa_SM_WB_Md_Re3")], genes=rawdata_Tpa_RBBH[,1:1])
y_WB_SB_RBBH_Tps_UF <- DGEList(counts=rawdata_Tps_RBBH[,c("Tps_SF_WB_Md_Re1","Tps_SF_WB_Md_Re2","Tps_SF_WB_Md_Re3","Tps_SM_WB_Md_Re1","Tps_SM_WB_Md_Re2","Tps_SM_WB_Md_Re3")], genes=rawdata_Tps_RBBH[,1:1])

## for sex-asex
y_WB_SA_RBBH_Tbi_UF <- DGEList(counts=rawdata_Tbi_RBBH[,c("Tbi_SF_WB_Md_Re1","Tbi_SF_WB_Md_Re2","Tbi_SF_WB_Md_Re3","Tte_AF_WB_Vi_Re1","Tte_AF_WB_Vi_Re2","Tte_AF_WB_Vi_Re3")], genes=rawdata_Tbi_RBBH[,1:1])
y_WB_SA_RBBH_Tce_UF <- DGEList(counts=rawdata_Tce_RBBH[,c("Tce_SF_WB_Md_Re1","Tce_SF_WB_Md_Re2","Tce_SF_WB_Md_Re3","Tms_AF_WB_Vi_Re1","Tms_AF_WB_Vi_Re2","Tms_AF_WB_Vi_Re3")], genes=rawdata_Tce_RBBH[,1:1])
y_WB_SA_RBBH_Tcm_UF <- DGEList(counts=rawdata_Tcm_RBBH[,c("Tcm_SF_WB_Md_Re1","Tcm_SF_WB_Md_Re2","Tcm_SF_WB_Md_Re3","Tsi_AF_WB_Vi_Re1","Tsi_AF_WB_Vi_Re2","Tsi_AF_WB_Vi_Re3")], genes=rawdata_Tcm_RBBH[,1:1])
y_WB_SA_RBBH_Tpa_UF <- DGEList(counts=rawdata_Tpa_RBBH[,c("Tpa_SF_WB_Md_Re1","Tpa_SF_WB_Md_Re2","Tpa_SF_WB_Md_Re3","Tge_AF_WB_Vi_Re1","Tge_AF_WB_Vi_Re2","Tge_AF_WB_Vi_Re3")], genes=rawdata_Tpa_RBBH[,1:1])
y_WB_SA_RBBH_Tps_UF <- DGEList(counts=rawdata_Tps_RBBH[,c("Tps_SF_WB_Md_Re1","Tps_SF_WB_Md_Re2","Tps_SF_WB_Md_Re3","Tdi_AF_WB_Vi_Re1","Tdi_AF_WB_Vi_Re2","Tdi_AF_WB_Vi_Re3")], genes=rawdata_Tps_RBBH[,1:1])

#### RT

## for SB
y_RT_SB_RBBH_Tbi_UF <- DGEList(counts=rawdata_Tbi_RBBH[,c("Tbi_SF_RT_Md_Re1","Tbi_SF_RT_Md_Re2","Tbi_SF_RT_Md_Re3","Tbi_SM_RT_Md_Re1","Tbi_SM_RT_Md_Re2","Tbi_SM_RT_Md_Re3")], genes=rawdata_Tbi_RBBH[,1:1])
y_RT_SB_RBBH_Tce_UF <- DGEList(counts=rawdata_Tce_RBBH[,c("Tce_SF_RT_Md_Re1","Tce_SF_RT_Md_Re2","Tce_SF_RT_Md_Re3","Tce_SM_RT_Md_Re1","Tce_SM_RT_Md_Re2","Tce_SM_RT_Md_Re3")], genes=rawdata_Tce_RBBH[,1:1])
y_RT_SB_RBBH_Tcm_UF <- DGEList(counts=rawdata_Tcm_RBBH[,c("Tcm_SF_RT_Md_Re1","Tcm_SF_RT_Md_Re2","Tcm_SF_RT_Md_Re3","Tcm_SM_RT_Md_Re1","Tcm_SM_RT_Md_Re2","Tcm_SM_RT_Md_Re3")], genes=rawdata_Tcm_RBBH[,1:1])
y_RT_SB_RBBH_Tpa_UF <- DGEList(counts=rawdata_Tpa_RBBH[,c("Tpa_SF_RT_Md_Re1","Tpa_SF_RT_Md_Re2","Tpa_SF_RT_Md_Re3","Tpa_SM_RT_Md_Re1","Tpa_SM_RT_Md_Re2","Tpa_SM_RT_Md_Re3")], genes=rawdata_Tpa_RBBH[,1:1])
y_RT_SB_RBBH_Tps_UF <- DGEList(counts=rawdata_Tps_RBBH[,c("Tps_SF_RT_Md_Re1","Tps_SF_RT_Md_Re2","Tps_SF_RT_Md_Re3","Tps_SM_RT_Md_Re1","Tps_SM_RT_Md_Re2","Tps_SM_RT_Md_Re3")], genes=rawdata_Tps_RBBH[,1:1])

## for sex-asex
y_RT_SA_RBBH_Tbi_UF <- DGEList(counts=rawdata_Tbi_RBBH[,c("Tbi_SF_RT_Md_Re1","Tbi_SF_RT_Md_Re2","Tbi_SF_RT_Md_Re3","Tte_AF_RT_Vi_Re1","Tte_AF_RT_Vi_Re2","Tte_AF_RT_Vi_Re3")], genes=rawdata_Tbi_RBBH[,1:1])
y_RT_SA_RBBH_Tce_UF <- DGEList(counts=rawdata_Tce_RBBH[,c("Tce_SF_RT_Md_Re1","Tce_SF_RT_Md_Re2","Tce_SF_RT_Md_Re3","Tms_AF_RT_Vi_Re1","Tms_AF_RT_Vi_Re2","Tms_AF_RT_Vi_Re3")], genes=rawdata_Tce_RBBH[,1:1])
y_RT_SA_RBBH_Tcm_UF <- DGEList(counts=rawdata_Tcm_RBBH[,c("Tcm_SF_RT_Md_Re1","Tcm_SF_RT_Md_Re2","Tcm_SF_RT_Md_Re3","Tsi_AF_RT_Vi_Re1","Tsi_AF_RT_Vi_Re2","Tsi_AF_RT_Vi_Re3")], genes=rawdata_Tcm_RBBH[,1:1])
y_RT_SA_RBBH_Tpa_UF <- DGEList(counts=rawdata_Tpa_RBBH[,c("Tpa_SF_RT_Md_Re1","Tpa_SF_RT_Md_Re2","Tpa_SF_RT_Md_Re3","Tge_AF_RT_Vi_Re1","Tge_AF_RT_Vi_Re2","Tge_AF_RT_Vi_Re3")], genes=rawdata_Tpa_RBBH[,1:1])
y_RT_SA_RBBH_Tps_UF <- DGEList(counts=rawdata_Tps_RBBH[,c("Tps_SF_RT_Md_Re1","Tps_SF_RT_Md_Re2","Tps_SF_RT_Md_Re3","Tdi_AF_RT_Vi_Re1","Tdi_AF_RT_Vi_Re2","Tdi_AF_RT_Vi_Re3")], genes=rawdata_Tps_RBBH[,1:1])


#### LG

## for SB
y_LG_SB_RBBH_Tbi_UF <- DGEList(counts=rawdata_Tbi_RBBH[,c("Tbi_SF_LG_Md_Re1","Tbi_SF_LG_Md_Re2","Tbi_SF_LG_Md_Re3","Tbi_SM_LG_Md_Re1","Tbi_SM_LG_Md_Re2","Tbi_SM_LG_Md_Re3")], genes=rawdata_Tbi_RBBH[,1:1])
y_LG_SB_RBBH_Tce_UF <- DGEList(counts=rawdata_Tce_RBBH[,c("Tce_SF_LG_Md_Re1","Tce_SF_LG_Md_Re2","Tce_SF_LG_Md_Re3","Tce_SM_LG_Md_Re1","Tce_SM_LG_Md_Re2","Tce_SM_LG_Md_Re3")], genes=rawdata_Tce_RBBH[,1:1])
y_LG_SB_RBBH_Tcm_UF <- DGEList(counts=rawdata_Tcm_RBBH[,c("Tcm_SF_LG_Md_Re1","Tcm_SF_LG_Md_Re2","Tcm_SF_LG_Md_Re3","Tcm_SM_LG_Md_Re1","Tcm_SM_LG_Md_Re2","Tcm_SM_LG_Md_Re3")], genes=rawdata_Tcm_RBBH[,1:1])
y_LG_SB_RBBH_Tpa_UF <- DGEList(counts=rawdata_Tpa_RBBH[,c("Tpa_SF_LG_Md_Re1","Tpa_SF_LG_Md_Re2","Tpa_SF_LG_Md_Re3","Tpa_SM_LG_Md_Re1","Tpa_SM_LG_Md_Re2","Tpa_SM_LG_Md_Re3")], genes=rawdata_Tpa_RBBH[,1:1])
y_LG_SB_RBBH_Tps_UF <- DGEList(counts=rawdata_Tps_RBBH[,c("Tps_SF_LG_Md_Re1","Tps_SF_LG_Md_Re2","Tps_SF_LG_Md_Re3","Tps_SM_LG_Md_Re1","Tps_SM_LG_Md_Re2","Tps_SM_LG_Md_Re3")], genes=rawdata_Tps_RBBH[,1:1])

## for sex-asex
y_LG_SA_RBBH_Tbi_UF <- DGEList(counts=rawdata_Tbi_RBBH[,c("Tbi_SF_LG_Md_Re1","Tbi_SF_LG_Md_Re2","Tbi_SF_LG_Md_Re3","Tte_AF_LG_Vi_Re1","Tte_AF_LG_Vi_Re2","Tte_AF_LG_Vi_Re3")], genes=rawdata_Tbi_RBBH[,1:1])
y_LG_SA_RBBH_Tce_UF <- DGEList(counts=rawdata_Tce_RBBH[,c("Tce_SF_LG_Md_Re1","Tce_SF_LG_Md_Re2","Tce_SF_LG_Md_Re3","Tms_AF_LG_Vi_Re1","Tms_AF_LG_Vi_Re2","Tms_AF_LG_Vi_Re3")], genes=rawdata_Tce_RBBH[,1:1])
y_LG_SA_RBBH_Tcm_UF <- DGEList(counts=rawdata_Tcm_RBBH[,c("Tcm_SF_LG_Md_Re1","Tcm_SF_LG_Md_Re2","Tcm_SF_LG_Md_Re3","Tsi_AF_LG_Vi_Re1","Tsi_AF_LG_Vi_Re2","Tsi_AF_LG_Vi_Re3")], genes=rawdata_Tcm_RBBH[,1:1])
y_LG_SA_RBBH_Tpa_UF <- DGEList(counts=rawdata_Tpa_RBBH[,c("Tpa_SF_LG_Md_Re1","Tpa_SF_LG_Md_Re2","Tpa_SF_LG_Md_Re3","Tge_AF_LG_Vi_Re1","Tge_AF_LG_Vi_Re2","Tge_AF_LG_Vi_Re3")], genes=rawdata_Tpa_RBBH[,1:1])
y_LG_SA_RBBH_Tps_UF <- DGEList(counts=rawdata_Tps_RBBH[,c("Tps_SF_LG_Md_Re1","Tps_SF_LG_Md_Re2","Tps_SF_LG_Md_Re3","Tdi_AF_LG_Vi_Re1","Tdi_AF_LG_Vi_Re2","Tdi_AF_LG_Vi_Re3")], genes=rawdata_Tps_RBBH[,1:1])


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

y_WB_SB_RBBH_Tbi <- filt_and_norm_maj(y_WB_SB_RBBH_Tbi_UF,0.5, 2)
y_WB_SB_RBBH_Tce <- filt_and_norm_maj(y_WB_SB_RBBH_Tce_UF,0.5, 2)
y_WB_SB_RBBH_Tcm <- filt_and_norm_maj(y_WB_SB_RBBH_Tcm_UF,0.5, 2)
y_WB_SB_RBBH_Tpa <- filt_and_norm_maj(y_WB_SB_RBBH_Tpa_UF,0.5, 2)
y_WB_SB_RBBH_Tps <- filt_and_norm_maj(y_WB_SB_RBBH_Tps_UF,0.5, 2)


# Number number of genes / samples in orig data
# [1] 17583     6
# [1] 17170     6
# [1] 15500     6
# [1] 16113     6
# [1] 15932     6

# Number number of genes / samples after filtering
# [1] 15401     6
# [1] 15082     6
# [1] 14154     6
# [1] 14111     6
# [1] 13994     6

y_WB_SA_RBBH_Tbi <- filt_and_norm_maj(y_WB_SA_RBBH_Tbi_UF,0.5, 2)
y_WB_SA_RBBH_Tce <- filt_and_norm_maj(y_WB_SA_RBBH_Tce_UF,0.5, 2)
y_WB_SA_RBBH_Tcm <- filt_and_norm_maj(y_WB_SA_RBBH_Tcm_UF,0.5, 2)
y_WB_SA_RBBH_Tpa <- filt_and_norm_maj(y_WB_SA_RBBH_Tpa_UF,0.5, 2)
y_WB_SA_RBBH_Tps <- filt_and_norm_maj(y_WB_SA_RBBH_Tps_UF,0.5, 2)


# Number number of genes / samples in orig data
# [1] 17583     6
# [1] 17170     6
# [1] 15500     6
# [1] 16113     6
# [1] 15932     

# Number number of genes / samples after filtering
# [1] 15793     6
# [1] 15308     6
# [1] 14354     6
# [1] 14205     6
# [1] 14257     6

#### RT

y_RT_SB_RBBH_Tbi <- filt_and_norm_maj(y_RT_SB_RBBH_Tbi_UF,0.5, 2)
y_RT_SB_RBBH_Tce <- filt_and_norm_maj(y_RT_SB_RBBH_Tce_UF,0.5, 2)
y_RT_SB_RBBH_Tcm <- filt_and_norm_maj(y_RT_SB_RBBH_Tcm_UF,0.5, 2)
y_RT_SB_RBBH_Tpa <- filt_and_norm_maj(y_RT_SB_RBBH_Tpa_UF,0.5, 2)
y_RT_SB_RBBH_Tps <- filt_and_norm_maj(y_RT_SB_RBBH_Tps_UF,0.5, 2)

# Number number of genes / samples in orig data
# [1] 17583     6
# [1] 17170     6
# [1] 15500     6
# [1] 16113     6
# [1] 15932     6

# Number number of genes / samples after filtering
# [1] 13068     6
# [1] 12484     6
# [1] 11656     6
# [1] 11404     6
# [1] 11977     6

y_RT_SA_RBBH_Tbi <- filt_and_norm_maj(y_RT_SA_RBBH_Tbi_UF,0.5, 2)
y_RT_SA_RBBH_Tce <- filt_and_norm_maj(y_RT_SA_RBBH_Tce_UF,0.5, 2)
y_RT_SA_RBBH_Tcm <- filt_and_norm_maj(y_RT_SA_RBBH_Tcm_UF,0.5, 2)
y_RT_SA_RBBH_Tpa <- filt_and_norm_maj(y_RT_SA_RBBH_Tpa_UF,0.5, 2)
y_RT_SA_RBBH_Tps <- filt_and_norm_maj(y_RT_SA_RBBH_Tps_UF,0.5, 2)

# # Number number of genes / samples in orig data
# [1] 17583     6
# [1] 17170     6
# [1] 15500     6
# [1] 16113     6
# [1] 15932     6

# Number number of genes / samples after filtering
# [1] 13196     6
# [1] 12599     6
# [1] 11654     6
# [1] 11322     6
# [1] 11990     6


#### LG

y_LG_SB_RBBH_Tbi <- filt_and_norm_maj(y_LG_SB_RBBH_Tbi_UF,0.5, 2)
y_LG_SB_RBBH_Tce <- filt_and_norm_maj(y_LG_SB_RBBH_Tce_UF,0.5, 2)
y_LG_SB_RBBH_Tcm <- filt_and_norm_maj(y_LG_SB_RBBH_Tcm_UF,0.5, 2)
y_LG_SB_RBBH_Tpa <- filt_and_norm_maj(y_LG_SB_RBBH_Tpa_UF,0.5, 2)
y_LG_SB_RBBH_Tps <- filt_and_norm_maj(y_LG_SB_RBBH_Tps_UF,0.5, 2)


# Number number of genes / samples in orig data
# [1] 17583     6
# [1] 17170     6
# [1] 15500     6
# [1] 16113     6
# [1] 15932     6

# Number number of genes / samples after filtering
# [1] 12335     6
# [1] 11120     6
# [1] 11349     6
# [1] 10565     6
# [1] 11483     6


y_LG_SA_RBBH_Tbi <- filt_and_norm_maj(y_LG_SA_RBBH_Tbi_UF,0.5, 2)
y_LG_SA_RBBH_Tce <- filt_and_norm_maj(y_LG_SA_RBBH_Tce_UF,0.5, 2)
y_LG_SA_RBBH_Tcm <- filt_and_norm_maj(y_LG_SA_RBBH_Tcm_UF,0.5, 2)
y_LG_SA_RBBH_Tpa <- filt_and_norm_maj(y_LG_SA_RBBH_Tpa_UF,0.5, 2)
y_LG_SA_RBBH_Tps <- filt_and_norm_maj(y_LG_SA_RBBH_Tps_UF,0.5, 2)


# Number number of genes / samples in orig data
# [1] 17583     6
# [1] 17170     6
# [1] 15500     6
# [1] 16113     6
# [1] 15932     6

# Number number of genes / samples after filtering
# [1] 12302     6
# [1] 11068     6
# [1] 11223     6
# [1] 10057     6
# [1] 11313     6


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

data.frame(Sample=colnames(y_WB_SB_RBBH_Tbi),sex_SB)


### design matrix WB 
design_WB_SB_RBBH_Tbi <- model.matrix(~sex_SB) 
rownames(design_WB_SB_RBBH_Tbi) <- colnames(y_WB_SB_RBBH_Tbi)
design_WB_SB_RBBH_Tbi 

### design matrix RT 
design_RT_SB_RBBH_Tbi <- model.matrix(~sex_SB) 
rownames(design_RT_SB_RBBH_Tbi ) <- colnames(y_RT_SB_RBBH_Tbi)
design_RT_SB_RBBH_Tbi 

### design matrix LG 
design_LG_SB_RBBH_Tbi  <- model.matrix(~sex_SB ) 
rownames(design_LG_SB_RBBH_Tbi ) <- colnames(y_LG_SB_RBBH_Tbi )
design_LG_SB_RBBH_Tbi 

###### Est dispersion (samples)

y_WB_SB_RBBH_Tbi  <- estimateDisp(y_WB_SB_RBBH_Tbi, design_WB_SB_RBBH_Tbi)
y_RT_SB_RBBH_Tbi  <- estimateDisp(y_RT_SB_RBBH_Tbi, design_RT_SB_RBBH_Tbi)
y_LG_SB_RBBH_Tbi  <- estimateDisp(y_LG_SB_RBBH_Tbi, design_LG_SB_RBBH_Tbi)

### design matrix WB 
design_WB_SB_RBBH_Tce <- model.matrix(~sex_SB) 
rownames(design_WB_SB_RBBH_Tce) <- colnames(y_WB_SB_RBBH_Tce)
design_WB_SB_RBBH_Tce 

### design matrix RT 
design_RT_SB_RBBH_Tce <- model.matrix(~sex_SB) 
rownames(design_RT_SB_RBBH_Tce ) <- colnames(y_RT_SB_RBBH_Tce)
design_RT_SB_RBBH_Tce 

### design matrix LG 
design_LG_SB_RBBH_Tce  <- model.matrix(~sex_SB ) 
rownames(design_LG_SB_RBBH_Tce ) <- colnames(y_LG_SB_RBBH_Tce )
design_LG_SB_RBBH_Tce 

###### Est dispersion (samples)

y_WB_SB_RBBH_Tce  <- estimateDisp(y_WB_SB_RBBH_Tce, design_WB_SB_RBBH_Tce)
y_RT_SB_RBBH_Tce  <- estimateDisp(y_RT_SB_RBBH_Tce, design_RT_SB_RBBH_Tce)
y_LG_SB_RBBH_Tce  <- estimateDisp(y_LG_SB_RBBH_Tce, design_LG_SB_RBBH_Tce)


### design matrix WB 
design_WB_SB_RBBH_Tcm <- model.matrix(~sex_SB) 
rownames(design_WB_SB_RBBH_Tcm) <- colnames(y_WB_SB_RBBH_Tcm)
design_WB_SB_RBBH_Tcm 

### design matrix RT 
design_RT_SB_RBBH_Tcm <- model.matrix(~sex_SB) 
rownames(design_RT_SB_RBBH_Tcm ) <- colnames(y_RT_SB_RBBH_Tcm)
design_RT_SB_RBBH_Tcm 

### design matrix LG 
design_LG_SB_RBBH_Tcm  <- model.matrix(~sex_SB ) 
rownames(design_LG_SB_RBBH_Tcm ) <- colnames(y_LG_SB_RBBH_Tcm )
design_LG_SB_RBBH_Tcm 

###### Est dispersion (samples)

y_WB_SB_RBBH_Tcm  <- estimateDisp(y_WB_SB_RBBH_Tcm, design_WB_SB_RBBH_Tcm)
y_RT_SB_RBBH_Tcm  <- estimateDisp(y_RT_SB_RBBH_Tcm, design_RT_SB_RBBH_Tcm)
y_LG_SB_RBBH_Tcm  <- estimateDisp(y_LG_SB_RBBH_Tcm, design_LG_SB_RBBH_Tcm)


### design matrix WB 
design_WB_SB_RBBH_Tpa <- model.matrix(~sex_SB) 
rownames(design_WB_SB_RBBH_Tpa) <- colnames(y_WB_SB_RBBH_Tpa)
design_WB_SB_RBBH_Tpa 

### design matrix RT 
design_RT_SB_RBBH_Tpa <- model.matrix(~sex_SB) 
rownames(design_RT_SB_RBBH_Tpa ) <- colnames(y_RT_SB_RBBH_Tpa)
design_RT_SB_RBBH_Tpa 

### design matrix LG 
design_LG_SB_RBBH_Tpa  <- model.matrix(~sex_SB ) 
rownames(design_LG_SB_RBBH_Tpa ) <- colnames(y_LG_SB_RBBH_Tpa )
design_LG_SB_RBBH_Tpa 

###### Est dispersion (samples)

y_WB_SB_RBBH_Tpa  <- estimateDisp(y_WB_SB_RBBH_Tpa, design_WB_SB_RBBH_Tpa)
y_RT_SB_RBBH_Tpa  <- estimateDisp(y_RT_SB_RBBH_Tpa, design_RT_SB_RBBH_Tpa)
y_LG_SB_RBBH_Tpa  <- estimateDisp(y_LG_SB_RBBH_Tpa, design_LG_SB_RBBH_Tpa)


### design matrix WB 
design_WB_SB_RBBH_Tps <- model.matrix(~sex_SB) 
rownames(design_WB_SB_RBBH_Tps) <- colnames(y_WB_SB_RBBH_Tps)
design_WB_SB_RBBH_Tps 

### design matrix RT 
design_RT_SB_RBBH_Tps <- model.matrix(~sex_SB) 
rownames(design_RT_SB_RBBH_Tps ) <- colnames(y_RT_SB_RBBH_Tps)
design_RT_SB_RBBH_Tps 

### design matrix LG 
design_LG_SB_RBBH_Tps  <- model.matrix(~sex_SB ) 
rownames(design_LG_SB_RBBH_Tps ) <- colnames(y_LG_SB_RBBH_Tps )
design_LG_SB_RBBH_Tps 

###### Est dispersion (samples)

y_WB_SB_RBBH_Tps  <- estimateDisp(y_WB_SB_RBBH_Tps, design_WB_SB_RBBH_Tps)
y_RT_SB_RBBH_Tps  <- estimateDisp(y_RT_SB_RBBH_Tps, design_RT_SB_RBBH_Tps)
y_LG_SB_RBBH_Tps  <- estimateDisp(y_LG_SB_RBBH_Tps, design_LG_SB_RBBH_Tps)



#####################################################################################################################################################################
###### GET BCV ## this is the (common) BCV


BCV_SB_df <- as.data.frame(cbind(
c("Tbi","Tce", "Tcm", "Tpa", "Tps"),
c(
sqrt(y_WB_SB_RBBH_Tbi$common.dispersion),
sqrt(y_WB_SB_RBBH_Tce$common.dispersion),
sqrt(y_WB_SB_RBBH_Tcm$common.dispersion),
sqrt(y_WB_SB_RBBH_Tpa$common.dispersion),
sqrt(y_WB_SB_RBBH_Tps$common.dispersion)),
c(
sqrt(y_RT_SB_RBBH_Tbi$common.dispersion),
sqrt(y_RT_SB_RBBH_Tce$common.dispersion),
sqrt(y_RT_SB_RBBH_Tcm$common.dispersion),
sqrt(y_RT_SB_RBBH_Tpa$common.dispersion),
sqrt(y_RT_SB_RBBH_Tps$common.dispersion)),

c(
sqrt(y_LG_SB_RBBH_Tbi$common.dispersion),
sqrt(y_LG_SB_RBBH_Tce$common.dispersion),
sqrt(y_LG_SB_RBBH_Tcm$common.dispersion),
sqrt(y_LG_SB_RBBH_Tpa$common.dispersion),
sqrt(y_LG_SB_RBBH_Tps$common.dispersion))))

colnames(BCV_SB_df) <- c("sp","WB","RT","LG")

write.csv(BCV_SB_df, "BCV_SB_df.csv", row.names = F)


#####################################################################################################################################################################
###### fit glm 

fit_WB_SB_RBBH_Tbi <- glmFit(y_WB_SB_RBBH_Tbi, design_WB_SB_RBBH_Tbi,robust=TRUE)
fit_RT_SB_RBBH_Tbi <- glmFit(y_RT_SB_RBBH_Tbi, design_RT_SB_RBBH_Tbi,robust=TRUE)
fit_LG_SB_RBBH_Tbi <- glmFit(y_LG_SB_RBBH_Tbi, design_LG_SB_RBBH_Tbi,robust=TRUE)

fit_WB_SB_RBBH_Tce <- glmFit(y_WB_SB_RBBH_Tce, design_WB_SB_RBBH_Tce,robust=TRUE)
fit_RT_SB_RBBH_Tce <- glmFit(y_RT_SB_RBBH_Tce, design_RT_SB_RBBH_Tce,robust=TRUE)
fit_LG_SB_RBBH_Tce <- glmFit(y_LG_SB_RBBH_Tce, design_LG_SB_RBBH_Tce,robust=TRUE)

fit_WB_SB_RBBH_Tcm <- glmFit(y_WB_SB_RBBH_Tcm, design_WB_SB_RBBH_Tcm,robust=TRUE)
fit_RT_SB_RBBH_Tcm <- glmFit(y_RT_SB_RBBH_Tcm, design_RT_SB_RBBH_Tcm,robust=TRUE)
fit_LG_SB_RBBH_Tcm <- glmFit(y_LG_SB_RBBH_Tcm, design_LG_SB_RBBH_Tcm,robust=TRUE)

fit_WB_SB_RBBH_Tpa <- glmFit(y_WB_SB_RBBH_Tpa, design_WB_SB_RBBH_Tpa,robust=TRUE)
fit_RT_SB_RBBH_Tpa <- glmFit(y_RT_SB_RBBH_Tpa, design_RT_SB_RBBH_Tpa,robust=TRUE)
fit_LG_SB_RBBH_Tpa <- glmFit(y_LG_SB_RBBH_Tpa, design_LG_SB_RBBH_Tpa,robust=TRUE)

fit_WB_SB_RBBH_Tps <- glmFit(y_WB_SB_RBBH_Tps, design_WB_SB_RBBH_Tps,robust=TRUE)
fit_RT_SB_RBBH_Tps <- glmFit(y_RT_SB_RBBH_Tps, design_RT_SB_RBBH_Tps,robust=TRUE)
fit_LG_SB_RBBH_Tps <- glmFit(y_LG_SB_RBBH_Tps, design_LG_SB_RBBH_Tps,robust=TRUE)

colnames(fit_LG_SB_RBBH_Tps)
colnames(fit_WB_SB_RBBH_Tbi)

# [1] "(Intercept)" "sex_SB1F"   


## get SB genes ### NOTE +ve FC = higher experssion in females (female-biased)

fit_WB_SB_RBBH_Tbi_c <- glmLRT(fit_WB_SB_RBBH_Tbi,coef=2)
fit_RT_SB_RBBH_Tbi_c <- glmLRT(fit_RT_SB_RBBH_Tbi,coef=2)
fit_LG_SB_RBBH_Tbi_c <- glmLRT(fit_LG_SB_RBBH_Tbi,coef=2)

TTT_WB_SB_RBBH_Tbi_sex_bias <- get_DE_genes(fit_WB_SB_RBBH_Tbi_c , 0.05)
TTT_RT_SB_RBBH_Tbi_sex_bias <- get_DE_genes(fit_RT_SB_RBBH_Tbi_c , 0.05)
TTT_LG_SB_RBBH_Tbi_sex_bias <- get_DE_genes(fit_LG_SB_RBBH_Tbi_c , 0.05)

fit_WB_SB_RBBH_Tce_c <- glmLRT(fit_WB_SB_RBBH_Tce,coef=2)
fit_RT_SB_RBBH_Tce_c <- glmLRT(fit_RT_SB_RBBH_Tce,coef=2)
fit_LG_SB_RBBH_Tce_c <- glmLRT(fit_LG_SB_RBBH_Tce,coef=2)

TTT_WB_SB_RBBH_Tce_sex_bias <- get_DE_genes(fit_WB_SB_RBBH_Tce_c , 0.05)
TTT_RT_SB_RBBH_Tce_sex_bias <- get_DE_genes(fit_RT_SB_RBBH_Tce_c , 0.05)
TTT_LG_SB_RBBH_Tce_sex_bias <- get_DE_genes(fit_LG_SB_RBBH_Tce_c , 0.05)

fit_WB_SB_RBBH_Tcm_c <- glmLRT(fit_WB_SB_RBBH_Tcm,coef=2)
fit_RT_SB_RBBH_Tcm_c <- glmLRT(fit_RT_SB_RBBH_Tcm,coef=2)
fit_LG_SB_RBBH_Tcm_c <- glmLRT(fit_LG_SB_RBBH_Tcm,coef=2)

TTT_WB_SB_RBBH_Tcm_sex_bias <- get_DE_genes(fit_WB_SB_RBBH_Tcm_c , 0.05)
TTT_RT_SB_RBBH_Tcm_sex_bias <- get_DE_genes(fit_RT_SB_RBBH_Tcm_c , 0.05)
TTT_LG_SB_RBBH_Tcm_sex_bias <- get_DE_genes(fit_LG_SB_RBBH_Tcm_c , 0.05)

fit_WB_SB_RBBH_Tpa_c <- glmLRT(fit_WB_SB_RBBH_Tpa,coef=2)
fit_RT_SB_RBBH_Tpa_c <- glmLRT(fit_RT_SB_RBBH_Tpa,coef=2)
fit_LG_SB_RBBH_Tpa_c <- glmLRT(fit_LG_SB_RBBH_Tpa,coef=2)

TTT_WB_SB_RBBH_Tpa_sex_bias <- get_DE_genes(fit_WB_SB_RBBH_Tpa_c , 0.05)
TTT_RT_SB_RBBH_Tpa_sex_bias <- get_DE_genes(fit_RT_SB_RBBH_Tpa_c , 0.05)
TTT_LG_SB_RBBH_Tpa_sex_bias <- get_DE_genes(fit_LG_SB_RBBH_Tpa_c , 0.05)

fit_WB_SB_RBBH_Tps_c <- glmLRT(fit_WB_SB_RBBH_Tps,coef=2)
fit_RT_SB_RBBH_Tps_c <- glmLRT(fit_RT_SB_RBBH_Tps,coef=2)
fit_LG_SB_RBBH_Tps_c <- glmLRT(fit_LG_SB_RBBH_Tps,coef=2)

TTT_WB_SB_RBBH_Tps_sex_bias <- get_DE_genes(fit_WB_SB_RBBH_Tps_c , 0.05)
TTT_RT_SB_RBBH_Tps_sex_bias <- get_DE_genes(fit_RT_SB_RBBH_Tps_c , 0.05)
TTT_LG_SB_RBBH_Tps_sex_bias <- get_DE_genes(fit_LG_SB_RBBH_Tps_c , 0.05)


####### export full tables

write.table(TTT_WB_SB_RBBH_Tbi_sex_bias$table, "TTT_lrt_Tbi_sex_bias_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SB_RBBH_Tbi_sex_bias$table, "TTT_lrt_Tbi_sex_bias_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SB_RBBH_Tbi_sex_bias$table, "TTT_lrt_Tbi_sex_bias_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)

write.table(TTT_WB_SB_RBBH_Tce_sex_bias$table, "TTT_lrt_Tce_sex_bias_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SB_RBBH_Tce_sex_bias$table, "TTT_lrt_Tce_sex_bias_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SB_RBBH_Tce_sex_bias$table, "TTT_lrt_Tce_sex_bias_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)

write.table(TTT_WB_SB_RBBH_Tcm_sex_bias$table, "TTT_lrt_Tcm_sex_bias_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SB_RBBH_Tcm_sex_bias$table, "TTT_lrt_Tcm_sex_bias_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SB_RBBH_Tcm_sex_bias$table, "TTT_lrt_Tcm_sex_bias_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)

write.table(TTT_WB_SB_RBBH_Tpa_sex_bias$table, "TTT_lrt_Tpa_sex_bias_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SB_RBBH_Tpa_sex_bias$table, "TTT_lrt_Tpa_sex_bias_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SB_RBBH_Tpa_sex_bias$table, "TTT_lrt_Tpa_sex_bias_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)

write.table(TTT_WB_SB_RBBH_Tps_sex_bias$table, "TTT_lrt_Tps_sex_bias_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SB_RBBH_Tps_sex_bias$table, "TTT_lrt_Tps_sex_bias_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SB_RBBH_Tps_sex_bias$table, "TTT_lrt_Tps_sex_bias_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)




##################################################################################################
####### get (sex-asex) SA_genes

####### design matrix

## using 0SF for sexual female and 1AF for asexual female samples

RM_SA = factor(c(
"0SF","0SF","0SF","1AF","1AF","1AF"
))

data.frame(Sample=colnames(y_WB_SA_RBBH_Tbi),RM_SA)


### design matrix WB 
design_WB_SA_RBBH_Tbi <- model.matrix(~RM_SA) 
rownames(design_WB_SA_RBBH_Tbi) <- colnames(y_WB_SA_RBBH_Tbi)
design_WB_SA_RBBH_Tbi 

### design matrix RT 
design_RT_SA_RBBH_Tbi <- model.matrix(~RM_SA) 
rownames(design_RT_SA_RBBH_Tbi ) <- colnames(y_RT_SA_RBBH_Tbi)
design_RT_SA_RBBH_Tbi 

### design matrix LG 
design_LG_SA_RBBH_Tbi  <- model.matrix(~RM_SA ) 
rownames(design_LG_SA_RBBH_Tbi ) <- colnames(y_LG_SA_RBBH_Tbi )
design_LG_SA_RBBH_Tbi 

###### Est dispersion (samples)

y_WB_SA_RBBH_Tbi  <- estimateDisp(y_WB_SA_RBBH_Tbi, design_WB_SA_RBBH_Tbi)
y_RT_SA_RBBH_Tbi  <- estimateDisp(y_RT_SA_RBBH_Tbi, design_RT_SA_RBBH_Tbi)
y_LG_SA_RBBH_Tbi  <- estimateDisp(y_LG_SA_RBBH_Tbi, design_LG_SA_RBBH_Tbi)

### design matrix WB 
design_WB_SA_RBBH_Tce <- model.matrix(~RM_SA) 
rownames(design_WB_SA_RBBH_Tce) <- colnames(y_WB_SA_RBBH_Tce)
design_WB_SA_RBBH_Tce 

### design matrix RT 
design_RT_SA_RBBH_Tce <- model.matrix(~RM_SA) 
rownames(design_RT_SA_RBBH_Tce ) <- colnames(y_RT_SA_RBBH_Tce)
design_RT_SA_RBBH_Tce 

### design matrix LG 
design_LG_SA_RBBH_Tce  <- model.matrix(~RM_SA ) 
rownames(design_LG_SA_RBBH_Tce ) <- colnames(y_LG_SA_RBBH_Tce )
design_LG_SA_RBBH_Tce 

###### Est dispersion (samples)

y_WB_SA_RBBH_Tce  <- estimateDisp(y_WB_SA_RBBH_Tce, design_WB_SA_RBBH_Tce)
y_RT_SA_RBBH_Tce  <- estimateDisp(y_RT_SA_RBBH_Tce, design_RT_SA_RBBH_Tce)
y_LG_SA_RBBH_Tce  <- estimateDisp(y_LG_SA_RBBH_Tce, design_LG_SA_RBBH_Tce)

### design matrix WB 
design_WB_SA_RBBH_Tcm <- model.matrix(~RM_SA) 
rownames(design_WB_SA_RBBH_Tcm) <- colnames(y_WB_SA_RBBH_Tcm)
design_WB_SA_RBBH_Tcm 

### design matrix RT 
design_RT_SA_RBBH_Tcm <- model.matrix(~RM_SA) 
rownames(design_RT_SA_RBBH_Tcm ) <- colnames(y_RT_SA_RBBH_Tcm)
design_RT_SA_RBBH_Tcm 

### design matrix LG 
design_LG_SA_RBBH_Tcm  <- model.matrix(~RM_SA ) 
rownames(design_LG_SA_RBBH_Tcm ) <- colnames(y_LG_SA_RBBH_Tcm )
design_LG_SA_RBBH_Tcm 

###### Est dispersion (samples)

y_WB_SA_RBBH_Tcm  <- estimateDisp(y_WB_SA_RBBH_Tcm, design_WB_SA_RBBH_Tcm)
y_RT_SA_RBBH_Tcm  <- estimateDisp(y_RT_SA_RBBH_Tcm, design_RT_SA_RBBH_Tcm)
y_LG_SA_RBBH_Tcm  <- estimateDisp(y_LG_SA_RBBH_Tcm, design_LG_SA_RBBH_Tcm)


### design matrix WB 
design_WB_SA_RBBH_Tpa <- model.matrix(~RM_SA) 
rownames(design_WB_SA_RBBH_Tpa) <- colnames(y_WB_SA_RBBH_Tpa)
design_WB_SA_RBBH_Tpa 

### design matrix RT 
design_RT_SA_RBBH_Tpa <- model.matrix(~RM_SA) 
rownames(design_RT_SA_RBBH_Tpa ) <- colnames(y_RT_SA_RBBH_Tpa)
design_RT_SA_RBBH_Tpa 

### design matrix LG 
design_LG_SA_RBBH_Tpa  <- model.matrix(~RM_SA ) 
rownames(design_LG_SA_RBBH_Tpa ) <- colnames(y_LG_SA_RBBH_Tpa )
design_LG_SA_RBBH_Tpa 

###### Est dispersion (samples)

y_WB_SA_RBBH_Tpa  <- estimateDisp(y_WB_SA_RBBH_Tpa, design_WB_SA_RBBH_Tpa)
y_RT_SA_RBBH_Tpa  <- estimateDisp(y_RT_SA_RBBH_Tpa, design_RT_SA_RBBH_Tpa)
y_LG_SA_RBBH_Tpa  <- estimateDisp(y_LG_SA_RBBH_Tpa, design_LG_SA_RBBH_Tpa)


### design matrix WB 
design_WB_SA_RBBH_Tps <- model.matrix(~RM_SA) 
rownames(design_WB_SA_RBBH_Tps) <- colnames(y_WB_SA_RBBH_Tps)
design_WB_SA_RBBH_Tps 

### design matrix RT 
design_RT_SA_RBBH_Tps <- model.matrix(~RM_SA) 
rownames(design_RT_SA_RBBH_Tps ) <- colnames(y_RT_SA_RBBH_Tps)
design_RT_SA_RBBH_Tps 

### design matrix LG 
design_LG_SA_RBBH_Tps  <- model.matrix(~RM_SA ) 
rownames(design_LG_SA_RBBH_Tps ) <- colnames(y_LG_SA_RBBH_Tps )
design_LG_SA_RBBH_Tps 

###### Est dispersion (samples)

y_WB_SA_RBBH_Tps  <- estimateDisp(y_WB_SA_RBBH_Tps, design_WB_SA_RBBH_Tps)
y_RT_SA_RBBH_Tps  <- estimateDisp(y_RT_SA_RBBH_Tps, design_RT_SA_RBBH_Tps)
y_LG_SA_RBBH_Tps  <- estimateDisp(y_LG_SA_RBBH_Tps, design_LG_SA_RBBH_Tps)



#####################################################################################################################################################################
###### fit glm 

fit_WB_SA_RBBH_Tbi <- glmFit(y_WB_SA_RBBH_Tbi, design_WB_SA_RBBH_Tbi,robust=TRUE)
fit_RT_SA_RBBH_Tbi <- glmFit(y_RT_SA_RBBH_Tbi, design_RT_SA_RBBH_Tbi,robust=TRUE)
fit_LG_SA_RBBH_Tbi <- glmFit(y_LG_SA_RBBH_Tbi, design_LG_SA_RBBH_Tbi,robust=TRUE)

fit_WB_SA_RBBH_Tce <- glmFit(y_WB_SA_RBBH_Tce, design_WB_SA_RBBH_Tce,robust=TRUE)
fit_RT_SA_RBBH_Tce <- glmFit(y_RT_SA_RBBH_Tce, design_RT_SA_RBBH_Tce,robust=TRUE)
fit_LG_SA_RBBH_Tce <- glmFit(y_LG_SA_RBBH_Tce, design_LG_SA_RBBH_Tce,robust=TRUE)

fit_WB_SA_RBBH_Tcm <- glmFit(y_WB_SA_RBBH_Tcm, design_WB_SA_RBBH_Tcm,robust=TRUE)
fit_RT_SA_RBBH_Tcm <- glmFit(y_RT_SA_RBBH_Tcm, design_RT_SA_RBBH_Tcm,robust=TRUE)
fit_LG_SA_RBBH_Tcm <- glmFit(y_LG_SA_RBBH_Tcm, design_LG_SA_RBBH_Tcm,robust=TRUE)

fit_WB_SA_RBBH_Tpa <- glmFit(y_WB_SA_RBBH_Tpa, design_WB_SA_RBBH_Tpa,robust=TRUE)
fit_RT_SA_RBBH_Tpa <- glmFit(y_RT_SA_RBBH_Tpa, design_RT_SA_RBBH_Tpa,robust=TRUE)
fit_LG_SA_RBBH_Tpa <- glmFit(y_LG_SA_RBBH_Tpa, design_LG_SA_RBBH_Tpa,robust=TRUE)

fit_WB_SA_RBBH_Tps <- glmFit(y_WB_SA_RBBH_Tps, design_WB_SA_RBBH_Tps,robust=TRUE)
fit_RT_SA_RBBH_Tps <- glmFit(y_RT_SA_RBBH_Tps, design_RT_SA_RBBH_Tps,robust=TRUE)
fit_LG_SA_RBBH_Tps <- glmFit(y_LG_SA_RBBH_Tps, design_LG_SA_RBBH_Tps,robust=TRUE)

colnames(fit_LG_SA_RBBH_Tps)
colnames(fit_WB_SA_RBBH_Tbi)

# [1] "(Intercept)" "RM_SA1AF"     


## get Sex-asex genes ### NOTE +ve FC = higher experssion in ASEX

fit_WB_SA_RBBH_Tbi_c <- glmLRT(fit_WB_SA_RBBH_Tbi,coef=2)
fit_RT_SA_RBBH_Tbi_c <- glmLRT(fit_RT_SA_RBBH_Tbi,coef=2)
fit_LG_SA_RBBH_Tbi_c <- glmLRT(fit_LG_SA_RBBH_Tbi,coef=2)

TTT_WB_SA_RBBH_Tbi_sex_asex <- get_DE_genes(fit_WB_SA_RBBH_Tbi_c , 0.05)
TTT_RT_SA_RBBH_Tbi_sex_asex <- get_DE_genes(fit_RT_SA_RBBH_Tbi_c , 0.05)
TTT_LG_SA_RBBH_Tbi_sex_asex <- get_DE_genes(fit_LG_SA_RBBH_Tbi_c , 0.05)

fit_WB_SA_RBBH_Tce_c <- glmLRT(fit_WB_SA_RBBH_Tce,coef=2)
fit_RT_SA_RBBH_Tce_c <- glmLRT(fit_RT_SA_RBBH_Tce,coef=2)
fit_LG_SA_RBBH_Tce_c <- glmLRT(fit_LG_SA_RBBH_Tce,coef=2)

TTT_WB_SA_RBBH_Tce_sex_asex <- get_DE_genes(fit_WB_SA_RBBH_Tce_c , 0.05)
TTT_RT_SA_RBBH_Tce_sex_asex <- get_DE_genes(fit_RT_SA_RBBH_Tce_c , 0.05)
TTT_LG_SA_RBBH_Tce_sex_asex <- get_DE_genes(fit_LG_SA_RBBH_Tce_c , 0.05)

fit_WB_SA_RBBH_Tcm_c <- glmLRT(fit_WB_SA_RBBH_Tcm,coef=2)
fit_RT_SA_RBBH_Tcm_c <- glmLRT(fit_RT_SA_RBBH_Tcm,coef=2)
fit_LG_SA_RBBH_Tcm_c <- glmLRT(fit_LG_SA_RBBH_Tcm,coef=2)

TTT_WB_SA_RBBH_Tcm_sex_asex <- get_DE_genes(fit_WB_SA_RBBH_Tcm_c , 0.05)
TTT_RT_SA_RBBH_Tcm_sex_asex <- get_DE_genes(fit_RT_SA_RBBH_Tcm_c , 0.05)
TTT_LG_SA_RBBH_Tcm_sex_asex <- get_DE_genes(fit_LG_SA_RBBH_Tcm_c , 0.05)

fit_WB_SA_RBBH_Tpa_c <- glmLRT(fit_WB_SA_RBBH_Tpa,coef=2)
fit_RT_SA_RBBH_Tpa_c <- glmLRT(fit_RT_SA_RBBH_Tpa,coef=2)
fit_LG_SA_RBBH_Tpa_c <- glmLRT(fit_LG_SA_RBBH_Tpa,coef=2)

TTT_WB_SA_RBBH_Tpa_sex_asex <- get_DE_genes(fit_WB_SA_RBBH_Tpa_c , 0.05)
TTT_RT_SA_RBBH_Tpa_sex_asex <- get_DE_genes(fit_RT_SA_RBBH_Tpa_c , 0.05)
TTT_LG_SA_RBBH_Tpa_sex_asex <- get_DE_genes(fit_LG_SA_RBBH_Tpa_c , 0.05)

fit_WB_SA_RBBH_Tps_c <- glmLRT(fit_WB_SA_RBBH_Tps,coef=2)
fit_RT_SA_RBBH_Tps_c <- glmLRT(fit_RT_SA_RBBH_Tps,coef=2)
fit_LG_SA_RBBH_Tps_c <- glmLRT(fit_LG_SA_RBBH_Tps,coef=2)

TTT_WB_SA_RBBH_Tps_sex_asex <- get_DE_genes(fit_WB_SA_RBBH_Tps_c , 0.05)
TTT_RT_SA_RBBH_Tps_sex_asex <- get_DE_genes(fit_RT_SA_RBBH_Tps_c , 0.05)
TTT_LG_SA_RBBH_Tps_sex_asex <- get_DE_genes(fit_LG_SA_RBBH_Tps_c , 0.05)

####### export full tables

write.table(TTT_WB_SA_RBBH_Tbi_sex_asex$table, "TTT_lrt_Tbi_sex_asex_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SA_RBBH_Tbi_sex_asex$table, "TTT_lrt_Tbi_sex_asex_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SA_RBBH_Tbi_sex_asex$table, "TTT_lrt_Tbi_sex_asex_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)

write.table(TTT_WB_SA_RBBH_Tce_sex_asex$table, "TTT_lrt_Tce_sex_asex_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SA_RBBH_Tce_sex_asex$table, "TTT_lrt_Tce_sex_asex_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SA_RBBH_Tce_sex_asex$table, "TTT_lrt_Tce_sex_asex_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)

write.table(TTT_WB_SA_RBBH_Tcm_sex_asex$table, "TTT_lrt_Tcm_sex_asex_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SA_RBBH_Tcm_sex_asex$table, "TTT_lrt_Tcm_sex_asex_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SA_RBBH_Tcm_sex_asex$table, "TTT_lrt_Tcm_sex_asex_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)

write.table(TTT_WB_SA_RBBH_Tpa_sex_asex$table, "TTT_lrt_Tpa_sex_asex_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SA_RBBH_Tpa_sex_asex$table, "TTT_lrt_Tpa_sex_asex_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SA_RBBH_Tpa_sex_asex$table, "TTT_lrt_Tpa_sex_asex_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)

write.table(TTT_WB_SA_RBBH_Tps_sex_asex$table, "TTT_lrt_Tps_sex_asex_WB.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_SA_RBBH_Tps_sex_asex$table, "TTT_lrt_Tps_sex_asex_RT.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_SA_RBBH_Tps_sex_asex$table, "TTT_lrt_Tps_sex_asex_LG.csv", sep = ',', quote = FALSE, row.names = FALSE)




#######################################################################################################################
##### summarise the data a bit

## N sex-biased genes

N_sex_bias_genes <- as.data.frame(c(
length(TTT_WB_SB_RBBH_Tbi_sex_bias$S_gene_list),
length(TTT_WB_SB_RBBH_Tce_sex_bias$S_gene_list),
length(TTT_WB_SB_RBBH_Tcm_sex_bias$S_gene_list),
length(TTT_WB_SB_RBBH_Tpa_sex_bias$S_gene_list),
length(TTT_WB_SB_RBBH_Tps_sex_bias$S_gene_list),
length(TTT_RT_SB_RBBH_Tbi_sex_bias$S_gene_list),
length(TTT_RT_SB_RBBH_Tce_sex_bias$S_gene_list),
length(TTT_RT_SB_RBBH_Tcm_sex_bias$S_gene_list),
length(TTT_RT_SB_RBBH_Tpa_sex_bias$S_gene_list),
length(TTT_RT_SB_RBBH_Tps_sex_bias$S_gene_list),
length(TTT_LG_SB_RBBH_Tbi_sex_bias$S_gene_list),
length(TTT_LG_SB_RBBH_Tce_sex_bias$S_gene_list),
length(TTT_LG_SB_RBBH_Tcm_sex_bias$S_gene_list),
length(TTT_LG_SB_RBBH_Tpa_sex_bias$S_gene_list),
length(TTT_LG_SB_RBBH_Tps_sex_bias$S_gene_list)
))


N_sex_bias_genes <- as.data.frame(cbind(N_sex_bias_genes,c(
length(TTT_WB_SB_RBBH_Tbi_sex_bias$table[,1]),
length(TTT_WB_SB_RBBH_Tce_sex_bias$table[,1]),
length(TTT_WB_SB_RBBH_Tcm_sex_bias$table[,1]),
length(TTT_WB_SB_RBBH_Tpa_sex_bias$table[,1]),
length(TTT_WB_SB_RBBH_Tps_sex_bias$table[,1]),
length(TTT_RT_SB_RBBH_Tbi_sex_bias$table[,1]),
length(TTT_RT_SB_RBBH_Tce_sex_bias$table[,1]),
length(TTT_RT_SB_RBBH_Tcm_sex_bias$table[,1]),
length(TTT_RT_SB_RBBH_Tpa_sex_bias$table[,1]),
length(TTT_RT_SB_RBBH_Tps_sex_bias$table[,1]),
length(TTT_LG_SB_RBBH_Tbi_sex_bias$table[,1]),
length(TTT_LG_SB_RBBH_Tce_sex_bias$table[,1]),
length(TTT_LG_SB_RBBH_Tcm_sex_bias$table[,1]),
length(TTT_LG_SB_RBBH_Tpa_sex_bias$table[,1]),
length(TTT_LG_SB_RBBH_Tps_sex_bias$table[,1])
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


write.csv(N_sex_bias_genes, "N_sex_bias_genes_RBBH.csv")



p3 <- ggplot(N_sex_bias_genes, aes(factor(group_ordered), N_sex_bias_genes, fill = tiss )) + 
	geom_bar(stat="identity", position = "dodge") + 
	theme_bw() +
	xlab ("Species pair") + 
	ylab ("Number of sex biased genes, FDR < 0.05")  +
	scale_y_continuous(expand = c(0,0)) + 
	coord_cartesian(ylim=c(-0,8000)) 
	

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
writeLines(capture.output(sessionInfo()), "sex_bias_edgeR_sessionInfo.txt")

