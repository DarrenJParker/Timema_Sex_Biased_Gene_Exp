### divergence_and_change_in_SB

library(ggplot2)
library(pheatmap)
library(cowplot)
library(grid)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(cowplot)
library(gtable)
library(RColorBrewer)
library(pvclust)
library("VennDiagram")
library("SuperExactTest")
library(raster)

print (sessionInfo())

# R version 3.5.1 (2018-07-02)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.14.3

# Matrix products: default
# BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
 # [1] raster_2.8-4         sp_1.3-1             SuperExactTest_1.0.4 VennDiagram_1.6.20   futile.logger_1.4.3  pvclust_2.0-0        RColorBrewer_1.1-2  
 # [8] gtable_0.2.0         lattice_0.20-38      gridExtra_2.3        cowplot_0.9.3        pheatmap_1.0.10      ggplot2_3.1.0       

# loaded via a namespace (and not attached):
 # [1] Rcpp_1.0.0           pillar_1.3.0         compiler_3.5.1       formatR_1.5          plyr_1.8.4           bindr_0.1.1          futile.options_1.0.1
 # [8] digest_0.6.18        tibble_1.4.2         pkgconfig_2.0.2      rlang_0.3.0.1        bindrcpp_0.2.2       withr_2.1.2          dplyr_0.7.8         
# [15] tidyselect_0.2.5     glue_1.3.0           R6_2.3.0             purrr_0.2.5          lambda.r_1.2.3       magrittr_1.5         scales_1.0.0        
# [22] codetools_0.2-15     assertthat_0.2.0     colorspace_1.3-2     labeling_0.3         lazyeval_0.2.1       munsell_0.5.0        crayon_1.3.4        
# > 

###### fun

five_sp_DE_venn <- function(Tbi,Tce,Tcm,Tpa,Tps,title){
	venny.plot <- venn.diagram(
	list("Tbi" = Tbi, "Tce" = Tce, "Tcm" = Tcm, "Tpa" = Tpa, "Tps" = Tps ), filename = NULL,
                            fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                            cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                            margin = 0.6, cat.dist = 0.23, main = title, main.pos = c(0.5,0.8), main.cex = 2, main.fontface = "bold", cat.cex = 2)
	return(venny.plot)
}

theme_update(plot.title = element_text(hjust = 0.5))



################################################################################################################################################
##### read in data


setwd("output/DE_joined/")

dat_Tbi_WB_RBBH_SB_asex = read.csv("Tbi_WB_RBBH_disp_allsepar_SB_asex.csv")
dat_Tbi_RT_RBBH_SB_asex = read.csv("Tbi_RT_RBBH_disp_allsepar_SB_asex.csv")
dat_Tbi_LG_RBBH_SB_asex = read.csv("Tbi_LG_RBBH_disp_allsepar_SB_asex.csv")

dat_Tce_WB_RBBH_SB_asex = read.csv("Tce_WB_RBBH_disp_allsepar_SB_asex.csv")
dat_Tce_RT_RBBH_SB_asex = read.csv("Tce_RT_RBBH_disp_allsepar_SB_asex.csv")
dat_Tce_LG_RBBH_SB_asex = read.csv("Tce_LG_RBBH_disp_allsepar_SB_asex.csv")

dat_Tcm_WB_RBBH_SB_asex = read.csv("Tcm_WB_RBBH_disp_allsepar_SB_asex.csv")
dat_Tcm_RT_RBBH_SB_asex = read.csv("Tcm_RT_RBBH_disp_allsepar_SB_asex.csv")
dat_Tcm_LG_RBBH_SB_asex = read.csv("Tcm_LG_RBBH_disp_allsepar_SB_asex.csv")

dat_Tpa_WB_RBBH_SB_asex = read.csv("Tpa_WB_RBBH_disp_allsepar_SB_asex.csv")
dat_Tpa_RT_RBBH_SB_asex = read.csv("Tpa_RT_RBBH_disp_allsepar_SB_asex.csv")
dat_Tpa_LG_RBBH_SB_asex = read.csv("Tpa_LG_RBBH_disp_allsepar_SB_asex.csv")

dat_Tps_WB_RBBH_SB_asex = read.csv("Tps_WB_RBBH_disp_allsepar_SB_asex.csv")
dat_Tps_RT_RBBH_SB_asex = read.csv("Tps_RT_RBBH_disp_allsepar_SB_asex.csv")
dat_Tps_LG_RBBH_SB_asex = read.csv("Tps_LG_RBBH_disp_allsepar_SB_asex.csv")


head(dat_Tbi_WB_RBBH_SB_asex, n = 20)

#################################################################################################################################################
###### refine sex bias by FC (not log FC)
	
FC_cutoff = 2

dat_Tbi_WB_RBBH_SB_asex$Tbi_WB_sexbias2 <- ifelse(2^(sqrt(dat_Tbi_WB_RBBH_SB_asex$Tbi_WB_log2FC_SB * dat_Tbi_WB_RBBH_SB_asex$Tbi_WB_log2FC_SB)) > FC_cutoff , as.character(dat_Tbi_WB_RBBH_SB_asex$Tbi_WB_sexbias), "Unbiased")
dat_Tce_WB_RBBH_SB_asex$Tce_WB_sexbias2 <- ifelse(2^(sqrt(dat_Tce_WB_RBBH_SB_asex$Tce_WB_log2FC_SB * dat_Tce_WB_RBBH_SB_asex$Tce_WB_log2FC_SB)) > FC_cutoff , as.character(dat_Tce_WB_RBBH_SB_asex$Tce_WB_sexbias), "Unbiased")
dat_Tcm_WB_RBBH_SB_asex$Tcm_WB_sexbias2 <- ifelse(2^(sqrt(dat_Tcm_WB_RBBH_SB_asex$Tcm_WB_log2FC_SB * dat_Tcm_WB_RBBH_SB_asex$Tcm_WB_log2FC_SB)) > FC_cutoff , as.character(dat_Tcm_WB_RBBH_SB_asex$Tcm_WB_sexbias), "Unbiased")
dat_Tpa_WB_RBBH_SB_asex$Tpa_WB_sexbias2 <- ifelse(2^(sqrt(dat_Tpa_WB_RBBH_SB_asex$Tpa_WB_log2FC_SB * dat_Tpa_WB_RBBH_SB_asex$Tpa_WB_log2FC_SB)) > FC_cutoff , as.character(dat_Tpa_WB_RBBH_SB_asex$Tpa_WB_sexbias), "Unbiased")
dat_Tps_WB_RBBH_SB_asex$Tps_WB_sexbias2 <- ifelse(2^(sqrt(dat_Tps_WB_RBBH_SB_asex$Tps_WB_log2FC_SB * dat_Tps_WB_RBBH_SB_asex$Tps_WB_log2FC_SB)) > FC_cutoff , as.character(dat_Tps_WB_RBBH_SB_asex$Tps_WB_sexbias), "Unbiased")

dat_Tbi_RT_RBBH_SB_asex$Tbi_RT_sexbias2 <- ifelse(2^(sqrt(dat_Tbi_RT_RBBH_SB_asex$Tbi_RT_log2FC_SB * dat_Tbi_RT_RBBH_SB_asex$Tbi_RT_log2FC_SB)) > FC_cutoff , as.character(dat_Tbi_RT_RBBH_SB_asex$Tbi_RT_sexbias), "Unbiased")
dat_Tce_RT_RBBH_SB_asex$Tce_RT_sexbias2 <- ifelse(2^(sqrt(dat_Tce_RT_RBBH_SB_asex$Tce_RT_log2FC_SB * dat_Tce_RT_RBBH_SB_asex$Tce_RT_log2FC_SB)) > FC_cutoff , as.character(dat_Tce_RT_RBBH_SB_asex$Tce_RT_sexbias), "Unbiased")
dat_Tcm_RT_RBBH_SB_asex$Tcm_RT_sexbias2 <- ifelse(2^(sqrt(dat_Tcm_RT_RBBH_SB_asex$Tcm_RT_log2FC_SB * dat_Tcm_RT_RBBH_SB_asex$Tcm_RT_log2FC_SB)) > FC_cutoff , as.character(dat_Tcm_RT_RBBH_SB_asex$Tcm_RT_sexbias), "Unbiased")
dat_Tpa_RT_RBBH_SB_asex$Tpa_RT_sexbias2 <- ifelse(2^(sqrt(dat_Tpa_RT_RBBH_SB_asex$Tpa_RT_log2FC_SB * dat_Tpa_RT_RBBH_SB_asex$Tpa_RT_log2FC_SB)) > FC_cutoff , as.character(dat_Tpa_RT_RBBH_SB_asex$Tpa_RT_sexbias), "Unbiased")
dat_Tps_RT_RBBH_SB_asex$Tps_RT_sexbias2 <- ifelse(2^(sqrt(dat_Tps_RT_RBBH_SB_asex$Tps_RT_log2FC_SB * dat_Tps_RT_RBBH_SB_asex$Tps_RT_log2FC_SB)) > FC_cutoff , as.character(dat_Tps_RT_RBBH_SB_asex$Tps_RT_sexbias), "Unbiased")


dat_Tbi_LG_RBBH_SB_asex$Tbi_LG_sexbias2 <- ifelse(2^(sqrt(dat_Tbi_LG_RBBH_SB_asex$Tbi_LG_log2FC_SB * dat_Tbi_LG_RBBH_SB_asex$Tbi_LG_log2FC_SB)) > FC_cutoff , as.character(dat_Tbi_LG_RBBH_SB_asex$Tbi_LG_sexbias), "Unbiased")
dat_Tce_LG_RBBH_SB_asex$Tce_LG_sexbias2 <- ifelse(2^(sqrt(dat_Tce_LG_RBBH_SB_asex$Tce_LG_log2FC_SB * dat_Tce_LG_RBBH_SB_asex$Tce_LG_log2FC_SB)) > FC_cutoff , as.character(dat_Tce_LG_RBBH_SB_asex$Tce_LG_sexbias), "Unbiased")
dat_Tcm_LG_RBBH_SB_asex$Tcm_LG_sexbias2 <- ifelse(2^(sqrt(dat_Tcm_LG_RBBH_SB_asex$Tcm_LG_log2FC_SB * dat_Tcm_LG_RBBH_SB_asex$Tcm_LG_log2FC_SB)) > FC_cutoff , as.character(dat_Tcm_LG_RBBH_SB_asex$Tcm_LG_sexbias), "Unbiased")
dat_Tpa_LG_RBBH_SB_asex$Tpa_LG_sexbias2 <- ifelse(2^(sqrt(dat_Tpa_LG_RBBH_SB_asex$Tpa_LG_log2FC_SB * dat_Tpa_LG_RBBH_SB_asex$Tpa_LG_log2FC_SB)) > FC_cutoff , as.character(dat_Tpa_LG_RBBH_SB_asex$Tpa_LG_sexbias), "Unbiased")
dat_Tps_LG_RBBH_SB_asex$Tps_LG_sexbias2 <- ifelse(2^(sqrt(dat_Tps_LG_RBBH_SB_asex$Tps_LG_log2FC_SB * dat_Tps_LG_RBBH_SB_asex$Tps_LG_log2FC_SB)) > FC_cutoff , as.character(dat_Tps_LG_RBBH_SB_asex$Tps_LG_sexbias), "Unbiased")

head(dat_Tpa_WB_RBBH_SB_asex)

#### add divergence time (values from Bast J, Parker DJ, Dumas Z, Jalvingh KM, Tran Van P, Jaron KS, et al. (2018). Consequences of asexuality in natural populations: insights from stick insects. Mol Biol Evol 35: 1668–1677.)

dat_Tbi_WB_RBBH_SB_asex$Tbi_WB_div	<- rep(0.791674, length(dat_Tbi_WB_RBBH_SB_asex[,1]))
dat_Tce_WB_RBBH_SB_asex$Tce_WB_div	<- rep(1.10925,  length(dat_Tce_WB_RBBH_SB_asex[,1]))
dat_Tcm_WB_RBBH_SB_asex$Tcm_WB_div	<- rep(1.40629,  length(dat_Tcm_WB_RBBH_SB_asex[,1]))
dat_Tpa_WB_RBBH_SB_asex$Tpa_WB_div	<- rep(1.86055,  length(dat_Tpa_WB_RBBH_SB_asex[,1]))
dat_Tps_WB_RBBH_SB_asex$Tps_WB_div	<- rep(1.22685,  length(dat_Tps_WB_RBBH_SB_asex[,1]))
	
dat_Tbi_RT_RBBH_SB_asex$Tbi_RT_div	<- rep(0.791674, length(dat_Tbi_RT_RBBH_SB_asex[,1]))
dat_Tce_RT_RBBH_SB_asex$Tce_RT_div	<- rep(1.10925,  length(dat_Tce_RT_RBBH_SB_asex[,1]))
dat_Tcm_RT_RBBH_SB_asex$Tcm_RT_div	<- rep(1.40629,  length(dat_Tcm_RT_RBBH_SB_asex[,1]))
dat_Tpa_RT_RBBH_SB_asex$Tpa_RT_div	<- rep(1.86055,  length(dat_Tpa_RT_RBBH_SB_asex[,1]))
dat_Tps_RT_RBBH_SB_asex$Tps_RT_div	<- rep(1.22685,  length(dat_Tps_RT_RBBH_SB_asex[,1]))
	
dat_Tbi_LG_RBBH_SB_asex$Tbi_LG_div	<- rep(0.791674, length(dat_Tbi_LG_RBBH_SB_asex[,1]))
dat_Tce_LG_RBBH_SB_asex$Tce_LG_div	<- rep(1.10925,  length(dat_Tce_LG_RBBH_SB_asex[,1]))
dat_Tcm_LG_RBBH_SB_asex$Tcm_LG_div	<- rep(1.40629,  length(dat_Tcm_LG_RBBH_SB_asex[,1]))
dat_Tpa_LG_RBBH_SB_asex$Tpa_LG_div	<- rep(1.86055,  length(dat_Tpa_LG_RBBH_SB_asex[,1]))
dat_Tps_LG_RBBH_SB_asex$Tps_LG_div	<- rep(1.22685,  length(dat_Tps_LG_RBBH_SB_asex[,1]))

head(dat_Tps_LG_RBBH_SB_asex)

dat_Tbi_WB_RBBH_SB_asex$genename 	<- as.character(dat_Tbi_WB_RBBH_SB_asex$genename)
dat_Tce_WB_RBBH_SB_asex$genename 	<- as.character(dat_Tce_WB_RBBH_SB_asex$genename)
dat_Tcm_WB_RBBH_SB_asex$genename 	<- as.character(dat_Tcm_WB_RBBH_SB_asex$genename)
dat_Tpa_WB_RBBH_SB_asex$genename 	<- as.character(dat_Tpa_WB_RBBH_SB_asex$genename)
dat_Tps_WB_RBBH_SB_asex$genename 	<- as.character(dat_Tps_WB_RBBH_SB_asex$genename)

dat_Tbi_RT_RBBH_SB_asex$genename 	<- as.character(dat_Tbi_RT_RBBH_SB_asex$genename)
dat_Tce_RT_RBBH_SB_asex$genename 	<- as.character(dat_Tce_RT_RBBH_SB_asex$genename)
dat_Tcm_RT_RBBH_SB_asex$genename 	<- as.character(dat_Tcm_RT_RBBH_SB_asex$genename)
dat_Tpa_RT_RBBH_SB_asex$genename 	<- as.character(dat_Tpa_RT_RBBH_SB_asex$genename)
dat_Tps_RT_RBBH_SB_asex$genename 	<- as.character(dat_Tps_RT_RBBH_SB_asex$genename)

dat_Tbi_LG_RBBH_SB_asex$genename 	<- as.character(dat_Tbi_LG_RBBH_SB_asex$genename)
dat_Tce_LG_RBBH_SB_asex$genename 	<- as.character(dat_Tce_LG_RBBH_SB_asex$genename)
dat_Tcm_LG_RBBH_SB_asex$genename 	<- as.character(dat_Tcm_LG_RBBH_SB_asex$genename)
dat_Tpa_LG_RBBH_SB_asex$genename 	<- as.character(dat_Tpa_LG_RBBH_SB_asex$genename)
dat_Tps_LG_RBBH_SB_asex$genename 	<- as.character(dat_Tps_LG_RBBH_SB_asex$genename)

dat_Tbi_WB_RBBH_SB_asex$Tbi_WB_sexbias 	<- as.character(dat_Tbi_WB_RBBH_SB_asex$Tbi_WB_sexbias)
dat_Tce_WB_RBBH_SB_asex$Tce_WB_sexbias 	<- as.character(dat_Tce_WB_RBBH_SB_asex$Tce_WB_sexbias)
dat_Tcm_WB_RBBH_SB_asex$Tcm_WB_sexbias 	<- as.character(dat_Tcm_WB_RBBH_SB_asex$Tcm_WB_sexbias)
dat_Tpa_WB_RBBH_SB_asex$Tpa_WB_sexbias 	<- as.character(dat_Tpa_WB_RBBH_SB_asex$Tpa_WB_sexbias)
dat_Tps_WB_RBBH_SB_asex$Tps_WB_sexbias 	<- as.character(dat_Tps_WB_RBBH_SB_asex$Tps_WB_sexbias)

dat_Tbi_RT_RBBH_SB_asex$Tbi_RT_sexbias 	<- as.character(dat_Tbi_RT_RBBH_SB_asex$Tbi_RT_sexbias)
dat_Tce_RT_RBBH_SB_asex$Tce_RT_sexbias 	<- as.character(dat_Tce_RT_RBBH_SB_asex$Tce_RT_sexbias)
dat_Tcm_RT_RBBH_SB_asex$Tcm_RT_sexbias 	<- as.character(dat_Tcm_RT_RBBH_SB_asex$Tcm_RT_sexbias)
dat_Tpa_RT_RBBH_SB_asex$Tpa_RT_sexbias 	<- as.character(dat_Tpa_RT_RBBH_SB_asex$Tpa_RT_sexbias)
dat_Tps_RT_RBBH_SB_asex$Tps_RT_sexbias 	<- as.character(dat_Tps_RT_RBBH_SB_asex$Tps_RT_sexbias)

dat_Tbi_LG_RBBH_SB_asex$Tbi_LG_sexbias 	<- as.character(dat_Tbi_LG_RBBH_SB_asex$Tbi_LG_sexbias)
dat_Tce_LG_RBBH_SB_asex$Tce_LG_sexbias 	<- as.character(dat_Tce_LG_RBBH_SB_asex$Tce_LG_sexbias)
dat_Tcm_LG_RBBH_SB_asex$Tcm_LG_sexbias 	<- as.character(dat_Tcm_LG_RBBH_SB_asex$Tcm_LG_sexbias)
dat_Tpa_LG_RBBH_SB_asex$Tpa_LG_sexbias 	<- as.character(dat_Tpa_LG_RBBH_SB_asex$Tpa_LG_sexbias)
dat_Tps_LG_RBBH_SB_asex$Tps_LG_sexbias 	<- as.character(dat_Tps_LG_RBBH_SB_asex$Tps_LG_sexbias)

#### cat together 

head(dat_Tbi_WB_RBBH_SB_asex)
str(dat_Tbi_WB_RBBH_SB_asex)

dat_allsp_WB_RBBH_SB_asex <- as.data.frame(mapply(c,
dat_Tbi_WB_RBBH_SB_asex,
dat_Tce_WB_RBBH_SB_asex,
dat_Tcm_WB_RBBH_SB_asex,
dat_Tpa_WB_RBBH_SB_asex,
dat_Tps_WB_RBBH_SB_asex
))

dat_allsp_RT_RBBH_SB_asex <- as.data.frame(mapply(c,
dat_Tbi_RT_RBBH_SB_asex,
dat_Tce_RT_RBBH_SB_asex,
dat_Tcm_RT_RBBH_SB_asex,
dat_Tpa_RT_RBBH_SB_asex,
dat_Tps_RT_RBBH_SB_asex
))

dat_allsp_LG_RBBH_SB_asex <- as.data.frame(mapply(c,
dat_Tbi_LG_RBBH_SB_asex,
dat_Tce_LG_RBBH_SB_asex,
dat_Tcm_LG_RBBH_SB_asex,
dat_Tpa_LG_RBBH_SB_asex,
dat_Tps_LG_RBBH_SB_asex
))

dat_allsp_WBRTLG_RBBH_SB_asex <- as.data.frame(mapply(c,
dat_Tbi_WB_RBBH_SB_asex,
dat_Tce_WB_RBBH_SB_asex,
dat_Tcm_WB_RBBH_SB_asex,
dat_Tpa_WB_RBBH_SB_asex,
dat_Tps_WB_RBBH_SB_asex,
dat_Tbi_RT_RBBH_SB_asex,
dat_Tce_RT_RBBH_SB_asex,
dat_Tcm_RT_RBBH_SB_asex,
dat_Tpa_RT_RBBH_SB_asex,
dat_Tps_RT_RBBH_SB_asex,
dat_Tbi_LG_RBBH_SB_asex,
dat_Tce_LG_RBBH_SB_asex,
dat_Tcm_LG_RBBH_SB_asex,
dat_Tpa_LG_RBBH_SB_asex,
dat_Tps_LG_RBBH_SB_asex
))


colnames(dat_allsp_WB_RBBH_SB_asex) <- c("genename" ,"log2FC_SB", "FDR_SB", "sexbias", "log2FC_SA" ,"FDR_SA", "FPKM_SF_WB_Md_Re1", "FPKM_SF_WB_Md_Re2", "FPKM_SF_WB_Md_Re3", "FPKM_SM_WB_Md_Re1", "FPKM_SM_WB_Md_Re2", "FPKM_SM_WB_Md_Re3", "FPKM_AF_WB_Vi_Re1", "FPKM_AF_WB_Vi_Re2", "FPKM_AF_WB_Vi_Re3", "sexbias2", "div" )
colnames(dat_allsp_RT_RBBH_SB_asex) <- c("genename" ,"log2FC_SB", "FDR_SB", "sexbias", "log2FC_SA" ,"FDR_SA", "FPKM_SF_RT_Md_Re1", "FPKM_SF_RT_Md_Re2", "FPKM_SF_RT_Md_Re3", "FPKM_SM_RT_Md_Re1", "FPKM_SM_RT_Md_Re2", "FPKM_SM_RT_Md_Re3", "FPKM_AF_RT_Vi_Re1", "FPKM_AF_RT_Vi_Re2", "FPKM_AF_RT_Vi_Re3","sexbias2", "div" )
colnames(dat_allsp_LG_RBBH_SB_asex) <- c("genename" ,"log2FC_SB", "FDR_SB", "sexbias", "log2FC_SA" ,"FDR_SA", "FPKM_SF_LG_Md_Re1", "FPKM_SF_LG_Md_Re2", "FPKM_SF_LG_Md_Re3", "FPKM_SM_LG_Md_Re1", "FPKM_SM_LG_Md_Re2", "FPKM_SM_LG_Md_Re3", "FPKM_AF_LG_Vi_Re1", "FPKM_AF_LG_Vi_Re2", "FPKM_AF_LG_Vi_Re3","sexbias2", "div" )
colnames(dat_allsp_WBRTLG_RBBH_SB_asex) <- c("genename" ,"log2FC_SB", "FDR_SB", "sexbias", "log2FC_SA" ,"FDR_SA", "FPKM_SF_Md_Re1", "FPKM_SF_Md_Re2", "FPKM_SF_Md_Re3", "FPKM_SM_Md_Re1", "FPKM_SM_Md_Re2", "FPKM_SM_Md_Re3", "FPKM_AF_Vi_Re1", "FPKM_AF_Vi_Re2", "FPKM_AF_Vi_Re3","sexbias2", "div" )


###
dat_allsp_WBRTLG_RBBH_SB_asex$tiss <- c(rep("WB", length(dat_allsp_WB_RBBH_SB_asex[,1])), rep("RT", length(dat_allsp_RT_RBBH_SB_asex[,1])), rep("LG", length(dat_allsp_LG_RBBH_SB_asex[,1])))
dat_allsp_WBRTLG_RBBH_SB_asex$tiss <- as.factor(dat_allsp_WBRTLG_RBBH_SB_asex$tiss )


dat_allsp_WB_RBBH_SB_asex$log2FC_SA <- as.numeric(as.character(dat_allsp_WB_RBBH_SB_asex$log2FC_SA))
dat_allsp_WB_RBBH_SB_asex$div <- as.numeric(as.character(dat_allsp_WB_RBBH_SB_asex$div))
dat_allsp_RT_RBBH_SB_asex$log2FC_SA <- as.numeric(as.character(dat_allsp_RT_RBBH_SB_asex$log2FC_SA))
dat_allsp_RT_RBBH_SB_asex$div <- as.numeric(as.character(dat_allsp_RT_RBBH_SB_asex$div))
dat_allsp_LG_RBBH_SB_asex$log2FC_SA <- as.numeric(as.character(dat_allsp_LG_RBBH_SB_asex$log2FC_SA))
dat_allsp_LG_RBBH_SB_asex$div <- as.numeric(as.character(dat_allsp_LG_RBBH_SB_asex$div))
dat_allsp_WBRTLG_RBBH_SB_asex$log2FC_SA <- as.numeric(as.character(dat_allsp_WBRTLG_RBBH_SB_asex$log2FC_SA))
dat_allsp_WBRTLG_RBBH_SB_asex$div <- as.numeric(as.character(dat_allsp_WBRTLG_RBBH_SB_asex$div))


### ANOVA
str(dat_allsp_WBRTLG_RBBH_SB_asex)

dat_allsp_WBRTLG_RBBH_SB_asex_MB <- subset(dat_allsp_WBRTLG_RBBH_SB_asex, dat_allsp_WBRTLG_RBBH_SB_asex$sexbias2 == "male_biased")
dat_allsp_WBRTLG_RBBH_SB_asex_FB <- subset(dat_allsp_WBRTLG_RBBH_SB_asex, dat_allsp_WBRTLG_RBBH_SB_asex$sexbias2 == "female_biased")

str(dat_allsp_WBRTLG_RBBH_SB_asex_MB )
str(dat_allsp_WBRTLG_RBBH_SB_asex_FB )

m1_MB <- glm(dat_allsp_WBRTLG_RBBH_SB_asex_MB$log2FC_SA ~ dat_allsp_WBRTLG_RBBH_SB_asex_MB$div * dat_allsp_WBRTLG_RBBH_SB_asex_MB$tiss)
m1_FB <- glm(dat_allsp_WBRTLG_RBBH_SB_asex_FB$log2FC_SA ~ dat_allsp_WBRTLG_RBBH_SB_asex_FB$div * dat_allsp_WBRTLG_RBBH_SB_asex_FB$tiss)

m1_FB_out <- drop1(m1_FB,~.,test="F") 
m1_MB_out <- drop1(m1_MB,~.,test="F") 


########

### overall test
## permutes all vals
perm_all_glm <- function(df){
	
	df$tiss   <-   as.character(df$tiss)
	
	dd    <- as.data.frame(cbind(
		sample(as.character(df$log2FC_SA)),
		df$tiss,
		as.character(df$div)
		))
	
	colnames(dd) <- c("log2FC_SA", "tiss", "div")
	
	dd$log2FC_SA <- as.numeric(as.character(dd$log2FC_SA))
	dd$div       <- as.numeric(as.character(dd$div))
			
	model_1 <- glm(dd$log2FC_SA ~ dd$div * dd$tiss)
	out <- drop1(model_1,~.,test="F") 
	
	F_div = out$F[2]
	F_tis = out$F[3]
	F_int = out$F[4]


	P_div = out$P[2]
	P_tis = out$P[3]
	P_int = out$P[4]

	out_vals <- c(F_div, F_tis, F_int, P_div, P_tis, P_int)
	return(out_vals)
		

}


############## ALSO test main effects more specifically
### permutes div within tisses type
perm_div_glm <- function(df){
	
	df_WB <- subset(df, df$tiss == "WB")
	df_RT <- subset(df, df$tiss == "RT")
	df_LG <- subset(df, df$tiss == "LG")	
	
	dd    <- as.data.frame(cbind(
		as.character(df$log2FC_SA),
		as.character(df$tiss),
		c(sample(df_WB$div), sample(df_RT$div), sample(df_LG$div))
		
		))
	
	colnames(dd) <- c("log2FC_SA", "tiss", "div")
	
	dd$log2FC_SA <- as.numeric(as.character(dd$log2FC_SA))
	dd$div       <- as.numeric(as.character(dd$div))
		
	model_1 <- glm(dd$log2FC_SA ~ dd$div * dd$tiss)
	out <- drop1(model_1,~.,test="F") 
	
	F_div = out$F[2]
	F_tis = out$F[3]
	F_int = out$F[4]


	P_div = out$P[2]
	P_tis = out$P[3]
	P_int = out$P[4]
		
	out_vals <- c(F_div, F_tis, F_int, P_div, P_tis, P_int)
	return(out_vals)
	
}

### permutes tissue within div classes
perm_tis_glm <- function(df){
	
	df$tiss   <-   as.character(df$tiss)
	
	df_1 <- subset(df, df$div == 0.791674)
	df_2 <- subset(df, df$div == 1.10925)
	df_3 <- subset(df, df$div == 1.40629)
	df_4 <- subset(df, df$div == 1.86055)
	df_5 <- subset(df, df$div == 1.22685)	
	
	
	
	dd    <- as.data.frame(cbind(
		as.character(df$log2FC_SA),
		c(sample(df_1$tiss), sample(df_2$tiss), sample(df_3$tiss), sample(df_4$tiss), sample(df_5$tiss)),
		as.character(df$div)

		))
	
	colnames(dd) <- c("log2FC_SA", "tiss", "div")
	
	dd$log2FC_SA <- as.numeric(as.character(dd$log2FC_SA))
	dd$div       <- as.numeric(as.character(dd$div))
	#dd$tiss   <-   as.character(dd$tiss)
	
	print(str(dd))
			
	model_1 <- glm(dd$log2FC_SA ~ dd$div * dd$tiss)
	out <- drop1(model_1,~.,test="F") 
	
	F_div = out$F[2]
	F_tis = out$F[3]
	F_int = out$F[4]


	P_div = out$P[2]
	P_tis = out$P[3]
	P_int = out$P[4]

	out_vals <- c(F_div, F_tis, F_int, P_div, P_tis, P_int)
	return(out_vals)
		

}


#### run perms

all_perm_out_MB <- data.frame()
all_perm_out_FB <- data.frame()
div_perm_out_MB <- data.frame()
div_perm_out_FB <- data.frame()
tis_perm_out_MB <- data.frame()
tis_perm_out_FB <- data.frame()


N_perm = 10000

for(i in seq(1:N_perm)){
	run_all_MB <- perm_all_glm(dat_allsp_WBRTLG_RBBH_SB_asex_MB)	
	run_all_FB <- perm_all_glm(dat_allsp_WBRTLG_RBBH_SB_asex_FB)	
	run_div_MB <- perm_div_glm(dat_allsp_WBRTLG_RBBH_SB_asex_MB)	
	run_div_FB <- perm_div_glm(dat_allsp_WBRTLG_RBBH_SB_asex_FB)	
	run_tis_MB <- perm_div_glm(dat_allsp_WBRTLG_RBBH_SB_asex_MB)	
	run_tis_FB <- perm_div_glm(dat_allsp_WBRTLG_RBBH_SB_asex_FB)	

	all_perm_out_MB <- rbind(all_perm_out_MB, run_all_MB)	
	all_perm_out_FB <- rbind(all_perm_out_FB, run_all_FB)	
	div_perm_out_MB <- rbind(div_perm_out_MB, run_div_MB)	
	div_perm_out_FB <- rbind(div_perm_out_FB, run_div_FB)	
	tis_perm_out_MB <- rbind(tis_perm_out_MB, run_tis_MB)	
	tis_perm_out_FB <- rbind(tis_perm_out_FB, run_tis_FB)	
}


colnames(all_perm_out_MB) <- c("F_div", "F_tis", "F_int", "P_div", "P_tis", "P_int")
colnames(all_perm_out_FB) <- c("F_div", "F_tis", "F_int", "P_div", "P_tis", "P_int")
colnames(div_perm_out_MB) <- c("F_div", "F_tis", "F_int", "P_div", "P_tis", "P_int")
colnames(div_perm_out_FB) <- c("F_div", "F_tis", "F_int", "P_div", "P_tis", "P_int")
colnames(tis_perm_out_MB) <- c("F_div", "F_tis", "F_int", "P_div", "P_tis", "P_int")
colnames(tis_perm_out_FB) <- c("F_div", "F_tis", "F_int", "P_div", "P_tis", "P_int")



#### write perms out to csv

write.csv(all_perm_out_MB, file=paste("all_perm_out_MB_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)
write.csv(all_perm_out_FB, file=paste("all_perm_out_FB_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)

write.csv(div_perm_out_MB, file=paste("div_perm_out_MB_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)
write.csv(div_perm_out_FB, file=paste("div_perm_out_FB_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)

write.csv(tis_perm_out_MB, file=paste("tis_perm_out_MB_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)
write.csv(tis_perm_out_FB, file=paste("tis_perm_out_FB_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)




### what is the chance of obs my values by chance?


str(m1_MB_out)
m1_MB_out

MB_div_m1 <- m1_MB_out$F[2] 
MB_tis_m1 <- m1_MB_out$F[3] 
MB_int_m1 <- m1_MB_out$F[4] 

FB_div_m1 <- m1_FB_out$F[2] 
FB_tis_m1 <- m1_FB_out$F[3] 
FB_int_m1 <- m1_FB_out$F[4] 



head(all_perm_out_MB)
hist(all_perm_out_MB$F_div)

get_P_val <- function(perm_vector, orig_TS){
	v1 <- ifelse(perm_vector > orig_TS , perm_vector, NA)
	v2 <- v1[!is.na(v1)]
	N_over = length(v2)
	P = N_over / length(perm_vector)
	print(N_over)
	return(P)	
}


perm_output_df <- as.data.frame(
rbind(
c(
get_P_val(all_perm_out_MB$F_div, MB_div_m1),
get_P_val(all_perm_out_FB$F_div, FB_div_m1),
get_P_val(all_perm_out_MB$F_tis, MB_tis_m1),
get_P_val(all_perm_out_FB$F_tis, FB_tis_m1),
get_P_val(all_perm_out_MB$F_int, MB_int_m1),
get_P_val(all_perm_out_FB$F_int, FB_int_m1)),


c(
get_P_val(div_perm_out_MB$F_div, MB_div_m1),
get_P_val(div_perm_out_FB$F_div, FB_div_m1),
get_P_val(div_perm_out_MB$F_tis, MB_tis_m1),
get_P_val(div_perm_out_FB$F_tis, FB_tis_m1),
get_P_val(div_perm_out_MB$F_int, MB_int_m1),
get_P_val(div_perm_out_FB$F_int, FB_int_m1)),

c(
get_P_val(tis_perm_out_MB$F_div, MB_div_m1),
get_P_val(tis_perm_out_FB$F_div, FB_div_m1),
get_P_val(tis_perm_out_MB$F_tis, MB_tis_m1),
get_P_val(tis_perm_out_FB$F_tis, FB_tis_m1),
get_P_val(tis_perm_out_MB$F_int, MB_int_m1),
get_P_val(tis_perm_out_FB$F_int, FB_int_m1))

))

rownames(perm_output_df) <- c("all_perm", "div_perm", "tis_perm")
colnames(perm_output_df) <- c("MB_div", "FB_div", "MB_tis", "FB_tis",  "MB_int", "FB_int")  
perm_output_df

# output

write.csv(perm_output_df, file=paste("perm_output_df_nperm=", N_perm, ".csv", sep = ""), row.names=TRUE)

### plots

library(plotrix)

head(dat_Tbi_WB_RBBH_SB_asex)

### means and ST error


divs = c(0.791674, 1.10925, 1.40629, 1.86055, 1.22685)
sp   = c("Tbi", "Tce", "Tcm", "Tpa", "Tps")
ylim_1 = -1.7
ylim_2 = 1.7

FB_WB_means_df <- as.data.frame(cbind(
c(
mean(subset(dat_Tbi_WB_RBBH_SB_asex, dat_Tbi_WB_RBBH_SB_asex$Tbi_WB_sexbias2 == "female_biased")$Tbi_WB_log2FC_SA, na.rm=T),
mean(subset(dat_Tce_WB_RBBH_SB_asex, dat_Tce_WB_RBBH_SB_asex$Tce_WB_sexbias2 == "female_biased")$Tce_WB_log2FC_SA, na.rm=T),
mean(subset(dat_Tcm_WB_RBBH_SB_asex, dat_Tcm_WB_RBBH_SB_asex$Tcm_WB_sexbias2 == "female_biased")$Tcm_WB_log2FC_SA, na.rm=T),
mean(subset(dat_Tpa_WB_RBBH_SB_asex, dat_Tpa_WB_RBBH_SB_asex$Tpa_WB_sexbias2 == "female_biased")$Tpa_WB_log2FC_SA, na.rm=T),
mean(subset(dat_Tps_WB_RBBH_SB_asex, dat_Tps_WB_RBBH_SB_asex$Tps_WB_sexbias2 == "female_biased")$Tps_WB_log2FC_SA, na.rm=T)
),
c(
std.error(subset(dat_Tbi_WB_RBBH_SB_asex, dat_Tbi_WB_RBBH_SB_asex$Tbi_WB_sexbias2 == "female_biased")$Tbi_WB_log2FC_SA, na.rm=T),
std.error(subset(dat_Tce_WB_RBBH_SB_asex, dat_Tce_WB_RBBH_SB_asex$Tce_WB_sexbias2 == "female_biased")$Tce_WB_log2FC_SA, na.rm=T),
std.error(subset(dat_Tcm_WB_RBBH_SB_asex, dat_Tcm_WB_RBBH_SB_asex$Tcm_WB_sexbias2 == "female_biased")$Tcm_WB_log2FC_SA, na.rm=T),
std.error(subset(dat_Tpa_WB_RBBH_SB_asex, dat_Tpa_WB_RBBH_SB_asex$Tpa_WB_sexbias2 == "female_biased")$Tpa_WB_log2FC_SA, na.rm=T),
std.error(subset(dat_Tps_WB_RBBH_SB_asex, dat_Tps_WB_RBBH_SB_asex$Tps_WB_sexbias2 == "female_biased")$Tps_WB_log2FC_SA, na.rm=T)
),
divs, sp
))

colnames(FB_WB_means_df) <- c("mean_change", "STerror", "div" ,"sp")

FB_WB_means_df$mean_change <- as.numeric(as.character(FB_WB_means_df$mean_change ))
FB_WB_means_df$STerror     <- as.numeric(as.character(FB_WB_means_df$STerror))
FB_WB_means_df$div         <- as.numeric(as.character(FB_WB_means_df$div ))

FB_WB_means_df$lower <- FB_WB_means_df$mean_change - FB_WB_means_df$STerror
FB_WB_means_df$upper <- FB_WB_means_df$mean_change + FB_WB_means_df$STerror

FB_WB_means_P1 <- ggplot(data = FB_WB_means_df,aes(x = div,y = mean_change, col = "dd")) + 
    theme_bw() +
    geom_point() + 
    geom_errorbar(aes(ymin = lower,ymax = upper)) + ylim(ylim_1,ylim_2) + 
    scale_color_manual(values=c("darkred")) + 
    xlab("Jukes–Cantor corrected divergence") +
	ylab("Mean log2 fold change in asexual females") + theme(legend.position="none",plot.title = element_text(hjust = 0.5)) +     
	ggtitle("Whole-body")  



FB_RT_means_df <- as.data.frame(cbind(
c(
mean(subset(dat_Tbi_RT_RBBH_SB_asex, dat_Tbi_RT_RBBH_SB_asex$Tbi_RT_sexbias2 == "female_biased")$Tbi_RT_log2FC_SA, na.rm=T),
mean(subset(dat_Tce_RT_RBBH_SB_asex, dat_Tce_RT_RBBH_SB_asex$Tce_RT_sexbias2 == "female_biased")$Tce_RT_log2FC_SA, na.rm=T),
mean(subset(dat_Tcm_RT_RBBH_SB_asex, dat_Tcm_RT_RBBH_SB_asex$Tcm_RT_sexbias2 == "female_biased")$Tcm_RT_log2FC_SA, na.rm=T),
mean(subset(dat_Tpa_RT_RBBH_SB_asex, dat_Tpa_RT_RBBH_SB_asex$Tpa_RT_sexbias2 == "female_biased")$Tpa_RT_log2FC_SA, na.rm=T),
mean(subset(dat_Tps_RT_RBBH_SB_asex, dat_Tps_RT_RBBH_SB_asex$Tps_RT_sexbias2 == "female_biased")$Tps_RT_log2FC_SA, na.rm=T)
),
c(
std.error(subset(dat_Tbi_RT_RBBH_SB_asex, dat_Tbi_RT_RBBH_SB_asex$Tbi_RT_sexbias2 == "female_biased")$Tbi_RT_log2FC_SA, na.rm=T),
std.error(subset(dat_Tce_RT_RBBH_SB_asex, dat_Tce_RT_RBBH_SB_asex$Tce_RT_sexbias2 == "female_biased")$Tce_RT_log2FC_SA, na.rm=T),
std.error(subset(dat_Tcm_RT_RBBH_SB_asex, dat_Tcm_RT_RBBH_SB_asex$Tcm_RT_sexbias2 == "female_biased")$Tcm_RT_log2FC_SA, na.rm=T),
std.error(subset(dat_Tpa_RT_RBBH_SB_asex, dat_Tpa_RT_RBBH_SB_asex$Tpa_RT_sexbias2 == "female_biased")$Tpa_RT_log2FC_SA, na.rm=T),
std.error(subset(dat_Tps_RT_RBBH_SB_asex, dat_Tps_RT_RBBH_SB_asex$Tps_RT_sexbias2 == "female_biased")$Tps_RT_log2FC_SA, na.rm=T)
),
divs, sp
))

colnames(FB_RT_means_df) <- c("mean_change", "STerror", "div" ,"sp")

FB_RT_means_df$mean_change <- as.numeric(as.character(FB_RT_means_df$mean_change ))
FB_RT_means_df$STerror     <- as.numeric(as.character(FB_RT_means_df$STerror))
FB_RT_means_df$div         <- as.numeric(as.character(FB_RT_means_df$div ))

FB_RT_means_df$lower <- FB_RT_means_df$mean_change - FB_RT_means_df$STerror
FB_RT_means_df$upper <- FB_RT_means_df$mean_change + FB_RT_means_df$STerror

FB_RT_means_P1 <- ggplot(data = FB_RT_means_df,aes(x = div,y = mean_change, col = "dd")) + 
    theme_bw() +
    geom_point() + 
    geom_errorbar(aes(ymin = lower,ymax = upper)) + ylim(ylim_1,ylim_2) + 
    scale_color_manual(values=c("darkred")) + 
    xlab("Jukes–Cantor corrected divergence") +
	ylab("Mean log2 fold change in asexual females") + theme(legend.position="none",plot.title = element_text(hjust = 0.5)) +     
	ggtitle("Reproductive tract")  


FB_LG_means_df <- as.data.frame(cbind(
c(
mean(subset(dat_Tbi_LG_RBBH_SB_asex, dat_Tbi_LG_RBBH_SB_asex$Tbi_LG_sexbias2 == "female_biased")$Tbi_LG_log2FC_SA, na.rm=T),
mean(subset(dat_Tce_LG_RBBH_SB_asex, dat_Tce_LG_RBBH_SB_asex$Tce_LG_sexbias2 == "female_biased")$Tce_LG_log2FC_SA, na.rm=T),
mean(subset(dat_Tcm_LG_RBBH_SB_asex, dat_Tcm_LG_RBBH_SB_asex$Tcm_LG_sexbias2 == "female_biased")$Tcm_LG_log2FC_SA, na.rm=T),
mean(subset(dat_Tpa_LG_RBBH_SB_asex, dat_Tpa_LG_RBBH_SB_asex$Tpa_LG_sexbias2 == "female_biased")$Tpa_LG_log2FC_SA, na.rm=T),
mean(subset(dat_Tps_LG_RBBH_SB_asex, dat_Tps_LG_RBBH_SB_asex$Tps_LG_sexbias2 == "female_biased")$Tps_LG_log2FC_SA, na.rm=T)
),
c(
std.error(subset(dat_Tbi_LG_RBBH_SB_asex, dat_Tbi_LG_RBBH_SB_asex$Tbi_LG_sexbias2 == "female_biased")$Tbi_LG_log2FC_SA, na.rm=T),
std.error(subset(dat_Tce_LG_RBBH_SB_asex, dat_Tce_LG_RBBH_SB_asex$Tce_LG_sexbias2 == "female_biased")$Tce_LG_log2FC_SA, na.rm=T),
std.error(subset(dat_Tcm_LG_RBBH_SB_asex, dat_Tcm_LG_RBBH_SB_asex$Tcm_LG_sexbias2 == "female_biased")$Tcm_LG_log2FC_SA, na.rm=T),
std.error(subset(dat_Tpa_LG_RBBH_SB_asex, dat_Tpa_LG_RBBH_SB_asex$Tpa_LG_sexbias2 == "female_biased")$Tpa_LG_log2FC_SA, na.rm=T),
std.error(subset(dat_Tps_LG_RBBH_SB_asex, dat_Tps_LG_RBBH_SB_asex$Tps_LG_sexbias2 == "female_biased")$Tps_LG_log2FC_SA, na.rm=T)
),
divs, sp
))

colnames(FB_LG_means_df) <- c("mean_change", "STerror", "div" ,"sp")

FB_LG_means_df$mean_change <- as.numeric(as.character(FB_LG_means_df$mean_change ))
FB_LG_means_df$STerror     <- as.numeric(as.character(FB_LG_means_df$STerror))
FB_LG_means_df$div         <- as.numeric(as.character(FB_LG_means_df$div ))

FB_LG_means_df$lower <- FB_LG_means_df$mean_change - FB_LG_means_df$STerror
FB_LG_means_df$upper <- FB_LG_means_df$mean_change + FB_LG_means_df$STerror

FB_LG_means_P1 <- ggplot(data = FB_LG_means_df,aes(x = div,y = mean_change, col = "dd")) + 
    theme_bw() +
    geom_point() + 
    geom_errorbar(aes(ymin = lower,ymax = upper)) + ylim(ylim_1,ylim_2) + 
    scale_color_manual(values=c("darkred")) + 
    xlab("Jukes–Cantor corrected divergence") +
	ylab("Mean log2 fold change in asexual females") + theme(legend.position="none",plot.title = element_text(hjust = 0.5)) +     
	ggtitle("Legs")  





MB_WB_means_df <- as.data.frame(cbind(
c(
mean(subset(dat_Tbi_WB_RBBH_SB_asex, dat_Tbi_WB_RBBH_SB_asex$Tbi_WB_sexbias2 == "male_biased")$Tbi_WB_log2FC_SA, na.rm=T),
mean(subset(dat_Tce_WB_RBBH_SB_asex, dat_Tce_WB_RBBH_SB_asex$Tce_WB_sexbias2 == "male_biased")$Tce_WB_log2FC_SA, na.rm=T),
mean(subset(dat_Tcm_WB_RBBH_SB_asex, dat_Tcm_WB_RBBH_SB_asex$Tcm_WB_sexbias2 == "male_biased")$Tcm_WB_log2FC_SA, na.rm=T),
mean(subset(dat_Tpa_WB_RBBH_SB_asex, dat_Tpa_WB_RBBH_SB_asex$Tpa_WB_sexbias2 == "male_biased")$Tpa_WB_log2FC_SA, na.rm=T),
mean(subset(dat_Tps_WB_RBBH_SB_asex, dat_Tps_WB_RBBH_SB_asex$Tps_WB_sexbias2 == "male_biased")$Tps_WB_log2FC_SA, na.rm=T)
),
c(
std.error(subset(dat_Tbi_WB_RBBH_SB_asex, dat_Tbi_WB_RBBH_SB_asex$Tbi_WB_sexbias2 == "male_biased")$Tbi_WB_log2FC_SA, na.rm=T),
std.error(subset(dat_Tce_WB_RBBH_SB_asex, dat_Tce_WB_RBBH_SB_asex$Tce_WB_sexbias2 == "male_biased")$Tce_WB_log2FC_SA, na.rm=T),
std.error(subset(dat_Tcm_WB_RBBH_SB_asex, dat_Tcm_WB_RBBH_SB_asex$Tcm_WB_sexbias2 == "male_biased")$Tcm_WB_log2FC_SA, na.rm=T),
std.error(subset(dat_Tpa_WB_RBBH_SB_asex, dat_Tpa_WB_RBBH_SB_asex$Tpa_WB_sexbias2 == "male_biased")$Tpa_WB_log2FC_SA, na.rm=T),
std.error(subset(dat_Tps_WB_RBBH_SB_asex, dat_Tps_WB_RBBH_SB_asex$Tps_WB_sexbias2 == "male_biased")$Tps_WB_log2FC_SA, na.rm=T)
),
divs, sp
))

colnames(MB_WB_means_df) <- c("mean_change", "STerror", "div" ,"sp")

MB_WB_means_df$mean_change <- as.numeric(as.character(MB_WB_means_df$mean_change ))
MB_WB_means_df$STerror     <- as.numeric(as.character(MB_WB_means_df$STerror))
MB_WB_means_df$div         <- as.numeric(as.character(MB_WB_means_df$div ))

MB_WB_means_df$lower <- MB_WB_means_df$mean_change - MB_WB_means_df$STerror
MB_WB_means_df$upper <- MB_WB_means_df$mean_change + MB_WB_means_df$STerror

MB_WB_means_P1 <- ggplot(data = MB_WB_means_df,aes(x = div,y = mean_change, col = "dd")) + 
    theme_bw() +
    geom_point() + 
    geom_errorbar(aes(ymin = lower,ymax = upper)) + ylim(ylim_1,ylim_2) + 
    scale_color_manual(values=c("royalblue")) + 
    xlab("Jukes–Cantor corrected divergence") +
	ylab("Mean log2 fold change in asexual females") + theme(legend.position="none",plot.title = element_text(hjust = 0.5)) +     
	ggtitle("Whole-body")  



MB_RT_means_df <- as.data.frame(cbind(
c(
mean(subset(dat_Tbi_RT_RBBH_SB_asex, dat_Tbi_RT_RBBH_SB_asex$Tbi_RT_sexbias2 == "male_biased")$Tbi_RT_log2FC_SA, na.rm=T),
mean(subset(dat_Tce_RT_RBBH_SB_asex, dat_Tce_RT_RBBH_SB_asex$Tce_RT_sexbias2 == "male_biased")$Tce_RT_log2FC_SA, na.rm=T),
mean(subset(dat_Tcm_RT_RBBH_SB_asex, dat_Tcm_RT_RBBH_SB_asex$Tcm_RT_sexbias2 == "male_biased")$Tcm_RT_log2FC_SA, na.rm=T),
mean(subset(dat_Tpa_RT_RBBH_SB_asex, dat_Tpa_RT_RBBH_SB_asex$Tpa_RT_sexbias2 == "male_biased")$Tpa_RT_log2FC_SA, na.rm=T),
mean(subset(dat_Tps_RT_RBBH_SB_asex, dat_Tps_RT_RBBH_SB_asex$Tps_RT_sexbias2 == "male_biased")$Tps_RT_log2FC_SA, na.rm=T)
),
c(
std.error(subset(dat_Tbi_RT_RBBH_SB_asex, dat_Tbi_RT_RBBH_SB_asex$Tbi_RT_sexbias2 == "male_biased")$Tbi_RT_log2FC_SA, na.rm=T),
std.error(subset(dat_Tce_RT_RBBH_SB_asex, dat_Tce_RT_RBBH_SB_asex$Tce_RT_sexbias2 == "male_biased")$Tce_RT_log2FC_SA, na.rm=T),
std.error(subset(dat_Tcm_RT_RBBH_SB_asex, dat_Tcm_RT_RBBH_SB_asex$Tcm_RT_sexbias2 == "male_biased")$Tcm_RT_log2FC_SA, na.rm=T),
std.error(subset(dat_Tpa_RT_RBBH_SB_asex, dat_Tpa_RT_RBBH_SB_asex$Tpa_RT_sexbias2 == "male_biased")$Tpa_RT_log2FC_SA, na.rm=T),
std.error(subset(dat_Tps_RT_RBBH_SB_asex, dat_Tps_RT_RBBH_SB_asex$Tps_RT_sexbias2 == "male_biased")$Tps_RT_log2FC_SA, na.rm=T)
),
divs, sp
))

colnames(MB_RT_means_df) <- c("mean_change", "STerror", "div" ,"sp")

MB_RT_means_df$mean_change <- as.numeric(as.character(MB_RT_means_df$mean_change ))
MB_RT_means_df$STerror     <- as.numeric(as.character(MB_RT_means_df$STerror))
MB_RT_means_df$div         <- as.numeric(as.character(MB_RT_means_df$div ))

MB_RT_means_df$lower <- MB_RT_means_df$mean_change - MB_RT_means_df$STerror
MB_RT_means_df$upper <- MB_RT_means_df$mean_change + MB_RT_means_df$STerror

MB_RT_means_P1 <- ggplot(data = MB_RT_means_df,aes(x = div,y = mean_change, col = "dd")) + 
    theme_bw() +
    geom_point() + 
    geom_errorbar(aes(ymin = lower,ymax = upper)) + ylim(ylim_1,ylim_2) + 
    scale_color_manual(values=c("royalblue")) + 
    xlab("Jukes–Cantor corrected divergence") +
	ylab("Mean log2 fold change in asexual females") + theme(legend.position="none",plot.title = element_text(hjust = 0.5)) +     
	ggtitle("Reproductive tract")  


MB_LG_means_df <- as.data.frame(cbind(
c(
mean(subset(dat_Tbi_LG_RBBH_SB_asex, dat_Tbi_LG_RBBH_SB_asex$Tbi_LG_sexbias2 == "male_biased")$Tbi_LG_log2FC_SA, na.rm=T),
mean(subset(dat_Tce_LG_RBBH_SB_asex, dat_Tce_LG_RBBH_SB_asex$Tce_LG_sexbias2 == "male_biased")$Tce_LG_log2FC_SA, na.rm=T),
mean(subset(dat_Tcm_LG_RBBH_SB_asex, dat_Tcm_LG_RBBH_SB_asex$Tcm_LG_sexbias2 == "male_biased")$Tcm_LG_log2FC_SA, na.rm=T),
mean(subset(dat_Tpa_LG_RBBH_SB_asex, dat_Tpa_LG_RBBH_SB_asex$Tpa_LG_sexbias2 == "male_biased")$Tpa_LG_log2FC_SA, na.rm=T),
mean(subset(dat_Tps_LG_RBBH_SB_asex, dat_Tps_LG_RBBH_SB_asex$Tps_LG_sexbias2 == "male_biased")$Tps_LG_log2FC_SA, na.rm=T)
),
c(
std.error(subset(dat_Tbi_LG_RBBH_SB_asex, dat_Tbi_LG_RBBH_SB_asex$Tbi_LG_sexbias2 == "male_biased")$Tbi_LG_log2FC_SA, na.rm=T),
std.error(subset(dat_Tce_LG_RBBH_SB_asex, dat_Tce_LG_RBBH_SB_asex$Tce_LG_sexbias2 == "male_biased")$Tce_LG_log2FC_SA, na.rm=T),
std.error(subset(dat_Tcm_LG_RBBH_SB_asex, dat_Tcm_LG_RBBH_SB_asex$Tcm_LG_sexbias2 == "male_biased")$Tcm_LG_log2FC_SA, na.rm=T),
std.error(subset(dat_Tpa_LG_RBBH_SB_asex, dat_Tpa_LG_RBBH_SB_asex$Tpa_LG_sexbias2 == "male_biased")$Tpa_LG_log2FC_SA, na.rm=T),
std.error(subset(dat_Tps_LG_RBBH_SB_asex, dat_Tps_LG_RBBH_SB_asex$Tps_LG_sexbias2 == "male_biased")$Tps_LG_log2FC_SA, na.rm=T)
),
divs, sp
))

colnames(MB_LG_means_df) <- c("mean_change", "STerror", "div" ,"sp")

MB_LG_means_df$mean_change <- as.numeric(as.character(MB_LG_means_df$mean_change ))
MB_LG_means_df$STerror     <- as.numeric(as.character(MB_LG_means_df$STerror))
MB_LG_means_df$div         <- as.numeric(as.character(MB_LG_means_df$div ))

MB_LG_means_df$lower <- MB_LG_means_df$mean_change - MB_LG_means_df$STerror
MB_LG_means_df$upper <- MB_LG_means_df$mean_change + MB_LG_means_df$STerror

MB_LG_means_P1 <- ggplot(data = MB_LG_means_df,aes(x = div,y = mean_change, col = "dd")) + 
    theme_bw() +
    geom_point() + 
    geom_errorbar(aes(ymin = lower,ymax = upper)) + ylim(ylim_1,ylim_2) + 
    scale_color_manual(values=c("royalblue")) + 
    xlab("Jukes–Cantor corrected divergence") +
	ylab("Mean log2 fold change in asexual females") + theme(legend.position="none",plot.title = element_text(hjust = 0.5)) +     
	ggtitle("Legs")   




pdf("divtime_change_inSB_exp_dotplot.pdf", width = 7, height = 11)
plot_grid(
MB_WB_means_P1, FB_WB_means_P1,
MB_RT_means_P1, FB_RT_means_P1,
MB_LG_means_P1, FB_LG_means_P1, ncol = 2, nrow = 3, label_size=12)
dev.off()
getwd() ## where has my plot gone....?

########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), "divergence_and_change_in_SB.R.txt")


