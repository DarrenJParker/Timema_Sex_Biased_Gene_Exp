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


rawdata_10sp <- read.csv("Data/Read_counts/10sp_orth_counts_Arthropoda+Mixed+NOBLASTHIT+cont_K2edge.counts", check.names=FALSE, stringsAsFactors=FALSE)
head(rawdata_10sp)
colnames(rawdata_10sp)
length(colnames(rawdata_10sp))

dir.create("Output")
dir.create("Output/DE_10sp")
setwd("Output/DE_10sp")

###### DE analysis DOING for each tissue seperatly 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 


#### WB all females

y_WBallfemales_10sp <- DGEList(counts=rawdata_10sp[,c(
"Tbi_SF_WB_Vi_Re1","Tbi_SF_WB_Vi_Re2","Tbi_SF_WB_Vi_Re3","Tbi_SF_WB_Md_Re1","Tbi_SF_WB_Md_Re2","Tbi_SF_WB_Md_Re3","Tte_AF_WB_Vi_Re1","Tte_AF_WB_Vi_Re2","Tte_AF_WB_Vi_Re3",
"Tce_SF_WB_Vi_Re1","Tce_SF_WB_Vi_Re2","Tce_SF_WB_Vi_Re3","Tce_SF_WB_Md_Re1","Tce_SF_WB_Md_Re2","Tce_SF_WB_Md_Re3","Tms_AF_WB_Vi_Re1","Tms_AF_WB_Vi_Re2","Tms_AF_WB_Vi_Re3",
"Tcm_SF_WB_Vi_Re1","Tcm_SF_WB_Vi_Re2","Tcm_SF_WB_Vi_Re3","Tcm_SF_WB_Md_Re1","Tcm_SF_WB_Md_Re2","Tcm_SF_WB_Md_Re3","Tsi_AF_WB_Vi_Re1","Tsi_AF_WB_Vi_Re2","Tsi_AF_WB_Vi_Re3",
"Tpa_SF_WB_Vi_Re1","Tpa_SF_WB_Vi_Re2","Tpa_SF_WB_Vi_Re3","Tpa_SF_WB_Md_Re1","Tpa_SF_WB_Md_Re2","Tpa_SF_WB_Md_Re3","Tge_AF_WB_Vi_Re1","Tge_AF_WB_Vi_Re2","Tge_AF_WB_Vi_Re3",
"Tps_SF_WB_Vi_Re1","Tps_SF_WB_Vi_Re2","Tps_SF_WB_Vi_Re3","Tps_SF_WB_Md_Re1","Tps_SF_WB_Md_Re2","Tps_SF_WB_Md_Re3","Tdi_AF_WB_Vi_Re1","Tdi_AF_WB_Vi_Re2","Tdi_AF_WB_Vi_Re3"
)], genes=rawdata_10sp[,1:1])


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### filter and normalisation (code)

filt_and_norm_maj <- function(y,cpm_cut,cut_in_Nsams){
	
	cat("\nNumber number of genes / samples in orig data\n")
	print(dim(y)) ### number of genes / samples
	head(cpm(y)) 
	keep <- 
	  rowSums(cpm(y[,1:3])> cpm_cut) >= cut_in_Nsams & 
      rowSums(cpm(y[,4:6])> cpm_cut) >= cut_in_Nsams &
      rowSums(cpm(y[,7:9])> cpm_cut) >= cut_in_Nsams &    
      rowSums(cpm(y[,10:12])> cpm_cut) >= cut_in_Nsams &     
      rowSums(cpm(y[,13:15])> cpm_cut) >= cut_in_Nsams &    
      rowSums(cpm(y[,16:18])> cpm_cut) >= cut_in_Nsams &      
      rowSums(cpm(y[,19:21])> cpm_cut) >= cut_in_Nsams &    
      rowSums(cpm(y[,22:24])> cpm_cut) >= cut_in_Nsams &     
      rowSums(cpm(y[,25:27])> cpm_cut) >= cut_in_Nsams &
      rowSums(cpm(y[,28:30])> cpm_cut) >= cut_in_Nsams &
      rowSums(cpm(y[,31:33])> cpm_cut) >= cut_in_Nsams &
      rowSums(cpm(y[,34:36])> cpm_cut) >= cut_in_Nsams &
      rowSums(cpm(y[,37:39])> cpm_cut) >= cut_in_Nsams &
      rowSums(cpm(y[,40:42])> cpm_cut) >= cut_in_Nsams &
      rowSums(cpm(y[,43:45])> cpm_cut) >= cut_in_Nsams
            

	y <- y[keep,]

	y$samples$lib.size <- colSums(y$counts) # if filter recalc lib size
	y <- calcNormFactors(y)
	cat("\nLib norm factors\n")
	print(y$samples)	
	cat("\nNumber number of genes / samples after filtering\n")
	print(dim(y))
	
	return(y)
}


########################################################################################################################################################################
####### looking at WB vi females vs md females

y_WBallfemales_10sp_F <- filt_and_norm_maj(y_WBallfemales_10sp,0.5, 2)
cpm_df_y_WBallfemales_10sp_F_log = as.data.frame(cpm(y_WBallfemales_10sp_F, log = T))

Tbi_SF_WB_Vi_mlcpm <- rowMeans(cbind(cpm_df_y_WBallfemales_10sp_F_log$Tbi_SF_WB_Vi_Re1, cpm_df_y_WBallfemales_10sp_F_log$Tbi_SF_WB_Vi_Re2, cpm_df_y_WBallfemales_10sp_F_log$Tbi_SF_WB_Vi_Re3))
Tbi_SF_WB_Md_mlcpm <- rowMeans(cbind(cpm_df_y_WBallfemales_10sp_F_log$Tbi_SF_WB_Md_Re1, cpm_df_y_WBallfemales_10sp_F_log$Tbi_SF_WB_Md_Re2, cpm_df_y_WBallfemales_10sp_F_log$Tbi_SF_WB_Md_Re3))
Tte_AF_WB_Vi_mlcpm <- rowMeans(cbind(cpm_df_y_WBallfemales_10sp_F_log$Tte_AF_WB_Vi_Re1, cpm_df_y_WBallfemales_10sp_F_log$Tte_AF_WB_Vi_Re2, cpm_df_y_WBallfemales_10sp_F_log$Tte_AF_WB_Vi_Re3))
Tce_SF_WB_Vi_mlcpm <- rowMeans(cbind(cpm_df_y_WBallfemales_10sp_F_log$Tce_SF_WB_Vi_Re1, cpm_df_y_WBallfemales_10sp_F_log$Tce_SF_WB_Vi_Re2, cpm_df_y_WBallfemales_10sp_F_log$Tce_SF_WB_Vi_Re3))
Tce_SF_WB_Md_mlcpm <- rowMeans(cbind(cpm_df_y_WBallfemales_10sp_F_log$Tce_SF_WB_Md_Re1, cpm_df_y_WBallfemales_10sp_F_log$Tce_SF_WB_Md_Re2, cpm_df_y_WBallfemales_10sp_F_log$Tce_SF_WB_Md_Re3))
Tms_AF_WB_Vi_mlcpm <- rowMeans(cbind(cpm_df_y_WBallfemales_10sp_F_log$Tms_AF_WB_Vi_Re1, cpm_df_y_WBallfemales_10sp_F_log$Tms_AF_WB_Vi_Re2, cpm_df_y_WBallfemales_10sp_F_log$Tms_AF_WB_Vi_Re3))
Tcm_SF_WB_Vi_mlcpm <- rowMeans(cbind(cpm_df_y_WBallfemales_10sp_F_log$Tcm_SF_WB_Vi_Re1, cpm_df_y_WBallfemales_10sp_F_log$Tcm_SF_WB_Vi_Re2, cpm_df_y_WBallfemales_10sp_F_log$Tcm_SF_WB_Vi_Re3))
Tcm_SF_WB_Md_mlcpm <- rowMeans(cbind(cpm_df_y_WBallfemales_10sp_F_log$Tcm_SF_WB_Md_Re1, cpm_df_y_WBallfemales_10sp_F_log$Tcm_SF_WB_Md_Re2, cpm_df_y_WBallfemales_10sp_F_log$Tcm_SF_WB_Md_Re3))
Tsi_AF_WB_Vi_mlcpm <- rowMeans(cbind(cpm_df_y_WBallfemales_10sp_F_log$Tsi_AF_WB_Vi_Re1, cpm_df_y_WBallfemales_10sp_F_log$Tsi_AF_WB_Vi_Re2, cpm_df_y_WBallfemales_10sp_F_log$Tsi_AF_WB_Vi_Re3))
Tpa_SF_WB_Vi_mlcpm <- rowMeans(cbind(cpm_df_y_WBallfemales_10sp_F_log$Tpa_SF_WB_Vi_Re1, cpm_df_y_WBallfemales_10sp_F_log$Tpa_SF_WB_Vi_Re2, cpm_df_y_WBallfemales_10sp_F_log$Tpa_SF_WB_Vi_Re3))
Tpa_SF_WB_Md_mlcpm <- rowMeans(cbind(cpm_df_y_WBallfemales_10sp_F_log$Tpa_SF_WB_Md_Re1, cpm_df_y_WBallfemales_10sp_F_log$Tpa_SF_WB_Md_Re2, cpm_df_y_WBallfemales_10sp_F_log$Tpa_SF_WB_Md_Re3))
Tge_AF_WB_Vi_mlcpm <- rowMeans(cbind(cpm_df_y_WBallfemales_10sp_F_log$Tge_AF_WB_Vi_Re1, cpm_df_y_WBallfemales_10sp_F_log$Tge_AF_WB_Vi_Re2, cpm_df_y_WBallfemales_10sp_F_log$Tge_AF_WB_Vi_Re3))
Tps_SF_WB_Vi_mlcpm <- rowMeans(cbind(cpm_df_y_WBallfemales_10sp_F_log$Tps_SF_WB_Vi_Re1, cpm_df_y_WBallfemales_10sp_F_log$Tps_SF_WB_Vi_Re2, cpm_df_y_WBallfemales_10sp_F_log$Tps_SF_WB_Vi_Re3))
Tps_SF_WB_Md_mlcpm <- rowMeans(cbind(cpm_df_y_WBallfemales_10sp_F_log$Tps_SF_WB_Md_Re1, cpm_df_y_WBallfemales_10sp_F_log$Tps_SF_WB_Md_Re2, cpm_df_y_WBallfemales_10sp_F_log$Tps_SF_WB_Md_Re3))
Tdi_AF_WB_Vi_mlcpm <- rowMeans(cbind(cpm_df_y_WBallfemales_10sp_F_log$Tdi_AF_WB_Vi_Re1, cpm_df_y_WBallfemales_10sp_F_log$Tdi_AF_WB_Vi_Re2, cpm_df_y_WBallfemales_10sp_F_log$Tdi_AF_WB_Vi_Re3))

mean_cpm_df_y_WBallfemales_10sp_F_log <- as.data.frame(cbind(
Tbi_SF_WB_Vi_mlcpm, Tbi_SF_WB_Md_mlcpm, Tte_AF_WB_Vi_mlcpm,
Tce_SF_WB_Vi_mlcpm, Tce_SF_WB_Md_mlcpm, Tms_AF_WB_Vi_mlcpm,
Tcm_SF_WB_Vi_mlcpm, Tcm_SF_WB_Md_mlcpm, Tsi_AF_WB_Vi_mlcpm,
Tpa_SF_WB_Vi_mlcpm, Tpa_SF_WB_Md_mlcpm, Tge_AF_WB_Vi_mlcpm,
Tps_SF_WB_Vi_mlcpm, Tps_SF_WB_Md_mlcpm, Tdi_AF_WB_Vi_mlcpm
))

rownames(mean_cpm_df_y_WBallfemales_10sp_F_log) <- y_WBallfemales_10sp_F$genes[,1]

head(mean_cpm_df_y_WBallfemales_10sp_F_log)

### plot heatmap
png(filename = "heatmap_Md_vs_Vi_all_females_meanlogcpm.png", width = 6, height = 7, units = "in", pointsize = 12, bg = "white", res = 300)
pheatmap(mean_cpm_df_y_WBallfemales_10sp_F_log,clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete")
dev.off()

pdf("heatmap_Md_vs_Vi_all_females_meanlogcpm.pdf", width = 6, height = 7)
pheatmap(mean_cpm_df_y_WBallfemales_10sp_F_log,clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete")
dev.off()

## support vals

Pv_clust_fit_mean_cpm_df_y_WBallfemales_10sp_F_log <- pvclust(mean_cpm_df_y_WBallfemales_10sp_F_log, method.hclust="ward.D2", method.dist="euclidean", nboot=10000) 
plot(Pv_clust_fit_mean_cpm_df_y_WBallfemales_10sp_F_log)

png(filename = "PPv_clust_fit_mean_cpm_df_y_WBallfemales_10sp_F_log_n=10000.png", width = 5, height = 7, units = "in", pointsize = 12, bg = "white", res = 300)
plot(Pv_clust_fit_mean_cpm_df_y_WBallfemales_10sp_F_log)
dev.off()

########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), "Sex_bias_edgeR_10sp_orths_virgin_mated_asexual_heatmap.R_sessionInfo.txt")






















