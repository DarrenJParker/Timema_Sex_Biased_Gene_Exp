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
 # [8] tibble_1.4.2         pkgconfig_2.0.2      rlang_0.3.0.1        bindrcpp_0.2.2       withr_2.1.2          dplyr_0.7.8          tidyselect_0.2.5    
# [15] glue_1.3.0           R6_2.3.0             purrr_0.2.5          lambda.r_1.2.3       magrittr_1.5         scales_1.0.0         codetools_0.2-15    
# [22] assertthat_0.2.0     colorspace_1.3-2     lazyeval_0.2.1       munsell_0.5.0        crayon_1.3.4        
# > 


################################################################################################################################################
##### read in data

setwd("output/DE_10sp_joined_Virgin")

dat_Tbi_WB_10sp_SB_asex = read.csv("Tbi_WB_10sp_disp_allsepar_SB_asex.csv")
dat_Tce_WB_10sp_SB_asex = read.csv("Tce_WB_10sp_disp_allsepar_SB_asex.csv")
dat_Tcm_WB_10sp_SB_asex = read.csv("Tcm_WB_10sp_disp_allsepar_SB_asex.csv")
dat_Tpa_WB_10sp_SB_asex = read.csv("Tpa_WB_10sp_disp_allsepar_SB_asex.csv")
dat_Tps_WB_10sp_SB_asex = read.csv("Tps_WB_10sp_disp_allsepar_SB_asex.csv")



head(dat_Tbi_WB_10sp_SB_asex, n = 20)

#################################################################################################################################################
###### refine sex bias by FC (not log FC)
	
FC_cutoff = 2

dat_Tbi_WB_10sp_SB_asex$Tbi_WB_sexbias2 <- ifelse(2^(sqrt(dat_Tbi_WB_10sp_SB_asex$Tbi_WB_log2FC_SB * dat_Tbi_WB_10sp_SB_asex$Tbi_WB_log2FC_SB)) > FC_cutoff , as.character(dat_Tbi_WB_10sp_SB_asex$Tbi_WB_sexbias), "Unbiased")
dat_Tce_WB_10sp_SB_asex$Tce_WB_sexbias2 <- ifelse(2^(sqrt(dat_Tce_WB_10sp_SB_asex$Tce_WB_log2FC_SB * dat_Tce_WB_10sp_SB_asex$Tce_WB_log2FC_SB)) > FC_cutoff , as.character(dat_Tce_WB_10sp_SB_asex$Tce_WB_sexbias), "Unbiased")
dat_Tcm_WB_10sp_SB_asex$Tcm_WB_sexbias2 <- ifelse(2^(sqrt(dat_Tcm_WB_10sp_SB_asex$Tcm_WB_log2FC_SB * dat_Tcm_WB_10sp_SB_asex$Tcm_WB_log2FC_SB)) > FC_cutoff , as.character(dat_Tcm_WB_10sp_SB_asex$Tcm_WB_sexbias), "Unbiased")
dat_Tpa_WB_10sp_SB_asex$Tpa_WB_sexbias2 <- ifelse(2^(sqrt(dat_Tpa_WB_10sp_SB_asex$Tpa_WB_log2FC_SB * dat_Tpa_WB_10sp_SB_asex$Tpa_WB_log2FC_SB)) > FC_cutoff , as.character(dat_Tpa_WB_10sp_SB_asex$Tpa_WB_sexbias), "Unbiased")
dat_Tps_WB_10sp_SB_asex$Tps_WB_sexbias2 <- ifelse(2^(sqrt(dat_Tps_WB_10sp_SB_asex$Tps_WB_log2FC_SB * dat_Tps_WB_10sp_SB_asex$Tps_WB_log2FC_SB)) > FC_cutoff , as.character(dat_Tps_WB_10sp_SB_asex$Tps_WB_sexbias), "Unbiased")



################################################################################################################################################
############### boxplot SB 
################################################################################################################################################

make_long_table = function(dat_Tbi,dat_Tce,dat_Tcm,dat_Tpa,dat_Tps,tiss){
	
	sp_u = "Tbi"
	df_t1 = as.data.frame(cbind(
		eval(parse(text=paste(dat_Tbi,'$', sp_u, '_' ,tiss,'_log2FC_SA',sep=''))), 
		as.character(eval(parse(text=paste(dat_Tbi,'$', sp_u, '_'  ,tiss,'_sexbias',sep='')))), 
		as.character(eval(parse(text=paste(dat_Tbi,'$', sp_u, '_'  ,tiss,'_sexbias2',sep='')))))) 
	df_t1$sp     = rep(sp_u, length(df_t1[,1]))
	df_t1$group  = paste(df_t1$sp, df_t1$V2, sep = "_")
	df_t1$group2 = paste(df_t1$sp, df_t1$V3, sep = "_")
	
	sp_u = "Tce"
	df_t2 = as.data.frame(cbind(
		eval(parse(text=paste(dat_Tce,'$', sp_u, '_' ,tiss,'_log2FC_SA',sep=''))), 
		as.character(eval(parse(text=paste(dat_Tce,'$', sp_u, '_' ,tiss,'_sexbias',sep='')))), 
		as.character(eval(parse(text=paste(dat_Tce,'$', sp_u, '_' ,tiss,'_sexbias2',sep='')))))) 
	df_t2$sp     = rep(sp_u, length(df_t2[,1]))
	df_t2$group  = paste(df_t2$sp, df_t2$V2, sep = "_")
	df_t2$group2 = paste(df_t2$sp, df_t2$V3, sep = "_")
	
	sp_u = "Tcm"
	df_t3 = as.data.frame(cbind(
		eval(parse(text=paste(dat_Tcm,'$', sp_u, '_' ,tiss,'_log2FC_SA',sep=''))), 
		as.character(eval(parse(text=paste(dat_Tcm,'$', sp_u, '_' ,tiss,'_sexbias',sep='')))), 
		as.character(eval(parse(text=paste(dat_Tcm,'$', sp_u, '_' ,tiss,'_sexbias2',sep='')))))) 
	df_t3$sp     = rep(sp_u, length(df_t3[,1]))
	df_t3$group  = paste(df_t3$sp, df_t3$V2, sep = "_")
	df_t3$group2 = paste(df_t3$sp, df_t3$V3, sep = "_")
	
	sp_u = "Tpa"
	df_t4 = as.data.frame(cbind(
		eval(parse(text=paste(dat_Tpa,'$', sp_u, '_' ,tiss,'_log2FC_SA',sep=''))), 
		as.character(eval(parse(text=paste(dat_Tpa,'$', sp_u, '_' ,tiss,'_sexbias',sep='')))), 
		as.character(eval(parse(text=paste(dat_Tpa,'$', sp_u, '_' ,tiss,'_sexbias2',sep='')))))) 
	df_t4$sp     = rep(sp_u, length(df_t4[,1]))
	df_t4$group  = paste(df_t4$sp, df_t4$V2, sep = "_")
	df_t4$group2 = paste(df_t4$sp, df_t4$V3, sep = "_")
	
	sp_u = "Tps"
	df_t5 = as.data.frame(cbind(
		eval(parse(text=paste(dat_Tps,'$', sp_u, '_' ,tiss,'_log2FC_SA',sep=''))), 
		as.character(eval(parse(text=paste(dat_Tps,'$', sp_u, '_' ,tiss,'_sexbias',sep='')))), 
		as.character(eval(parse(text=paste(dat_Tps,'$', sp_u, '_' ,tiss,'_sexbias2',sep='')))))) 
	df_t5$sp     = rep(sp_u, length(df_t5[,1]))
	df_t5$group  = paste(df_t5$sp, df_t5$V2, sep = "_")
	df_t5$group2 = paste(df_t5$sp, df_t5$V3, sep = "_")

		

	df_t = rbind(df_t1,df_t2,df_t3,df_t4,df_t5)
	colnames(df_t) <- c("log2FC_SA", "sexbias", "sexbias2", "sp", "SBgroup", "SBgroup2")
	
	df_t$log2FC_SA <- as.numeric(as.character(df_t$log2FC_SA))
	df_t$SBgroup <- as.factor(df_t$SBgroup)
	df_t$SBgroup2 <- as.factor(df_t$SBgroup2)
	df_t$SBgroup_ord <- ordered(df_t$SBgroup, levels = c(
	"Tbi_female_biased", "Tbi_male_biased",   "Tbi_Unbiased", 
	"Tce_female_biased", "Tce_male_biased",   "Tce_Unbiased",
	"Tps_female_biased", "Tps_male_biased",  "Tps_Unbiased",
	"Tcm_female_biased", "Tcm_male_biased",   "Tcm_Unbiased",  
	"Tpa_female_biased", "Tpa_male_biased",   "Tpa_Unbiased"
	))
	
	df_t$SBgroup2_ord <- ordered(df_t$SBgroup2, levels = c(
	"Tbi_female_biased", "Tbi_male_biased",   "Tbi_Unbiased", 
	"Tce_female_biased", "Tce_male_biased",   "Tce_Unbiased",
	"Tps_female_biased", "Tps_male_biased",  "Tps_Unbiased",
	"Tcm_female_biased", "Tcm_male_biased",   "Tcm_Unbiased",  
	"Tpa_female_biased", "Tpa_male_biased",   "Tpa_Unbiased"
	))	
	

	df_t$sp_ord <- ordered(df_t$sp, levels = c("Tbi", "Tce", "Tps", "Tcm", "Tpa"	))	

	return(df_t)	
}

dat_ALL_WB_10sp_SB_asex_long <- make_long_table("dat_Tbi_WB_10sp_SB_asex", "dat_Tce_WB_10sp_SB_asex", "dat_Tcm_WB_10sp_SB_asex", "dat_Tpa_WB_10sp_SB_asex", "dat_Tps_WB_10sp_SB_asex",  "WB")


######### remove NAs

dat_ALL_WB_10sp_SB_asex_long_c <- na.omit(dat_ALL_WB_10sp_SB_asex_long)

##### plot boxplots

head(dat_ALL_WB_10sp_SB_asex_long_c)

dat_ALL_WB_10sp_SB_asex_long_c$sexbias2_ord <- ordered(dat_ALL_WB_10sp_SB_asex_long_c$sexbias2, levels = c("female_biased",  "Unbiased", "male_biased"))	

WB_SA_box_s2 <- ggplot(dat_ALL_WB_10sp_SB_asex_long_c, aes(sp_ord, log2FC_SA)) + 
	theme_classic() +
	geom_boxplot(aes(fill = factor(sexbias2_ord)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, lwd = 0.5, fatten = 1, outlier.shape = NA) +
	coord_cartesian(ylim=c(-3.5,3.5)) +
	ylab ("log2 fold change in asexual females") +
	xlab ("Gene-class by species") + 
	scale_fill_manual(values=c("firebrick2", "grey", "royalblue2")) + geom_hline(yintercept = 0) + labs(fill='') +
	ggtitle("Whole-body - Virgin females| SB class = FDR <0.05 FC > 2")
	

## pdf

pdf(paste("Sex_asex_10sp_WB_FC_boxplot_s2_FDR005FC2_Vi.pdf",sep = ""), width = 6, height = 8)
WB_SA_box_s2
dev.off()
getwd() ## wh


########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), "sex_bias_plotsetc_withVIfemales_sessionInfo.txt")




