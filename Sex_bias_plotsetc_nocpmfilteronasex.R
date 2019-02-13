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

setwd("Output/DE_joined_nocpmfilteronasex")

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

dat_ALL_WB_RBBH_SB_asex_long <- make_long_table("dat_Tbi_WB_RBBH_SB_asex", "dat_Tce_WB_RBBH_SB_asex", "dat_Tcm_WB_RBBH_SB_asex", "dat_Tpa_WB_RBBH_SB_asex", "dat_Tps_WB_RBBH_SB_asex",  "WB")
dat_ALL_RT_RBBH_SB_asex_long <- make_long_table("dat_Tbi_RT_RBBH_SB_asex", "dat_Tce_RT_RBBH_SB_asex", "dat_Tcm_RT_RBBH_SB_asex", "dat_Tpa_RT_RBBH_SB_asex", "dat_Tps_RT_RBBH_SB_asex",  "RT")
dat_ALL_LG_RBBH_SB_asex_long <- make_long_table("dat_Tbi_LG_RBBH_SB_asex", "dat_Tce_LG_RBBH_SB_asex", "dat_Tcm_LG_RBBH_SB_asex", "dat_Tpa_LG_RBBH_SB_asex", "dat_Tps_LG_RBBH_SB_asex",  "LG")


#### join tissues together

dat_ALL_WB_RBBH_SB_asex_long$tiss <- rep("WB", length(dat_ALL_WB_RBBH_SB_asex_long[,1]))
dat_ALL_WB_RBBH_SB_asex_long$sp_group_tiss <- as.factor(paste(dat_ALL_WB_RBBH_SB_asex_long$SBgroup, dat_ALL_WB_RBBH_SB_asex_long$tiss, sep = "_"))
dat_ALL_WB_RBBH_SB_asex_long$sp_tiss <- as.factor(paste(dat_ALL_WB_RBBH_SB_asex_long$sp, dat_ALL_WB_RBBH_SB_asex_long$tiss, sep = "_"))

dat_ALL_RT_RBBH_SB_asex_long$tiss <- rep("RT", length(dat_ALL_RT_RBBH_SB_asex_long[,1]))
dat_ALL_RT_RBBH_SB_asex_long$sp_group_tiss <- as.factor(paste(dat_ALL_RT_RBBH_SB_asex_long$SBgroup, dat_ALL_RT_RBBH_SB_asex_long$tiss, sep = "_"))
dat_ALL_RT_RBBH_SB_asex_long$sp_tiss <- as.factor(paste(dat_ALL_RT_RBBH_SB_asex_long$sp, dat_ALL_RT_RBBH_SB_asex_long$tiss, sep = "_"))

dat_ALL_LG_RBBH_SB_asex_long$tiss <- rep("LG", length(dat_ALL_LG_RBBH_SB_asex_long[,1]))
dat_ALL_LG_RBBH_SB_asex_long$sp_group_tiss <- as.factor(paste(dat_ALL_LG_RBBH_SB_asex_long$SBgroup, dat_ALL_LG_RBBH_SB_asex_long$tiss, sep = "_"))
dat_ALL_LG_RBBH_SB_asex_long$sp_tiss <- as.factor(paste(dat_ALL_LG_RBBH_SB_asex_long$sp, dat_ALL_LG_RBBH_SB_asex_long$tiss, sep = "_"))


dat_ALL_WBRTLG_RBBH_SB_asex_long <- rbind(dat_ALL_WB_RBBH_SB_asex_long, dat_ALL_RT_RBBH_SB_asex_long, dat_ALL_LG_RBBH_SB_asex_long)
dat_ALL_WBRTLG_RBBH_SB_asex_long$sp_tiss_ord <- ordered(dat_ALL_WBRTLG_RBBH_SB_asex_long$sp_tiss, levels = c("Tbi_WB", "Tce_WB", "Tps_WB", "Tcm_WB", "Tpa_WB", "Tbi_RT", "Tce_RT", "Tps_RT", "Tcm_RT", "Tpa_RT", "Tbi_LG", "Tce_LG", "Tps_LG", "Tcm_LG", "Tpa_LG"))	
dat_ALL_WBRTLG_RBBH_SB_asex_long$sexbias_ord <- ordered(dat_ALL_WBRTLG_RBBH_SB_asex_long$sexbias, levels = c("female_biased",  "Unbiased", "male_biased"))	
dat_ALL_WBRTLG_RBBH_SB_asex_long$sexbias2_ord <- ordered(dat_ALL_WBRTLG_RBBH_SB_asex_long$sexbias2, levels = c("female_biased",  "Unbiased", "male_biased"))	

str(dat_ALL_WBRTLG_RBBH_SB_asex_long)

######### remove NAs

dat_ALL_WB_RBBH_SB_asex_long_c <- na.omit(dat_ALL_WB_RBBH_SB_asex_long)
dat_ALL_RT_RBBH_SB_asex_long_c <- na.omit(dat_ALL_RT_RBBH_SB_asex_long)
dat_ALL_LG_RBBH_SB_asex_long_c <- na.omit(dat_ALL_LG_RBBH_SB_asex_long)

dat_ALL_WBRTLG_RBBH_SB_asex_long_c <- na.omit(dat_ALL_WBRTLG_RBBH_SB_asex_long)

###### add some empty vals to space the boxplot


dat_ALL_WBRTLG_RBBH_SB_asex_long_c$sp_tiss2 <- as.character(dat_ALL_WBRTLG_RBBH_SB_asex_long_c$sp_tiss)
str(dat_ALL_WBRTLG_RBBH_SB_asex_long_c)

dat_ALL_WBRTLG_RBBH_SB_asex_long_c_2 <- rbind(dat_ALL_WBRTLG_RBBH_SB_asex_long_c, c(NA,"Unbiased","Unbiased",NA,NA,NA,NA,NA,NA,"WB", NA,NA,NA,NA,NA, "skip_WB"))
dat_ALL_WBRTLG_RBBH_SB_asex_long_c_2 <- rbind(dat_ALL_WBRTLG_RBBH_SB_asex_long_c_2, c(NA,"Unbiased","Unbiased",NA,NA,NA,NA,NA,NA,"RT", NA,NA,NA,NA,NA, "skip_RT"))
dat_ALL_WBRTLG_RBBH_SB_asex_long_c_2$log2FC_SA <- as.numeric(as.character(dat_ALL_WBRTLG_RBBH_SB_asex_long_c_2$log2FC_SA ))

dat_ALL_WBRTLG_RBBH_SB_asex_long_c_2$sp_tiss2 <- as.factor(dat_ALL_WBRTLG_RBBH_SB_asex_long_c_2$sp_tiss2)
dat_ALL_WBRTLG_RBBH_SB_asex_long_c_2$sp_tiss2_ord <- ordered(dat_ALL_WBRTLG_RBBH_SB_asex_long_c_2$sp_tiss2, levels = c("Tbi_WB", "Tce_WB", "Tps_WB", "Tcm_WB", "Tpa_WB", "skip_WB", "Tbi_RT", "Tce_RT", "Tps_RT", "Tcm_RT", "Tpa_RT", "skip_RT", "Tbi_LG", "Tce_LG", "Tps_LG", "Tcm_LG", "Tpa_LG"))	

str(dat_ALL_WBRTLG_RBBH_SB_asex_long_c_2)

##### plot boxplot

WBRTLG_SA_box_s2_noout <- ggplot(dat_ALL_WBRTLG_RBBH_SB_asex_long_c_2, aes(sp_tiss2_ord, log2FC_SA)) + 
	theme_classic() +
	geom_boxplot(aes(fill = factor(sexbias2_ord)),position=position_dodge(0.65), width = 0.45, outlier.size = 0, lwd = 0.5, fatten = 1, outlier.shape = NA) +
	coord_cartesian(ylim=c(-3,3)) +
	ylab ("log2 fold change in asexual females") +
	xlab ("Gene-class by species") + 
	scale_fill_manual(values=c("firebrick2", "grey", "royalblue2")) + geom_hline(yintercept = 0) + labs(fill='') +
	ggtitle("SB class = FDR <0.05 FC > 2")


pdf(paste("Sex_asex_RBBH_LG_FC_boxplot_s2_FDR005FC2_nocpmfiltonasex.pdf",sep = ""), width = 8.2, height = 5.5)
WBRTLG_SA_box_s2_noout
dev.off()
getwd() ## wh


#########################################################################################################################
##### wilcox

wilk_test <- function(tissue,dataframe){

	##### subset
	dat1_Tbi_female_biased_s1 <-  subset(dataframe, dataframe$SBgroup == "Tbi_female_biased") 
	dat1_Tbi_male_biased_s1 <-  subset(dataframe, dataframe$SBgroup == "Tbi_male_biased")   
	dat1_Tbi_Unbiased_s1 <-  subset(dataframe, dataframe$SBgroup == "Tbi_Unbiased")      
	dat1_Tce_female_biased_s1 <-  subset(dataframe, dataframe$SBgroup == "Tce_female_biased") 
	dat1_Tce_male_biased_s1 <-  subset(dataframe, dataframe$SBgroup == "Tce_male_biased")   
	dat1_Tce_Unbiased_s1 <-  subset(dataframe, dataframe$SBgroup == "Tce_Unbiased")      
	dat1_Tcm_female_biased_s1 <-  subset(dataframe, dataframe$SBgroup == "Tcm_female_biased") 
	dat1_Tcm_male_biased_s1 <-  subset(dataframe, dataframe$SBgroup == "Tcm_male_biased")   
	dat1_Tcm_Unbiased_s1 <-  subset(dataframe, dataframe$SBgroup == "Tcm_Unbiased")      
	dat1_Tpa_female_biased_s1 <-  subset(dataframe, dataframe$SBgroup == "Tpa_female_biased") 
	dat1_Tpa_male_biased_s1 <-  subset(dataframe, dataframe$SBgroup == "Tpa_male_biased")  
	dat1_Tpa_Unbiased_s1 <-  subset(dataframe, dataframe$SBgroup == "Tpa_Unbiased")      
	dat1_Tps_female_biased_s1 <-  subset(dataframe, dataframe$SBgroup == "Tps_female_biased") 
	dat1_Tps_male_biased_s1 <-  subset(dataframe, dataframe$SBgroup == "Tps_male_biased")   
	dat1_Tps_Unbiased_s1 <-  subset(dataframe, dataframe$SBgroup == "Tps_Unbiased")  
	
	wil_Tbi_female_biased_s1 <- wilcox.test(dat1_Tbi_female_biased_s1$log2FC_SA, dat1_Tbi_Unbiased_s1$log2FC_SA, paired = FALSE)		
	wil_Tce_female_biased_s1 <- wilcox.test(dat1_Tce_female_biased_s1$log2FC_SA, dat1_Tce_Unbiased_s1$log2FC_SA, paired = FALSE)
	wil_Tcm_female_biased_s1 <- wilcox.test(dat1_Tcm_female_biased_s1$log2FC_SA, dat1_Tcm_Unbiased_s1$log2FC_SA, paired = FALSE)
	wil_Tpa_female_biased_s1 <- wilcox.test(dat1_Tpa_female_biased_s1$log2FC_SA, dat1_Tpa_Unbiased_s1$log2FC_SA, paired = FALSE)
	wil_Tps_female_biased_s1 <- wilcox.test(dat1_Tps_female_biased_s1$log2FC_SA, dat1_Tps_Unbiased_s1$log2FC_SA, paired = FALSE)
	
	wil_Tbi_male_biased_s1 <- wilcox.test(dat1_Tbi_male_biased_s1$log2FC_SA, dat1_Tbi_Unbiased_s1$log2FC_SA, paired = FALSE)		
	wil_Tce_male_biased_s1 <- wilcox.test(dat1_Tce_male_biased_s1$log2FC_SA, dat1_Tce_Unbiased_s1$log2FC_SA, paired = FALSE)
	wil_Tcm_male_biased_s1 <- wilcox.test(dat1_Tcm_male_biased_s1$log2FC_SA, dat1_Tcm_Unbiased_s1$log2FC_SA, paired = FALSE)
	wil_Tpa_male_biased_s1 <- wilcox.test(dat1_Tpa_male_biased_s1$log2FC_SA, dat1_Tpa_Unbiased_s1$log2FC_SA, paired = FALSE)
	wil_Tps_male_biased_s1 <- wilcox.test(dat1_Tps_male_biased_s1$log2FC_SA, dat1_Tps_Unbiased_s1$log2FC_SA, paired = FALSE)

	med_Tbi_female_biased_s1 <- median(dat1_Tbi_female_biased_s1$log2FC_SA)
	med_Tce_female_biased_s1 <- median(dat1_Tce_female_biased_s1$log2FC_SA)
	med_Tcm_female_biased_s1 <- median(dat1_Tcm_female_biased_s1$log2FC_SA)
	med_Tpa_female_biased_s1 <- median(dat1_Tpa_female_biased_s1$log2FC_SA)
	med_Tps_female_biased_s1 <- median(dat1_Tps_female_biased_s1$log2FC_SA)
	
	med_Tbi_male_biased_s1 <- median(dat1_Tbi_male_biased_s1$log2FC_SA)	
	med_Tce_male_biased_s1 <- median(dat1_Tce_male_biased_s1$log2FC_SA)
	med_Tcm_male_biased_s1 <- median(dat1_Tcm_male_biased_s1$log2FC_SA)
	med_Tpa_male_biased_s1 <- median(dat1_Tpa_male_biased_s1$log2FC_SA)
	med_Tps_male_biased_s1 <- median(dat1_Tps_male_biased_s1$log2FC_SA)

	med_Tbi_Unbiased_s1 <- median(dat1_Tbi_Unbiased_s1$log2FC_SA)	
	med_Tce_Unbiased_s1 <- median(dat1_Tce_Unbiased_s1$log2FC_SA)
	med_Tcm_Unbiased_s1 <- median(dat1_Tcm_Unbiased_s1$log2FC_SA)
	med_Tpa_Unbiased_s1 <- median(dat1_Tpa_Unbiased_s1$log2FC_SA)
	med_Tps_Unbiased_s1 <- median(dat1_Tps_Unbiased_s1$log2FC_SA)


	FB_N_s1 <- c(length(dat1_Tbi_female_biased_s1[,1]), length(dat1_Tce_female_biased_s1[,1]), length(dat1_Tps_female_biased_s1[,1]), length(dat1_Tcm_female_biased_s1[,1]), length(dat1_Tpa_female_biased_s1[,1]))
	MB_N_s1 <- c(length(dat1_Tbi_male_biased_s1[,1]), length(dat1_Tce_male_biased_s1[,1]), length(dat1_Tps_male_biased_s1[,1]), length(dat1_Tcm_male_biased_s1[,1]), length(dat1_Tpa_male_biased_s1[,1]))
	UB_N_s1 <- c(length(dat1_Tbi_Unbiased_s1[,1]), length(dat1_Tce_Unbiased_s1[,1]), length(dat1_Tps_Unbiased_s1[,1]), length(dat1_Tcm_Unbiased_s1[,1]), length(dat1_Tpa_Unbiased_s1[,1]))


	sps = c("Tbi", "Tce", "Tps", "Tcm", "Tpa")	
	tiss = rep(tissue, 5)
	all_med_female_biased_s1 <- c(med_Tbi_female_biased_s1,med_Tce_female_biased_s1,med_Tps_female_biased_s1,med_Tcm_female_biased_s1,med_Tpa_female_biased_s1)
	all_med_male_biased_s1 <- c(med_Tbi_male_biased_s1,med_Tce_male_biased_s1,med_Tps_male_biased_s1,med_Tcm_male_biased_s1,med_Tpa_male_biased_s1)	
	all_med_Unbiased_s1 <- c(med_Tbi_Unbiased_s1,med_Tce_Unbiased_s1,med_Tps_Unbiased_s1,med_Tcm_Unbiased_s1,med_Tpa_Unbiased_s1)
	
	FB_wil_s1 <- as.vector(c(wil_Tbi_female_biased_s1$statistic,wil_Tce_female_biased_s1$statistic,wil_Tps_female_biased_s1$statistic,wil_Tcm_female_biased_s1$statistic,wil_Tpa_female_biased_s1$statistic))
	MB_wil_s1 <- as.vector(c(wil_Tbi_male_biased_s1$statistic,wil_Tce_male_biased_s1$statistic,wil_Tps_male_biased_s1$statistic,wil_Tcm_male_biased_s1$statistic,wil_Tpa_male_biased_s1$statistic))
		
	FB_wilp_s1 <- as.vector(c(wil_Tbi_female_biased_s1$p.value,wil_Tce_female_biased_s1$p.value,wil_Tps_female_biased_s1$p.value,wil_Tcm_female_biased_s1$p.value,wil_Tpa_female_biased_s1$p.value))
	MB_wilp_s1 <- as.vector(c(wil_Tbi_male_biased_s1$p.value,wil_Tce_male_biased_s1$p.value,wil_Tps_male_biased_s1$p.value,wil_Tcm_male_biased_s1$p.value,wil_Tpa_male_biased_s1$p.value))
	

			
	s1_out_table <- as.data.frame(cbind(sps, tiss, FB_N_s1, MB_N_s1, UB_N_s1, all_med_female_biased_s1, all_med_male_biased_s1, all_med_Unbiased_s1, FB_wil_s1, MB_wil_s1, FB_wilp_s1, MB_wilp_s1) )

	##### subset
	dat1_Tbi_female_biased_s2 <-  subset(dataframe, dataframe$SBgroup2 == "Tbi_female_biased") 
	dat1_Tbi_male_biased_s2 <-  subset(dataframe, dataframe$SBgroup2 == "Tbi_male_biased")   
	dat1_Tbi_Unbiased_s2 <-  subset(dataframe, dataframe$SBgroup2 == "Tbi_Unbiased")      
	dat1_Tce_female_biased_s2 <-  subset(dataframe, dataframe$SBgroup2 == "Tce_female_biased") 
	dat1_Tce_male_biased_s2 <-  subset(dataframe, dataframe$SBgroup2 == "Tce_male_biased")   
	dat1_Tce_Unbiased_s2 <-  subset(dataframe, dataframe$SBgroup2 == "Tce_Unbiased")      
	dat1_Tcm_female_biased_s2 <-  subset(dataframe, dataframe$SBgroup2 == "Tcm_female_biased") 
	dat1_Tcm_male_biased_s2 <-  subset(dataframe, dataframe$SBgroup2 == "Tcm_male_biased")   
	dat1_Tcm_Unbiased_s2 <-  subset(dataframe, dataframe$SBgroup2 == "Tcm_Unbiased")      
	dat1_Tpa_female_biased_s2 <-  subset(dataframe, dataframe$SBgroup2 == "Tpa_female_biased") 
	dat1_Tpa_male_biased_s2 <-  subset(dataframe, dataframe$SBgroup2 == "Tpa_male_biased")  
	dat1_Tpa_Unbiased_s2 <-  subset(dataframe, dataframe$SBgroup2 == "Tpa_Unbiased")      
	dat1_Tps_female_biased_s2 <-  subset(dataframe, dataframe$SBgroup2 == "Tps_female_biased") 
	dat1_Tps_male_biased_s2 <-  subset(dataframe, dataframe$SBgroup2 == "Tps_male_biased")   
	dat1_Tps_Unbiased_s2 <-  subset(dataframe, dataframe$SBgroup2 == "Tps_Unbiased")  
	
	wil_Tbi_female_biased_s2 <- wilcox.test(dat1_Tbi_female_biased_s2$log2FC_SA, dat1_Tbi_Unbiased_s2$log2FC_SA, paired = FALSE)		
	wil_Tce_female_biased_s2 <- wilcox.test(dat1_Tce_female_biased_s2$log2FC_SA, dat1_Tce_Unbiased_s2$log2FC_SA, paired = FALSE)
	wil_Tcm_female_biased_s2 <- wilcox.test(dat1_Tcm_female_biased_s2$log2FC_SA, dat1_Tcm_Unbiased_s2$log2FC_SA, paired = FALSE)
	wil_Tpa_female_biased_s2 <- wilcox.test(dat1_Tpa_female_biased_s2$log2FC_SA, dat1_Tpa_Unbiased_s2$log2FC_SA, paired = FALSE)
	wil_Tps_female_biased_s2 <- wilcox.test(dat1_Tps_female_biased_s2$log2FC_SA, dat1_Tps_Unbiased_s2$log2FC_SA, paired = FALSE)
	
	wil_Tbi_male_biased_s2 <- wilcox.test(dat1_Tbi_male_biased_s2$log2FC_SA, dat1_Tbi_Unbiased_s2$log2FC_SA, paired = FALSE)		
	wil_Tce_male_biased_s2 <- wilcox.test(dat1_Tce_male_biased_s2$log2FC_SA, dat1_Tce_Unbiased_s2$log2FC_SA, paired = FALSE)
	wil_Tcm_male_biased_s2 <- wilcox.test(dat1_Tcm_male_biased_s2$log2FC_SA, dat1_Tcm_Unbiased_s2$log2FC_SA, paired = FALSE)
	wil_Tpa_male_biased_s2 <- wilcox.test(dat1_Tpa_male_biased_s2$log2FC_SA, dat1_Tpa_Unbiased_s2$log2FC_SA, paired = FALSE)
	wil_Tps_male_biased_s2 <- wilcox.test(dat1_Tps_male_biased_s2$log2FC_SA, dat1_Tps_Unbiased_s2$log2FC_SA, paired = FALSE)

	med_Tbi_female_biased_s2 <- median(dat1_Tbi_female_biased_s2$log2FC_SA)
	med_Tce_female_biased_s2 <- median(dat1_Tce_female_biased_s2$log2FC_SA)
	med_Tcm_female_biased_s2 <- median(dat1_Tcm_female_biased_s2$log2FC_SA)
	med_Tpa_female_biased_s2 <- median(dat1_Tpa_female_biased_s2$log2FC_SA)
	med_Tps_female_biased_s2 <- median(dat1_Tps_female_biased_s2$log2FC_SA)
	
	med_Tbi_male_biased_s2 <- median(dat1_Tbi_male_biased_s2$log2FC_SA)	
	med_Tce_male_biased_s2 <- median(dat1_Tce_male_biased_s2$log2FC_SA)
	med_Tcm_male_biased_s2 <- median(dat1_Tcm_male_biased_s2$log2FC_SA)
	med_Tpa_male_biased_s2 <- median(dat1_Tpa_male_biased_s2$log2FC_SA)
	med_Tps_male_biased_s2 <- median(dat1_Tps_male_biased_s2$log2FC_SA)

	med_Tbi_Unbiased_s2 <- median(dat1_Tbi_Unbiased_s2$log2FC_SA)	
	med_Tce_Unbiased_s2 <- median(dat1_Tce_Unbiased_s2$log2FC_SA)
	med_Tcm_Unbiased_s2 <- median(dat1_Tcm_Unbiased_s2$log2FC_SA)
	med_Tpa_Unbiased_s2 <- median(dat1_Tpa_Unbiased_s2$log2FC_SA)
	med_Tps_Unbiased_s2 <- median(dat1_Tps_Unbiased_s2$log2FC_SA)


	FB_N_s2 <- c(length(dat1_Tbi_female_biased_s2[,1]), length(dat1_Tce_female_biased_s2[,1]), length(dat1_Tps_female_biased_s2[,1]), length(dat1_Tcm_female_biased_s2[,1]), length(dat1_Tpa_female_biased_s2[,1]))
	MB_N_s2 <- c(length(dat1_Tbi_male_biased_s2[,1]), length(dat1_Tce_male_biased_s2[,1]), length(dat1_Tps_male_biased_s2[,1]), length(dat1_Tcm_male_biased_s2[,1]), length(dat1_Tpa_male_biased_s2[,1]))
	UB_N_s2 <- c(length(dat1_Tbi_Unbiased_s2[,1]), length(dat1_Tce_Unbiased_s2[,1]), length(dat1_Tps_Unbiased_s2[,1]), length(dat1_Tcm_Unbiased_s2[,1]), length(dat1_Tpa_Unbiased_s2[,1]))

	all_med_female_biased_s2 <- c(med_Tbi_female_biased_s2,med_Tce_female_biased_s2,med_Tps_female_biased_s2,med_Tcm_female_biased_s2,med_Tpa_female_biased_s2)
	all_med_male_biased_s2 <- c(med_Tbi_male_biased_s2,med_Tce_male_biased_s2,med_Tps_male_biased_s2,med_Tcm_male_biased_s2,med_Tpa_male_biased_s2)	
	all_med_Unbiased_s2 <- c(med_Tbi_Unbiased_s2,med_Tce_Unbiased_s2,med_Tps_Unbiased_s2,med_Tcm_Unbiased_s2,med_Tpa_Unbiased_s2)
	
	FB_wil_s2 <- as.vector(c(wil_Tbi_female_biased_s2$statistic,wil_Tce_female_biased_s2$statistic,wil_Tps_female_biased_s2$statistic,wil_Tcm_female_biased_s2$statistic,wil_Tpa_female_biased_s2$statistic))
	MB_wil_s2 <- as.vector(c(wil_Tbi_male_biased_s2$statistic,wil_Tce_male_biased_s2$statistic,wil_Tps_male_biased_s2$statistic,wil_Tcm_male_biased_s2$statistic,wil_Tpa_male_biased_s2$statistic))
		
	FB_wilp_s2 <- as.vector(c(wil_Tbi_female_biased_s2$p.value,wil_Tce_female_biased_s2$p.value,wil_Tps_female_biased_s2$p.value,wil_Tcm_female_biased_s2$p.value,wil_Tpa_female_biased_s2$p.value))
	MB_wilp_s2 <- as.vector(c(wil_Tbi_male_biased_s2$p.value,wil_Tce_male_biased_s2$p.value,wil_Tps_male_biased_s2$p.value,wil_Tcm_male_biased_s2$p.value,wil_Tpa_male_biased_s2$p.value))
	
	print(FB_wil_s2 )
			
	s2_out_table <- as.data.frame(cbind(sps, tiss, FB_N_s2, MB_N_s2, UB_N_s2, all_med_female_biased_s2, all_med_male_biased_s2, all_med_Unbiased_s2, FB_wil_s2, MB_wil_s2, FB_wilp_s2, MB_wilp_s2) )
	

	out_list <- list("s1_out_table" = s1_out_table, "s2_out_table" = s2_out_table)
	return(out_list)

}	
	
wilcox_S1_FDR005 <- as.data.frame(rbind(
wilk_test("WB", dat_ALL_WB_RBBH_SB_asex_long_c)$s1_out_table,	
wilk_test("RT", dat_ALL_RT_RBBH_SB_asex_long_c)$s1_out_table,
wilk_test("LG", dat_ALL_LG_RBBH_SB_asex_long_c)$s1_out_table))		

wilcox_S1_FDR005$FB_wilp_s1 <- as.numeric(as.character(wilcox_S1_FDR005$FB_wilp_s1))
wilcox_S1_FDR005$MB_wilp_s1 <- as.numeric(as.character(wilcox_S1_FDR005$MB_wilp_s1))

wilcox_S1_FDR005$FB_wilFDR_s1 <- p.adjust(wilcox_S1_FDR005$FB_wilp_s1, method = "BH")
wilcox_S1_FDR005$MB_wilFDR_s1 <- p.adjust(wilcox_S1_FDR005$MB_wilp_s1, method = "BH")



wilcox_S2_FDR005_FC2 <- as.data.frame(rbind(
wilk_test("WB", dat_ALL_WB_RBBH_SB_asex_long_c)$s2_out_table,	
wilk_test("RT", dat_ALL_RT_RBBH_SB_asex_long_c)$s2_out_table,	
wilk_test("LG", dat_ALL_LG_RBBH_SB_asex_long_c)$s2_out_table))		

wilcox_S2_FDR005_FC2$FB_wilp_s2 <- as.numeric(as.character(wilcox_S2_FDR005_FC2$FB_wilp_s2))
wilcox_S2_FDR005_FC2$MB_wilp_s2 <- as.numeric(as.character(wilcox_S2_FDR005_FC2$MB_wilp_s2))

wilcox_S2_FDR005_FC2$FB_wilFDR_s1 <- p.adjust(wilcox_S2_FDR005_FC2$FB_wilp_s2, method = "BH")
wilcox_S2_FDR005_FC2$MB_wilFDR_s1 <- p.adjust(wilcox_S2_FDR005_FC2$MB_wilp_s2, method = "BH")




write.csv(wilcox_S1_FDR005, "wilcox_S1_FDR005.csv")
write.csv(wilcox_S2_FDR005_FC2, "wilcox_S2_FDR005_FC2.csv")

########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), "sex_bias_plotsetc_nocpmfiltonasex_sessionInfo.txt")



