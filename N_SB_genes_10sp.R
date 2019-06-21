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
library(fitdistrplus)
library(lme4)

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
 # [1] lme4_1.1-19          Matrix_1.2-15        fitdistrplus_1.0-11  npsurv_0.4-0         lsei_1.2-0           survival_2.43-1      MASS_7.3-51.1        raster_2.8-4         sp_1.3-1             SuperExactTest_1.0.4 VennDiagram_1.6.20  
# [12] futile.logger_1.4.3  pvclust_2.0-0        RColorBrewer_1.1-2   gtable_0.2.0         lattice_0.20-38      gridExtra_2.3        cowplot_0.9.3        pheatmap_1.0.10      ggplot2_3.1.0       

# loaded via a namespace (and not attached):
 # [1] Rcpp_1.0.0           nloptr_1.2.1         pillar_1.3.0         compiler_3.5.1       formatR_1.5          plyr_1.8.4           bindr_0.1.1          tools_3.5.1          futile.options_1.0.1 nlme_3.1-137         tibble_1.4.2        
# [12] pkgconfig_2.0.2      rlang_0.3.0.1        bindrcpp_0.2.2       withr_2.1.2          dplyr_0.7.8          tidyselect_0.2.5     glue_1.3.0           R6_2.3.0             minqa_1.2.4          purrr_0.2.5          lambda.r_1.2.3      
# [23] magrittr_1.5         splines_3.5.1        scales_1.0.0         codetools_0.2-15     assertthat_0.2.0     colorspace_1.3-2     labeling_0.3         lazyeval_0.2.1       munsell_0.5.0        crayon_1.3.4        

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

setwd("Output/DE_joined_10sp")

dat_ALL_WB_10sp_SB_asex2 = read.csv("ALL_WB_10sp_disp_allsepar_wCPM_SB_asex_2.csv")
dat_ALL_RT_10sp_SB_asex2 = read.csv("ALL_RT_10sp_disp_allsepar_wCPM_SB_asex_2.csv")
dat_ALL_LG_10sp_SB_asex2 = read.csv("ALL_LG_10sp_disp_allsepar_wCPM_SB_asex_2.csv")

head(dat_ALL_WB_10sp_SB_asex2, n = 20)

################################################################################################################################################
############### boxplot SB 
################################################################################################################################################

make_long_table_wNsp = function(dat,tiss){
	
	sp_u = "Tbi"
	df_t1 = as.data.frame(cbind(
		eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_log2FC_SA',sep=''))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_'  ,tiss,'_sexbias',sep='')))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_'  ,tiss,'_sexbias2',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_FB',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_MB',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_FB_wFC',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_MB_wFC',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','genename',sep=''))))		
		))
		
		 
	df_t1$sp     = rep(sp_u, length(df_t1[,1]))
	df_t1$group  = paste(df_t1$sp, df_t1$V2, sep = "_")
	df_t1$group2 = paste(df_t1$sp, df_t1$V3, sep = "_")
	df_t1$N_sp_FB
	
	sp_u = "Tce"
	df_t2 = as.data.frame(cbind(
		eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_log2FC_SA',sep=''))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_sexbias',sep='')))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_sexbias2',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_FB',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_MB',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_FB_wFC',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_MB_wFC',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','genename',sep=''))))		
		))
		
	df_t2$sp     = rep(sp_u, length(df_t2[,1]))
	df_t2$group  = paste(df_t2$sp, df_t2$V2, sep = "_")
	df_t2$group2 = paste(df_t2$sp, df_t2$V3, sep = "_")
	
	sp_u = "Tcm"
	df_t3 = as.data.frame(cbind(
		eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_log2FC_SA',sep=''))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_sexbias',sep='')))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_sexbias2',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_FB',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_MB',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_FB_wFC',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_MB_wFC',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','genename',sep=''))))		
		))
		
		
	df_t3$sp     = rep(sp_u, length(df_t3[,1]))
	df_t3$group  = paste(df_t3$sp, df_t3$V2, sep = "_")
	df_t3$group2 = paste(df_t3$sp, df_t3$V3, sep = "_")
	
	sp_u = "Tpa"
	df_t4 = as.data.frame(cbind(
		eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_log2FC_SA',sep=''))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_sexbias',sep='')))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_sexbias2',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_FB',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_MB',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_FB_wFC',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_MB_wFC',sep=''))))	,
		as.character(eval(parse(text=paste(dat,'$','genename',sep=''))))	
		))
		
		
	df_t4$sp     = rep(sp_u, length(df_t4[,1]))
	df_t4$group  = paste(df_t4$sp, df_t4$V2, sep = "_")
	df_t4$group2 = paste(df_t4$sp, df_t4$V3, sep = "_")
	
	sp_u = "Tps"
	df_t5 = as.data.frame(cbind(
		eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_log2FC_SA',sep=''))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_sexbias',sep='')))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_sexbias2',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_FB',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_MB',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_FB_wFC',sep='')))),
		as.character(eval(parse(text=paste(dat,'$','N_sp_MB_wFC',sep=''))))	,
		as.character(eval(parse(text=paste(dat,'$','genename',sep=''))))	
		))
		
		
	df_t5$sp     = rep(sp_u, length(df_t5[,1]))
	df_t5$group  = paste(df_t5$sp, df_t5$V2, sep = "_")
	df_t5$group2 = paste(df_t5$sp, df_t5$V3, sep = "_")
	


	df_t = rbind(df_t1,df_t2,df_t3,df_t4,df_t5)
	colnames(df_t) <- c("log2FC_SA", "sexbias", "sexbias2", "N_sp_FB", "N_sp_MB", "N_sp_FB_wFC", "N_sp_MB_wFC",  "genename", "sp", "SBgroup", "SBgroup2")
	
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

	df_t$sexbias2_Nsp = ifelse(df_t$sexbias2 == "male_biased",   paste(df_t$sexbias2, df_t$N_sp_MB_wFC, sep = "_"), 
                        ifelse(df_t$sexbias2 == "female_biased", paste(df_t$sexbias2, df_t$N_sp_FB_wFC, sep = "_"), paste(df_t$sexbias2, sep = "")))

    df_t$sexbias2_Nsp_ord <-  ordered(df_t$sexbias2_Nsp, levels = c(
    "female_biased_1", "female_biased_2", "female_biased_3", "female_biased_4", "female_biased_5",
    "Unbiased", 
    "male_biased_1", "male_biased_2", "male_biased_3", "male_biased_4", "male_biased_5"
    ))	
	return(df_t)	
}

dat_ALL_WB_10sp_SB_asex_long2 <- make_long_table_wNsp("dat_ALL_WB_10sp_SB_asex2", "WB")
dat_ALL_RT_10sp_SB_asex_long2 <- make_long_table_wNsp("dat_ALL_RT_10sp_SB_asex2", "RT")
dat_ALL_LG_10sp_SB_asex_long2 <- make_long_table_wNsp("dat_ALL_LG_10sp_SB_asex2", "LG")

head(dat_ALL_WB_10sp_SB_asex_long2, n = 10)

######### remove NAs

dat_ALL_WB_10sp_SB_asex_long2_c <- na.omit(dat_ALL_WB_10sp_SB_asex_long2)
dat_ALL_RT_10sp_SB_asex_long2_c <- na.omit(dat_ALL_RT_10sp_SB_asex_long2)
dat_ALL_LG_10sp_SB_asex_long2_c <- na.omit(dat_ALL_LG_10sp_SB_asex_long2)


### add tiss


dat_ALL_WB_10sp_SB_asex_long2_c$tiss  <- rep("WB", length(dat_ALL_WB_10sp_SB_asex_long2_c[,1]))
dat_ALL_RT_10sp_SB_asex_long2_c$tiss  <- rep("RT", length(dat_ALL_RT_10sp_SB_asex_long2_c[,1]))
dat_ALL_LG_10sp_SB_asex_long2_c$tiss  <- rep("LG", length(dat_ALL_LG_10sp_SB_asex_long2_c[,1]))


### join and add some spacing

dat_ALL_WBRTLG_10sp_SB_asex_long2_c <- rbind(dat_ALL_WB_10sp_SB_asex_long2_c, dat_ALL_RT_10sp_SB_asex_long2_c, dat_ALL_LG_10sp_SB_asex_long2_c)
head(dat_ALL_WBRTLG_10sp_SB_asex_long2_c)
dat_ALL_WBRTLG_10sp_SB_asex_long2_c$sexbias2_Nsp_tiss <- paste(dat_ALL_WBRTLG_10sp_SB_asex_long2_c$sexbias2_Nsp, dat_ALL_WBRTLG_10sp_SB_asex_long2_c$tiss, sep = "_")
head(dat_ALL_WBRTLG_10sp_SB_asex_long2_c)

dat_ALL_WBRTLG_10sp_SB_asex_long2_c <- rbind(dat_ALL_WBRTLG_10sp_SB_asex_long2_c, c(NA,"Unbiased","Unbiased",NA,NA,NA,NA,NA,NA,NA, NA,NA,NA,NA,NA,NA,  "skip_WB"))
dat_ALL_WBRTLG_10sp_SB_asex_long2_c <- rbind(dat_ALL_WBRTLG_10sp_SB_asex_long2_c, c(NA,"Unbiased","Unbiased",NA,NA,NA,NA,NA,NA,NA, NA,NA,NA,NA,NA, NA, "skip_RT"))
dat_ALL_WBRTLG_10sp_SB_asex_long2_c <- rbind(dat_ALL_WBRTLG_10sp_SB_asex_long2_c, c(NA,"Unbiased","Unbiased",NA,NA,NA,NA,NA,NA,NA, NA,NA,NA,NA,NA, NA, "male_biased_5_WB"))
dat_ALL_WBRTLG_10sp_SB_asex_long2_c$log2FC_SA <- as.numeric(as.character(dat_ALL_WBRTLG_10sp_SB_asex_long2_c$log2FC_SA ))
dat_ALL_WBRTLG_10sp_SB_asex_long2_c$sexbias2_Nsp_tiss <- as.factor(dat_ALL_WBRTLG_10sp_SB_asex_long2_c$sexbias2_Nsp_tiss)
head(dat_ALL_WBRTLG_10sp_SB_asex_long2_c)
levels(dat_ALL_WBRTLG_10sp_SB_asex_long2_c$sexbias2_Nsp_tiss)

dat_ALL_WBRTLG_10sp_SB_asex_long2_c$sexbias2_Nsp_tiss_ord <- ordered(dat_ALL_WBRTLG_10sp_SB_asex_long2_c$sexbias2_Nsp_tiss, levels = c(
"female_biased_1_WB",
"female_biased_2_WB",
"female_biased_3_WB",
"female_biased_4_WB",
"female_biased_5_WB",
"Unbiased_WB",
"male_biased_1_WB",
"male_biased_2_WB",
"male_biased_3_WB",
"male_biased_4_WB",
"male_biased_5_WB",
"skip_WB",
"female_biased_1_RT",
"female_biased_2_RT",
"female_biased_3_RT",
"female_biased_4_RT",
"female_biased_5_RT",
"Unbiased_RT",
"male_biased_1_RT",
"male_biased_2_RT",
"male_biased_3_RT",
"male_biased_4_RT",
"male_biased_5_RT",
"skip_RT",
"female_biased_1_LG",
"female_biased_2_LG",
"female_biased_3_LG",
"female_biased_4_LG",
"female_biased_5_LG",
"Unbiased_LG",
"male_biased_1_LG",
"male_biased_2_LG",
"male_biased_3_LG",
"male_biased_4_LG",
"male_biased_5_LG"))


str(dat_ALL_WBRTLG_10sp_SB_asex_long2_c)

##### plot boxplots



WBRTLG_SA_box_s2_noout_Nsp <- ggplot(dat_ALL_WBRTLG_10sp_SB_asex_long2_c, aes(sexbias2_Nsp_tiss_ord, log2FC_SA)) + 
	theme_classic() +
	geom_boxplot(aes(fill = factor(sexbias2_Nsp_tiss_ord)),position=position_dodge(0.65), width = 0.45, outlier.size = 0, lwd = 0.5, fatten = 1, outlier.shape = NA) +
	coord_cartesian(ylim=c(-3,3)) +
	ylab ("log2 fold change in asexual females") +
	xlab ("Gene-class by species") + 
	scale_fill_manual(values=c(
	"firebrick2", "firebrick2", "firebrick2", "firebrick2", "firebrick2", "grey", "royalblue2", "royalblue2", "royalblue2", "royalblue2",
	"firebrick2", "firebrick2", "firebrick2", "firebrick2", "firebrick2", "grey", "royalblue2", "royalblue2", "royalblue2", "royalblue2", "royalblue2",
	"firebrick2", "firebrick2", "firebrick2", "firebrick2", "firebrick2", "grey", "royalblue2", "royalblue2", "royalblue2", "royalblue2", "royalblue2"
	)) + geom_hline(yintercept = 0) + labs(fill='') +
	ggtitle("SB class = FDR <0.05 FC > 2")


WBRTLG_SA_box_s2_noout_Nsp2 <- WBRTLG_SA_box_s2_noout_Nsp  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.position="none")


pdf(paste("Sex_asex_10sp_WBRTLG_FC_boxplot_s2_FDR005FC2_noout_Nsp.pdf",sep = ""), width = 5, height = 6.5)
WBRTLG_SA_box_s2_noout_Nsp2
dev.off()
getwd() 



##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
#### test the relationship between the number of species a gene is sex-biased in and the shift in expression in asexual females 

##### subset the data

dat_ALL_WB_10sp_FB_asex_long2_c <- subset(dat_ALL_WB_10sp_SB_asex_long2_c, dat_ALL_WB_10sp_SB_asex_long2_c$sexbias2 == "female_biased")
dat_ALL_WB_10sp_MB_asex_long2_c <- subset(dat_ALL_WB_10sp_SB_asex_long2_c, dat_ALL_WB_10sp_SB_asex_long2_c$sexbias2 == "male_biased")
dat_ALL_RT_10sp_FB_asex_long2_c <- subset(dat_ALL_RT_10sp_SB_asex_long2_c, dat_ALL_RT_10sp_SB_asex_long2_c$sexbias2 == "female_biased")
dat_ALL_RT_10sp_MB_asex_long2_c <- subset(dat_ALL_RT_10sp_SB_asex_long2_c, dat_ALL_RT_10sp_SB_asex_long2_c$sexbias2 == "male_biased")
dat_ALL_LG_10sp_FB_asex_long2_c <- subset(dat_ALL_LG_10sp_SB_asex_long2_c, dat_ALL_LG_10sp_SB_asex_long2_c$sexbias2 == "female_biased")
dat_ALL_LG_10sp_MB_asex_long2_c <- subset(dat_ALL_LG_10sp_SB_asex_long2_c, dat_ALL_LG_10sp_SB_asex_long2_c$sexbias2 == "male_biased")

dat_ALL_WB_10sp_FB_asex_long2_c$N_sp_FB_wFC <- as.numeric(as.character(dat_ALL_WB_10sp_FB_asex_long2_c$N_sp_FB_wFC))
dat_ALL_RT_10sp_FB_asex_long2_c$N_sp_FB_wFC <- as.numeric(as.character(dat_ALL_RT_10sp_FB_asex_long2_c$N_sp_FB_wFC))
dat_ALL_LG_10sp_FB_asex_long2_c$N_sp_FB_wFC <- as.numeric(as.character(dat_ALL_LG_10sp_FB_asex_long2_c$N_sp_FB_wFC))
dat_ALL_WB_10sp_MB_asex_long2_c$N_sp_MB_wFC <- as.numeric(as.character(dat_ALL_WB_10sp_MB_asex_long2_c$N_sp_MB_wFC))
dat_ALL_RT_10sp_MB_asex_long2_c$N_sp_MB_wFC <- as.numeric(as.character(dat_ALL_RT_10sp_MB_asex_long2_c$N_sp_MB_wFC))
dat_ALL_LG_10sp_MB_asex_long2_c$N_sp_MB_wFC <- as.numeric(as.character(dat_ALL_LG_10sp_MB_asex_long2_c$N_sp_MB_wFC))


##### GLMM

fit_WB_FB     <- fitdist(as.numeric(na.omit(dat_ALL_WB_10sp_FB_asex_long2_c$log2FC_SA + 10)), "norm")
fit_WB_FB_LN  <- fitdist(as.numeric(na.omit(dat_ALL_WB_10sp_FB_asex_long2_c$log2FC_SA + 10)), "lnorm")
fit_WB_FB_G   <- fitdist(as.numeric(na.omit(dat_ALL_WB_10sp_FB_asex_long2_c$log2FC_SA + 10)), "gamma")

fit_RT_FB     <- fitdist(as.numeric(na.omit(dat_ALL_RT_10sp_FB_asex_long2_c$log2FC_SA + 10)), "norm")
fit_RT_FB_LN  <- fitdist(as.numeric(na.omit(dat_ALL_RT_10sp_FB_asex_long2_c$log2FC_SA + 10)), "lnorm")
fit_RT_FB_G   <- fitdist(as.numeric(na.omit(dat_ALL_RT_10sp_FB_asex_long2_c$log2FC_SA + 10)), "gamma")

fit_LG_FB     <- fitdist(as.numeric(na.omit(dat_ALL_LG_10sp_FB_asex_long2_c$log2FC_SA + 10)), "norm")
fit_LG_FB_LN  <- fitdist(as.numeric(na.omit(dat_ALL_LG_10sp_FB_asex_long2_c$log2FC_SA + 10)), "lnorm")
fit_LG_FB_G   <- fitdist(as.numeric(na.omit(dat_ALL_LG_10sp_FB_asex_long2_c$log2FC_SA + 10)), "gamma")


fit_WB_MB     <- fitdist(as.numeric(na.omit(dat_ALL_WB_10sp_MB_asex_long2_c$log2FC_SA + 10)), "norm")
fit_WB_MB_LN  <- fitdist(as.numeric(na.omit(dat_ALL_WB_10sp_MB_asex_long2_c$log2FC_SA + 10)), "lnorm")
fit_WB_MB_G   <- fitdist(as.numeric(na.omit(dat_ALL_WB_10sp_MB_asex_long2_c$log2FC_SA + 10)), "gamma")

fit_RT_MB     <- fitdist(as.numeric(na.omit(dat_ALL_RT_10sp_MB_asex_long2_c$log2FC_SA + 10)), "norm")
fit_RT_MB_LN  <- fitdist(as.numeric(na.omit(dat_ALL_RT_10sp_MB_asex_long2_c$log2FC_SA + 10)), "lnorm")
fit_RT_MB_G   <- fitdist(as.numeric(na.omit(dat_ALL_RT_10sp_MB_asex_long2_c$log2FC_SA + 10)), "gamma")

fit_LG_MB     <- fitdist(as.numeric(na.omit(dat_ALL_LG_10sp_MB_asex_long2_c$log2FC_SA + 10)), "norm")
fit_LG_MB_LN  <- fitdist(as.numeric(na.omit(dat_ALL_LG_10sp_MB_asex_long2_c$log2FC_SA + 10)), "lnorm")
fit_LG_MB_G   <- fitdist(as.numeric(na.omit(dat_ALL_LG_10sp_MB_asex_long2_c$log2FC_SA + 10)), "gamma")


### normal fits the best
denscomp(list(fit_WB_FB, fit_WB_FB_LN, fit_WB_FB_G), legendtext = c("norm", "log-norm", "gamma"))
denscomp(list(fit_RT_FB, fit_RT_FB_LN, fit_RT_FB_G), legendtext = c("norm", "log-norm", "gamma"))
denscomp(list(fit_LG_FB, fit_LG_FB_LN, fit_LG_FB_G), legendtext = c("norm", "log-norm", "gamma"))
denscomp(list(fit_WB_MB, fit_WB_MB_LN, fit_WB_MB_G), legendtext = c("norm", "log-norm", "gamma"))
denscomp(list(fit_RT_MB, fit_RT_MB_LN, fit_RT_MB_G), legendtext = c("norm", "log-norm", "gamma"))
denscomp(list(fit_LG_MB, fit_LG_MB_LN, fit_LG_MB_G), legendtext = c("norm", "log-norm", "gamma"))

qqcomp(list(fit_WB_FB, fit_WB_FB_LN, fit_WB_FB_G), legendtext = c("norm", "log-norm", "gamma"))
qqcomp(list(fit_RT_FB, fit_RT_FB_LN, fit_RT_FB_G), legendtext = c("norm", "log-norm", "gamma"))
qqcomp(list(fit_LG_FB, fit_LG_FB_LN, fit_LG_FB_G), legendtext = c("norm", "log-norm", "gamma"))
qqcomp(list(fit_WB_MB, fit_WB_MB_LN, fit_WB_MB_G), legendtext = c("norm", "log-norm", "gamma"))
qqcomp(list(fit_RT_MB, fit_RT_MB_LN, fit_RT_MB_G), legendtext = c("norm", "log-norm", "gamma"))
qqcomp(list(fit_LG_MB, fit_LG_MB_LN, fit_LG_MB_G), legendtext = c("norm", "log-norm", "gamma"))


WB_FB_m1 = lmer(log2FC_SA ~ N_sp_FB_wFC +  (1|genename), data = dat_ALL_WB_10sp_FB_asex_long2_c) 
RT_FB_m1 = lmer(log2FC_SA ~ N_sp_FB_wFC +  (1|genename), data = dat_ALL_RT_10sp_FB_asex_long2_c) 
LG_FB_m1 = lmer(log2FC_SA ~ N_sp_FB_wFC +  (1|genename), data = dat_ALL_LG_10sp_FB_asex_long2_c) 

WB_MB_m1 = lmer(log2FC_SA ~ N_sp_MB_wFC +  (1|genename), data = dat_ALL_WB_10sp_MB_asex_long2_c) 
RT_MB_m1 = lmer(log2FC_SA ~ N_sp_MB_wFC +  (1|genename), data = dat_ALL_RT_10sp_MB_asex_long2_c) 
LG_MB_m1 = lmer(log2FC_SA ~ N_sp_MB_wFC +  (1|genename), data = dat_ALL_LG_10sp_MB_asex_long2_c) 

#### FB all have +ve slope
summary(WB_FB_m1) 
summary(RT_FB_m1)
summary(LG_FB_m1)

#### MB all have -ve slope
summary(WB_MB_m1)
summary(RT_MB_m1)
summary(LG_MB_m1)

### get FDRs

GLMM_FDRs <- p.adjust(c(drop1(WB_FB_m1,test="Chisq")$P[2],
drop1(RT_FB_m1,test="Chisq")$P[2],
drop1(LG_FB_m1,test="Chisq")$P[2],
drop1(WB_MB_m1,test="Chisq")$P[2],
drop1(RT_MB_m1,test="Chisq")$P[2],
drop1(LG_MB_m1,test="Chisq")$P[2]))

max(GLMM_FDRs)

##### some of the GLMMs have a singlar fit - this is because in some models there are very few shared genes making it very difficult to fit the random effect
#### check I get the same result with glms

WB_FB_m0 = lm(log2FC_SA ~ N_sp_FB_wFC, data = dat_ALL_WB_10sp_FB_asex_long2_c) 
RT_FB_m0 = lm(log2FC_SA ~ N_sp_FB_wFC, data = dat_ALL_RT_10sp_FB_asex_long2_c) 
LG_FB_m0 = lm(log2FC_SA ~ N_sp_FB_wFC, data = dat_ALL_LG_10sp_FB_asex_long2_c) 

WB_MB_m0 = lm(log2FC_SA ~ N_sp_MB_wFC, data = dat_ALL_WB_10sp_MB_asex_long2_c) 
RT_MB_m0 = lm(log2FC_SA ~ N_sp_MB_wFC, data = dat_ALL_RT_10sp_MB_asex_long2_c) 
LG_MB_m0 = lm(log2FC_SA ~ N_sp_MB_wFC, data = dat_ALL_LG_10sp_MB_asex_long2_c) 

#### FB all have +ve slope
summary(WB_FB_m0) 
summary(RT_FB_m0)
summary(LG_FB_m0)

#### MB all have -ve slope
summary(WB_MB_m0)
summary(RT_MB_m0)
summary(LG_MB_m0)

### get FDRs

GLM_FDRs <- p.adjust(c(drop1(WB_FB_m0,test="Chisq")$P[2],
drop1(RT_FB_m0,test="Chisq")$P[2],
drop1(LG_FB_m0,test="Chisq")$P[2],
drop1(WB_MB_m0,test="Chisq")$P[2],
drop1(RT_MB_m0,test="Chisq")$P[2],
drop1(LG_MB_m0,test="Chisq")$P[2]))

max(GLM_FDRs)

### output


#### FB all have +ve slope
str(summary(WB_FB_m1) )
summary(RT_FB_m1)
summary(LG_FB_m1)

#### MB all have -ve slope
summary(WB_MB_m1)
summary(RT_MB_m1)
summary(LG_MB_m1)


GLMM_GLM_output <- as.data.frame(cbind(
c("FB", "FB", "FB", "MB", "MB", "MB"),
c("WB", "RT", "LG", "WB", "RT", "LG"),
c(
summary(WB_FB_m1)$coefficients[2],
summary(RT_FB_m1)$coefficients[2],
summary(LG_FB_m1)$coefficients[2],
summary(WB_MB_m1)$coefficients[2],
summary(RT_MB_m1)$coefficients[2],
summary(LG_MB_m1)$coefficients[2]),

GLM_FDRs,
GLMM_FDRs
))

colnames(GLMM_GLM_output) <- c("sex-bias", "tissue-type", "gradient", "GLM_FDR", "GLMM_FDR")
write.csv(GLMM_GLM_output, file="N_SB_GLMM_GLM_output.csv", row.names=F)


#### perm anova (gives the ~same results as the glmm )

## permutes all vals


perm_all_glm_FB <- function(df){
	
	#df$tiss   <-   as.character(df$tiss)
	
	dd    <- as.data.frame(cbind(
		sample(as.character(df$log2FC_SA)),
		df$N_sp_FB_wFC
		))
	
	colnames(dd) <- c("log2FC_SA", "N_sp_FB_wFC")
	
	dd$log2FC_SA     <- as.numeric(as.character(dd$log2FC_SA))
	dd$N_sp_FB_wFC   <- as.numeric(as.character(dd$N_sp_FB_wFC))
	
	#print(str(dd))
			
	model_1 <- glm(dd$log2FC_SA ~ dd$N_sp_FB_wFC)
	out <- drop1(model_1,~.,test="F") 
	
	F_N_SB = out$F[2]
	P_N_SB = out$P[2]


	#print(head(dd, n = 20))	
	#print(out)	
	out_vals <- c(F_N_SB, P_N_SB)
	return(out_vals)
		
}


perm_all_glm_MB <- function(df){
	
	#df$tiss   <-   as.character(df$tiss)
	
	dd    <- as.data.frame(cbind(
		sample(as.character(df$log2FC_SA)),
		df$N_sp_MB_wFC
		))
	
	colnames(dd) <- c("log2FC_SA", "N_sp_MB_wFC")
	
	dd$log2FC_SA     <- as.numeric(as.character(dd$log2FC_SA))
	dd$N_sp_MB_wFC   <- as.numeric(as.character(dd$N_sp_MB_wFC))
	
	#print(str(dd))
			
	model_1 <- glm(dd$log2FC_SA ~ dd$N_sp_MB_wFC)
	out <- drop1(model_1,~.,test="F") 
	
	F_N_SB = out$F[2]
	P_N_SB = out$P[2]


	#print(head(dd, n = 20))	
	#print(out)	
	out_vals <- c(F_N_SB, P_N_SB)
	return(out_vals)
		
}

# perm_all_glm_FB(dat_ALL_WB_10sp_FB_asex_long2_c)
# perm_all_glm_MB(dat_ALL_WB_10sp_MB_asex_long2_c)
# perm_all_glm_MB(dat_ALL_RT_10sp_MB_asex_long2_c)
# perm_all_glm_MB(dat_ALL_RT_10sp_FB_asex_long2_c)
# perm_all_glm_MB(dat_ALL_LG_10sp_MB_asex_long2_c)
# perm_all_glm_MB(dat_ALL_LG_10sp_FB_asex_long2_c)

#### run perms

all_perm_out_FB_WB <- data.frame()
all_perm_out_FB_RT <- data.frame()
all_perm_out_FB_LG <- data.frame()

all_perm_out_MB_WB <- data.frame()
all_perm_out_MB_RT <- data.frame()
all_perm_out_MB_LG <- data.frame()

N_perm = 10000

for(i in seq(1:N_perm)){
	
	run_all_FB_WB <- perm_all_glm_FB(dat_ALL_WB_10sp_FB_asex_long2_c)
	run_all_FB_RT <- perm_all_glm_FB(dat_ALL_RT_10sp_FB_asex_long2_c)
	run_all_FB_LG <- perm_all_glm_FB(dat_ALL_LG_10sp_FB_asex_long2_c)	
	run_all_MB_WB <- perm_all_glm_MB(dat_ALL_WB_10sp_MB_asex_long2_c)
	run_all_MB_RT <- perm_all_glm_MB(dat_ALL_RT_10sp_MB_asex_long2_c)
	run_all_MB_LG <- perm_all_glm_MB(dat_ALL_LG_10sp_MB_asex_long2_c)	

	all_perm_out_FB_WB <- rbind(all_perm_out_FB_WB, run_all_FB_WB)	
	all_perm_out_FB_RT <- rbind(all_perm_out_FB_RT, run_all_FB_RT)	
	all_perm_out_FB_LG <- rbind(all_perm_out_FB_LG, run_all_FB_LG)	
	all_perm_out_MB_WB <- rbind(all_perm_out_MB_WB, run_all_MB_WB)	
	all_perm_out_MB_RT <- rbind(all_perm_out_MB_RT, run_all_MB_RT)	
	all_perm_out_MB_LG <- rbind(all_perm_out_MB_LG, run_all_MB_LG)	

}



colnames(all_perm_out_FB_WB) <- c("F_N_SB","P_N_SB")
colnames(all_perm_out_FB_RT) <- c("F_N_SB","P_N_SB")
colnames(all_perm_out_FB_LG) <- c("F_N_SB","P_N_SB")
colnames(all_perm_out_MB_WB) <- c("F_N_SB","P_N_SB")
colnames(all_perm_out_MB_RT) <- c("F_N_SB","P_N_SB")
colnames(all_perm_out_MB_LG) <- c("F_N_SB","P_N_SB")


#### write perms out to csv

write.csv(all_perm_out_FB_WB, file=paste("all_perm_out_FB_WB_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)
write.csv(all_perm_out_FB_RT, file=paste("all_perm_out_FB_RT_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)
write.csv(all_perm_out_FB_LG, file=paste("all_perm_out_FB_LG_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)
write.csv(all_perm_out_MB_WB, file=paste("all_perm_out_MB_WB_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)
write.csv(all_perm_out_MB_RT, file=paste("all_perm_out_MB_RT_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)
write.csv(all_perm_out_MB_LG, file=paste("all_perm_out_MB_LG_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)



### what is the chance of obs my values by chance?


str(m1_MB_out)
m1_MB_out

m1_WB_FB_F <- m1_WB_FB_out$F[2] 
m1_RT_FB_F <- m1_RT_FB_out$F[2] 
m1_LG_FB_F <- m1_LG_FB_out$F[2] 
m1_WB_MB_F <- m1_WB_MB_out$F[2] 
m1_RT_MB_F <- m1_RT_MB_out$F[2] 
m1_LG_MB_F <- m1_LG_MB_out$F[2] 


get_P_val <- function(perm_vector, orig_TS){
	v1 <- ifelse(perm_vector > orig_TS , perm_vector, NA)
	v2 <- v1[!is.na(v1)]
	N_over = length(v2)
	P = N_over / length(perm_vector)
	print(N_over)
	return(P)	
}

### FDR correct p-values


padj_df <- as.data.frame(cbind(
c("WB_FB_adj_p",
"RT_FB_adj_p",
"LG_FB_adj_p",
"WB_MB_adj_p",
"RT_MB_adj_p",
"LG_MB_adj_p"),
p.adjust(c(
get_P_val(all_perm_out_FB_WB$F_N_SB, m1_WB_FB_F),
get_P_val(all_perm_out_FB_RT$F_N_SB, m1_WB_FB_F),
get_P_val(all_perm_out_FB_LG$F_N_SB, m1_WB_FB_F),
get_P_val(all_perm_out_MB_WB$F_N_SB, m1_WB_MB_F),
get_P_val(all_perm_out_MB_RT$F_N_SB, m1_WB_MB_F),
get_P_val(all_perm_out_MB_LG$F_N_SB, m1_WB_MB_F)))))

# output
colnames(padj_df) <- c("group", "FDR")
write.csv(padj_df, file=paste("N_SB_perm_output_df_nperm=", N_perm, ".csv", sep = ""), row.names=F)




######################################################################################################################################################################## 
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), "N_SB_genes_10sp.R_sessionInfo.txt")



