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
		as.character(eval(parse(text=paste(dat,'$','N_sp_MB_wFC',sep=''))))	
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
		as.character(eval(parse(text=paste(dat,'$','N_sp_MB_wFC',sep=''))))	
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
		as.character(eval(parse(text=paste(dat,'$','N_sp_MB_wFC',sep=''))))	
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
		as.character(eval(parse(text=paste(dat,'$','N_sp_MB_wFC',sep=''))))	
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
		as.character(eval(parse(text=paste(dat,'$','N_sp_MB_wFC',sep=''))))	
		))
		
		
	df_t5$sp     = rep(sp_u, length(df_t5[,1]))
	df_t5$group  = paste(df_t5$sp, df_t5$V2, sep = "_")
	df_t5$group2 = paste(df_t5$sp, df_t5$V3, sep = "_")
	


	df_t = rbind(df_t1,df_t2,df_t3,df_t4,df_t5)
	colnames(df_t) <- c("log2FC_SA", "sexbias", "sexbias2", "N_sp_FB", "N_sp_MB", "N_sp_FB_wFC", "N_sp_MB_wFC", "sp", "SBgroup", "SBgroup2")
	
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

head(dat_ALL_RT_10sp_SB_asex_long2, n = 10)

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




########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), "N_SB_genes_10sp.R_sessionInfo.txt")
















