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

dat_Tce_WB_RBBH_SB_asex = read.csv("Data/With_chr_info/Tce_WB_with_chrinfo.csv")
dat_Tce_RT_RBBH_SB_asex = read.csv("Data/With_chr_info/Tce_RT_with_chrinfo.csv")
dat_Tce_LG_RBBH_SB_asex = read.csv("Data/With_chr_info/Tce_LG_with_chrinfo.csv")


head(dat_Tce_WB_RBBH_SB_asex , n = 20)

#################################################################################################################################################
###### refine sex bias by FC (not log FC)
	
FC_cutoff = 2

dat_Tce_WB_RBBH_SB_asex$Tce_WB_sexbias2 <- ifelse(2^(sqrt(dat_Tce_WB_RBBH_SB_asex$Tce_WB_log2FC_SB * dat_Tce_WB_RBBH_SB_asex$Tce_WB_log2FC_SB)) > FC_cutoff , as.character(dat_Tce_WB_RBBH_SB_asex$Tce_WB_sexbias), "Unbiased")
dat_Tce_RT_RBBH_SB_asex$Tce_RT_sexbias2 <- ifelse(2^(sqrt(dat_Tce_RT_RBBH_SB_asex$Tce_RT_log2FC_SB * dat_Tce_RT_RBBH_SB_asex$Tce_RT_log2FC_SB)) > FC_cutoff , as.character(dat_Tce_RT_RBBH_SB_asex$Tce_RT_sexbias), "Unbiased")
dat_Tce_LG_RBBH_SB_asex$Tce_LG_sexbias2 <- ifelse(2^(sqrt(dat_Tce_LG_RBBH_SB_asex$Tce_LG_log2FC_SB * dat_Tce_LG_RBBH_SB_asex$Tce_LG_log2FC_SB)) > FC_cutoff , as.character(dat_Tce_LG_RBBH_SB_asex$Tce_LG_sexbias), "Unbiased")

head(dat_Tce_WB_RBBH_SB_asex)

#### get autosomal transcripts

dat_Tce_WB_RBBH_SB_asex_A = subset(dat_Tce_WB_RBBH_SB_asex, dat_Tce_WB_RBBH_SB_asex$chr_type == "A")
dat_Tce_RT_RBBH_SB_asex_A = subset(dat_Tce_RT_RBBH_SB_asex, dat_Tce_RT_RBBH_SB_asex$chr_type == "A")
dat_Tce_LG_RBBH_SB_asex_A = subset(dat_Tce_LG_RBBH_SB_asex, dat_Tce_LG_RBBH_SB_asex$chr_type == "A")

################################################################################################################################################
############### boxplot SB 
################################################################################################################################################

make_long_table = function(dat_Tce,tiss){
	
	sp_u = "Tce"
	df_t2 = as.data.frame(cbind(
		eval(parse(text=paste(dat_Tce,'$', sp_u, '_' ,tiss,'_log2FC_SA',sep=''))), 
		as.character(eval(parse(text=paste(dat_Tce,'$', sp_u, '_' ,tiss,'_sexbias',sep='')))), 
		as.character(eval(parse(text=paste(dat_Tce,'$', sp_u, '_' ,tiss,'_sexbias2',sep='')))))) 
	df_t2$sp     = rep(sp_u, length(df_t2[,1]))
	df_t2$group  = paste(df_t2$sp, df_t2$V2, sep = "_")
	df_t2$group2 = paste(df_t2$sp, df_t2$V3, sep = "_")
	

	df_t = rbind(df_t2)
	colnames(df_t) <- c("log2FC_SA", "sexbias", "sexbias2", "sp", "SBgroup", "SBgroup2")
	
	df_t$log2FC_SA <- as.numeric(as.character(df_t$log2FC_SA))
	df_t$SBgroup <- as.factor(df_t$SBgroup)
	df_t$SBgroup2 <- as.factor(df_t$SBgroup2)
	df_t$SBgroup_ord <- ordered(df_t$SBgroup, levels = c(
	"Tce_female_biased", "Tce_male_biased",   "Tce_Unbiased"
	))
	
	df_t$SBgroup2_ord <- ordered(df_t$SBgroup2, levels = c(
	"Tce_female_biased", "Tce_male_biased",   "Tce_Unbiased"
	))	
	

	df_t$sp_ord <- ordered(df_t$sp, levels = c("Tce"))	

	return(df_t)	
}

dat_ALL_WB_RBBH_SB_asex_long_A <- make_long_table("dat_Tce_WB_RBBH_SB_asex_A", "WB")
dat_ALL_RT_RBBH_SB_asex_long_A <- make_long_table("dat_Tce_RT_RBBH_SB_asex_A", "RT")
dat_ALL_LG_RBBH_SB_asex_long_A <- make_long_table("dat_Tce_LG_RBBH_SB_asex_A", "LG")

dat_ALL_WB_RBBH_SB_asex_long_A$tiss <- rep("WB", length(dat_ALL_WB_RBBH_SB_asex_long_A[,1]))
dat_ALL_WB_RBBH_SB_asex_long_A$sp_group_tiss <- as.factor(paste(dat_ALL_WB_RBBH_SB_asex_long_A$SBgroup, dat_ALL_WB_RBBH_SB_asex_long_A$tiss, sep = "_"))
dat_ALL_WB_RBBH_SB_asex_long_A$sp_tiss <- as.factor(paste(dat_ALL_WB_RBBH_SB_asex_long_A$sp, dat_ALL_WB_RBBH_SB_asex_long_A$tiss, sep = "_"))

dat_ALL_RT_RBBH_SB_asex_long_A$tiss <- rep("RT", length(dat_ALL_RT_RBBH_SB_asex_long_A[,1]))
dat_ALL_RT_RBBH_SB_asex_long_A$sp_group_tiss <- as.factor(paste(dat_ALL_RT_RBBH_SB_asex_long_A$SBgroup, dat_ALL_RT_RBBH_SB_asex_long_A$tiss, sep = "_"))
dat_ALL_RT_RBBH_SB_asex_long_A$sp_tiss <- as.factor(paste(dat_ALL_RT_RBBH_SB_asex_long_A$sp, dat_ALL_RT_RBBH_SB_asex_long_A$tiss, sep = "_"))

dat_ALL_LG_RBBH_SB_asex_long_A$tiss <- rep("LG", length(dat_ALL_LG_RBBH_SB_asex_long_A[,1]))
dat_ALL_LG_RBBH_SB_asex_long_A$sp_group_tiss <- as.factor(paste(dat_ALL_LG_RBBH_SB_asex_long_A$SBgroup, dat_ALL_LG_RBBH_SB_asex_long_A$tiss, sep = "_"))
dat_ALL_LG_RBBH_SB_asex_long_A$sp_tiss <- as.factor(paste(dat_ALL_LG_RBBH_SB_asex_long_A$sp, dat_ALL_LG_RBBH_SB_asex_long_A$tiss, sep = "_"))


dat_ALL_WBRTLG_RBBH_SB_asex_long_A <- rbind(dat_ALL_WB_RBBH_SB_asex_long_A, dat_ALL_RT_RBBH_SB_asex_long_A, dat_ALL_LG_RBBH_SB_asex_long_A)
dat_ALL_WBRTLG_RBBH_SB_asex_long_A$sp_tiss_ord <- ordered(dat_ALL_WBRTLG_RBBH_SB_asex_long_A$sp_tiss, levels = c("Tce_WB", "Tce_RT", "Tce_LG"))	
dat_ALL_WBRTLG_RBBH_SB_asex_long_A$sexbias_ord <- ordered(dat_ALL_WBRTLG_RBBH_SB_asex_long_A$sexbias, levels = c("female_biased",  "Unbiased", "male_biased"))	
dat_ALL_WBRTLG_RBBH_SB_asex_long_A$sexbias2_ord <- ordered(dat_ALL_WBRTLG_RBBH_SB_asex_long_A$sexbias2, levels = c("female_biased",  "Unbiased", "male_biased"))	

str(dat_ALL_WBRTLG_RBBH_SB_asex_long_A)

######### remove NAs

dat_ALL_WB_RBBH_SB_asex_long_A_c <- na.omit(dat_ALL_WB_RBBH_SB_asex_long_A)
dat_ALL_RT_RBBH_SB_asex_long_A_c <- na.omit(dat_ALL_RT_RBBH_SB_asex_long_A)
dat_ALL_LG_RBBH_SB_asex_long_A_c <- na.omit(dat_ALL_LG_RBBH_SB_asex_long_A)

dat_ALL_WBRTLG_RBBH_SB_asex_long_A_c <- na.omit(dat_ALL_WBRTLG_RBBH_SB_asex_long_A)

head(dat_ALL_WBRTLG_RBBH_SB_asex_long_A)

##### plot boxplots


WBRTLG_SA_box_s2_noout_A <- ggplot(dat_ALL_WBRTLG_RBBH_SB_asex_long_A_c, aes(sp_tiss_ord, log2FC_SA)) + 
	theme_classic() +
	geom_boxplot(aes(fill = factor(sexbias2_ord)),position=position_dodge(0.65), width = 0.45, outlier.size = 0, lwd = 0.5, fatten = 1, outlier.shape = NA) +
	coord_cartesian(ylim=c(-3,3)) +
	ylab ("log2 fold change in asexual females") +
	xlab ("Gene-class by species") + 
	scale_fill_manual(values=c("firebrick2", "grey", "royalblue2")) + geom_hline(yintercept = 0) + labs(fill='') +
	ggtitle("SB class = FDR <0.05 | FC > 2 | Autosomes")

WBRTLG_SA_box_s2_noout_A2 <- WBRTLG_SA_box_s2_noout_A + theme(axis.text.x = element_text(angle = 90, hjust = 1))


### output

dir.create("Output/withchr")

## pdf

pdf("Output/withchr/boxplot_s2_FDR005FC2_Autosomesonly.pdf", width = 4, height = 8)
WBRTLG_SA_box_s2_noout_A2 
dev.off()
getwd() ## wh


########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), "Output/withchr/sex_bias_plotsetc_withchrinfo_sessionInfo.txt")

