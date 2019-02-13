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

# ###### fun

five_sp_DE_venn <- function(Tbi,Tce,Tcm,Tpa,Tps,title){
	venny.plot <- venn.diagram(
	list("Tbi" = Tbi, "Tce" = Tce, "Tcm" = Tcm, "Tpa" = Tpa, "Tps" = Tps ), filename = NULL,
                            fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                            cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                            margin = 0.6, cat.dist = 0.23, main = title, main.pos = c(0.5,0.8), main.cex = 2, main.fontface = "bold", cat.cex = 2)
	return(venny.plot)
}


five_sp_DE_venn_single_col_F <- function(Tbi,Tce,Tcm,Tpa,Tps,title){
	venny.plot <- venn.diagram(
	list("Tbi" = Tbi, "Tce" = Tce, "Tcm" = Tcm, "Tpa" = Tpa, "Tps" = Tps ), filename = NULL,
                            fill = c("firebrick2", "firebrick2", "firebrick2", "firebrick2", "firebrick2"),
                            cat.col = c("firebrick2", "firebrick2", "firebrick2", "firebrick2", "firebrick2"),
                            margin = 0.6, cat.dist = 0.23, main = title, main.pos = c(0.5,0.8), main.cex = 2, main.fontface = "bold", cat.cex = 2, alpha = 0.3)
	return(venny.plot)
}


five_sp_DE_venn_single_col_M <- function(Tbi,Tce,Tcm,Tpa,Tps,title){
	venny.plot <- venn.diagram(
	list("Tbi" = Tbi, "Tce" = Tce, "Tcm" = Tcm, "Tpa" = Tpa, "Tps" = Tps ), filename = NULL,
                            fill = c("royalblue2", "royalblue2", "royalblue2", "royalblue2", "royalblue2"),
                            cat.col = c("royalblue2", "royalblue2", "royalblue2", "royalblue2", "royalblue2"),
                            margin = 0.6, cat.dist = 0.23, main = title, main.pos = c(0.5,0.8), main.cex = 2, main.fontface = "bold", cat.cex = 2, alpha = 0.35)
	return(venny.plot)
}


theme_update(plot.title = element_text(hjust = 0.5))

################################################################################################################################################
##### read in data

setwd("output/DE_joined_10sp")

dat_ALL_WB_10sp_SB_asex = read.csv("ALL_WB_10sp_disp_allsepar_wCPM_SB_asex.csv")
dat_ALL_RT_10sp_SB_asex = read.csv("ALL_RT_10sp_disp_allsepar_wCPM_SB_asex.csv")
dat_ALL_LG_10sp_SB_asex = read.csv("ALL_LG_10sp_disp_allsepar_wCPM_SB_asex.csv")


head(dat_ALL_RT_10sp_SB_asex, n = 20)

#################################################################################################################################################
###### refine sex bias by FC (not log FC)
	
FC_cutoff = 2

dat_ALL_WB_10sp_SB_asex$Tbi_WB_sexbias2 <- ifelse(2^(sqrt(dat_ALL_WB_10sp_SB_asex$Tbi_WB_log2FC_SB * dat_ALL_WB_10sp_SB_asex$Tbi_WB_log2FC_SB)) > FC_cutoff , as.character(dat_ALL_WB_10sp_SB_asex$Tbi_WB_sexbias), "Unbiased")
dat_ALL_WB_10sp_SB_asex$Tce_WB_sexbias2 <- ifelse(2^(sqrt(dat_ALL_WB_10sp_SB_asex$Tce_WB_log2FC_SB * dat_ALL_WB_10sp_SB_asex$Tce_WB_log2FC_SB)) > FC_cutoff , as.character(dat_ALL_WB_10sp_SB_asex$Tce_WB_sexbias), "Unbiased")
dat_ALL_WB_10sp_SB_asex$Tcm_WB_sexbias2 <- ifelse(2^(sqrt(dat_ALL_WB_10sp_SB_asex$Tcm_WB_log2FC_SB * dat_ALL_WB_10sp_SB_asex$Tcm_WB_log2FC_SB)) > FC_cutoff , as.character(dat_ALL_WB_10sp_SB_asex$Tcm_WB_sexbias), "Unbiased")
dat_ALL_WB_10sp_SB_asex$Tpa_WB_sexbias2 <- ifelse(2^(sqrt(dat_ALL_WB_10sp_SB_asex$Tpa_WB_log2FC_SB * dat_ALL_WB_10sp_SB_asex$Tpa_WB_log2FC_SB)) > FC_cutoff , as.character(dat_ALL_WB_10sp_SB_asex$Tpa_WB_sexbias), "Unbiased")
dat_ALL_WB_10sp_SB_asex$Tps_WB_sexbias2 <- ifelse(2^(sqrt(dat_ALL_WB_10sp_SB_asex$Tps_WB_log2FC_SB * dat_ALL_WB_10sp_SB_asex$Tps_WB_log2FC_SB)) > FC_cutoff , as.character(dat_ALL_WB_10sp_SB_asex$Tps_WB_sexbias), "Unbiased")

dat_ALL_RT_10sp_SB_asex$Tbi_RT_sexbias2 <- ifelse(2^(sqrt(dat_ALL_RT_10sp_SB_asex$Tbi_RT_log2FC_SB * dat_ALL_RT_10sp_SB_asex$Tbi_RT_log2FC_SB)) > FC_cutoff , as.character(dat_ALL_RT_10sp_SB_asex$Tbi_RT_sexbias), "Unbiased")
dat_ALL_RT_10sp_SB_asex$Tce_RT_sexbias2 <- ifelse(2^(sqrt(dat_ALL_RT_10sp_SB_asex$Tce_RT_log2FC_SB * dat_ALL_RT_10sp_SB_asex$Tce_RT_log2FC_SB)) > FC_cutoff , as.character(dat_ALL_RT_10sp_SB_asex$Tce_RT_sexbias), "Unbiased")
dat_ALL_RT_10sp_SB_asex$Tcm_RT_sexbias2 <- ifelse(2^(sqrt(dat_ALL_RT_10sp_SB_asex$Tcm_RT_log2FC_SB * dat_ALL_RT_10sp_SB_asex$Tcm_RT_log2FC_SB)) > FC_cutoff , as.character(dat_ALL_RT_10sp_SB_asex$Tcm_RT_sexbias), "Unbiased")
dat_ALL_RT_10sp_SB_asex$Tpa_RT_sexbias2 <- ifelse(2^(sqrt(dat_ALL_RT_10sp_SB_asex$Tpa_RT_log2FC_SB * dat_ALL_RT_10sp_SB_asex$Tpa_RT_log2FC_SB)) > FC_cutoff , as.character(dat_ALL_RT_10sp_SB_asex$Tpa_RT_sexbias), "Unbiased")
dat_ALL_RT_10sp_SB_asex$Tps_RT_sexbias2 <- ifelse(2^(sqrt(dat_ALL_RT_10sp_SB_asex$Tps_RT_log2FC_SB * dat_ALL_RT_10sp_SB_asex$Tps_RT_log2FC_SB)) > FC_cutoff , as.character(dat_ALL_RT_10sp_SB_asex$Tps_RT_sexbias), "Unbiased")

dat_ALL_LG_10sp_SB_asex$Tbi_LG_sexbias2 <- ifelse(2^(sqrt(dat_ALL_LG_10sp_SB_asex$Tbi_LG_log2FC_SB * dat_ALL_LG_10sp_SB_asex$Tbi_LG_log2FC_SB)) > FC_cutoff , as.character(dat_ALL_LG_10sp_SB_asex$Tbi_LG_sexbias), "Unbiased")
dat_ALL_LG_10sp_SB_asex$Tce_LG_sexbias2 <- ifelse(2^(sqrt(dat_ALL_LG_10sp_SB_asex$Tce_LG_log2FC_SB * dat_ALL_LG_10sp_SB_asex$Tce_LG_log2FC_SB)) > FC_cutoff , as.character(dat_ALL_LG_10sp_SB_asex$Tce_LG_sexbias), "Unbiased")
dat_ALL_LG_10sp_SB_asex$Tcm_LG_sexbias2 <- ifelse(2^(sqrt(dat_ALL_LG_10sp_SB_asex$Tcm_LG_log2FC_SB * dat_ALL_LG_10sp_SB_asex$Tcm_LG_log2FC_SB)) > FC_cutoff , as.character(dat_ALL_LG_10sp_SB_asex$Tcm_LG_sexbias), "Unbiased")
dat_ALL_LG_10sp_SB_asex$Tpa_LG_sexbias2 <- ifelse(2^(sqrt(dat_ALL_LG_10sp_SB_asex$Tpa_LG_log2FC_SB * dat_ALL_LG_10sp_SB_asex$Tpa_LG_log2FC_SB)) > FC_cutoff , as.character(dat_ALL_LG_10sp_SB_asex$Tpa_LG_sexbias), "Unbiased")
dat_ALL_LG_10sp_SB_asex$Tps_LG_sexbias2 <- ifelse(2^(sqrt(dat_ALL_LG_10sp_SB_asex$Tps_LG_log2FC_SB * dat_ALL_LG_10sp_SB_asex$Tps_LG_log2FC_SB)) > FC_cutoff , as.character(dat_ALL_LG_10sp_SB_asex$Tps_LG_sexbias), "Unbiased")


##############################################################################################################################################################
###### find how many species a gene is sex-biased in

get_SB_in_any_comp = function(df){

SB_all = c()
SB_dir_all = c()
SB_all_wFC = c()
SB_dir_all_wFC = c()	

N_sp_FB = c()
N_sp_MB = c()
N_sp_SB = c()
N_sp_FB_wFC = c()
N_sp_MB_wFC = c()
N_sp_SB_wFC = c()

for(i in seq(1,length(df[,1]))){

	SB_v   = c(as.character(df[i,4]), 
	            as.character(df[i,9]), 
	            as.character(df[i,14]), 
	            as.character(df[i,19]), 
	            as.character(df[i,24]))
	
	a1 = "Not_SB"
	
	if("female_biased" %in% SB_v){ a1 = "SB"}
	if("male_biased" %in% SB_v){ a1 = "SB"}
	
	b1 = "Not_SB"
	b2 = "Not_SB"
	
	FB_N = sum(SB_v == "female_biased")
	MB_N = sum(SB_v == "male_biased")
	SB_N = FB_N + MB_N
	
	N_sp_FB = c(N_sp_FB, FB_N)
	N_sp_MB = c(N_sp_MB, MB_N)
	N_sp_SB = c(N_sp_SB, SB_N)
			
	if("female_biased" %in% SB_v){ b1 = "female_biased"}
	if("male_biased" %in% SB_v){ b2 = "male_biased"}
	
	b3 = paste(b1, b2, sep = "_")	
	
	SB_all = c(SB_all, a1)
	SB_dir_all = c(SB_dir_all, b3)
	
	
	SB_v_wFC = c(as.character(df[i,72]), 
	             as.character(df[i,73]), 
	             as.character(df[i,74]), 
	             as.character(df[i,75]), 
	             as.character(df[i,76]))	
	
	a1_wFC = "Not_SB"
	
	if("female_biased" %in% SB_v_wFC){ a1_wFC = "SB"}
	if("male_biased" %in% SB_v_wFC){ a1_wFC = "SB"}
	
	b1_wFC = "Not_SB"
	b2_wFC = "Not_SB"

	FB_N_wFC = sum(SB_v_wFC == "female_biased")
	MB_N_wFC = sum(SB_v_wFC == "male_biased")
	SB_N_wFC = FB_N_wFC + MB_N_wFC
	
	N_sp_FB_wFC = c(N_sp_FB_wFC, FB_N_wFC)
	N_sp_MB_wFC = c(N_sp_MB_wFC, MB_N_wFC)
	N_sp_SB_wFC = c(N_sp_SB_wFC, SB_N_wFC)
			
	if("female_biased" %in% SB_v_wFC){ b1_wFC = "female_biased"}
	if("male_biased" %in% SB_v_wFC){ b2_wFC = "male_biased"}
	
	b3_wFC = paste(b1_wFC, b2_wFC, sep = "_")	
		
	SB_all_wFC = c(SB_all_wFC, a1_wFC)
	SB_dir_all_wFC = c(SB_dir_all_wFC, b3_wFC)	
	
	
	
	
	# print(i)

	# print(SB_v)	
	# print(a1 )
	# print(b3)
	# print(SB_v_wFC)
	# print(a1_wFC)
	# print(b3_wFC)	
}

df$SB_all           <- SB_all 
df$SB_dir_all       <- SB_dir_all
df$SB_all_wFC       <- SB_all_wFC
df$SB_dir_all_wFC   <- SB_dir_all_wFC
df$N_sp_SB <- N_sp_SB
df$N_sp_FB <- N_sp_FB
df$N_sp_MB <- N_sp_MB
df$N_sp_SB_wFC <- N_sp_SB_wFC
df$N_sp_FB_wFC <- N_sp_FB_wFC
df$N_sp_MB_wFC <- N_sp_MB_wFC
#print(SB_all)
return(df)

}


# data

dat_ALL_WB_10sp_SB_asex_2 <- get_SB_in_any_comp(dat_ALL_WB_10sp_SB_asex)
dat_ALL_RT_10sp_SB_asex_2 <- get_SB_in_any_comp(dat_ALL_RT_10sp_SB_asex)
dat_ALL_LG_10sp_SB_asex_2 <- get_SB_in_any_comp(dat_ALL_LG_10sp_SB_asex)
head(dat_ALL_WB_10sp_SB_asex_2, n = 20)

########### EXPORT

write.table(dat_ALL_WB_10sp_SB_asex_2, "ALL_WB_10sp_disp_allsepar_wCPM_SB_asex_2.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(dat_ALL_RT_10sp_SB_asex_2, "ALL_RT_10sp_disp_allsepar_wCPM_SB_asex_2.csv", sep = ',', quote = FALSE, row.names = FALSE)
write.table(dat_ALL_LG_10sp_SB_asex_2, "ALL_LG_10sp_disp_allsepar_wCPM_SB_asex_2.csv", sep = ',', quote = FALSE, row.names = FALSE)


################################################################################################################################################
############### boxplot SB 
################################################################################################################################################

make_long_table = function(dat,tiss){
	
	sp_u = "Tbi"
	df_t1 = as.data.frame(cbind(
		eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_log2FC_SA',sep=''))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_'  ,tiss,'_sexbias',sep='')))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_'  ,tiss,'_sexbias2',sep='')))))) 
	df_t1$sp     = rep(sp_u, length(df_t1[,1]))
	df_t1$group  = paste(df_t1$sp, df_t1$V2, sep = "_")
	df_t1$group2 = paste(df_t1$sp, df_t1$V3, sep = "_")
	
	sp_u = "Tce"
	df_t2 = as.data.frame(cbind(
		eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_log2FC_SA',sep=''))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_sexbias',sep='')))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_sexbias2',sep='')))))) 
	df_t2$sp     = rep(sp_u, length(df_t2[,1]))
	df_t2$group  = paste(df_t2$sp, df_t2$V2, sep = "_")
	df_t2$group2 = paste(df_t2$sp, df_t2$V3, sep = "_")
	
	sp_u = "Tcm"
	df_t3 = as.data.frame(cbind(
		eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_log2FC_SA',sep=''))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_sexbias',sep='')))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_sexbias2',sep='')))))) 
	df_t3$sp     = rep(sp_u, length(df_t3[,1]))
	df_t3$group  = paste(df_t3$sp, df_t3$V2, sep = "_")
	df_t3$group2 = paste(df_t3$sp, df_t3$V3, sep = "_")
	
	sp_u = "Tpa"
	df_t4 = as.data.frame(cbind(
		eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_log2FC_SA',sep=''))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_sexbias',sep='')))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_sexbias2',sep='')))))) 
	df_t4$sp     = rep(sp_u, length(df_t4[,1]))
	df_t4$group  = paste(df_t4$sp, df_t4$V2, sep = "_")
	df_t4$group2 = paste(df_t4$sp, df_t4$V3, sep = "_")
	
	sp_u = "Tps"
	df_t5 = as.data.frame(cbind(
		eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_log2FC_SA',sep=''))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_sexbias',sep='')))), 
		as.character(eval(parse(text=paste(dat,'$', sp_u, '_' ,tiss,'_sexbias2',sep='')))))) 
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

dat_ALL_WB_10sp_SB_asex_long <- make_long_table("dat_ALL_WB_10sp_SB_asex", "WB")
dat_ALL_RT_10sp_SB_asex_long <- make_long_table("dat_ALL_RT_10sp_SB_asex", "RT")
dat_ALL_LG_10sp_SB_asex_long <- make_long_table("dat_ALL_LG_10sp_SB_asex", "LG")

str(dat_ALL_WB_10sp_SB_asex_long)



#### join tissues together

dat_ALL_WB_10sp_SB_asex_long$tiss <- rep("WB", length(dat_ALL_WB_10sp_SB_asex_long[,1]))
dat_ALL_WB_10sp_SB_asex_long$sp_group_tiss <- as.factor(paste(dat_ALL_WB_10sp_SB_asex_long$SBgroup, dat_ALL_WB_10sp_SB_asex_long$tiss, sep = "_"))
dat_ALL_WB_10sp_SB_asex_long$sp_tiss <- as.factor(paste(dat_ALL_WB_10sp_SB_asex_long$sp, dat_ALL_WB_10sp_SB_asex_long$tiss, sep = "_"))

dat_ALL_RT_10sp_SB_asex_long$tiss <- rep("RT", length(dat_ALL_RT_10sp_SB_asex_long[,1]))
dat_ALL_RT_10sp_SB_asex_long$sp_group_tiss <- as.factor(paste(dat_ALL_RT_10sp_SB_asex_long$SBgroup, dat_ALL_RT_10sp_SB_asex_long$tiss, sep = "_"))
dat_ALL_RT_10sp_SB_asex_long$sp_tiss <- as.factor(paste(dat_ALL_RT_10sp_SB_asex_long$sp, dat_ALL_RT_10sp_SB_asex_long$tiss, sep = "_"))

dat_ALL_LG_10sp_SB_asex_long$tiss <- rep("LG", length(dat_ALL_LG_10sp_SB_asex_long[,1]))
dat_ALL_LG_10sp_SB_asex_long$sp_group_tiss <- as.factor(paste(dat_ALL_LG_10sp_SB_asex_long$SBgroup, dat_ALL_LG_10sp_SB_asex_long$tiss, sep = "_"))
dat_ALL_LG_10sp_SB_asex_long$sp_tiss <- as.factor(paste(dat_ALL_LG_10sp_SB_asex_long$sp, dat_ALL_LG_10sp_SB_asex_long$tiss, sep = "_"))


dat_ALL_WBRTLG_10sp_SB_asex_long <- rbind(dat_ALL_WB_10sp_SB_asex_long, dat_ALL_RT_10sp_SB_asex_long, dat_ALL_LG_10sp_SB_asex_long)
dat_ALL_WBRTLG_10sp_SB_asex_long$sp_tiss_ord <- ordered(dat_ALL_WBRTLG_10sp_SB_asex_long$sp_tiss, levels = c("Tbi_WB", "Tce_WB", "Tps_WB", "Tcm_WB", "Tpa_WB", "Tbi_RT", "Tce_RT", "Tps_RT", "Tcm_RT", "Tpa_RT", "Tbi_LG", "Tce_LG", "Tps_LG", "Tcm_LG", "Tpa_LG"))	
dat_ALL_WBRTLG_10sp_SB_asex_long$sexbias_ord <- ordered(dat_ALL_WBRTLG_10sp_SB_asex_long$sexbias, levels = c("female_biased",  "Unbiased", "male_biased"))	
dat_ALL_WBRTLG_10sp_SB_asex_long$sexbias2_ord <- ordered(dat_ALL_WBRTLG_10sp_SB_asex_long$sexbias2, levels = c("female_biased",  "Unbiased", "male_biased"))	

str(dat_ALL_WBRTLG_10sp_SB_asex_long)

######### remove NAs

dat_ALL_WB_10sp_SB_asex_long_c <- na.omit(dat_ALL_WB_10sp_SB_asex_long)
dat_ALL_RT_10sp_SB_asex_long_c <- na.omit(dat_ALL_RT_10sp_SB_asex_long)
dat_ALL_LG_10sp_SB_asex_long_c <- na.omit(dat_ALL_LG_10sp_SB_asex_long)

dat_ALL_WBRTLG_10sp_SB_asex_long_c <- na.omit(dat_ALL_WBRTLG_10sp_SB_asex_long)

###### add some empty vals to space the boxplot



dat_ALL_WBRTLG_10sp_SB_asex_long_c$sp_tiss2 <- as.character(dat_ALL_WBRTLG_10sp_SB_asex_long_c$sp_tiss)
str(dat_ALL_WBRTLG_10sp_SB_asex_long_c)

dat_ALL_WBRTLG_10sp_SB_asex_long_c_2 <- rbind(dat_ALL_WBRTLG_10sp_SB_asex_long_c, c(NA,"Unbiased","Unbiased",NA,NA,NA,NA,NA,NA,"WB", NA,NA,NA,NA,NA, "skip_WB"))
dat_ALL_WBRTLG_10sp_SB_asex_long_c_2 <- rbind(dat_ALL_WBRTLG_10sp_SB_asex_long_c_2, c(NA,"Unbiased","Unbiased",NA,NA,NA,NA,NA,NA,"RT", NA,NA,NA,NA,NA, "skip_RT"))
dat_ALL_WBRTLG_10sp_SB_asex_long_c_2$log2FC_SA <- as.numeric(as.character(dat_ALL_WBRTLG_10sp_SB_asex_long_c_2$log2FC_SA ))

dat_ALL_WBRTLG_10sp_SB_asex_long_c_2$sp_tiss2 <- as.factor(dat_ALL_WBRTLG_10sp_SB_asex_long_c_2$sp_tiss2)
dat_ALL_WBRTLG_10sp_SB_asex_long_c_2$sp_tiss2_ord <- ordered(dat_ALL_WBRTLG_10sp_SB_asex_long_c_2$sp_tiss2, levels = c("Tbi_WB", "Tce_WB", "Tps_WB", "Tcm_WB", "Tpa_WB", "skip_WB", "Tbi_RT", "Tce_RT", "Tps_RT", "Tcm_RT", "Tpa_RT", "skip_RT", "Tbi_LG", "Tce_LG", "Tps_LG", "Tcm_LG", "Tpa_LG"))	

str(dat_ALL_WBRTLG_10sp_SB_asex_long_c_2)

##### plot boxplots

WBRTLG_SA_box_s1 <- ggplot(dat_ALL_WBRTLG_10sp_SB_asex_long_c_2, aes(sp_tiss2_ord, log2FC_SA)) + 
	theme_classic() +
	geom_boxplot(aes(fill = factor(sexbias_ord)),position=position_dodge(0.65), width = 0.45, outlier.size = 0, lwd = 0.5, fatten = 1, outlier.shape = NA) +
	coord_cartesian(ylim=c(-3,3)) +
	ylab ("log2 fold change in asexual females") +
	xlab ("Gene-class by species") + 
	scale_fill_manual(values=c("firebrick2", "grey", "royalblue2")) + geom_hline(yintercept = 0) + labs(fill='') +
	ggtitle("SB class = FDR <0.05")


WBRTLG_SA_box_s2 <- ggplot(dat_ALL_WBRTLG_10sp_SB_asex_long_c_2, aes(sp_tiss2_ord, log2FC_SA)) + 
	theme_classic() +
	geom_boxplot(aes(fill = factor(sexbias2_ord)),position=position_dodge(0.65), width = 0.45, outlier.size = 0, lwd = 0.5, fatten = 1, outlier.shape = NA) +
	coord_cartesian(ylim=c(-3,3)) +
	ylab ("log2 fold change in asexual females") +
	xlab ("Gene-class by species") + 
	scale_fill_manual(values=c("firebrick2", "grey", "royalblue2")) + geom_hline(yintercept = 0) + labs(fill='') +
	ggtitle("SB class = FDR <0.05 FC > 2")


pdf(paste("Sex_asex_10sp_LG_FC_boxplot_s2_FDR005.pdf",sep = ""), width = 8.2, height = 8)
WBRTLG_SA_box_s1
dev.off()
getwd() ## wh


pdf(paste("Sex_asex_10sp_LG_FC_boxplot_s2_FDR005FC2.pdf",sep = ""), width = 8.2, height = 8)
WBRTLG_SA_box_s2
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
wilk_test("WB", dat_ALL_WB_10sp_SB_asex_long_c)$s1_out_table,	
wilk_test("RT", dat_ALL_RT_10sp_SB_asex_long_c)$s1_out_table,
wilk_test("LG", dat_ALL_LG_10sp_SB_asex_long_c)$s1_out_table))		

wilcox_S1_FDR005$FB_wilp_s1 <- as.numeric(as.character(wilcox_S1_FDR005$FB_wilp_s1))
wilcox_S1_FDR005$MB_wilp_s1 <- as.numeric(as.character(wilcox_S1_FDR005$MB_wilp_s1))

wilcox_S1_FDR005$FB_wilFDR_s1 <- p.adjust(wilcox_S1_FDR005$FB_wilp_s1, method = "BH")
wilcox_S1_FDR005$MB_wilFDR_s1 <- p.adjust(wilcox_S1_FDR005$MB_wilp_s1, method = "BH")


wilcox_S2_FDR005_FC2 <- as.data.frame(rbind(
wilk_test("WB", dat_ALL_WB_10sp_SB_asex_long_c)$s2_out_table,	
wilk_test("RT", dat_ALL_RT_10sp_SB_asex_long_c)$s2_out_table,
wilk_test("LG", dat_ALL_LG_10sp_SB_asex_long_c)$s2_out_table))	

wilcox_S2_FDR005_FC2$FB_wilp_s2 <- as.numeric(as.character(wilcox_S2_FDR005_FC2$FB_wilp_s2))
wilcox_S2_FDR005_FC2$MB_wilp_s2 <- as.numeric(as.character(wilcox_S2_FDR005_FC2$MB_wilp_s2))

wilcox_S2_FDR005_FC2$FB_wilFDR_s1 <- p.adjust(wilcox_S2_FDR005_FC2$FB_wilp_s2, method = "BH")
wilcox_S2_FDR005_FC2$MB_wilFDR_s1 <- p.adjust(wilcox_S2_FDR005_FC2$MB_wilp_s2, method = "BH")


write.csv(wilcox_S1_FDR005, "wilcox_S1_FDR005.csv")
write.csv(wilcox_S2_FDR005_FC2, "wilcox_S2_FDR005_FC2.csv")


################################################################################################################################################
############### Plot N SB,MB,FB genes 
################################################################################################################################################

## get N genes

N_sex_bias_genes = wilcox_S1_FDR005[,1:5]
N_sex_bias_genes$group = paste(N_sex_bias_genes$sps, N_sex_bias_genes$tiss, sep = "_")

N_sex_bias_genes$group_ordered <-  ordered(N_sex_bias_genes$group, levels=c(
"Tbi_WB","Tbi_RT","Tbi_LG",
"Tce_WB","Tce_RT","Tce_LG",
"Tps_WB","Tps_RT","Tps_LG",
"Tcm_WB","Tcm_RT","Tcm_LG",
"Tpa_WB","Tpa_RT","Tpa_LG"
))

N_sex_bias_genes$MB_N_s1 = as.numeric(as.character(N_sex_bias_genes$MB_N_s1))
N_sex_bias_genes$FB_N_s1 = as.numeric(as.character(N_sex_bias_genes$FB_N_s1))
N_sex_bias_genes$UB_N_s1 = as.numeric(as.character(N_sex_bias_genes$UB_N_s1))


N_sex_bias_genes$SB_N_s1 =  N_sex_bias_genes$MB_N_s1 + N_sex_bias_genes$FB_N_s1
write.csv(N_sex_bias_genes, "N_sex_bias_genes_10sp.csv")

N_sex_bias_genes_wFC = wilcox_S2_FDR005_FC2[,1:5]
N_sex_bias_genes_wFC$group = paste(N_sex_bias_genes_wFC$sps, N_sex_bias_genes_wFC$tiss, sep = "_")

N_sex_bias_genes_wFC$group_ordered <-  ordered(N_sex_bias_genes_wFC$group, levels=c(
"Tbi_WB","Tbi_RT","Tbi_LG",
"Tce_WB","Tce_RT","Tce_LG",
"Tps_WB","Tps_RT","Tps_LG",
"Tcm_WB","Tcm_RT","Tcm_LG",
"Tpa_WB","Tpa_RT","Tpa_LG"
))

N_sex_bias_genes_wFC$MB_N_s2 = as.numeric(as.character(N_sex_bias_genes_wFC$MB_N_s2))
N_sex_bias_genes_wFC$FB_N_s2 = as.numeric(as.character(N_sex_bias_genes_wFC$FB_N_s2))
N_sex_bias_genes_wFC$UB_N_s2 = as.numeric(as.character(N_sex_bias_genes_wFC$UB_N_s2))

N_sex_bias_genes_wFC$SB_N_s2 =  N_sex_bias_genes_wFC$MB_N_s2 + N_sex_bias_genes_wFC$FB_N_s2
write.csv(N_sex_bias_genes_wFC, "N_sex_bias_genes_10sp_wFC.csv")


p1_NSB <- ggplot(N_sex_bias_genes, aes(factor(group_ordered), SB_N_s1, fill = tiss )) + 
	geom_bar(stat="identity", position = "dodge") + 
	theme_bw() +
	xlab ("Species pair") + 
	ylab ("Number of sex-biased genes, FDR < 0.05")  +
	scale_y_continuous(expand = c(0,0)) + 
	coord_cartesian(ylim=c(-0,2000)) +
	ggtitle("Number of sex-biased genes | FDR < 0.05")
	
p2_NSB = p1_NSB + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
   scale_fill_manual(values=c("black", "firebrick3", "#56B4E9"))


p1_NSB_wFC <- ggplot(N_sex_bias_genes_wFC, aes(factor(group_ordered), SB_N_s2, fill = tiss )) + 
	geom_bar(stat="identity", position = "dodge") + 
	theme_bw() +
	xlab ("Species pair") + 
	ylab ("Number of sex-biased genes, FDR < 0.05 | FC > 2")  +
	scale_y_continuous(expand = c(0,0)) + 
	coord_cartesian(ylim=c(-0,1200)) +
	ggtitle("Number of sex-biased genes | FDR < 0.05 | FC > 2")
	
p2_NSB_wFC = p1_NSB_wFC + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
   scale_fill_manual(values=c("black", "firebrick3", "#56B4E9"))

### export

pdf("N_SB_genes_10sp.pdf", width = 6, height = 7)
p2_NSB
dev.off()

pdf("N_SB_genes_10sp_wFC.pdf", width = 6, height = 7)
p2_NSB_wFC
dev.off()



### full table


N_sex_bias_genes_t2 <- cbind(N_sex_bias_genes, rep("FDR", length( N_sex_bias_genes[,1])))
N_sex_bias_genes_wFC_t2 <- cbind(N_sex_bias_genes_wFC, rep("FDR + FC", length( N_sex_bias_genes_wFC[,1])))


colnames(N_sex_bias_genes_t2) = c("sp", "tiss", "FB_N", "MB_N", "UB_N", "group", "group_ordered", "SB_N", "Threshold")
colnames(N_sex_bias_genes_wFC_t2) = c("sp", "tiss", "FB_N", "MB_N", "UB_N", "group", "group_ordered", "SB_N", "Threshold")

N_sex_bias_genes_all <- rbind(N_sex_bias_genes_t2,N_sex_bias_genes_wFC_t2)

N_sex_bias_genes_all$SB_dir = ifelse(N_sex_bias_genes_all$FB_N > N_sex_bias_genes_all$MB_N, "FB",
							  ifelse(N_sex_bias_genes_all$FB_N < N_sex_bias_genes_all$MB_N, "MB", "NB"))
							  							  
# test_bias


df_temp = as.data.frame(cbind(N_sex_bias_genes_all$FB_N,N_sex_bias_genes_all$MB_N))
df_temp1 <- cbind(df_temp, t(apply(df_temp, 1, function(x) {
             ch <- chisq.test(x, p = c(1/2, 1/2), correct = FALSE)
             c(unname(ch$statistic), ch$p.value)})))

colnames(df_temp1)[3:4] <- c('xsquared', 'pvalue')

df_temp1$FDR <- p.adjust(df_temp1$pvalue, method = "BH")
df_temp1


N_sex_bias_genes_all2 <- cbind(N_sex_bias_genes_all, df_temp1$xsquared, df_temp1$pvalue,df_temp1$FDR)

colnames(N_sex_bias_genes_all2) = c("sp", "tiss", "FB_N", "MB_N", "UB_N", "group", "group_ordered", "SB_N", "Threshold","SB_dir", "xsquared", "pvalue", "FDR")

# export

write.csv(N_sex_bias_genes_all2, "N_sex_bias_genes_10sp_ALL.csv")


################################################################################################
############ CPM heat maps

## get av CPM

head(dat_ALL_WB_10sp_SB_asex_2)

# all

## WB

dat_sexuals_cpm_WB_10sp_allg <- as.data.frame(
cbind(
((dat_ALL_WB_10sp_SB_asex_2$cpm_Tbi_SF_WB_Md_Re1 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tbi_SF_WB_Md_Re2 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tbi_SF_WB_Md_Re3) / 3),
((dat_ALL_WB_10sp_SB_asex_2$cpm_Tbi_SM_WB_Md_Re1 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tbi_SM_WB_Md_Re2 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tbi_SM_WB_Md_Re3) / 3),
((dat_ALL_WB_10sp_SB_asex_2$cpm_Tce_SF_WB_Md_Re1 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tce_SF_WB_Md_Re2 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tce_SF_WB_Md_Re3) / 3),
((dat_ALL_WB_10sp_SB_asex_2$cpm_Tce_SM_WB_Md_Re1 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tce_SM_WB_Md_Re2 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tce_SM_WB_Md_Re3) / 3),
((dat_ALL_WB_10sp_SB_asex_2$cpm_Tcm_SF_WB_Md_Re1 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tcm_SF_WB_Md_Re2 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tcm_SF_WB_Md_Re3) / 3),
((dat_ALL_WB_10sp_SB_asex_2$cpm_Tcm_SM_WB_Md_Re1 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tcm_SM_WB_Md_Re2 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tcm_SM_WB_Md_Re3) / 3),
((dat_ALL_WB_10sp_SB_asex_2$cpm_Tpa_SF_WB_Md_Re1 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tpa_SF_WB_Md_Re2 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tpa_SF_WB_Md_Re3) / 3),
((dat_ALL_WB_10sp_SB_asex_2$cpm_Tpa_SM_WB_Md_Re1 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tpa_SM_WB_Md_Re2 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tpa_SM_WB_Md_Re3) / 3),
((dat_ALL_WB_10sp_SB_asex_2$cpm_Tps_SF_WB_Md_Re1 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tps_SF_WB_Md_Re2 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tps_SF_WB_Md_Re3) / 3),
((dat_ALL_WB_10sp_SB_asex_2$cpm_Tps_SM_WB_Md_Re1 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tps_SM_WB_Md_Re2 + dat_ALL_WB_10sp_SB_asex_2$cpm_Tps_SM_WB_Md_Re3) / 3)))

colnames(dat_sexuals_cpm_WB_10sp_allg) <- c(
"Tbi_WB_SF_cpm", "Tbi_WB_SM_cpm", 
"Tce_WB_SF_cpm", "Tce_WB_SM_cpm",
"Tcm_WB_SF_cpm", "Tcm_WB_SM_cpm", 
"Tpa_WB_SF_cpm", "Tpa_WB_SM_cpm", 
"Tps_WB_SF_cpm", "Tps_WB_SM_cpm")
rownames(dat_sexuals_cpm_WB_10sp_allg) <- dat_ALL_WB_10sp_SB_asex_2$genename

head(dat_sexuals_cpm_WB_10sp_allg)

png(filename = "heatmap_WB_10sp_logcpm_allg.png", width = 4, height = 10.5, units = "in", bg = "white", res = 300)
pheatmap(dat_sexuals_cpm_WB_10sp_allg, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "ward.D2", show_rownames = F, 
main = paste("WB_logcpm - allg - N = ", length(dat_sexuals_cpm_WB_10sp_allg[,1])))
dev.off()
getwd()


## RT

dat_sexuals_cpm_RT_10sp_allg <- as.data.frame(
cbind(
((dat_ALL_RT_10sp_SB_asex_2$cpm_Tbi_SF_RT_Md_Re1 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tbi_SF_RT_Md_Re2 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tbi_SF_RT_Md_Re3) / 3),
((dat_ALL_RT_10sp_SB_asex_2$cpm_Tbi_SM_RT_Md_Re1 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tbi_SM_RT_Md_Re2 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tbi_SM_RT_Md_Re3) / 3),
((dat_ALL_RT_10sp_SB_asex_2$cpm_Tce_SF_RT_Md_Re1 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tce_SF_RT_Md_Re2 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tce_SF_RT_Md_Re3) / 3),
((dat_ALL_RT_10sp_SB_asex_2$cpm_Tce_SM_RT_Md_Re1 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tce_SM_RT_Md_Re2 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tce_SM_RT_Md_Re3) / 3),
((dat_ALL_RT_10sp_SB_asex_2$cpm_Tcm_SF_RT_Md_Re1 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tcm_SF_RT_Md_Re2 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tcm_SF_RT_Md_Re3) / 3),
((dat_ALL_RT_10sp_SB_asex_2$cpm_Tcm_SM_RT_Md_Re1 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tcm_SM_RT_Md_Re2 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tcm_SM_RT_Md_Re3) / 3),
((dat_ALL_RT_10sp_SB_asex_2$cpm_Tpa_SF_RT_Md_Re1 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tpa_SF_RT_Md_Re2 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tpa_SF_RT_Md_Re3) / 3),
((dat_ALL_RT_10sp_SB_asex_2$cpm_Tpa_SM_RT_Md_Re1 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tpa_SM_RT_Md_Re2 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tpa_SM_RT_Md_Re3) / 3),
((dat_ALL_RT_10sp_SB_asex_2$cpm_Tps_SF_RT_Md_Re1 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tps_SF_RT_Md_Re2 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tps_SF_RT_Md_Re3) / 3),
((dat_ALL_RT_10sp_SB_asex_2$cpm_Tps_SM_RT_Md_Re1 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tps_SM_RT_Md_Re2 + dat_ALL_RT_10sp_SB_asex_2$cpm_Tps_SM_RT_Md_Re3) / 3)))

colnames(dat_sexuals_cpm_RT_10sp_allg) <- c(
"Tbi_RT_SF_cpm", "Tbi_RT_SM_cpm", 
"Tce_RT_SF_cpm", "Tce_RT_SM_cpm",
"Tcm_RT_SF_cpm", "Tcm_RT_SM_cpm", 
"Tpa_RT_SF_cpm", "Tpa_RT_SM_cpm", 
"Tps_RT_SF_cpm", "Tps_RT_SM_cpm")
rownames(dat_sexuals_cpm_RT_10sp_allg) <- dat_ALL_RT_10sp_SB_asex_2$genename

head(dat_sexuals_cpm_RT_10sp_allg)

png(filename = "heatmap_RT_10sp_logcpm_allg.png", width = 4, height = 10.5, units = "in", bg = "white", res = 300)
pheatmap(dat_sexuals_cpm_RT_10sp_allg, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "ward.D2", show_rownames = F, 
main = paste("RT_logcpm - allg - N = ", length(dat_sexuals_cpm_RT_10sp_allg[,1])))
dev.off()
getwd()


## LG

dat_sexuals_cpm_LG_10sp_allg <- as.data.frame(
cbind(
((dat_ALL_LG_10sp_SB_asex_2$cpm_Tbi_SF_LG_Md_Re1 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tbi_SF_LG_Md_Re2 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tbi_SF_LG_Md_Re3) / 3),
((dat_ALL_LG_10sp_SB_asex_2$cpm_Tbi_SM_LG_Md_Re1 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tbi_SM_LG_Md_Re2 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tbi_SM_LG_Md_Re3) / 3),
((dat_ALL_LG_10sp_SB_asex_2$cpm_Tce_SF_LG_Md_Re1 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tce_SF_LG_Md_Re2 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tce_SF_LG_Md_Re3) / 3),
((dat_ALL_LG_10sp_SB_asex_2$cpm_Tce_SM_LG_Md_Re1 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tce_SM_LG_Md_Re2 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tce_SM_LG_Md_Re3) / 3),
((dat_ALL_LG_10sp_SB_asex_2$cpm_Tcm_SF_LG_Md_Re1 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tcm_SF_LG_Md_Re2 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tcm_SF_LG_Md_Re3) / 3),
((dat_ALL_LG_10sp_SB_asex_2$cpm_Tcm_SM_LG_Md_Re1 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tcm_SM_LG_Md_Re2 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tcm_SM_LG_Md_Re3) / 3),
((dat_ALL_LG_10sp_SB_asex_2$cpm_Tpa_SF_LG_Md_Re1 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tpa_SF_LG_Md_Re2 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tpa_SF_LG_Md_Re3) / 3),
((dat_ALL_LG_10sp_SB_asex_2$cpm_Tpa_SM_LG_Md_Re1 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tpa_SM_LG_Md_Re2 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tpa_SM_LG_Md_Re3) / 3),
((dat_ALL_LG_10sp_SB_asex_2$cpm_Tps_SF_LG_Md_Re1 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tps_SF_LG_Md_Re2 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tps_SF_LG_Md_Re3) / 3),
((dat_ALL_LG_10sp_SB_asex_2$cpm_Tps_SM_LG_Md_Re1 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tps_SM_LG_Md_Re2 + dat_ALL_LG_10sp_SB_asex_2$cpm_Tps_SM_LG_Md_Re3) / 3)))

colnames(dat_sexuals_cpm_LG_10sp_allg) <- c(
"Tbi_LG_SF_cpm", "Tbi_LG_SM_cpm", 
"Tce_LG_SF_cpm", "Tce_LG_SM_cpm",
"Tcm_LG_SF_cpm", "Tcm_LG_SM_cpm", 
"Tpa_LG_SF_cpm", "Tpa_LG_SM_cpm", 
"Tps_LG_SF_cpm", "Tps_LG_SM_cpm")
rownames(dat_sexuals_cpm_LG_10sp_allg) <- dat_ALL_LG_10sp_SB_asex_2$genename

head(dat_sexuals_cpm_LG_10sp_allg)

png(filename = "heatmap_LG_10sp_logcpm_allg.png", width = 4, height = 10.5, units = "in", bg = "white", res = 300)
pheatmap(dat_sexuals_cpm_LG_10sp_allg, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "ward.D2", show_rownames = F, 
main = paste("LG_logcpm - allg - N = ", length(dat_sexuals_cpm_LG_10sp_allg[,1])))
dev.off()
getwd()


########## bootstrap

# allg

N_boot = 10000

png(filename = paste("heatmap_bootstrap_WB_10sp_logcpm_allg_N=", N_boot, ".png", sep = ""), width = 8, height = 10.5, units = "in", bg = "white", res = 300)
Pv_clust_all_WB <- pvclust(dat_sexuals_cpm_WB_10sp_allg, method.hclust="ward.D2", method.dist="euclidean", nboot= N_boot) ## put nboot to 10000 for serious work
plot(Pv_clust_all_WB)
dev.off()
getwd()

png(filename = paste("heatmap_bootstrap_RT_10sp_logcpm_allg_N=", N_boot, ".png", sep = ""), width = 8, height = 10.5, units = "in", bg = "white", res = 300)
Pv_clust_all_RT <- pvclust(dat_sexuals_cpm_RT_10sp_allg, method.hclust="ward.D2", method.dist="euclidean", nboot= N_boot) ## put nboot to 10000 for serious work
plot(Pv_clust_all_RT)
dev.off()
getwd()

png(filename = paste("heatmap_bootstrap_LG_10sp_logcpm_allg_N=", N_boot, ".png", sep = ""), width = 8, height = 10.5, units = "in", bg = "white", res = 300)
Pv_clust_all_LG <- pvclust(dat_sexuals_cpm_LG_10sp_allg, method.hclust="ward.D2", method.dist="euclidean", nboot= N_boot) ## put nboot to 10000 for serious work
plot(Pv_clust_all_LG)
dev.off()
getwd()




###############################################################################################################################################
# are sex biased genes the same?
###############################################################################################################################################


head(dat_ALL_WB_10sp_SB_asex_2)


Tbi_WB_SB_FB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tbi_WB_sexbias == "female_biased")[,1])
Tce_WB_SB_FB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tce_WB_sexbias == "female_biased")[,1])
Tcm_WB_SB_FB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tcm_WB_sexbias == "female_biased")[,1])
Tpa_WB_SB_FB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tpa_WB_sexbias == "female_biased")[,1])
Tps_WB_SB_FB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tps_WB_sexbias == "female_biased")[,1])

Tbi_WB_SB_MB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tbi_WB_sexbias == "male_biased")[,1])
Tce_WB_SB_MB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tce_WB_sexbias == "male_biased")[,1])
Tcm_WB_SB_MB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tcm_WB_sexbias == "male_biased")[,1])
Tpa_WB_SB_MB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tpa_WB_sexbias == "male_biased")[,1])
Tps_WB_SB_MB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tps_WB_sexbias == "male_biased")[,1])

Tbi_WB_SB_SB_genes <- c(Tbi_WB_SB_FB_genes, Tbi_WB_SB_MB_genes )
Tce_WB_SB_SB_genes <- c(Tce_WB_SB_FB_genes, Tce_WB_SB_MB_genes )
Tcm_WB_SB_SB_genes <- c(Tcm_WB_SB_FB_genes, Tcm_WB_SB_MB_genes )
Tpa_WB_SB_SB_genes <- c(Tpa_WB_SB_FB_genes, Tpa_WB_SB_MB_genes )
Tps_WB_SB_SB_genes <- c(Tps_WB_SB_FB_genes, Tps_WB_SB_MB_genes )


Tbi_RT_SB_FB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tbi_RT_sexbias == "female_biased")[,1])
Tce_RT_SB_FB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tce_RT_sexbias == "female_biased")[,1])
Tcm_RT_SB_FB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tcm_RT_sexbias == "female_biased")[,1])
Tpa_RT_SB_FB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tpa_RT_sexbias == "female_biased")[,1])
Tps_RT_SB_FB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tps_RT_sexbias == "female_biased")[,1])

Tbi_RT_SB_MB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tbi_RT_sexbias == "male_biased")[,1])
Tce_RT_SB_MB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tce_RT_sexbias == "male_biased")[,1])
Tcm_RT_SB_MB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tcm_RT_sexbias == "male_biased")[,1])
Tpa_RT_SB_MB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tpa_RT_sexbias == "male_biased")[,1])
Tps_RT_SB_MB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tps_RT_sexbias == "male_biased")[,1])

Tbi_RT_SB_SB_genes <- c(Tbi_RT_SB_FB_genes, Tbi_RT_SB_MB_genes )
Tce_RT_SB_SB_genes <- c(Tce_RT_SB_FB_genes, Tce_RT_SB_MB_genes )
Tcm_RT_SB_SB_genes <- c(Tcm_RT_SB_FB_genes, Tcm_RT_SB_MB_genes )
Tpa_RT_SB_SB_genes <- c(Tpa_RT_SB_FB_genes, Tpa_RT_SB_MB_genes )
Tps_RT_SB_SB_genes <- c(Tps_RT_SB_FB_genes, Tps_RT_SB_MB_genes )

Tbi_LG_SB_FB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tbi_LG_sexbias == "female_biased")[,1])
Tce_LG_SB_FB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tce_LG_sexbias == "female_biased")[,1])
Tcm_LG_SB_FB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tcm_LG_sexbias == "female_biased")[,1])
Tpa_LG_SB_FB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tpa_LG_sexbias == "female_biased")[,1])
Tps_LG_SB_FB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tps_LG_sexbias == "female_biased")[,1])

Tbi_LG_SB_MB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tbi_LG_sexbias == "male_biased")[,1])
Tce_LG_SB_MB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tce_LG_sexbias == "male_biased")[,1])
Tcm_LG_SB_MB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tcm_LG_sexbias == "male_biased")[,1])
Tpa_LG_SB_MB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tpa_LG_sexbias == "male_biased")[,1])
Tps_LG_SB_MB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tps_LG_sexbias == "male_biased")[,1])

Tbi_LG_SB_SB_genes <- c(Tbi_LG_SB_FB_genes, Tbi_LG_SB_MB_genes )
Tce_LG_SB_SB_genes <- c(Tce_LG_SB_FB_genes, Tce_LG_SB_MB_genes )
Tcm_LG_SB_SB_genes <- c(Tcm_LG_SB_FB_genes, Tcm_LG_SB_MB_genes )
Tpa_LG_SB_SB_genes <- c(Tpa_LG_SB_FB_genes, Tpa_LG_SB_MB_genes )
Tps_LG_SB_SB_genes <- c(Tps_LG_SB_FB_genes, Tps_LG_SB_MB_genes )



#### output ### between species

venn_WB_SB_SB <- five_sp_DE_venn(Tbi_WB_SB_SB_genes, Tce_WB_SB_SB_genes,Tcm_WB_SB_SB_genes, Tpa_WB_SB_SB_genes, Tps_WB_SB_SB_genes , "WB_SB_SB")
venn_RT_SB_SB <- five_sp_DE_venn(Tbi_RT_SB_SB_genes, Tce_RT_SB_SB_genes,Tcm_RT_SB_SB_genes, Tpa_RT_SB_SB_genes, Tps_RT_SB_SB_genes , "RT_SB_SB")
venn_LG_SB_SB <- five_sp_DE_venn(Tbi_LG_SB_SB_genes, Tce_LG_SB_SB_genes,Tcm_LG_SB_SB_genes, Tpa_LG_SB_SB_genes, Tps_LG_SB_SB_genes , "LG_SB_SB")

venn_WB_SB_FB <- five_sp_DE_venn(Tbi_WB_SB_FB_genes, Tce_WB_SB_FB_genes,Tcm_WB_SB_FB_genes, Tpa_WB_SB_FB_genes, Tps_WB_SB_FB_genes , "WB_SB_FB")
venn_RT_SB_FB <- five_sp_DE_venn(Tbi_RT_SB_FB_genes, Tce_RT_SB_FB_genes,Tcm_RT_SB_FB_genes, Tpa_RT_SB_FB_genes, Tps_RT_SB_FB_genes , "RT_SB_FB")
venn_LG_SB_FB <- five_sp_DE_venn(Tbi_LG_SB_FB_genes, Tce_LG_SB_FB_genes,Tcm_LG_SB_FB_genes, Tpa_LG_SB_FB_genes, Tps_LG_SB_FB_genes , "LG_SB_FB")

venn_WB_SB_MB <- five_sp_DE_venn(Tbi_WB_SB_MB_genes, Tce_WB_SB_MB_genes,Tcm_WB_SB_MB_genes, Tpa_WB_SB_MB_genes, Tps_WB_SB_MB_genes , "WB_SB_MB")
venn_RT_SB_MB <- five_sp_DE_venn(Tbi_RT_SB_MB_genes, Tce_RT_SB_MB_genes,Tcm_RT_SB_MB_genes, Tpa_RT_SB_MB_genes, Tps_RT_SB_MB_genes , "RT_SB_MB")
venn_LG_SB_MB <- five_sp_DE_venn(Tbi_LG_SB_MB_genes, Tce_LG_SB_MB_genes,Tcm_LG_SB_MB_genes, Tpa_LG_SB_MB_genes, Tps_LG_SB_MB_genes , "LG_SB_MB")


venn_WBRTLG_SB_SB <- arrangeGrob(
gTree(children=venn_WB_SB_SB), 
gTree(children=venn_RT_SB_SB),
gTree(children=venn_LG_SB_SB),
ncol = 1 )

ggsave(file="venn_WBRTLG_SB_SB.png", venn_WBRTLG_SB_SB, width = 6, height = 18)


venn_WBRTLG_SB_FB <- arrangeGrob(
gTree(children=venn_WB_SB_FB), 
gTree(children=venn_RT_SB_FB),
gTree(children=venn_LG_SB_FB),
ncol = 1 )

ggsave(file="venn_WBRTLG_SB_FB.png", venn_WBRTLG_SB_FB, width = 6, height = 18)

venn_WBRTLG_SB_MB <- arrangeGrob(
gTree(children=venn_WB_SB_MB), 
gTree(children=venn_RT_SB_MB),
gTree(children=venn_LG_SB_MB),
ncol = 1 )

ggsave(file="venn_WBRTLG_SB_MB.png", venn_WBRTLG_SB_MB, width = 6, height = 18)



## SBwFC

Tbi_WB_SBwFC_FB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tbi_WB_sexbias2 == "female_biased")[,1])
Tce_WB_SBwFC_FB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tce_WB_sexbias2 == "female_biased")[,1])
Tcm_WB_SBwFC_FB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tcm_WB_sexbias2 == "female_biased")[,1])
Tpa_WB_SBwFC_FB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tpa_WB_sexbias2 == "female_biased")[,1])
Tps_WB_SBwFC_FB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tps_WB_sexbias2 == "female_biased")[,1])

Tbi_WB_SBwFC_MB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tbi_WB_sexbias2 == "male_biased")[,1])
Tce_WB_SBwFC_MB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tce_WB_sexbias2 == "male_biased")[,1])
Tcm_WB_SBwFC_MB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tcm_WB_sexbias2 == "male_biased")[,1])
Tpa_WB_SBwFC_MB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tpa_WB_sexbias2 == "male_biased")[,1])
Tps_WB_SBwFC_MB_genes <- as.character(subset(dat_ALL_WB_10sp_SB_asex_2, dat_ALL_WB_10sp_SB_asex_2$Tps_WB_sexbias2 == "male_biased")[,1])

Tbi_WB_SBwFC_SB_genes <- c(Tbi_WB_SBwFC_FB_genes, Tbi_WB_SBwFC_MB_genes )
Tce_WB_SBwFC_SB_genes <- c(Tce_WB_SBwFC_FB_genes, Tce_WB_SBwFC_MB_genes )
Tcm_WB_SBwFC_SB_genes <- c(Tcm_WB_SBwFC_FB_genes, Tcm_WB_SBwFC_MB_genes )
Tpa_WB_SBwFC_SB_genes <- c(Tpa_WB_SBwFC_FB_genes, Tpa_WB_SBwFC_MB_genes )
Tps_WB_SBwFC_SB_genes <- c(Tps_WB_SBwFC_FB_genes, Tps_WB_SBwFC_MB_genes )


Tbi_RT_SBwFC_FB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tbi_RT_sexbias2 == "female_biased")[,1])
Tce_RT_SBwFC_FB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tce_RT_sexbias2 == "female_biased")[,1])
Tcm_RT_SBwFC_FB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tcm_RT_sexbias2 == "female_biased")[,1])
Tpa_RT_SBwFC_FB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tpa_RT_sexbias2 == "female_biased")[,1])
Tps_RT_SBwFC_FB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tps_RT_sexbias2 == "female_biased")[,1])

Tbi_RT_SBwFC_MB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tbi_RT_sexbias2 == "male_biased")[,1])
Tce_RT_SBwFC_MB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tce_RT_sexbias2 == "male_biased")[,1])
Tcm_RT_SBwFC_MB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tcm_RT_sexbias2 == "male_biased")[,1])
Tpa_RT_SBwFC_MB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tpa_RT_sexbias2 == "male_biased")[,1])
Tps_RT_SBwFC_MB_genes <- as.character(subset(dat_ALL_RT_10sp_SB_asex_2, dat_ALL_RT_10sp_SB_asex_2$Tps_RT_sexbias2 == "male_biased")[,1])

Tbi_RT_SBwFC_SB_genes <- c(Tbi_RT_SBwFC_FB_genes, Tbi_RT_SBwFC_MB_genes )
Tce_RT_SBwFC_SB_genes <- c(Tce_RT_SBwFC_FB_genes, Tce_RT_SBwFC_MB_genes )
Tcm_RT_SBwFC_SB_genes <- c(Tcm_RT_SBwFC_FB_genes, Tcm_RT_SBwFC_MB_genes )
Tpa_RT_SBwFC_SB_genes <- c(Tpa_RT_SBwFC_FB_genes, Tpa_RT_SBwFC_MB_genes )
Tps_RT_SBwFC_SB_genes <- c(Tps_RT_SBwFC_FB_genes, Tps_RT_SBwFC_MB_genes )


Tbi_LG_SBwFC_FB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tbi_LG_sexbias2 == "female_biased")[,1])
Tce_LG_SBwFC_FB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tce_LG_sexbias2 == "female_biased")[,1])
Tcm_LG_SBwFC_FB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tcm_LG_sexbias2 == "female_biased")[,1])
Tpa_LG_SBwFC_FB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tpa_LG_sexbias2 == "female_biased")[,1])
Tps_LG_SBwFC_FB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tps_LG_sexbias2 == "female_biased")[,1])

Tbi_LG_SBwFC_MB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tbi_LG_sexbias2 == "male_biased")[,1])
Tce_LG_SBwFC_MB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tce_LG_sexbias2 == "male_biased")[,1])
Tcm_LG_SBwFC_MB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tcm_LG_sexbias2 == "male_biased")[,1])
Tpa_LG_SBwFC_MB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tpa_LG_sexbias2 == "male_biased")[,1])
Tps_LG_SBwFC_MB_genes <- as.character(subset(dat_ALL_LG_10sp_SB_asex_2, dat_ALL_LG_10sp_SB_asex_2$Tps_LG_sexbias2 == "male_biased")[,1])

Tbi_LG_SBwFC_SB_genes <- c(Tbi_LG_SBwFC_FB_genes, Tbi_LG_SBwFC_MB_genes )
Tce_LG_SBwFC_SB_genes <- c(Tce_LG_SBwFC_FB_genes, Tce_LG_SBwFC_MB_genes )
Tcm_LG_SBwFC_SB_genes <- c(Tcm_LG_SBwFC_FB_genes, Tcm_LG_SBwFC_MB_genes )
Tpa_LG_SBwFC_SB_genes <- c(Tpa_LG_SBwFC_FB_genes, Tpa_LG_SBwFC_MB_genes )
Tps_LG_SBwFC_SB_genes <- c(Tps_LG_SBwFC_FB_genes, Tps_LG_SBwFC_MB_genes )




#### output

venn_WB_SBwFC_SB <- five_sp_DE_venn(Tbi_WB_SBwFC_SB_genes, Tce_WB_SBwFC_SB_genes,Tcm_WB_SBwFC_SB_genes, Tpa_WB_SBwFC_SB_genes, Tps_WB_SBwFC_SB_genes , "WB_SBwFC_SB")
venn_RT_SBwFC_SB <- five_sp_DE_venn(Tbi_RT_SBwFC_SB_genes, Tce_RT_SBwFC_SB_genes,Tcm_RT_SBwFC_SB_genes, Tpa_RT_SBwFC_SB_genes, Tps_RT_SBwFC_SB_genes , "RT_SBwFC_SB")
venn_LG_SBwFC_SB <- five_sp_DE_venn(Tbi_LG_SBwFC_SB_genes, Tce_LG_SBwFC_SB_genes,Tcm_LG_SBwFC_SB_genes, Tpa_LG_SBwFC_SB_genes, Tps_LG_SBwFC_SB_genes , "LG_SBwFC_SB")

venn_WB_SBwFC_FB <- five_sp_DE_venn(Tbi_WB_SBwFC_FB_genes, Tce_WB_SBwFC_FB_genes,Tcm_WB_SBwFC_FB_genes, Tpa_WB_SBwFC_FB_genes, Tps_WB_SBwFC_FB_genes , "WB_SBwFC_FB")
venn_RT_SBwFC_FB <- five_sp_DE_venn(Tbi_RT_SBwFC_FB_genes, Tce_RT_SBwFC_FB_genes,Tcm_RT_SBwFC_FB_genes, Tpa_RT_SBwFC_FB_genes, Tps_RT_SBwFC_FB_genes , "RT_SBwFC_FB")
venn_LG_SBwFC_FB <- five_sp_DE_venn(Tbi_LG_SBwFC_FB_genes, Tce_LG_SBwFC_FB_genes,Tcm_LG_SBwFC_FB_genes, Tpa_LG_SBwFC_FB_genes, Tps_LG_SBwFC_FB_genes , "LG_SBwFC_FB")

venn_WB_SBwFC_MB <- five_sp_DE_venn(Tbi_WB_SBwFC_MB_genes, Tce_WB_SBwFC_MB_genes,Tcm_WB_SBwFC_MB_genes, Tpa_WB_SBwFC_MB_genes, Tps_WB_SBwFC_MB_genes , "WB_SBwFC_MB")
venn_RT_SBwFC_MB <- five_sp_DE_venn(Tbi_RT_SBwFC_MB_genes, Tce_RT_SBwFC_MB_genes,Tcm_RT_SBwFC_MB_genes, Tpa_RT_SBwFC_MB_genes, Tps_RT_SBwFC_MB_genes , "RT_SBwFC_MB")
venn_LG_SBwFC_MB <- five_sp_DE_venn(Tbi_LG_SBwFC_MB_genes, Tce_LG_SBwFC_MB_genes,Tcm_LG_SBwFC_MB_genes, Tpa_LG_SBwFC_MB_genes, Tps_LG_SBwFC_MB_genes , "LG_SBwFC_MB")

venn_WB_SBwFC_FB_sc <- five_sp_DE_venn_single_col_F(Tbi_WB_SBwFC_FB_genes, Tce_WB_SBwFC_FB_genes,Tcm_WB_SBwFC_FB_genes, Tpa_WB_SBwFC_FB_genes, Tps_WB_SBwFC_FB_genes , "WB_SBwFC_FB")
venn_RT_SBwFC_FB_sc <- five_sp_DE_venn_single_col_F(Tbi_RT_SBwFC_FB_genes, Tce_RT_SBwFC_FB_genes,Tcm_RT_SBwFC_FB_genes, Tpa_RT_SBwFC_FB_genes, Tps_RT_SBwFC_FB_genes , "RT_SBwFC_FB")
venn_LG_SBwFC_FB_sc <- five_sp_DE_venn_single_col_F(Tbi_LG_SBwFC_FB_genes, Tce_LG_SBwFC_FB_genes,Tcm_LG_SBwFC_FB_genes, Tpa_LG_SBwFC_FB_genes, Tps_LG_SBwFC_FB_genes , "LG_SBwFC_FB")

venn_WB_SBwFC_MB_sc <- five_sp_DE_venn_single_col_M(Tbi_WB_SBwFC_MB_genes, Tce_WB_SBwFC_MB_genes,Tcm_WB_SBwFC_MB_genes, Tpa_WB_SBwFC_MB_genes, Tps_WB_SBwFC_MB_genes , "WB_SBwFC_MB")
venn_RT_SBwFC_MB_sc <- five_sp_DE_venn_single_col_M(Tbi_RT_SBwFC_MB_genes, Tce_RT_SBwFC_MB_genes,Tcm_RT_SBwFC_MB_genes, Tpa_RT_SBwFC_MB_genes, Tps_RT_SBwFC_MB_genes , "RT_SBwFC_MB")
venn_LG_SBwFC_MB_sc <- five_sp_DE_venn_single_col_M(Tbi_LG_SBwFC_MB_genes, Tce_LG_SBwFC_MB_genes,Tcm_LG_SBwFC_MB_genes, Tpa_LG_SBwFC_MB_genes, Tps_LG_SBwFC_MB_genes , "LG_SBwFC_MB")




########

venn_WBRTLG_SBwFC_MB_sc <- arrangeGrob(
gTree(children=venn_WB_SBwFC_MB_sc), 
gTree(children=venn_RT_SBwFC_MB_sc),
gTree(children=venn_LG_SBwFC_MB_sc),
ncol = 3 , nrow = 1)

ggsave(file="venn_WBRTLG_SBwFC_MB_sc.png", venn_WBRTLG_SBwFC_MB_sc, width = 18, height = 6)


venn_WBRTLG_SBwFC_FB_sc <- arrangeGrob(
gTree(children=venn_WB_SBwFC_FB_sc), 
gTree(children=venn_RT_SBwFC_FB_sc),
gTree(children=venn_LG_SBwFC_FB_sc),
ncol = 3 , nrow = 1)

ggsave(file="venn_WBRTLG_SBwFC_FB_sc.png", venn_WBRTLG_SBwFC_FB_sc, width = 18, height = 6)


##################### overlaps

WB_SB_SB_genes_list <- list("Tbi" = Tbi_WB_SB_SB_genes, "Tce" = Tce_WB_SB_SB_genes, "Tcm" = Tcm_WB_SB_SB_genes, "Tpa" = Tpa_WB_SB_SB_genes, "Tps" = Tps_WB_SB_SB_genes )
RT_SB_SB_genes_list <- list("Tbi" = Tbi_RT_SB_SB_genes, "Tce" = Tce_RT_SB_SB_genes, "Tcm" = Tcm_RT_SB_SB_genes, "Tpa" = Tpa_RT_SB_SB_genes, "Tps" = Tps_RT_SB_SB_genes )
LG_SB_SB_genes_list <- list("Tbi" = Tbi_LG_SB_SB_genes, "Tce" = Tce_LG_SB_SB_genes, "Tcm" = Tcm_LG_SB_SB_genes, "Tpa" = Tpa_LG_SB_SB_genes, "Tps" = Tps_LG_SB_SB_genes )

WB_SB_SB_genes_FET <- supertest(WB_SB_SB_genes_list , n= length(dat_ALL_WB_10sp_SB_asex_2[,1]))
write.csv(summary(WB_SB_SB_genes_FET)$Table, file="WB_SB_SB_genes_FET_summary.table.csv", row.names=FALSE)
RT_SB_SB_genes_FET <- supertest(RT_SB_SB_genes_list , n= length(dat_ALL_RT_10sp_SB_asex_2[,1]))
write.csv(summary(RT_SB_SB_genes_FET)$Table, file="RT_SB_SB_genes_FET_summary.table.csv", row.names=FALSE)
LG_SB_SB_genes_FET <- supertest(LG_SB_SB_genes_list , n= length(dat_ALL_LG_10sp_SB_asex_2[,1]))
write.csv(summary(LG_SB_SB_genes_FET)$Table, file="LG_SB_SB_genes_FET_summary.table.csv", row.names=FALSE)

WB_SBwFC_SB_genes_list <- list("Tbi" = Tbi_WB_SBwFC_SB_genes, "Tce" = Tce_WB_SBwFC_SB_genes, "Tcm" = Tcm_WB_SBwFC_SB_genes, "Tpa" = Tpa_WB_SBwFC_SB_genes, "Tps" = Tps_WB_SBwFC_SB_genes )
RT_SBwFC_SB_genes_list <- list("Tbi" = Tbi_RT_SBwFC_SB_genes, "Tce" = Tce_RT_SBwFC_SB_genes, "Tcm" = Tcm_RT_SBwFC_SB_genes, "Tpa" = Tpa_RT_SBwFC_SB_genes, "Tps" = Tps_RT_SBwFC_SB_genes )
LG_SBwFC_SB_genes_list <- list("Tbi" = Tbi_LG_SBwFC_SB_genes, "Tce" = Tce_LG_SBwFC_SB_genes, "Tcm" = Tcm_LG_SBwFC_SB_genes, "Tpa" = Tpa_LG_SBwFC_SB_genes, "Tps" = Tps_LG_SBwFC_SB_genes )

WB_SBwFC_SB_genes_FET <- supertest(WB_SBwFC_SB_genes_list , n= length(dat_ALL_WB_10sp_SB_asex_2[,1]))
write.csv(summary(WB_SBwFC_SB_genes_FET)$Table, file="WB_SBwFC_SB_genes_FET_summary.table.csv", row.names=FALSE)
RT_SBwFC_SB_genes_FET <- supertest(RT_SBwFC_SB_genes_list , n= length(dat_ALL_RT_10sp_SB_asex_2[,1]))
write.csv(summary(RT_SBwFC_SB_genes_FET)$Table, file="RT_SBwFC_SB_genes_FET_summary.table.csv", row.names=FALSE)
LG_SBwFC_SB_genes_FET <- supertest(LG_SBwFC_SB_genes_list , n= length(dat_ALL_LG_10sp_SB_asex_2[,1]))
write.csv(summary(LG_SBwFC_SB_genes_FET)$Table, file="LG_SBwFC_SB_genes_FET_summary.table.csv", row.names=FALSE)



WB_SB_FB_genes_list <- list("Tbi" = Tbi_WB_SB_FB_genes, "Tce" = Tce_WB_SB_FB_genes, "Tcm" = Tcm_WB_SB_FB_genes, "Tpa" = Tpa_WB_SB_FB_genes, "Tps" = Tps_WB_SB_FB_genes )
RT_SB_FB_genes_list <- list("Tbi" = Tbi_RT_SB_FB_genes, "Tce" = Tce_RT_SB_FB_genes, "Tcm" = Tcm_RT_SB_FB_genes, "Tpa" = Tpa_RT_SB_FB_genes, "Tps" = Tps_RT_SB_FB_genes )
LG_SB_FB_genes_list <- list("Tbi" = Tbi_LG_SB_FB_genes, "Tce" = Tce_LG_SB_FB_genes, "Tcm" = Tcm_LG_SB_FB_genes, "Tpa" = Tpa_LG_SB_FB_genes, "Tps" = Tps_LG_SB_FB_genes )

WB_SB_FB_genes_FET <- supertest(WB_SB_FB_genes_list , n= length(dat_ALL_WB_10sp_SB_asex_2[,1]))
write.csv(summary(WB_SB_FB_genes_FET)$Table, file="WB_SB_FB_genes_FET_summary.table.csv", row.names=FALSE)
RT_SB_FB_genes_FET <- supertest(RT_SB_FB_genes_list , n= length(dat_ALL_RT_10sp_SB_asex_2[,1]))
write.csv(summary(RT_SB_FB_genes_FET)$Table, file="RT_SB_FB_genes_FET_summary.table.csv", row.names=FALSE)
LG_SB_FB_genes_FET <- supertest(LG_SB_FB_genes_list , n= length(dat_ALL_LG_10sp_SB_asex_2[,1]))
write.csv(summary(LG_SB_FB_genes_FET)$Table, file="LG_SB_FB_genes_FET_summary.table.csv", row.names=FALSE)

WB_SBwFC_FB_genes_list <- list("Tbi" = Tbi_WB_SBwFC_FB_genes, "Tce" = Tce_WB_SBwFC_FB_genes, "Tcm" = Tcm_WB_SBwFC_FB_genes, "Tpa" = Tpa_WB_SBwFC_FB_genes, "Tps" = Tps_WB_SBwFC_FB_genes )
RT_SBwFC_FB_genes_list <- list("Tbi" = Tbi_RT_SBwFC_FB_genes, "Tce" = Tce_RT_SBwFC_FB_genes, "Tcm" = Tcm_RT_SBwFC_FB_genes, "Tpa" = Tpa_RT_SBwFC_FB_genes, "Tps" = Tps_RT_SBwFC_FB_genes )
LG_SBwFC_FB_genes_list <- list("Tbi" = Tbi_LG_SBwFC_FB_genes, "Tce" = Tce_LG_SBwFC_FB_genes, "Tcm" = Tcm_LG_SBwFC_FB_genes, "Tpa" = Tpa_LG_SBwFC_FB_genes, "Tps" = Tps_LG_SBwFC_FB_genes )

WB_SBwFC_FB_genes_FET <- supertest(WB_SBwFC_FB_genes_list , n= length(dat_ALL_WB_10sp_SB_asex_2[,1]))
write.csv(summary(WB_SBwFC_FB_genes_FET)$Table, file="WB_SBwFC_FB_genes_FET_summary.table.csv", row.names=FALSE)
RT_SBwFC_FB_genes_FET <- supertest(RT_SBwFC_FB_genes_list , n= length(dat_ALL_RT_10sp_SB_asex_2[,1]))
write.csv(summary(RT_SBwFC_FB_genes_FET)$Table, file="RT_SBwFC_FB_genes_FET_summary.table.csv", row.names=FALSE)
LG_SBwFC_FB_genes_FET <- supertest(LG_SBwFC_FB_genes_list , n= length(dat_ALL_LG_10sp_SB_asex_2[,1]))
write.csv(summary(LG_SBwFC_FB_genes_FET)$Table, file="LG_SBwFC_FB_genes_FET_summary.table.csv", row.names=FALSE)


WB_SB_MB_genes_list <- list("Tbi" = Tbi_WB_SB_MB_genes, "Tce" = Tce_WB_SB_MB_genes, "Tcm" = Tcm_WB_SB_MB_genes, "Tpa" = Tpa_WB_SB_MB_genes, "Tps" = Tps_WB_SB_MB_genes )
RT_SB_MB_genes_list <- list("Tbi" = Tbi_RT_SB_MB_genes, "Tce" = Tce_RT_SB_MB_genes, "Tcm" = Tcm_RT_SB_MB_genes, "Tpa" = Tpa_RT_SB_MB_genes, "Tps" = Tps_RT_SB_MB_genes )
LG_SB_MB_genes_list <- list("Tbi" = Tbi_LG_SB_MB_genes, "Tce" = Tce_LG_SB_MB_genes, "Tcm" = Tcm_LG_SB_MB_genes, "Tpa" = Tpa_LG_SB_MB_genes, "Tps" = Tps_LG_SB_MB_genes )

WB_SB_MB_genes_FET <- supertest(WB_SB_MB_genes_list , n= length(dat_ALL_WB_10sp_SB_asex_2[,1]))
write.csv(summary(WB_SB_MB_genes_FET)$Table, file="WB_SB_MB_genes_FET_summary.table.csv", row.names=FALSE)
RT_SB_MB_genes_FET <- supertest(RT_SB_MB_genes_list , n= length(dat_ALL_RT_10sp_SB_asex_2[,1]))
write.csv(summary(RT_SB_MB_genes_FET)$Table, file="RT_SB_MB_genes_FET_summary.table.csv", row.names=FALSE)
LG_SB_MB_genes_FET <- supertest(LG_SB_MB_genes_list , n= length(dat_ALL_LG_10sp_SB_asex_2[,1]))
write.csv(summary(LG_SB_MB_genes_FET)$Table, file="LG_SB_MB_genes_FET_summary.table.csv", row.names=FALSE)

WB_SBwFC_MB_genes_list <- list("Tbi" = Tbi_WB_SBwFC_MB_genes, "Tce" = Tce_WB_SBwFC_MB_genes, "Tcm" = Tcm_WB_SBwFC_MB_genes, "Tpa" = Tpa_WB_SBwFC_MB_genes, "Tps" = Tps_WB_SBwFC_MB_genes )
RT_SBwFC_MB_genes_list <- list("Tbi" = Tbi_RT_SBwFC_MB_genes, "Tce" = Tce_RT_SBwFC_MB_genes, "Tcm" = Tcm_RT_SBwFC_MB_genes, "Tpa" = Tpa_RT_SBwFC_MB_genes, "Tps" = Tps_RT_SBwFC_MB_genes )
LG_SBwFC_MB_genes_list <- list("Tbi" = Tbi_LG_SBwFC_MB_genes, "Tce" = Tce_LG_SBwFC_MB_genes, "Tcm" = Tcm_LG_SBwFC_MB_genes, "Tpa" = Tpa_LG_SBwFC_MB_genes, "Tps" = Tps_LG_SBwFC_MB_genes )

WB_SBwFC_MB_genes_FET <- supertest(WB_SBwFC_MB_genes_list , n= length(dat_ALL_WB_10sp_SB_asex_2[,1]))
write.csv(summary(WB_SBwFC_MB_genes_FET)$Table, file="WB_SBwFC_MB_genes_FET_summary.table.csv", row.names=FALSE)
RT_SBwFC_MB_genes_FET <- supertest(RT_SBwFC_MB_genes_list , n= length(dat_ALL_RT_10sp_SB_asex_2[,1]))
write.csv(summary(RT_SBwFC_MB_genes_FET)$Table, file="RT_SBwFC_MB_genes_FET_summary.table.csv", row.names=FALSE)
LG_SBwFC_MB_genes_FET <- supertest(LG_SBwFC_MB_genes_list , n= length(dat_ALL_LG_10sp_SB_asex_2[,1]))
write.csv(summary(LG_SBwFC_MB_genes_FET)$Table, file="LG_SBwFC_MB_genes_FET_summary.table.csv", row.names=FALSE)





########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), "sex_bias_plotsetc_10sp_sessionInfo.txt")

