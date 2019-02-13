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

setwd("Output/DE_joined_sexual_ref")

dat_Tbi_WB_longest_iso_SB_asex = read.csv("Tbi_WB_longest_iso_disp_allsepar_SB_asex.csv")
dat_Tbi_RT_longest_iso_SB_asex = read.csv("Tbi_RT_longest_iso_disp_allsepar_SB_asex.csv")
dat_Tbi_LG_longest_iso_SB_asex = read.csv("Tbi_LG_longest_iso_disp_allsepar_SB_asex.csv")

dat_Tce_WB_longest_iso_SB_asex = read.csv("Tce_WB_longest_iso_disp_allsepar_SB_asex.csv")
dat_Tce_RT_longest_iso_SB_asex = read.csv("Tce_RT_longest_iso_disp_allsepar_SB_asex.csv")
dat_Tce_LG_longest_iso_SB_asex = read.csv("Tce_LG_longest_iso_disp_allsepar_SB_asex.csv")

dat_Tcm_WB_longest_iso_SB_asex = read.csv("Tcm_WB_longest_iso_disp_allsepar_SB_asex.csv")
dat_Tcm_RT_longest_iso_SB_asex = read.csv("Tcm_RT_longest_iso_disp_allsepar_SB_asex.csv")
dat_Tcm_LG_longest_iso_SB_asex = read.csv("Tcm_LG_longest_iso_disp_allsepar_SB_asex.csv")

dat_Tpa_WB_longest_iso_SB_asex = read.csv("Tpa_WB_longest_iso_disp_allsepar_SB_asex.csv")
dat_Tpa_RT_longest_iso_SB_asex = read.csv("Tpa_RT_longest_iso_disp_allsepar_SB_asex.csv")
dat_Tpa_LG_longest_iso_SB_asex = read.csv("Tpa_LG_longest_iso_disp_allsepar_SB_asex.csv")

dat_Tps_WB_longest_iso_SB_asex = read.csv("Tps_WB_longest_iso_disp_allsepar_SB_asex.csv")
dat_Tps_RT_longest_iso_SB_asex = read.csv("Tps_RT_longest_iso_disp_allsepar_SB_asex.csv")
dat_Tps_LG_longest_iso_SB_asex = read.csv("Tps_LG_longest_iso_disp_allsepar_SB_asex.csv")


head(dat_Tbi_WB_longest_iso_SB_asex, n = 20)

#################################################################################################################################################
###### refine sex bias by FC (not log FC)
	
FC_cutoff = 2

dat_Tbi_WB_longest_iso_SB_asex$Tbi_WB_sexbias2 <- ifelse(2^(sqrt(dat_Tbi_WB_longest_iso_SB_asex$Tbi_WB_log2FC_SB * dat_Tbi_WB_longest_iso_SB_asex$Tbi_WB_log2FC_SB)) > FC_cutoff , as.character(dat_Tbi_WB_longest_iso_SB_asex$Tbi_WB_sexbias), "Unbiased")
dat_Tce_WB_longest_iso_SB_asex$Tce_WB_sexbias2 <- ifelse(2^(sqrt(dat_Tce_WB_longest_iso_SB_asex$Tce_WB_log2FC_SB * dat_Tce_WB_longest_iso_SB_asex$Tce_WB_log2FC_SB)) > FC_cutoff , as.character(dat_Tce_WB_longest_iso_SB_asex$Tce_WB_sexbias), "Unbiased")
dat_Tcm_WB_longest_iso_SB_asex$Tcm_WB_sexbias2 <- ifelse(2^(sqrt(dat_Tcm_WB_longest_iso_SB_asex$Tcm_WB_log2FC_SB * dat_Tcm_WB_longest_iso_SB_asex$Tcm_WB_log2FC_SB)) > FC_cutoff , as.character(dat_Tcm_WB_longest_iso_SB_asex$Tcm_WB_sexbias), "Unbiased")
dat_Tpa_WB_longest_iso_SB_asex$Tpa_WB_sexbias2 <- ifelse(2^(sqrt(dat_Tpa_WB_longest_iso_SB_asex$Tpa_WB_log2FC_SB * dat_Tpa_WB_longest_iso_SB_asex$Tpa_WB_log2FC_SB)) > FC_cutoff , as.character(dat_Tpa_WB_longest_iso_SB_asex$Tpa_WB_sexbias), "Unbiased")
dat_Tps_WB_longest_iso_SB_asex$Tps_WB_sexbias2 <- ifelse(2^(sqrt(dat_Tps_WB_longest_iso_SB_asex$Tps_WB_log2FC_SB * dat_Tps_WB_longest_iso_SB_asex$Tps_WB_log2FC_SB)) > FC_cutoff , as.character(dat_Tps_WB_longest_iso_SB_asex$Tps_WB_sexbias), "Unbiased")

dat_Tbi_RT_longest_iso_SB_asex$Tbi_RT_sexbias2 <- ifelse(2^(sqrt(dat_Tbi_RT_longest_iso_SB_asex$Tbi_RT_log2FC_SB * dat_Tbi_RT_longest_iso_SB_asex$Tbi_RT_log2FC_SB)) > FC_cutoff , as.character(dat_Tbi_RT_longest_iso_SB_asex$Tbi_RT_sexbias), "Unbiased")
dat_Tce_RT_longest_iso_SB_asex$Tce_RT_sexbias2 <- ifelse(2^(sqrt(dat_Tce_RT_longest_iso_SB_asex$Tce_RT_log2FC_SB * dat_Tce_RT_longest_iso_SB_asex$Tce_RT_log2FC_SB)) > FC_cutoff , as.character(dat_Tce_RT_longest_iso_SB_asex$Tce_RT_sexbias), "Unbiased")
dat_Tcm_RT_longest_iso_SB_asex$Tcm_RT_sexbias2 <- ifelse(2^(sqrt(dat_Tcm_RT_longest_iso_SB_asex$Tcm_RT_log2FC_SB * dat_Tcm_RT_longest_iso_SB_asex$Tcm_RT_log2FC_SB)) > FC_cutoff , as.character(dat_Tcm_RT_longest_iso_SB_asex$Tcm_RT_sexbias), "Unbiased")
dat_Tpa_RT_longest_iso_SB_asex$Tpa_RT_sexbias2 <- ifelse(2^(sqrt(dat_Tpa_RT_longest_iso_SB_asex$Tpa_RT_log2FC_SB * dat_Tpa_RT_longest_iso_SB_asex$Tpa_RT_log2FC_SB)) > FC_cutoff , as.character(dat_Tpa_RT_longest_iso_SB_asex$Tpa_RT_sexbias), "Unbiased")
dat_Tps_RT_longest_iso_SB_asex$Tps_RT_sexbias2 <- ifelse(2^(sqrt(dat_Tps_RT_longest_iso_SB_asex$Tps_RT_log2FC_SB * dat_Tps_RT_longest_iso_SB_asex$Tps_RT_log2FC_SB)) > FC_cutoff , as.character(dat_Tps_RT_longest_iso_SB_asex$Tps_RT_sexbias), "Unbiased")

dat_Tbi_LG_longest_iso_SB_asex$Tbi_LG_sexbias2 <- ifelse(2^(sqrt(dat_Tbi_LG_longest_iso_SB_asex$Tbi_LG_log2FC_SB * dat_Tbi_LG_longest_iso_SB_asex$Tbi_LG_log2FC_SB)) > FC_cutoff , as.character(dat_Tbi_LG_longest_iso_SB_asex$Tbi_LG_sexbias), "Unbiased")
dat_Tce_LG_longest_iso_SB_asex$Tce_LG_sexbias2 <- ifelse(2^(sqrt(dat_Tce_LG_longest_iso_SB_asex$Tce_LG_log2FC_SB * dat_Tce_LG_longest_iso_SB_asex$Tce_LG_log2FC_SB)) > FC_cutoff , as.character(dat_Tce_LG_longest_iso_SB_asex$Tce_LG_sexbias), "Unbiased")
dat_Tcm_LG_longest_iso_SB_asex$Tcm_LG_sexbias2 <- ifelse(2^(sqrt(dat_Tcm_LG_longest_iso_SB_asex$Tcm_LG_log2FC_SB * dat_Tcm_LG_longest_iso_SB_asex$Tcm_LG_log2FC_SB)) > FC_cutoff , as.character(dat_Tcm_LG_longest_iso_SB_asex$Tcm_LG_sexbias), "Unbiased")
dat_Tpa_LG_longest_iso_SB_asex$Tpa_LG_sexbias2 <- ifelse(2^(sqrt(dat_Tpa_LG_longest_iso_SB_asex$Tpa_LG_log2FC_SB * dat_Tpa_LG_longest_iso_SB_asex$Tpa_LG_log2FC_SB)) > FC_cutoff , as.character(dat_Tpa_LG_longest_iso_SB_asex$Tpa_LG_sexbias), "Unbiased")
dat_Tps_LG_longest_iso_SB_asex$Tps_LG_sexbias2 <- ifelse(2^(sqrt(dat_Tps_LG_longest_iso_SB_asex$Tps_LG_log2FC_SB * dat_Tps_LG_longest_iso_SB_asex$Tps_LG_log2FC_SB)) > FC_cutoff , as.character(dat_Tps_LG_longest_iso_SB_asex$Tps_LG_sexbias), "Unbiased")

head(dat_Tps_RT_longest_iso_SB_asex)





######################## get excluded genes

N_SB_Genes_kept_excluded <- function(df, df_name, sp, tiss){
	
	asex_sp <- ifelse(sp == "Tbi", "Tte", 
	           ifelse(sp == "Tce", "Tms", 
	           ifelse(sp == "Tcm", "Tsi", 
	           ifelse(sp == "Tpa", "Tge", 
               ifelse(sp == "Tps", "Tdi","FAIL")))))
	
	print(asex_sp)
	
	df_excl <- subset(df, is.na(eval(parse(text=paste(df_name,"$",sp,"_",tiss,"_log2FC_SA", sep=""))) == TRUE))

	
	N_MB_excl = length(subset(df_excl, eval(parse(text=paste("df_excl","$",sp,"_",tiss,"_sexbias2", sep=""))) == "male_biased")[,1])
	N_FB_excl = length(subset(df_excl, eval(parse(text=paste("df_excl","$",sp,"_",tiss,"_sexbias2", sep=""))) == "female_biased")[,1])
	N_UB_excl = length(subset(df_excl, eval(parse(text=paste("df_excl","$",sp,"_",tiss,"_sexbias2", sep=""))) == "Unbiased")[,1])
	total_excl = N_MB_excl + N_FB_excl + N_UB_excl
	
	print(N_MB_excl )
	print(N_FB_excl )	
	print(N_UB_excl )
	
	df_kept = na.omit(df)
	
	N_MB_kept = length(subset(df_kept, eval(parse(text=paste("df_kept","$",sp,"_",tiss,"_sexbias2", sep=""))) == "male_biased")[,1])
	N_FB_kept = length(subset(df_kept, eval(parse(text=paste("df_kept","$",sp,"_",tiss,"_sexbias2", sep=""))) == "female_biased")[,1])
	N_UB_kept = length(subset(df_kept, eval(parse(text=paste("df_kept","$",sp,"_",tiss,"_sexbias2", sep=""))) == "Unbiased")[,1])	
	total_kept = N_MB_kept + N_FB_kept + N_UB_kept
		
	print(N_MB_kept )
	print(N_FB_kept )	
	print(N_UB_kept )	
	
	
	FB_FT_mat <- matrix(c(N_FB_excl,(total_excl-N_FB_excl),N_FB_kept,(total_kept-N_FB_kept)), nrow = 2)
	MB_FT_mat <- matrix(c(N_MB_excl,(total_excl-N_MB_excl),N_MB_kept,(total_kept-N_MB_kept)), nrow = 2)
	FB_fisher <- fisher.test(FB_FT_mat,alternative="two.sided") ### enrichment (over or under-rep)
	MB_fisher <- fisher.test(MB_FT_mat,alternative="two.sided") ### enrichment (over or under-rep)
	
	
	print(FB_FT_mat)
	print(MB_FT_mat)	
	print(FB_fisher)
	print(MB_fisher)
			
	out_df = as.data.frame(cbind(
	c(N_MB_excl, N_FB_excl, N_UB_excl, N_MB_kept, N_FB_kept, N_UB_kept),
	c(total_excl,total_excl,total_excl,total_kept,total_kept,total_kept),
	c("MB_excl", "FB_excl", "UB_excl", "MB_kept", "FB_kept", "UB_kept"),
	c("MB", "FB", "UB", "MB", "FB", "UB"),
	rep(sp, 6),
	rep(tiss,6)
	))
	
	colnames(out_df) <- c("Ngenes", "total_genes", "gene_class_2", "bias","sp", "tiss")	
	
	
	out_df_2 = as.data.frame(cbind(
	c(total_excl,total_kept),
	c("excl", "kept"),
	rep(sp, 2),
	rep(tiss,2)
	))



	out_df3 <- as.data.frame(cbind(sp, tiss, N_MB_excl, N_FB_excl, total_excl, N_MB_kept, N_FB_kept, total_kept, FB_fisher$estimate,  MB_fisher$estimate,  FB_fisher$p.value,  MB_fisher$p.value))
	colnames(out_df3 ) <- c("sp", "tiss", "MB_excl", "FB_excl", "total_excl", "MB_kept", "FB_kept", "total_kept", "FB_fisher_OR",  "MB_fisher_OR",  "FB_fisher_p",  "MB_fisher_p")


	colnames(out_df_2) <- c("N_genes", "kept_or_excl","sp", "tiss")
	
	
	
	### calc mean FKPMS
	
	
	df_excl$SM_avgFPKM = (eval(parse(text=paste("df_excl","$FPKM_",sp,"_SM_",tiss,"_Md_Re1", sep=""))) + 
	                      eval(parse(text=paste("df_excl","$FPKM_",sp,"_SM_",tiss,"_Md_Re2", sep=""))) + 
	                      eval(parse(text=paste("df_excl","$FPKM_",sp,"_SM_",tiss,"_Md_Re3", sep="")))) / 3 
	
	df_excl$SF_avgFPKM = (eval(parse(text=paste("df_excl","$FPKM_",sp,"_SF_",tiss,"_Md_Re1", sep=""))) + 
	                      eval(parse(text=paste("df_excl","$FPKM_",sp,"_SF_",tiss,"_Md_Re2", sep=""))) + 
	                      eval(parse(text=paste("df_excl","$FPKM_",sp,"_SF_",tiss,"_Md_Re3", sep="")))) / 3 	
	
	df_excl$AF_avgFPKM = (eval(parse(text=paste("df_excl","$FPKM_",asex_sp,"_AF_",tiss,"_Vi_Re1", sep=""))) + 
	                      eval(parse(text=paste("df_excl","$FPKM_",asex_sp,"_AF_",tiss,"_Vi_Re2", sep=""))) + 
	                      eval(parse(text=paste("df_excl","$FPKM_",asex_sp,"_AF_",tiss,"_Vi_Re3", sep="")))) / 3 	
	     
	     
	
	## stick together
	
	
	excl_meanFPKM_df = as.data.frame(cbind(
	c(as.character(df_excl$genename), as.character(df_excl$genename), as.character(df_excl$genename)),
	rep(sp, length(df_excl[,1]) * 3),
	rep(tiss, length(df_excl[,1]) * 3),
	c(rep("SM", length(df_excl[,1])), rep("SF", length(df_excl[,1])),rep("AF", length(df_excl[,1]))),  
	c(eval(parse(text=paste("df_excl","$",sp,"_",tiss,"_sexbias2", sep=""))), eval(parse(text=paste("df_excl","$",sp,"_",tiss,"_sexbias2", sep=""))), eval(parse(text=paste("df_excl","$",sp,"_",tiss,"_sexbias2", sep="")))),
	c(df_excl$SM_avgFPKM, df_excl$SF_avgFPKM, df_excl$AF_avgFPKM)     
	))     
	     
	colnames(excl_meanFPKM_df) <- c("genename", "sp", "tiss", "sex", "sexbias2","avgFPKM")                 
	return(list("out_df" = out_df, "out_df_2" = out_df_2, "out_df3" = out_df3, "excl_df" = df_excl, "excl_meanFPKM_df" = excl_meanFPKM_df ))
		
}


### N kept

N_kept_excl_df <- rbind(
N_SB_Genes_kept_excluded(dat_Tbi_WB_longest_iso_SB_asex, "dat_Tbi_WB_longest_iso_SB_asex", "Tbi", "WB")$out_df_2,
N_SB_Genes_kept_excluded(dat_Tce_WB_longest_iso_SB_asex, "dat_Tce_WB_longest_iso_SB_asex", "Tce", "WB")$out_df_2,
N_SB_Genes_kept_excluded(dat_Tcm_WB_longest_iso_SB_asex, "dat_Tcm_WB_longest_iso_SB_asex", "Tcm", "WB")$out_df_2,
N_SB_Genes_kept_excluded(dat_Tpa_WB_longest_iso_SB_asex, "dat_Tpa_WB_longest_iso_SB_asex", "Tpa", "WB")$out_df_2,
N_SB_Genes_kept_excluded(dat_Tps_WB_longest_iso_SB_asex, "dat_Tps_WB_longest_iso_SB_asex", "Tps", "WB")$out_df_2,
N_SB_Genes_kept_excluded(dat_Tbi_RT_longest_iso_SB_asex, "dat_Tbi_RT_longest_iso_SB_asex", "Tbi", "RT")$out_df_2,
N_SB_Genes_kept_excluded(dat_Tce_RT_longest_iso_SB_asex, "dat_Tce_RT_longest_iso_SB_asex", "Tce", "RT")$out_df_2,
N_SB_Genes_kept_excluded(dat_Tcm_RT_longest_iso_SB_asex, "dat_Tcm_RT_longest_iso_SB_asex", "Tcm", "RT")$out_df_2,
N_SB_Genes_kept_excluded(dat_Tpa_RT_longest_iso_SB_asex, "dat_Tpa_RT_longest_iso_SB_asex", "Tpa", "RT")$out_df_2,
N_SB_Genes_kept_excluded(dat_Tps_RT_longest_iso_SB_asex, "dat_Tps_RT_longest_iso_SB_asex", "Tps", "RT")$out_df_2,
N_SB_Genes_kept_excluded(dat_Tbi_LG_longest_iso_SB_asex, "dat_Tbi_LG_longest_iso_SB_asex", "Tbi", "LG")$out_df_2,
N_SB_Genes_kept_excluded(dat_Tce_LG_longest_iso_SB_asex, "dat_Tce_LG_longest_iso_SB_asex", "Tce", "LG")$out_df_2,
N_SB_Genes_kept_excluded(dat_Tcm_LG_longest_iso_SB_asex, "dat_Tcm_LG_longest_iso_SB_asex", "Tcm", "LG")$out_df_2,
N_SB_Genes_kept_excluded(dat_Tpa_LG_longest_iso_SB_asex, "dat_Tpa_LG_longest_iso_SB_asex", "Tpa", "LG")$out_df_2,
N_SB_Genes_kept_excluded(dat_Tps_LG_longest_iso_SB_asex, "dat_Tps_LG_longest_iso_SB_asex", "Tps", "LG")$out_df_2)


write.table(N_kept_excl_df, "N_kept_excl_df_longest_iso.csv", sep = ',', quote = FALSE, row.names = FALSE)


### N kept by SB 
N_kept_excl_SB_df <- rbind(
N_SB_Genes_kept_excluded(dat_Tbi_WB_longest_iso_SB_asex, "dat_Tbi_WB_longest_iso_SB_asex", "Tbi", "WB")$out_df,
N_SB_Genes_kept_excluded(dat_Tce_WB_longest_iso_SB_asex, "dat_Tce_WB_longest_iso_SB_asex", "Tce", "WB")$out_df,
N_SB_Genes_kept_excluded(dat_Tcm_WB_longest_iso_SB_asex, "dat_Tcm_WB_longest_iso_SB_asex", "Tcm", "WB")$out_df,
N_SB_Genes_kept_excluded(dat_Tpa_WB_longest_iso_SB_asex, "dat_Tpa_WB_longest_iso_SB_asex", "Tpa", "WB")$out_df,
N_SB_Genes_kept_excluded(dat_Tps_WB_longest_iso_SB_asex, "dat_Tps_WB_longest_iso_SB_asex", "Tps", "WB")$out_df,
N_SB_Genes_kept_excluded(dat_Tbi_RT_longest_iso_SB_asex, "dat_Tbi_RT_longest_iso_SB_asex", "Tbi", "RT")$out_df,
N_SB_Genes_kept_excluded(dat_Tce_RT_longest_iso_SB_asex, "dat_Tce_RT_longest_iso_SB_asex", "Tce", "RT")$out_df,
N_SB_Genes_kept_excluded(dat_Tcm_RT_longest_iso_SB_asex, "dat_Tcm_RT_longest_iso_SB_asex", "Tcm", "RT")$out_df,
N_SB_Genes_kept_excluded(dat_Tpa_RT_longest_iso_SB_asex, "dat_Tpa_RT_longest_iso_SB_asex", "Tpa", "RT")$out_df,
N_SB_Genes_kept_excluded(dat_Tps_RT_longest_iso_SB_asex, "dat_Tps_RT_longest_iso_SB_asex", "Tps", "RT")$out_df,
N_SB_Genes_kept_excluded(dat_Tbi_LG_longest_iso_SB_asex, "dat_Tbi_LG_longest_iso_SB_asex", "Tbi", "LG")$out_df,
N_SB_Genes_kept_excluded(dat_Tce_LG_longest_iso_SB_asex, "dat_Tce_LG_longest_iso_SB_asex", "Tce", "LG")$out_df,
N_SB_Genes_kept_excluded(dat_Tcm_LG_longest_iso_SB_asex, "dat_Tcm_LG_longest_iso_SB_asex", "Tcm", "LG")$out_df,
N_SB_Genes_kept_excluded(dat_Tpa_LG_longest_iso_SB_asex, "dat_Tpa_LG_longest_iso_SB_asex", "Tpa", "LG")$out_df,
N_SB_Genes_kept_excluded(dat_Tps_LG_longest_iso_SB_asex, "dat_Tps_LG_longest_iso_SB_asex", "Tps", "LG")$out_df)


N_kept_excl_SB_df$Ngenes = as.numeric(as.character(N_kept_excl_SB_df$Ngenes)) 
N_kept_excl_SB_df$total_genes = as.numeric(as.character(N_kept_excl_SB_df$total_genes)) 
N_kept_excl_SB_df$percent_genes = N_kept_excl_SB_df$Ngenes / N_kept_excl_SB_df$total_genes * 100

str(N_kept_excl_SB_df)

### fishers test

N_kept_excl_SB_FT_df <- rbind(
N_SB_Genes_kept_excluded(dat_Tbi_WB_longest_iso_SB_asex, "dat_Tbi_WB_longest_iso_SB_asex", "Tbi", "WB")$out_df3,
N_SB_Genes_kept_excluded(dat_Tce_WB_longest_iso_SB_asex, "dat_Tce_WB_longest_iso_SB_asex", "Tce", "WB")$out_df3,
N_SB_Genes_kept_excluded(dat_Tcm_WB_longest_iso_SB_asex, "dat_Tcm_WB_longest_iso_SB_asex", "Tcm", "WB")$out_df3,
N_SB_Genes_kept_excluded(dat_Tpa_WB_longest_iso_SB_asex, "dat_Tpa_WB_longest_iso_SB_asex", "Tpa", "WB")$out_df3,
N_SB_Genes_kept_excluded(dat_Tps_WB_longest_iso_SB_asex, "dat_Tps_WB_longest_iso_SB_asex", "Tps", "WB")$out_df3,
N_SB_Genes_kept_excluded(dat_Tbi_RT_longest_iso_SB_asex, "dat_Tbi_RT_longest_iso_SB_asex", "Tbi", "RT")$out_df3,
N_SB_Genes_kept_excluded(dat_Tce_RT_longest_iso_SB_asex, "dat_Tce_RT_longest_iso_SB_asex", "Tce", "RT")$out_df3,
N_SB_Genes_kept_excluded(dat_Tcm_RT_longest_iso_SB_asex, "dat_Tcm_RT_longest_iso_SB_asex", "Tcm", "RT")$out_df3,
N_SB_Genes_kept_excluded(dat_Tpa_RT_longest_iso_SB_asex, "dat_Tpa_RT_longest_iso_SB_asex", "Tpa", "RT")$out_df3,
N_SB_Genes_kept_excluded(dat_Tps_RT_longest_iso_SB_asex, "dat_Tps_RT_longest_iso_SB_asex", "Tps", "RT")$out_df3,
N_SB_Genes_kept_excluded(dat_Tbi_LG_longest_iso_SB_asex, "dat_Tbi_LG_longest_iso_SB_asex", "Tbi", "LG")$out_df3,
N_SB_Genes_kept_excluded(dat_Tce_LG_longest_iso_SB_asex, "dat_Tce_LG_longest_iso_SB_asex", "Tce", "LG")$out_df3,
N_SB_Genes_kept_excluded(dat_Tcm_LG_longest_iso_SB_asex, "dat_Tcm_LG_longest_iso_SB_asex", "Tcm", "LG")$out_df3,
N_SB_Genes_kept_excluded(dat_Tpa_LG_longest_iso_SB_asex, "dat_Tpa_LG_longest_iso_SB_asex", "Tpa", "LG")$out_df3,
N_SB_Genes_kept_excluded(dat_Tps_LG_longest_iso_SB_asex, "dat_Tps_LG_longest_iso_SB_asex", "Tps", "LG")$out_df3)

N_kept_excl_SB_FT_df$MB_excl      <-  as.numeric(as.character(N_kept_excl_SB_FT_df$MB_excl))
N_kept_excl_SB_FT_df$FB_excl      <-  as.numeric(as.character(N_kept_excl_SB_FT_df$FB_excl))
N_kept_excl_SB_FT_df$total_excl   <-  as.numeric(as.character(N_kept_excl_SB_FT_df$total_excl))
N_kept_excl_SB_FT_df$MB_kept      <-  as.numeric(as.character(N_kept_excl_SB_FT_df$MB_kept))
N_kept_excl_SB_FT_df$FB_kept      <- as.numeric(as.character(N_kept_excl_SB_FT_df$FB_kept))
N_kept_excl_SB_FT_df$total_kept   <- as.numeric(as.character(N_kept_excl_SB_FT_df$total_kept))
N_kept_excl_SB_FT_df$FB_fisher_OR <- as.numeric(as.character(N_kept_excl_SB_FT_df$FB_fisher_OR))
N_kept_excl_SB_FT_df$MB_fisher_OR <- as.numeric(as.character(N_kept_excl_SB_FT_df$MB_fisher_OR))
N_kept_excl_SB_FT_df$FB_fisher_p  <- as.numeric(as.character(N_kept_excl_SB_FT_df$FB_fisher_p))
N_kept_excl_SB_FT_df$MB_fisher_p  <- as.numeric(as.character(N_kept_excl_SB_FT_df$MB_fisher_p))
N_kept_excl_SB_FT_df$sp_ord <- ordered(N_kept_excl_SB_FT_df$sp, levels = c("Tbi", "Tce", "Tps", "Tcm", "Tpa" ))

####### FDR correct

N_kept_excl_SB_FT_df$FB_fisher_FDR <- p.adjust(N_kept_excl_SB_FT_df$FB_fisher_p, method = "BH")
N_kept_excl_SB_FT_df$MB_fisher_FDR <- p.adjust(N_kept_excl_SB_FT_df$MB_fisher_p, method = "BH")

write.csv(N_kept_excl_SB_FT_df, "N_kept_excl_SB_FT_df.csv", row.names = F)

subset(N_kept_excl_SB_FT_df, N_kept_excl_SB_FT_df$tiss == "WB")
subset(N_kept_excl_SB_FT_df, N_kept_excl_SB_FT_df$tiss == "RT")
subset(N_kept_excl_SB_FT_df, N_kept_excl_SB_FT_df$tiss == "LG")

head(N_kept_excl_SB_FT_df)


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

dat_ALL_WB_longest_iso_SB_asex_long <- make_long_table("dat_Tbi_WB_longest_iso_SB_asex", "dat_Tce_WB_longest_iso_SB_asex", "dat_Tcm_WB_longest_iso_SB_asex", "dat_Tpa_WB_longest_iso_SB_asex", "dat_Tps_WB_longest_iso_SB_asex",  "WB")
dat_ALL_RT_longest_iso_SB_asex_long <- make_long_table("dat_Tbi_RT_longest_iso_SB_asex", "dat_Tce_RT_longest_iso_SB_asex", "dat_Tcm_RT_longest_iso_SB_asex", "dat_Tpa_RT_longest_iso_SB_asex", "dat_Tps_RT_longest_iso_SB_asex",  "RT")
dat_ALL_LG_longest_iso_SB_asex_long <- make_long_table("dat_Tbi_LG_longest_iso_SB_asex", "dat_Tce_LG_longest_iso_SB_asex", "dat_Tcm_LG_longest_iso_SB_asex", "dat_Tpa_LG_longest_iso_SB_asex", "dat_Tps_LG_longest_iso_SB_asex",  "LG")


#### join tissues together

dat_ALL_WB_longest_iso_SB_asex_long$tiss <- rep("WB", length(dat_ALL_WB_longest_iso_SB_asex_long[,1]))
dat_ALL_WB_longest_iso_SB_asex_long$sp_group_tiss <- as.factor(paste(dat_ALL_WB_longest_iso_SB_asex_long$SBgroup, dat_ALL_WB_longest_iso_SB_asex_long$tiss, sep = "_"))
dat_ALL_WB_longest_iso_SB_asex_long$sp_tiss <- as.factor(paste(dat_ALL_WB_longest_iso_SB_asex_long$sp, dat_ALL_WB_longest_iso_SB_asex_long$tiss, sep = "_"))

dat_ALL_RT_longest_iso_SB_asex_long$tiss <- rep("RT", length(dat_ALL_RT_longest_iso_SB_asex_long[,1]))
dat_ALL_RT_longest_iso_SB_asex_long$sp_group_tiss <- as.factor(paste(dat_ALL_RT_longest_iso_SB_asex_long$SBgroup, dat_ALL_RT_longest_iso_SB_asex_long$tiss, sep = "_"))
dat_ALL_RT_longest_iso_SB_asex_long$sp_tiss <- as.factor(paste(dat_ALL_RT_longest_iso_SB_asex_long$sp, dat_ALL_RT_longest_iso_SB_asex_long$tiss, sep = "_"))

dat_ALL_LG_longest_iso_SB_asex_long$tiss <- rep("LG", length(dat_ALL_LG_longest_iso_SB_asex_long[,1]))
dat_ALL_LG_longest_iso_SB_asex_long$sp_group_tiss <- as.factor(paste(dat_ALL_LG_longest_iso_SB_asex_long$SBgroup, dat_ALL_LG_longest_iso_SB_asex_long$tiss, sep = "_"))
dat_ALL_LG_longest_iso_SB_asex_long$sp_tiss <- as.factor(paste(dat_ALL_LG_longest_iso_SB_asex_long$sp, dat_ALL_LG_longest_iso_SB_asex_long$tiss, sep = "_"))


dat_ALL_WBRTLG_longest_iso_SB_asex_long <- rbind(dat_ALL_WB_longest_iso_SB_asex_long, dat_ALL_RT_longest_iso_SB_asex_long, dat_ALL_LG_longest_iso_SB_asex_long)
dat_ALL_WBRTLG_longest_iso_SB_asex_long$sp_tiss_ord <- ordered(dat_ALL_WBRTLG_longest_iso_SB_asex_long$sp_tiss, levels = c("Tbi_WB", "Tce_WB", "Tps_WB", "Tcm_WB", "Tpa_WB", "Tbi_RT", "Tce_RT", "Tps_RT", "Tcm_RT", "Tpa_RT", "Tbi_LG", "Tce_LG", "Tps_LG", "Tcm_LG", "Tpa_LG"))	
dat_ALL_WBRTLG_longest_iso_SB_asex_long$sexbias_ord <- ordered(dat_ALL_WBRTLG_longest_iso_SB_asex_long$sexbias, levels = c("female_biased",  "Unbiased", "male_biased"))	
dat_ALL_WBRTLG_longest_iso_SB_asex_long$sexbias2_ord <- ordered(dat_ALL_WBRTLG_longest_iso_SB_asex_long$sexbias2, levels = c("female_biased",  "Unbiased", "male_biased"))	

str(dat_ALL_WBRTLG_longest_iso_SB_asex_long)

######### remove NAs

dat_ALL_WB_longest_iso_SB_asex_long_c <- na.omit(dat_ALL_WB_longest_iso_SB_asex_long)
dat_ALL_RT_longest_iso_SB_asex_long_c <- na.omit(dat_ALL_RT_longest_iso_SB_asex_long)
dat_ALL_LG_longest_iso_SB_asex_long_c <- na.omit(dat_ALL_LG_longest_iso_SB_asex_long)

dat_ALL_WBRTLG_longest_iso_SB_asex_long_c <- na.omit(dat_ALL_WBRTLG_longest_iso_SB_asex_long)

###### add some empty vals to space the boxplot


dat_ALL_WBRTLG_longest_iso_SB_asex_long_c$sp_tiss2 <- as.character(dat_ALL_WBRTLG_longest_iso_SB_asex_long_c$sp_tiss)
str(dat_ALL_WBRTLG_longest_iso_SB_asex_long_c)

dat_ALL_WBRTLG_longest_iso_SB_asex_long_c_2 <- rbind(dat_ALL_WBRTLG_longest_iso_SB_asex_long_c, c(NA,"Unbiased","Unbiased",NA,NA,NA,NA,NA,NA,"WB", NA,NA,NA,NA,NA, "skip_WB"))
dat_ALL_WBRTLG_longest_iso_SB_asex_long_c_2 <- rbind(dat_ALL_WBRTLG_longest_iso_SB_asex_long_c_2, c(NA,"Unbiased","Unbiased",NA,NA,NA,NA,NA,NA,"RT", NA,NA,NA,NA,NA, "skip_RT"))
dat_ALL_WBRTLG_longest_iso_SB_asex_long_c_2$log2FC_SA <- as.numeric(as.character(dat_ALL_WBRTLG_longest_iso_SB_asex_long_c_2$log2FC_SA ))

dat_ALL_WBRTLG_longest_iso_SB_asex_long_c_2$sp_tiss2 <- as.factor(dat_ALL_WBRTLG_longest_iso_SB_asex_long_c_2$sp_tiss2)
dat_ALL_WBRTLG_longest_iso_SB_asex_long_c_2$sp_tiss2_ord <- ordered(dat_ALL_WBRTLG_longest_iso_SB_asex_long_c_2$sp_tiss2, levels = c("Tbi_WB", "Tce_WB", "Tps_WB", "Tcm_WB", "Tpa_WB", "skip_WB", "Tbi_RT", "Tce_RT", "Tps_RT", "Tcm_RT", "Tpa_RT", "skip_RT", "Tbi_LG", "Tce_LG", "Tps_LG", "Tcm_LG", "Tpa_LG"))	

str(dat_ALL_WBRTLG_longest_iso_SB_asex_long_c_2)

##### plot boxplots

WBRTLG_SA_box_s2 <- ggplot(dat_ALL_WBRTLG_longest_iso_SB_asex_long_c_2, aes(sp_tiss2_ord, log2FC_SA)) + 
	theme_classic() +
	geom_boxplot(aes(fill = factor(sexbias2_ord)),position=position_dodge(0.65), width = 0.45, outlier.size = 0, lwd = 0.5, fatten = 1, outlier.shape = NA) +
	coord_cartesian(ylim=c(-3,3)) +
	ylab ("log2 fold change in asexual females") +
	xlab ("Gene-class by species") + 
	scale_fill_manual(values=c("firebrick2", "grey", "royalblue2")) + geom_hline(yintercept = 0) + labs(fill='') +
	ggtitle("SB class = FDR <0.05 FC > 2")

WBRTLG_SA_box_s2_2 <- WBRTLG_SA_box_s2 + theme(axis.text.x = element_text(angle = 90, hjust = 1))



pdf(paste("Sex_asex_longest_iso_LG_FC_boxplot_s2_FDR005FC2.pdf",sep = ""), width = 8.2, height = 5.5)
WBRTLG_SA_box_s2_2
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
wilk_test("WB", dat_ALL_WB_longest_iso_SB_asex_long_c)$s1_out_table,	
wilk_test("RT", dat_ALL_RT_longest_iso_SB_asex_long_c)$s1_out_table,	
wilk_test("LG", dat_ALL_LG_longest_iso_SB_asex_long_c)$s1_out_table))		

wilcox_S1_FDR005$FB_wilp_s1 <- as.numeric(as.character(wilcox_S1_FDR005$FB_wilp_s1))
wilcox_S1_FDR005$MB_wilp_s1 <- as.numeric(as.character(wilcox_S1_FDR005$MB_wilp_s1))

wilcox_S1_FDR005$FB_wilFDR_s1 <- p.adjust(wilcox_S1_FDR005$FB_wilp_s1, method = "BH")
wilcox_S1_FDR005$MB_wilFDR_s1 <- p.adjust(wilcox_S1_FDR005$MB_wilp_s1, method = "BH")



wilcox_S2_FDR005_FC2 <- as.data.frame(rbind(
wilk_test("WB", dat_ALL_WB_longest_iso_SB_asex_long_c)$s2_out_table,	
wilk_test("RT", dat_ALL_RT_longest_iso_SB_asex_long_c)$s2_out_table,	
wilk_test("LG", dat_ALL_LG_longest_iso_SB_asex_long_c)$s2_out_table))		

wilcox_S2_FDR005_FC2$FB_wilp_s2 <- as.numeric(as.character(wilcox_S2_FDR005_FC2$FB_wilp_s2))
wilcox_S2_FDR005_FC2$MB_wilp_s2 <- as.numeric(as.character(wilcox_S2_FDR005_FC2$MB_wilp_s2))

wilcox_S2_FDR005_FC2$FB_wilFDR_s1 <- p.adjust(wilcox_S2_FDR005_FC2$FB_wilp_s2, method = "BH")
wilcox_S2_FDR005_FC2$MB_wilFDR_s1 <- p.adjust(wilcox_S2_FDR005_FC2$MB_wilp_s2, method = "BH")


write.csv(wilcox_S1_FDR005, "wilcox_S1_FDR005.csv")
write.csv(wilcox_S2_FDR005_FC2, "wilcox_S2_FDR005_FC2.csv")


########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), "sex_bias_plotsetc_sexual_ref_sessionInfo.txt")




