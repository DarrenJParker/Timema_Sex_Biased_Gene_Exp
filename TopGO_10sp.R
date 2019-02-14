### topGO
# install
# source("http://bioconductor.org/biocLite.R") 
# biocLite() 
# source("http://bioconductor.org/biocLite.R")   
# biocLite("topGO")
# biocLite("ALL")
# biocLite("affyLib")

library(topGO)
library(ALL)
library("VennDiagram")
library(gridExtra)
library(grid)
library(ggplot2)
library("SuperExactTest")

print (sessionInfo())

# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS High Sierra 10.13.4

# Matrix products: default
# BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
 # [1] grid      stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
 # [1] SuperExactTest_0.99.4 ggplot2_2.2.1         gridExtra_2.3         VennDiagram_1.6.17    futile.logger_1.4.3   ALL_1.18.0            topGO_2.28.0          SparseM_1.77         
 # [9] GO.db_3.4.1           AnnotationDbi_1.38.2  IRanges_2.10.5        S4Vectors_0.14.7      Biobase_2.36.2        graph_1.54.0          BiocGenerics_0.22.1  

# loaded via a namespace (and not attached):
 # [1] Rcpp_0.12.13         munsell_0.4.3        bit_1.1-12           colorspace_1.3-2     lattice_0.20-35      rlang_0.1.6          plyr_1.8.4           blob_1.1.0          
 # [9] tools_3.4.1          gtable_0.2.0         DBI_0.7              lambda.r_1.2         matrixStats_0.52.2   lazyeval_0.2.1       bit64_0.9-7          digest_0.6.12       
# [17] tibble_1.3.4         futile.options_1.0.0 memoise_1.1.0        RSQLite_2.0          compiler_3.4.1       scales_0.5.0         pkgconfig_2.0.1     


Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function (x) {  
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's. 
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}



#### load annotation 

setwd("Data/for_GOterm_analyses")

## nr_annotated
geneID2GO_Tbi_nr <- readMappings(file = "tbi_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_nr.annot_fortopgo.txt")
geneID2GO_Tce_nr <- readMappings(file = "tce_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_nr.annot_fortopgo.txt")
geneID2GO_Tcm_nr <- readMappings(file = "tcm_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_nr.annot_fortopgo.txt")
geneID2GO_Tpa_nr <- readMappings(file = "tpa_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_nr.annot_fortopgo.txt")
geneID2GO_Tps_nr <- readMappings(file = "tps_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_nr.annot_fortopgo.txt")

## Droso annotated
geneID2GO_Tbi_DROSO <- readMappings(file = "tbi_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_DROSO.annot_fortopgo.txt")
geneID2GO_Tce_DROSO <- readMappings(file = "tce_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_DROSO.annot_fortopgo.txt")
geneID2GO_Tcm_DROSO <- readMappings(file = "tcm_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_DROSO.annot_fortopgo.txt")
geneID2GO_Tpa_DROSO <- readMappings(file = "tpa_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_DROSO.annot_fortopgo.txt")
geneID2GO_Tps_DROSO <- readMappings(file = "tps_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_DROSO.annot_fortopgo.txt")


str(head(geneID2GO_Tcm_DROSO))

###############################################################################################################################################
#### read in tables with genename and rank

setwd("10_sp")

### Ranked gene lists of FB genes [female-biased FC > 2 + FDR < 0.05, ranked by FDR] + [female-biased rest, ranked by FDR] +  [male-biased FC < 2 + FDR < 0.05, rev-ranked by FDR] + [male-biased FC rest, rev-ranked by FDR]

RGL_Tbi_WB_FB_rankedlist_wFC <- as.list(read.table("lrt_Tbi_sex_bias_WB_FB_rankedlist_asranks_wFC.txt"))
RGL_Tce_WB_FB_rankedlist_wFC <- as.list(read.table("lrt_Tce_sex_bias_WB_FB_rankedlist_asranks_wFC.txt"))
RGL_Tcm_WB_FB_rankedlist_wFC <- as.list(read.table("lrt_Tcm_sex_bias_WB_FB_rankedlist_asranks_wFC.txt"))
RGL_Tpa_WB_FB_rankedlist_wFC <- as.list(read.table("lrt_Tpa_sex_bias_WB_FB_rankedlist_asranks_wFC.txt"))
RGL_Tps_WB_FB_rankedlist_wFC <- as.list(read.table("lrt_Tps_sex_bias_WB_FB_rankedlist_asranks_wFC.txt"))
RGL_Tbi_RT_FB_rankedlist_wFC <- as.list(read.table("lrt_Tbi_sex_bias_RT_FB_rankedlist_asranks_wFC.txt"))
RGL_Tce_RT_FB_rankedlist_wFC <- as.list(read.table("lrt_Tce_sex_bias_RT_FB_rankedlist_asranks_wFC.txt"))
RGL_Tcm_RT_FB_rankedlist_wFC <- as.list(read.table("lrt_Tcm_sex_bias_RT_FB_rankedlist_asranks_wFC.txt"))
RGL_Tpa_RT_FB_rankedlist_wFC <- as.list(read.table("lrt_Tpa_sex_bias_RT_FB_rankedlist_asranks_wFC.txt"))
RGL_Tps_RT_FB_rankedlist_wFC <- as.list(read.table("lrt_Tps_sex_bias_RT_FB_rankedlist_asranks_wFC.txt"))
RGL_Tbi_LG_FB_rankedlist_wFC <- as.list(read.table("lrt_Tbi_sex_bias_LG_FB_rankedlist_asranks_wFC.txt"))
RGL_Tce_LG_FB_rankedlist_wFC <- as.list(read.table("lrt_Tce_sex_bias_LG_FB_rankedlist_asranks_wFC.txt"))
RGL_Tcm_LG_FB_rankedlist_wFC <- as.list(read.table("lrt_Tcm_sex_bias_LG_FB_rankedlist_asranks_wFC.txt"))
RGL_Tpa_LG_FB_rankedlist_wFC <- as.list(read.table("lrt_Tpa_sex_bias_LG_FB_rankedlist_asranks_wFC.txt"))
RGL_Tps_LG_FB_rankedlist_wFC <- as.list(read.table("lrt_Tps_sex_bias_LG_FB_rankedlist_asranks_wFC.txt"))


### Ranked gene lists of MB genes [male-biased FC > 2 + FDR < 0.05, ranked by FDR] + [male-biased rest, ranked by FDR] +  [female-biased FC < 2 + FDR < 0.05, rev-ranked by FDR] + [female-biased FC rest, rev-ranked by FDR]

RGL_Tbi_WB_MB_rankedlist_wFC <- as.list(read.table("lrt_Tbi_sex_bias_WB_MB_rankedlist_asranks_wFC.txt"))
RGL_Tce_WB_MB_rankedlist_wFC <- as.list(read.table("lrt_Tce_sex_bias_WB_MB_rankedlist_asranks_wFC.txt"))
RGL_Tcm_WB_MB_rankedlist_wFC <- as.list(read.table("lrt_Tcm_sex_bias_WB_MB_rankedlist_asranks_wFC.txt"))
RGL_Tpa_WB_MB_rankedlist_wFC <- as.list(read.table("lrt_Tpa_sex_bias_WB_MB_rankedlist_asranks_wFC.txt"))
RGL_Tps_WB_MB_rankedlist_wFC <- as.list(read.table("lrt_Tps_sex_bias_WB_MB_rankedlist_asranks_wFC.txt"))
RGL_Tbi_RT_MB_rankedlist_wFC <- as.list(read.table("lrt_Tbi_sex_bias_RT_MB_rankedlist_asranks_wFC.txt"))
RGL_Tce_RT_MB_rankedlist_wFC <- as.list(read.table("lrt_Tce_sex_bias_RT_MB_rankedlist_asranks_wFC.txt"))
RGL_Tcm_RT_MB_rankedlist_wFC <- as.list(read.table("lrt_Tcm_sex_bias_RT_MB_rankedlist_asranks_wFC.txt"))
RGL_Tpa_RT_MB_rankedlist_wFC <- as.list(read.table("lrt_Tpa_sex_bias_RT_MB_rankedlist_asranks_wFC.txt"))
RGL_Tps_RT_MB_rankedlist_wFC <- as.list(read.table("lrt_Tps_sex_bias_RT_MB_rankedlist_asranks_wFC.txt"))
RGL_Tbi_LG_MB_rankedlist_wFC <- as.list(read.table("lrt_Tbi_sex_bias_LG_MB_rankedlist_asranks_wFC.txt"))
RGL_Tce_LG_MB_rankedlist_wFC <- as.list(read.table("lrt_Tce_sex_bias_LG_MB_rankedlist_asranks_wFC.txt"))
RGL_Tcm_LG_MB_rankedlist_wFC <- as.list(read.table("lrt_Tcm_sex_bias_LG_MB_rankedlist_asranks_wFC.txt"))
RGL_Tpa_LG_MB_rankedlist_wFC <- as.list(read.table("lrt_Tpa_sex_bias_LG_MB_rankedlist_asranks_wFC.txt"))
RGL_Tps_LG_MB_rankedlist_wFC <- as.list(read.table("lrt_Tps_sex_bias_LG_MB_rankedlist_asranks_wFC.txt"))

### Sig_female biased genes (FDR + FC) ranked by FC change in asexuals (+ve fold changes at the top, -ve at the bottom)

RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tbi_sex_asex_WB_SA_sig_wFC_FB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tce_sex_asex_WB_SA_sig_wFC_FB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tcm_sex_asex_WB_SA_sig_wFC_FB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tpa_sex_asex_WB_SA_sig_wFC_FB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tps_sex_asex_WB_SA_sig_wFC_FB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tbi_sex_asex_RT_SA_sig_wFC_FB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tce_sex_asex_RT_SA_sig_wFC_FB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tcm_sex_asex_RT_SA_sig_wFC_FB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tpa_sex_asex_RT_SA_sig_wFC_FB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tps_sex_asex_RT_SA_sig_wFC_FB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tbi_sex_asex_LG_SA_sig_wFC_FB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tce_sex_asex_LG_SA_sig_wFC_FB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tcm_sex_asex_LG_SA_sig_wFC_FB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tpa_sex_asex_LG_SA_sig_wFC_FB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tps_sex_asex_LG_SA_sig_wFC_FB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))

### Sig_male biased genes (FDR + FC) ranked by FC change in asexuals (+ve fold changes at the top, -ve at the bottom)

RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tbi_sex_asex_WB_SA_sig_wFC_MB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tce_sex_asex_WB_SA_sig_wFC_MB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tcm_sex_asex_WB_SA_sig_wFC_MB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tpa_sex_asex_WB_SA_sig_wFC_MB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tps_sex_asex_WB_SA_sig_wFC_MB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tbi_sex_asex_RT_SA_sig_wFC_MB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tce_sex_asex_RT_SA_sig_wFC_MB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tcm_sex_asex_RT_SA_sig_wFC_MB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tpa_sex_asex_RT_SA_sig_wFC_MB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tps_sex_asex_RT_SA_sig_wFC_MB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tbi_sex_asex_LG_SA_sig_wFC_MB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tce_sex_asex_LG_SA_sig_wFC_MB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tcm_sex_asex_LG_SA_sig_wFC_MB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tpa_sex_asex_LG_SA_sig_wFC_MB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))
RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_first <- as.list(read.table("lrt_Tps_sex_asex_LG_SA_sig_wFC_MB_list_sorted_positive_FC_first_rankedlist_asranks.txt"))

### Sig_female biased genes (FDR + FC) ranked by FC change in asexuals (-ve fold changes at the top, +ve at the bottom)

RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tbi_sex_asex_WB_SA_sig_wFC_FB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tce_sex_asex_WB_SA_sig_wFC_FB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tcm_sex_asex_WB_SA_sig_wFC_FB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tpa_sex_asex_WB_SA_sig_wFC_FB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tps_sex_asex_WB_SA_sig_wFC_FB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tbi_sex_asex_RT_SA_sig_wFC_FB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tce_sex_asex_RT_SA_sig_wFC_FB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tcm_sex_asex_RT_SA_sig_wFC_FB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tpa_sex_asex_RT_SA_sig_wFC_FB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tps_sex_asex_RT_SA_sig_wFC_FB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tbi_sex_asex_LG_SA_sig_wFC_FB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tce_sex_asex_LG_SA_sig_wFC_FB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tcm_sex_asex_LG_SA_sig_wFC_FB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tpa_sex_asex_LG_SA_sig_wFC_FB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tps_sex_asex_LG_SA_sig_wFC_FB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))

### Sig_male biased genes (FDR + FC) ranked by FC change in asexuals (-ve fold changes at the top, +ve at the bottom)

RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tbi_sex_asex_WB_SA_sig_wFC_MB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tce_sex_asex_WB_SA_sig_wFC_MB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tcm_sex_asex_WB_SA_sig_wFC_MB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tpa_sex_asex_WB_SA_sig_wFC_MB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tps_sex_asex_WB_SA_sig_wFC_MB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tbi_sex_asex_RT_SA_sig_wFC_MB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tce_sex_asex_RT_SA_sig_wFC_MB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tcm_sex_asex_RT_SA_sig_wFC_MB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tpa_sex_asex_RT_SA_sig_wFC_MB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tps_sex_asex_RT_SA_sig_wFC_MB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tbi_sex_asex_LG_SA_sig_wFC_MB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tce_sex_asex_LG_SA_sig_wFC_MB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tcm_sex_asex_LG_SA_sig_wFC_MB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tpa_sex_asex_LG_SA_sig_wFC_MB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))
RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_last <- as.list(read.table("lrt_Tps_sex_asex_LG_SA_sig_wFC_MB_list_sorted_positive_FC_last_rankedlist_asranks.txt"))


#### FOR OUTPUT

setwd("../../..")
dir.create("Output")
dir.create("Output/TopGO_out")
dir.create("Output/TopGO_out/10_sp")
setwd("Output/TopGO_out/10_sp")


#### need to get the genes as a named numeric vector

RGL_Tbi_WB_FB_rankedlist_wFC_GL <- RGL_Tbi_WB_FB_rankedlist_wFC$V2
RGL_Tce_WB_FB_rankedlist_wFC_GL <- RGL_Tce_WB_FB_rankedlist_wFC$V2
RGL_Tcm_WB_FB_rankedlist_wFC_GL <- RGL_Tcm_WB_FB_rankedlist_wFC$V2
RGL_Tpa_WB_FB_rankedlist_wFC_GL <- RGL_Tpa_WB_FB_rankedlist_wFC$V2
RGL_Tps_WB_FB_rankedlist_wFC_GL <- RGL_Tps_WB_FB_rankedlist_wFC$V2
RGL_Tbi_RT_FB_rankedlist_wFC_GL <- RGL_Tbi_RT_FB_rankedlist_wFC$V2
RGL_Tce_RT_FB_rankedlist_wFC_GL <- RGL_Tce_RT_FB_rankedlist_wFC$V2
RGL_Tcm_RT_FB_rankedlist_wFC_GL <- RGL_Tcm_RT_FB_rankedlist_wFC$V2
RGL_Tpa_RT_FB_rankedlist_wFC_GL <- RGL_Tpa_RT_FB_rankedlist_wFC$V2
RGL_Tps_RT_FB_rankedlist_wFC_GL <- RGL_Tps_RT_FB_rankedlist_wFC$V2
RGL_Tbi_LG_FB_rankedlist_wFC_GL <- RGL_Tbi_LG_FB_rankedlist_wFC$V2
RGL_Tce_LG_FB_rankedlist_wFC_GL <- RGL_Tce_LG_FB_rankedlist_wFC$V2
RGL_Tcm_LG_FB_rankedlist_wFC_GL <- RGL_Tcm_LG_FB_rankedlist_wFC$V2
RGL_Tpa_LG_FB_rankedlist_wFC_GL <- RGL_Tpa_LG_FB_rankedlist_wFC$V2
RGL_Tps_LG_FB_rankedlist_wFC_GL <- RGL_Tps_LG_FB_rankedlist_wFC$V2

RGL_Tbi_WB_MB_rankedlist_wFC_GL <- RGL_Tbi_WB_MB_rankedlist_wFC$V2
RGL_Tce_WB_MB_rankedlist_wFC_GL <- RGL_Tce_WB_MB_rankedlist_wFC$V2
RGL_Tcm_WB_MB_rankedlist_wFC_GL <- RGL_Tcm_WB_MB_rankedlist_wFC$V2
RGL_Tpa_WB_MB_rankedlist_wFC_GL <- RGL_Tpa_WB_MB_rankedlist_wFC$V2
RGL_Tps_WB_MB_rankedlist_wFC_GL <- RGL_Tps_WB_MB_rankedlist_wFC$V2
RGL_Tbi_RT_MB_rankedlist_wFC_GL <- RGL_Tbi_RT_MB_rankedlist_wFC$V2
RGL_Tce_RT_MB_rankedlist_wFC_GL <- RGL_Tce_RT_MB_rankedlist_wFC$V2
RGL_Tcm_RT_MB_rankedlist_wFC_GL <- RGL_Tcm_RT_MB_rankedlist_wFC$V2
RGL_Tpa_RT_MB_rankedlist_wFC_GL <- RGL_Tpa_RT_MB_rankedlist_wFC$V2
RGL_Tps_RT_MB_rankedlist_wFC_GL <- RGL_Tps_RT_MB_rankedlist_wFC$V2
RGL_Tbi_LG_MB_rankedlist_wFC_GL <- RGL_Tbi_LG_MB_rankedlist_wFC$V2
RGL_Tce_LG_MB_rankedlist_wFC_GL <- RGL_Tce_LG_MB_rankedlist_wFC$V2
RGL_Tcm_LG_MB_rankedlist_wFC_GL <- RGL_Tcm_LG_MB_rankedlist_wFC$V2
RGL_Tpa_LG_MB_rankedlist_wFC_GL <- RGL_Tpa_LG_MB_rankedlist_wFC$V2
RGL_Tps_LG_MB_rankedlist_wFC_GL <- RGL_Tps_LG_MB_rankedlist_wFC$V2


RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_first_GL <- RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_first$V2
RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_first_GL <- RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_first$V2
RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_first_GL <- RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_first$V2
RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_first_GL <- RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_first$V2
RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_first_GL <- RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_first$V2
RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_first_GL <- RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_first$V2
RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_first_GL <- RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_first$V2
RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_first_GL <- RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_first$V2
RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_first_GL <- RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_first$V2
RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_first_GL <- RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_first$V2
RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_first_GL <- RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_first$V2
RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_first_GL <- RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_first$V2
RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_first_GL <- RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_first$V2
RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_first_GL <- RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_first$V2
RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_first_GL <- RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_first$V2

RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_first_GL <- RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_first$V2
RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_first_GL <- RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_first$V2
RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_first_GL <- RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_first$V2
RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_first_GL <- RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_first$V2
RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_first_GL <- RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_first$V2
RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_first_GL <- RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_first$V2
RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_first_GL <- RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_first$V2
RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_first_GL <- RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_first$V2
RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_first_GL <- RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_first$V2
RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_first_GL <- RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_first$V2
RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_first_GL <- RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_first$V2
RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_first_GL <- RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_first$V2
RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_first_GL <- RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_first$V2
RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_first_GL <- RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_first$V2
RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_first_GL <- RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_first$V2

RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_last_GL <- RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_last$V2
RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_last_GL <- RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_last$V2
RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_last_GL <- RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_last$V2
RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_last_GL <- RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_last$V2
RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_last_GL <- RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_last$V2
RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_last_GL <- RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_last$V2
RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_last_GL <- RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_last$V2
RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_last_GL <- RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_last$V2
RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_last_GL <- RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_last$V2
RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_last_GL <- RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_last$V2
RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_last_GL <- RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_last$V2
RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_last_GL <- RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_last$V2
RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_last_GL <- RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_last$V2
RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_last_GL <- RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_last$V2
RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_last_GL <- RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_last$V2

RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_last_GL <- RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_last$V2
RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_last_GL <- RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_last$V2
RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_last_GL <- RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_last$V2
RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_last_GL <- RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_last$V2
RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_last_GL <- RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_last$V2
RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_last_GL <- RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_last$V2
RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_last_GL <- RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_last$V2
RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_last_GL <- RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_last$V2
RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_last_GL <- RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_last$V2
RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_last_GL <- RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_last$V2
RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_last_GL <- RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_last$V2
RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_last_GL <- RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_last$V2
RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_last_GL <- RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_last$V2
RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_last_GL <- RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_last$V2
RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_last_GL <- RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_last$V2


names(RGL_Tbi_WB_FB_rankedlist_wFC_GL) <- RGL_Tbi_WB_FB_rankedlist_wFC$V1
names(RGL_Tce_WB_FB_rankedlist_wFC_GL) <- RGL_Tce_WB_FB_rankedlist_wFC$V1
names(RGL_Tcm_WB_FB_rankedlist_wFC_GL) <- RGL_Tcm_WB_FB_rankedlist_wFC$V1
names(RGL_Tpa_WB_FB_rankedlist_wFC_GL) <- RGL_Tpa_WB_FB_rankedlist_wFC$V1
names(RGL_Tps_WB_FB_rankedlist_wFC_GL) <- RGL_Tps_WB_FB_rankedlist_wFC$V1
names(RGL_Tbi_RT_FB_rankedlist_wFC_GL) <- RGL_Tbi_RT_FB_rankedlist_wFC$V1
names(RGL_Tce_RT_FB_rankedlist_wFC_GL) <- RGL_Tce_RT_FB_rankedlist_wFC$V1
names(RGL_Tcm_RT_FB_rankedlist_wFC_GL) <- RGL_Tcm_RT_FB_rankedlist_wFC$V1
names(RGL_Tpa_RT_FB_rankedlist_wFC_GL) <- RGL_Tpa_RT_FB_rankedlist_wFC$V1
names(RGL_Tps_RT_FB_rankedlist_wFC_GL) <- RGL_Tps_RT_FB_rankedlist_wFC$V1
names(RGL_Tbi_LG_FB_rankedlist_wFC_GL) <- RGL_Tbi_LG_FB_rankedlist_wFC$V1
names(RGL_Tce_LG_FB_rankedlist_wFC_GL) <- RGL_Tce_LG_FB_rankedlist_wFC$V1
names(RGL_Tcm_LG_FB_rankedlist_wFC_GL) <- RGL_Tcm_LG_FB_rankedlist_wFC$V1
names(RGL_Tpa_LG_FB_rankedlist_wFC_GL) <- RGL_Tpa_LG_FB_rankedlist_wFC$V1
names(RGL_Tps_LG_FB_rankedlist_wFC_GL) <- RGL_Tps_LG_FB_rankedlist_wFC$V1

names(RGL_Tbi_WB_MB_rankedlist_wFC_GL) <- RGL_Tbi_WB_MB_rankedlist_wFC$V1
names(RGL_Tce_WB_MB_rankedlist_wFC_GL) <- RGL_Tce_WB_MB_rankedlist_wFC$V1
names(RGL_Tcm_WB_MB_rankedlist_wFC_GL) <- RGL_Tcm_WB_MB_rankedlist_wFC$V1
names(RGL_Tpa_WB_MB_rankedlist_wFC_GL) <- RGL_Tpa_WB_MB_rankedlist_wFC$V1
names(RGL_Tps_WB_MB_rankedlist_wFC_GL) <- RGL_Tps_WB_MB_rankedlist_wFC$V1
names(RGL_Tbi_RT_MB_rankedlist_wFC_GL) <- RGL_Tbi_RT_MB_rankedlist_wFC$V1
names(RGL_Tce_RT_MB_rankedlist_wFC_GL) <- RGL_Tce_RT_MB_rankedlist_wFC$V1
names(RGL_Tcm_RT_MB_rankedlist_wFC_GL) <- RGL_Tcm_RT_MB_rankedlist_wFC$V1
names(RGL_Tpa_RT_MB_rankedlist_wFC_GL) <- RGL_Tpa_RT_MB_rankedlist_wFC$V1
names(RGL_Tps_RT_MB_rankedlist_wFC_GL) <- RGL_Tps_RT_MB_rankedlist_wFC$V1
names(RGL_Tbi_LG_MB_rankedlist_wFC_GL) <- RGL_Tbi_LG_MB_rankedlist_wFC$V1
names(RGL_Tce_LG_MB_rankedlist_wFC_GL) <- RGL_Tce_LG_MB_rankedlist_wFC$V1
names(RGL_Tcm_LG_MB_rankedlist_wFC_GL) <- RGL_Tcm_LG_MB_rankedlist_wFC$V1
names(RGL_Tpa_LG_MB_rankedlist_wFC_GL) <- RGL_Tpa_LG_MB_rankedlist_wFC$V1
names(RGL_Tps_LG_MB_rankedlist_wFC_GL) <- RGL_Tps_LG_MB_rankedlist_wFC$V1

names(RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_first_GL) <- RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_first$V1
names(RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_first_GL) <- RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_first$V1
names(RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_first_GL) <- RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_first$V1
names(RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_first_GL) <- RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_first$V1
names(RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_first_GL) <- RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_first$V1
names(RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_first_GL) <- RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_first$V1
names(RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_first_GL) <- RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_first$V1
names(RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_first_GL) <- RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_first$V1
names(RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_first_GL) <- RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_first$V1
names(RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_first_GL) <- RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_first$V1
names(RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_first_GL) <- RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_first$V1
names(RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_first_GL) <- RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_first$V1
names(RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_first_GL) <- RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_first$V1
names(RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_first_GL) <- RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_first$V1
names(RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_first_GL) <- RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_first$V1

names(RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_first_GL) <- RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_first$V1
names(RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_first_GL) <- RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_first$V1
names(RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_first_GL) <- RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_first$V1
names(RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_first_GL) <- RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_first$V1
names(RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_first_GL) <- RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_first$V1
names(RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_first_GL) <- RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_first$V1
names(RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_first_GL) <- RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_first$V1
names(RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_first_GL) <- RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_first$V1
names(RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_first_GL) <- RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_first$V1
names(RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_first_GL) <- RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_first$V1
names(RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_first_GL) <- RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_first$V1
names(RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_first_GL) <- RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_first$V1
names(RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_first_GL) <- RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_first$V1
names(RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_first_GL) <- RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_first$V1
names(RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_first_GL) <- RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_first$V1

names(RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_last_GL) <- RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_last$V1
names(RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_last_GL) <- RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_last$V1
names(RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_last_GL) <- RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_last$V1
names(RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_last_GL) <- RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_last$V1
names(RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_last_GL) <- RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_last$V1
names(RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_last_GL) <- RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_last$V1
names(RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_last_GL) <- RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_last$V1
names(RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_last_GL) <- RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_last$V1
names(RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_last_GL) <- RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_last$V1
names(RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_last_GL) <- RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_last$V1
names(RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_last_GL) <- RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_last$V1
names(RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_last_GL) <- RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_last$V1
names(RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_last_GL) <- RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_last$V1
names(RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_last_GL) <- RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_last$V1
names(RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_last_GL) <- RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_last$V1

names(RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_last_GL) <- RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_last$V1
names(RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_last_GL) <- RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_last$V1
names(RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_last_GL) <- RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_last$V1
names(RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_last_GL) <- RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_last$V1
names(RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_last_GL) <- RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_last$V1
names(RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_last_GL) <- RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_last$V1
names(RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_last_GL) <- RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_last$V1
names(RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_last_GL) <- RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_last$V1
names(RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_last_GL) <- RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_last$V1
names(RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_last_GL) <- RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_last$V1
names(RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_last_GL) <- RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_last$V1
names(RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_last_GL) <- RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_last$V1
names(RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_last_GL) <- RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_last$V1
names(RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_last_GL) <- RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_last$V1
names(RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_last_GL) <- RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_last$V1


## set the sig_for_go to set the threshold for significant GO terms from the output of the GSEA
run_enrichment <- function(genelist, ref, sig_for_GO){
	
	### make rule for classing sig / non-sig - note this rule is not used for the GSEA
	
	topDiffGenes <- function(allScore) {return(allScore < 0.05)}
	# topDiffGenes <- function(allScore) {return(allScore < 1)} ## as a check - setting to one gives the same pvalues for the GSEA
	
	#### make GOdata object
	#### setting node size as 5 so at least 5 genes must be annot per GO terms 
	#### do enrichment test
	
	GODATA_BP = new("topGOdata", ontology = "BP", allGenes = genelist, geneSel = topDiffGenes,  annot = annFUN.gene2GO, gene2GO = ref, nodeSize = 5)

	### get N GOs used

	GO_term_use_BP_list = GODATA_BP@graph@nodes
	N_GO_term_use_BP = length(GODATA_BP@graph@nodes)
	result_GSEA_BP     <- runTest(GODATA_BP, algorithm = "elim", statistic = "ks")

	### combined tables
	allRes1_BP <- GenTable(GODATA_BP, GSEA = result_GSEA_BP, ranksOf = "GSEA", topNodes = length(GODATA_BP@graph@nodes), numChar = 200)

	sig_GSEA_BP_GO     = subset(allRes1_BP, allRes1_BP$GSEA < sig_for_GO)$GO.ID
	
	## return everything!
	out_list = list("N_GO_term_use_BP" = N_GO_term_use_BP, 
	                "GO_term_use_BP_list" = GO_term_use_BP_list, 
	                "allRes1_BP" = allRes1_BP, 
	                "sig_GSEA_BP_GO" = sig_GSEA_BP_GO,
	                "GODATA_BP" = GODATA_BP) 
	return(out_list)

}

#### run the enrichment stuff (0.05)

### nr

RGL_Tbi_WB_FB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tbi_WB_FB_rankedlist_wFC_GL, geneID2GO_Tbi_nr, 0.05)
RGL_Tce_WB_FB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tce_WB_FB_rankedlist_wFC_GL, geneID2GO_Tce_nr, 0.05)
RGL_Tcm_WB_FB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tcm_WB_FB_rankedlist_wFC_GL, geneID2GO_Tcm_nr, 0.05)
RGL_Tpa_WB_FB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tpa_WB_FB_rankedlist_wFC_GL, geneID2GO_Tpa_nr, 0.05)
RGL_Tps_WB_FB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tps_WB_FB_rankedlist_wFC_GL, geneID2GO_Tps_nr, 0.05)
RGL_Tbi_RT_FB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tbi_RT_FB_rankedlist_wFC_GL, geneID2GO_Tbi_nr, 0.05)
RGL_Tce_RT_FB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tce_RT_FB_rankedlist_wFC_GL, geneID2GO_Tce_nr, 0.05)
RGL_Tcm_RT_FB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tcm_RT_FB_rankedlist_wFC_GL, geneID2GO_Tcm_nr, 0.05)
RGL_Tpa_RT_FB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tpa_RT_FB_rankedlist_wFC_GL, geneID2GO_Tpa_nr, 0.05)
RGL_Tps_RT_FB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tps_RT_FB_rankedlist_wFC_GL, geneID2GO_Tps_nr, 0.05)
RGL_Tbi_LG_FB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tbi_LG_FB_rankedlist_wFC_GL, geneID2GO_Tbi_nr, 0.05)
RGL_Tce_LG_FB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tce_LG_FB_rankedlist_wFC_GL, geneID2GO_Tce_nr, 0.05)
RGL_Tcm_LG_FB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tcm_LG_FB_rankedlist_wFC_GL, geneID2GO_Tcm_nr, 0.05)
RGL_Tpa_LG_FB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tpa_LG_FB_rankedlist_wFC_GL, geneID2GO_Tpa_nr, 0.05)
RGL_Tps_LG_FB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tps_LG_FB_rankedlist_wFC_GL, geneID2GO_Tps_nr, 0.05)

RGL_Tbi_WB_MB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tbi_WB_MB_rankedlist_wFC_GL, geneID2GO_Tbi_nr, 0.05)
RGL_Tce_WB_MB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tce_WB_MB_rankedlist_wFC_GL, geneID2GO_Tce_nr, 0.05)
RGL_Tcm_WB_MB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tcm_WB_MB_rankedlist_wFC_GL, geneID2GO_Tcm_nr, 0.05)
RGL_Tpa_WB_MB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tpa_WB_MB_rankedlist_wFC_GL, geneID2GO_Tpa_nr, 0.05)
RGL_Tps_WB_MB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tps_WB_MB_rankedlist_wFC_GL, geneID2GO_Tps_nr, 0.05)
RGL_Tbi_RT_MB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tbi_RT_MB_rankedlist_wFC_GL, geneID2GO_Tbi_nr, 0.05)
RGL_Tce_RT_MB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tce_RT_MB_rankedlist_wFC_GL, geneID2GO_Tce_nr, 0.05)
RGL_Tcm_RT_MB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tcm_RT_MB_rankedlist_wFC_GL, geneID2GO_Tcm_nr, 0.05)
RGL_Tpa_RT_MB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tpa_RT_MB_rankedlist_wFC_GL, geneID2GO_Tpa_nr, 0.05)
RGL_Tps_RT_MB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tps_RT_MB_rankedlist_wFC_GL, geneID2GO_Tps_nr, 0.05)
RGL_Tbi_LG_MB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tbi_LG_MB_rankedlist_wFC_GL, geneID2GO_Tbi_nr, 0.05)
RGL_Tce_LG_MB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tce_LG_MB_rankedlist_wFC_GL, geneID2GO_Tce_nr, 0.05)
RGL_Tcm_LG_MB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tcm_LG_MB_rankedlist_wFC_GL, geneID2GO_Tcm_nr, 0.05)
RGL_Tpa_LG_MB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tpa_LG_MB_rankedlist_wFC_GL, geneID2GO_Tpa_nr, 0.05)
RGL_Tps_LG_MB_rankedlist_wFC_nr_enrich  <- run_enrichment(RGL_Tps_LG_MB_rankedlist_wFC_GL, geneID2GO_Tps_nr, 0.05)

RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tbi_nr, 0.05)
RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tce_nr, 0.05)
RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tcm_nr, 0.05)
RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tpa_nr, 0.05)
RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tps_nr, 0.05)
RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tbi_nr, 0.05)
RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tce_nr, 0.05)
RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tcm_nr, 0.05)
RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tpa_nr, 0.05)
RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tps_nr, 0.05)
RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tbi_nr, 0.05)
RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tce_nr, 0.05)
RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tcm_nr, 0.05)
RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tpa_nr, 0.05)
RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tps_nr, 0.05)

RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tbi_nr, 0.05)
RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tce_nr, 0.05)
RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tcm_nr, 0.05)
RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tpa_nr, 0.05)
RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tps_nr, 0.05)
RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tbi_nr, 0.05)
RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tce_nr, 0.05)
RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tcm_nr, 0.05)
RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tpa_nr, 0.05)
RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tps_nr, 0.05)
RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tbi_nr, 0.05)
RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tce_nr, 0.05)
RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tcm_nr, 0.05)
RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tpa_nr, 0.05)
RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich  <- run_enrichment(RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tps_nr, 0.05)

RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tbi_nr, 0.05)
RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tce_nr, 0.05)
RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tcm_nr, 0.05)
RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tpa_nr, 0.05)
RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tps_nr, 0.05)
RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tbi_nr, 0.05)
RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tce_nr, 0.05)
RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tcm_nr, 0.05)
RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tpa_nr, 0.05)
RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tps_nr, 0.05)
RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tbi_nr, 0.05)
RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tce_nr, 0.05)
RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tcm_nr, 0.05)
RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tpa_nr, 0.05)
RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tps_nr, 0.05)

RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tbi_nr, 0.05)
RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tce_nr, 0.05)
RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tcm_nr, 0.05)
RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tpa_nr, 0.05)
RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tps_nr, 0.05)
RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tbi_nr, 0.05)
RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tce_nr, 0.05)
RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tcm_nr, 0.05)
RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tpa_nr, 0.05)
RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tps_nr, 0.05)
RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tbi_nr, 0.05)
RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tce_nr, 0.05)
RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tcm_nr, 0.05)
RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tpa_nr, 0.05)
RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich  <- run_enrichment(RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tps_nr, 0.05)



### DROSO


RGL_Tbi_WB_FB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tbi_WB_FB_rankedlist_wFC_GL, geneID2GO_Tbi_DROSO, 0.05)
RGL_Tce_WB_FB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tce_WB_FB_rankedlist_wFC_GL, geneID2GO_Tce_DROSO, 0.05)
RGL_Tcm_WB_FB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tcm_WB_FB_rankedlist_wFC_GL, geneID2GO_Tcm_DROSO, 0.05)
RGL_Tpa_WB_FB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tpa_WB_FB_rankedlist_wFC_GL, geneID2GO_Tpa_DROSO, 0.05)
RGL_Tps_WB_FB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tps_WB_FB_rankedlist_wFC_GL, geneID2GO_Tps_DROSO, 0.05)
RGL_Tbi_RT_FB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tbi_RT_FB_rankedlist_wFC_GL, geneID2GO_Tbi_DROSO, 0.05)
RGL_Tce_RT_FB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tce_RT_FB_rankedlist_wFC_GL, geneID2GO_Tce_DROSO, 0.05)
RGL_Tcm_RT_FB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tcm_RT_FB_rankedlist_wFC_GL, geneID2GO_Tcm_DROSO, 0.05)
RGL_Tpa_RT_FB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tpa_RT_FB_rankedlist_wFC_GL, geneID2GO_Tpa_DROSO, 0.05)
RGL_Tps_RT_FB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tps_RT_FB_rankedlist_wFC_GL, geneID2GO_Tps_DROSO, 0.05)
RGL_Tbi_LG_FB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tbi_LG_FB_rankedlist_wFC_GL, geneID2GO_Tbi_DROSO, 0.05)
RGL_Tce_LG_FB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tce_LG_FB_rankedlist_wFC_GL, geneID2GO_Tce_DROSO, 0.05)
RGL_Tcm_LG_FB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tcm_LG_FB_rankedlist_wFC_GL, geneID2GO_Tcm_DROSO, 0.05)
RGL_Tpa_LG_FB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tpa_LG_FB_rankedlist_wFC_GL, geneID2GO_Tpa_DROSO, 0.05)
RGL_Tps_LG_FB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tps_LG_FB_rankedlist_wFC_GL, geneID2GO_Tps_DROSO, 0.05)

RGL_Tbi_WB_MB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tbi_WB_MB_rankedlist_wFC_GL, geneID2GO_Tbi_DROSO, 0.05)
RGL_Tce_WB_MB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tce_WB_MB_rankedlist_wFC_GL, geneID2GO_Tce_DROSO, 0.05)
RGL_Tcm_WB_MB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tcm_WB_MB_rankedlist_wFC_GL, geneID2GO_Tcm_DROSO, 0.05)
RGL_Tpa_WB_MB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tpa_WB_MB_rankedlist_wFC_GL, geneID2GO_Tpa_DROSO, 0.05)
RGL_Tps_WB_MB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tps_WB_MB_rankedlist_wFC_GL, geneID2GO_Tps_DROSO, 0.05)
RGL_Tbi_RT_MB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tbi_RT_MB_rankedlist_wFC_GL, geneID2GO_Tbi_DROSO, 0.05)
RGL_Tce_RT_MB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tce_RT_MB_rankedlist_wFC_GL, geneID2GO_Tce_DROSO, 0.05)
RGL_Tcm_RT_MB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tcm_RT_MB_rankedlist_wFC_GL, geneID2GO_Tcm_DROSO, 0.05)
RGL_Tpa_RT_MB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tpa_RT_MB_rankedlist_wFC_GL, geneID2GO_Tpa_DROSO, 0.05)
RGL_Tps_RT_MB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tps_RT_MB_rankedlist_wFC_GL, geneID2GO_Tps_DROSO, 0.05)
RGL_Tbi_LG_MB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tbi_LG_MB_rankedlist_wFC_GL, geneID2GO_Tbi_DROSO, 0.05)
RGL_Tce_LG_MB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tce_LG_MB_rankedlist_wFC_GL, geneID2GO_Tce_DROSO, 0.05)
RGL_Tcm_LG_MB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tcm_LG_MB_rankedlist_wFC_GL, geneID2GO_Tcm_DROSO, 0.05)
RGL_Tpa_LG_MB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tpa_LG_MB_rankedlist_wFC_GL, geneID2GO_Tpa_DROSO, 0.05)
RGL_Tps_LG_MB_rankedlist_wFC_DROSO_enrich  <- run_enrichment(RGL_Tps_LG_MB_rankedlist_wFC_GL, geneID2GO_Tps_DROSO, 0.05)


RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tbi_DROSO, 0.05)
RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tce_DROSO, 0.05)
RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tcm_DROSO, 0.05)
RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tpa_DROSO, 0.05)
RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tps_DROSO, 0.05)
RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tbi_DROSO, 0.05)
RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tce_DROSO, 0.05)
RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tcm_DROSO, 0.05)
RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tpa_DROSO, 0.05)
RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tps_DROSO, 0.05)
RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tbi_DROSO, 0.05)
RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tce_DROSO, 0.05)
RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tcm_DROSO, 0.05)
RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tpa_DROSO, 0.05)
RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_first_GL, geneID2GO_Tps_DROSO, 0.05)

RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tbi_DROSO, 0.05)
RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tce_DROSO, 0.05)
RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tcm_DROSO, 0.05)
RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tpa_DROSO, 0.05)
RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tps_DROSO, 0.05)
RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tbi_DROSO, 0.05)
RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tce_DROSO, 0.05)
RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tcm_DROSO, 0.05)
RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tpa_DROSO, 0.05)
RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tps_DROSO, 0.05)
RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tbi_DROSO, 0.05)
RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tce_DROSO, 0.05)
RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tcm_DROSO, 0.05)
RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tpa_DROSO, 0.05)
RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich  <- run_enrichment(RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_first_GL, geneID2GO_Tps_DROSO, 0.05)

RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tbi_DROSO, 0.05)
RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tce_DROSO, 0.05)
RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tcm_DROSO, 0.05)
RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tpa_DROSO, 0.05)
RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tps_DROSO, 0.05)
RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tbi_DROSO, 0.05)
RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tce_DROSO, 0.05)
RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tcm_DROSO, 0.05)
RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tpa_DROSO, 0.05)
RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tps_DROSO, 0.05)
RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tbi_DROSO, 0.05)
RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tce_DROSO, 0.05)
RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tcm_DROSO, 0.05)
RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tpa_DROSO, 0.05)
RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_last_GL, geneID2GO_Tps_DROSO, 0.05)

RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tbi_DROSO, 0.05)
RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tce_DROSO, 0.05)
RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tcm_DROSO, 0.05)
RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tpa_DROSO, 0.05)
RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tps_DROSO, 0.05)
RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tbi_DROSO, 0.05)
RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tce_DROSO, 0.05)
RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tcm_DROSO, 0.05)
RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tpa_DROSO, 0.05)
RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tps_DROSO, 0.05)
RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tbi_DROSO, 0.05)
RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tce_DROSO, 0.05)
RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tcm_DROSO, 0.05)
RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tpa_DROSO, 0.05)
RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich  <- run_enrichment(RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_last_GL, geneID2GO_Tps_DROSO, 0.05)



###################################################################################################################################
#### Venns
############################


five_sp_DE_venn <- function(Tbi,Tce,Tcm,Tpa,Tps,title){
	venny.plot <- venn.diagram(
	list("Tbi" = Tbi, "Tce" = Tce, "Tcm" = Tcm, "Tpa" = Tpa, "Tps" = Tps ), filename = NULL,
                            fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                            cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                            margin = 0.6, cat.dist = 0.23, main = title, main.pos = c(0.5,0.8), main.cex = 2, main.fontface = "bold", cat.cex = 2)

	return(venny.plot)
}


## FB
RGL_WB_FB_wFC_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
"FB 10 sp wFC orths, BP, whole-body, GSEA"
)

RGL_RT_FB_wFC_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
"FB 10 sp wFC orths, BP, rep tract, GSEA"
)

RGL_LG_FB_wFC_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
"FB 10 sp wFC orths, BP, legs, GSEA"
)


## MB
RGL_WB_MB_wFC_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
"MB 10 sp wFC orths, BP, whole-body, GSEA"
)

RGL_RT_MB_wFC_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
"MB 10 sp wFC orths, BP, rep tract, GSEA"
)

RGL_LG_MB_wFC_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
"MB 10 sp wFC orths, BP, legs, GSEA"
)


## FB
RGL_WB_FB_wFC_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
"FB 10 sp wFC orths, BP, whole-body, GSEA"
)

RGL_RT_FB_wFC_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
"FB 10 sp wFC orths, BP, rep tract, GSEA"
)

RGL_LG_FB_wFC_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
"FB 10 sp wFC orths, BP, legs, GSEA"
)


## MB
RGL_WB_MB_wFC_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
"MB 10 sp wFC orths, BP, whole-body, GSEA"
)

RGL_RT_MB_wFC_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
"MB 10 sp wFC orths, BP, rep tract, GSEA"
)

RGL_LG_MB_wFC_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
"MB 10 sp wFC orths, BP, legs, GSEA"
)



############# SEX-asex


RGL_WB_sig_wFCFB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_first, 10sp, WB"
)

RGL_RT_sig_wFCFB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_first, 10sp, RT"
)

RGL_LG_sig_wFCFB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_first, 10sp, LG"
)

RGL_WB_sig_wFCMB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_first, 10sp, WB"
)

RGL_RT_sig_wFCMB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_first, 10sp, RT"
)

RGL_LG_sig_wFCMB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_first, 10sp, LG"
)


RGL_WB_sig_wFCFB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_last, 10sp, WB"
)

RGL_RT_sig_wFCFB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_last, 10sp, RT"
)

RGL_LG_sig_wFCFB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_last, 10sp, LG"
)

RGL_WB_sig_wFCMB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_last, 10sp, WB"
)

RGL_RT_sig_wFCMB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_last, 10sp, RT"
)

RGL_LG_sig_wFCMB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_last, 10sp, LG"
)

RGL_WB_sig_wFCFB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_first, 10sp, WB"
)

RGL_RT_sig_wFCFB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_first, 10sp, RT"
)

RGL_LG_sig_wFCFB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_first, 10sp, LG"
)

RGL_WB_sig_wFCMB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_first, 10sp, WB"
)

RGL_RT_sig_wFCMB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_first, 10sp, RT"
)

RGL_LG_sig_wFCMB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_first, 10sp, LG"
)



RGL_WB_sig_wFCFB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_last, 10sp, WB"
)

RGL_RT_sig_wFCFB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_last, 10sp, RT"
)

RGL_LG_sig_wFCFB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_last, 10sp, LG"
)

RGL_WB_sig_wFCMB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_last, 10sp, WB"
)

RGL_RT_sig_wFCMB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_last, 10sp, RT"
)

RGL_LG_sig_wFCMB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_last, 10sp, LG"
)


### output


venn_RGL_FB_wFC_DROSO_enrich_venn_GSEA_BP_GO <- arrangeGrob(
gTree(children=RGL_WB_FB_wFC_DROSO_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_RT_FB_wFC_DROSO_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_LG_FB_wFC_DROSO_enrich_venn_GSEA_BP_GO), ncol=1)

ggsave(file="venn_RGL_FB_wFC_DROSO_enrich_venn_GSEA_BP_GO.png", venn_RGL_FB_wFC_DROSO_enrich_venn_GSEA_BP_GO, width = 6, height = 18)

venn_RGL_MB_wFC_DROSO_enrich_venn_GSEA_BP_GO <- arrangeGrob(
gTree(children=RGL_WB_MB_wFC_DROSO_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_RT_MB_wFC_DROSO_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_LG_MB_wFC_DROSO_enrich_venn_GSEA_BP_GO), ncol=1)

ggsave(file="venn_RGL_MB_wFC_DROSO_enrich_venn_GSEA_BP_GO.png", venn_RGL_MB_wFC_DROSO_enrich_venn_GSEA_BP_GO, width = 6, height = 18)


venn_RGL_FB_wFC_nr_enrich_venn_GSEA_BP_GO <- arrangeGrob(
gTree(children=RGL_WB_FB_wFC_nr_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_RT_FB_wFC_nr_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_LG_FB_wFC_nr_enrich_venn_GSEA_BP_GO), ncol=1)

ggsave(file="venn_RGL_FB_wFC_nr_enrich_venn_GSEA_BP_GO.png", venn_RGL_FB_wFC_nr_enrich_venn_GSEA_BP_GO, width = 6, height = 18)

venn_RGL_MB_wFC_nr_enrich_venn_GSEA_BP_GO <- arrangeGrob(
gTree(children=RGL_WB_MB_wFC_nr_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_RT_MB_wFC_nr_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_LG_MB_wFC_nr_enrich_venn_GSEA_BP_GO), ncol=1)

ggsave(file="venn_RGL_MB_wFC_nr_enrich_venn_GSEA_BP_GO.png", venn_RGL_MB_wFC_nr_enrich_venn_GSEA_BP_GO, width = 6, height = 18)


venn_RGL_sig_wFCMB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO <- arrangeGrob(
gTree(children=RGL_WB_sig_wFCMB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_RT_sig_wFCMB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_LG_sig_wFCMB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO), ncol=1)

ggsave(file="venn_RGL_sig_wFCMB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO.png", venn_RGL_sig_wFCMB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO, width = 6, height = 18)


venn_RGL_sig_wFCFB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO <- arrangeGrob(
gTree(children=RGL_WB_sig_wFCFB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_RT_sig_wFCFB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_LG_sig_wFCFB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO), ncol=1)

ggsave(file="venn_RGL_sig_wFCFB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO.png", venn_RGL_sig_wFCFB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO, width = 6, height = 18)


venn_RGL_sig_wFCMB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO <- arrangeGrob(
gTree(children=RGL_WB_sig_wFCMB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_RT_sig_wFCMB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_LG_sig_wFCMB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO), ncol=1)

ggsave(file="venn_RGL_sig_wFCMB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO.png", venn_RGL_sig_wFCMB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO, width = 6, height = 18)


venn_RGL_sig_wFCFB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO <- arrangeGrob(
gTree(children=RGL_WB_sig_wFCFB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_RT_sig_wFCFB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_LG_sig_wFCFB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO), ncol=1)

ggsave(file="venn_RGL_sig_wFCFB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO.png", venn_RGL_sig_wFCFB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO, width = 6, height = 18)


venn_RGL_sig_wFCMB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO <- arrangeGrob(
gTree(children=RGL_WB_sig_wFCMB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_RT_sig_wFCMB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_LG_sig_wFCMB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO), ncol=1)

ggsave(file="venn_RGL_sig_wFCMB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO.png", venn_RGL_sig_wFCMB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO, width = 6, height = 18)

venn_RGL_sig_wFCFB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO <- arrangeGrob(
gTree(children=RGL_WB_sig_wFCFB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_RT_sig_wFCFB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_LG_sig_wFCFB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO), ncol=1)

ggsave(file="venn_RGL_sig_wFCFB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO.png", venn_RGL_sig_wFCFB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO, width = 6, height = 18)



venn_RGL_sig_wFCMB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO <- arrangeGrob(
gTree(children=RGL_WB_sig_wFCMB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_RT_sig_wFCMB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_LG_sig_wFCMB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO), ncol=1)

ggsave(file="venn_RGL_sig_wFCMB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO.png", venn_RGL_sig_wFCMB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO, width = 6, height = 18)


venn_RGL_sig_wFCFB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO <- arrangeGrob(
gTree(children=RGL_WB_sig_wFCFB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_RT_sig_wFCFB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_LG_sig_wFCFB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO), ncol=1)

ggsave(file="venn_RGL_sig_wFCFB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO.png", venn_RGL_sig_wFCFB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO, width = 6, height = 18)





########################################################################################################################################################################
####### output session info

writeLines(capture.output(sessionInfo()), "TopGO_10sp.R_sessionInfo.txt")





