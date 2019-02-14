
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



#### load annotation from blast2GO.
#### Note was parsed with:  B2G_to_topGO.py (see 5_GO_term_analyses.sh)
#### (from TopGO man) "It is sufficient for the mapping to contain only the most specific GO annotations."

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



###############################################################################################################################################
#### read in tables with genename and rank


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

setwd("../../")
dir.create("Output")
dir.create("Output/TopGO_out")
setwd("Output/TopGO_out/")


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
"FB RBBH wFC, BP, whole-body, GSEA"
)

RGL_RT_FB_wFC_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
"FB RBBH wFC, BP, rep tract, GSEA"
)

RGL_LG_FB_wFC_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_FB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
"FB RBBH wFC, BP, legs, GSEA"
)


## MB
RGL_WB_MB_wFC_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
"MB RBBH wFC, BP, whole-body, GSEA"
)

RGL_RT_MB_wFC_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
"MB RBBH wFC, BP, rep tract, GSEA"
)

RGL_LG_MB_wFC_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_MB_rankedlist_wFC_DROSO_enrich$sig_GSEA_BP_GO,
"MB RBBH wFC, BP, legs, GSEA"
)


## FB
RGL_WB_FB_wFC_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
"FB RBBH wFC, BP, whole-body, GSEA"
)

RGL_RT_FB_wFC_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
"FB RBBH wFC, BP, rep tract, GSEA"
)

RGL_LG_FB_wFC_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_FB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
"FB RBBH wFC, BP, legs, GSEA"
)


## MB
RGL_WB_MB_wFC_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
"MB RBBH wFC, BP, whole-body, GSEA"
)

RGL_RT_MB_wFC_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
"MB RBBH wFC, BP, rep tract, GSEA"
)

RGL_LG_MB_wFC_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_MB_rankedlist_wFC_nr_enrich$sig_GSEA_BP_GO,
"MB RBBH wFC, BP, legs, GSEA"
)


############# SEX-asex
### for Tpa MB WB there were only 5 genes - so the enrichment could not run. 
### replace with missing_GOs 

missing_GOs = character(0)



RGL_WB_sig_wFCFB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_first, RBBH, WB"
)

RGL_RT_sig_wFCFB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_first, RBBH, RT"
)

RGL_LG_sig_wFCFB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_first, RBBH, LG"
)

RGL_WB_sig_wFCMB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
#RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
missing_GOs,
RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_first, RBBH, WB"
)

RGL_RT_sig_wFCMB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_first, RBBH, RT"
)

RGL_LG_sig_wFCMB_positive_FC_first_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_first, RBBH, LG"
)



RGL_WB_sig_wFCFB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_last, RBBH, WB"
)

RGL_RT_sig_wFCFB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_last, RBBH, RT"
)

RGL_LG_sig_wFCFB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_last, RBBH, LG"
)

RGL_WB_sig_wFCMB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
#RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
missing_GOs,
RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_last, RBBH, WB"
)

RGL_RT_sig_wFCMB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_last, RBBH, RT"
)

RGL_LG_sig_wFCMB_positive_FC_last_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_last, RBBH, LG"
)



RGL_WB_sig_wFCFB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_first, RBBH, WB"
)

RGL_RT_sig_wFCFB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_first, RBBH, RT"
)

RGL_LG_sig_wFCFB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_first, RBBH, LG"
)

RGL_WB_sig_wFCMB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
#RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
missing_GOs,
RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_first, RBBH, WB"
)

RGL_RT_sig_wFCMB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_first, RBBH, RT"
)

RGL_LG_sig_wFCMB_positive_FC_first_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_first, RBBH, LG"
)



RGL_WB_sig_wFCFB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_last, RBBH, WB"
)

RGL_RT_sig_wFCFB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_last, RBBH, RT"
)

RGL_LG_sig_wFCFB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCFB, SA positive_FC_last, RBBH, LG"
)

RGL_WB_sig_wFCMB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
#RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
missing_GOs,
RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_last, RBBH, WB"
)

RGL_RT_sig_wFCMB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_last, RBBH, RT"
)

RGL_LG_sig_wFCMB_positive_FC_last_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$sig_GSEA_BP_GO,
"sig_wFCMB, SA positive_FC_last, RBBH, LG"
)


### output venns


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



#### superexact test


### runs superexact test

run_superexact_GOs_BP <- function(basename){

	
	# GSEA
	GSEA_Tbi_GO_sig = eval(parse(text=paste('RGL_Tbi_',basename,'$','sig_GSEA_BP_GO',sep='')))
	GSEA_Tce_GO_sig = eval(parse(text=paste('RGL_Tce_',basename,'$','sig_GSEA_BP_GO',sep='')))
	GSEA_Tcm_GO_sig = eval(parse(text=paste('RGL_Tcm_',basename,'$','sig_GSEA_BP_GO',sep='')))
	GSEA_Tpa_GO_sig = eval(parse(text=paste('RGL_Tpa_',basename,'$','sig_GSEA_BP_GO',sep='')))	
	GSEA_Tps_GO_sig = eval(parse(text=paste('RGL_Tps_',basename,'$','sig_GSEA_BP_GO',sep='')))

	GSEA_all_GO_sig <- list(GSEA_Tbi_GO_sig,GSEA_Tce_GO_sig,GSEA_Tcm_GO_sig,GSEA_Tpa_GO_sig,GSEA_Tps_GO_sig)

	GSEA_Tbi_GO_allused = eval(parse(text=paste('RGL_Tbi_',basename,'$','GO_term_use_BP_list',sep='')))	
	GSEA_Tce_GO_allused = eval(parse(text=paste('RGL_Tce_',basename,'$','GO_term_use_BP_list',sep='')))	
	GSEA_Tcm_GO_allused = eval(parse(text=paste('RGL_Tcm_',basename,'$','GO_term_use_BP_list',sep='')))	
	GSEA_Tpa_GO_allused = eval(parse(text=paste('RGL_Tpa_',basename,'$','GO_term_use_BP_list',sep='')))	
	GSEA_Tps_GO_allused = eval(parse(text=paste('RGL_Tps_',basename,'$','GO_term_use_BP_list',sep='')))	
	
	GSEA_all_GO_allused <- list(GSEA_Tbi_GO_allused,GSEA_Tce_GO_allused,GSEA_Tcm_GO_allused,GSEA_Tpa_GO_allused,GSEA_Tps_GO_allused)	

	## lens of GOs differ to make more conservative I will use the size of the GO intersection of BP GOs
	GSEA_interlength <- length(Intersect(GSEA_all_GO_allused))
	
	print(GSEA_interlength)
	
	sup_out_GSEA = summary(supertest(GSEA_all_GO_sig, n= GSEA_interlength))

	out_list = list("sup_out_GSEA" = sup_out_GSEA)
	return(out_list)

}

### run ## some don't work as sets are too small.



sup_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich <- run_superexact_GOs_BP("WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich")$sup_out_GSEA
sup_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich <- run_superexact_GOs_BP("RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich")$sup_out_GSEA
sup_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich <- run_superexact_GOs_BP("LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich")$sup_out_GSEA
sup_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich <- run_superexact_GOs_BP("WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich")$sup_out_GSEA
sup_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich <- run_superexact_GOs_BP("RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich")$sup_out_GSEA
sup_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich <- run_superexact_GOs_BP("LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich")$sup_out_GSEA
sup_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich <- run_superexact_GOs_BP("WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich")$sup_out_GSEA
sup_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich <- run_superexact_GOs_BP("RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich")$sup_out_GSEA
sup_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich <- run_superexact_GOs_BP("LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich")$sup_out_GSEA
sup_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich <- run_superexact_GOs_BP("WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich")$sup_out_GSEA
sup_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich <- run_superexact_GOs_BP("RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich")$sup_out_GSEA
sup_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich <- run_superexact_GOs_BP("LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich")$sup_out_GSEA

###

sup_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich <- run_superexact_GOs_BP("WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich")$sup_out_GSEA
sup_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich <- run_superexact_GOs_BP("RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich")$sup_out_GSEA
sup_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich <- run_superexact_GOs_BP("LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich")$sup_out_GSEA
sup_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich <- run_superexact_GOs_BP("WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich")$sup_out_GSEA
sup_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich <- run_superexact_GOs_BP("RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich")$sup_out_GSEA
sup_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich <- run_superexact_GOs_BP("LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich")$sup_out_GSEA
sup_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich <- run_superexact_GOs_BP("WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich")$sup_out_GSEA
sup_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich <- run_superexact_GOs_BP("RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich")$sup_out_GSEA
sup_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich <- run_superexact_GOs_BP("LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich")$sup_out_GSEA
sup_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich <- run_superexact_GOs_BP("WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich")$sup_out_GSEA
sup_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich <- run_superexact_GOs_BP("RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich")$sup_out_GSEA
sup_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich <- run_superexact_GOs_BP("LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich")$sup_out_GSEA




## output


write.csv(sup_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$Table, file="WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich.csv", row.names=FALSE)
write.csv(sup_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$Table, file="RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich.csv", row.names=FALSE)
write.csv(sup_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$Table, file="LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich.csv", row.names=FALSE)
write.csv(sup_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$Table, file="WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich.csv", row.names=FALSE)
write.csv(sup_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$Table, file="RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich.csv", row.names=FALSE)
write.csv(sup_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$Table, file="LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich.csv", row.names=FALSE)
write.csv(sup_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$Table, file="WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich.csv", row.names=FALSE)
write.csv(sup_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$Table, file="RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich.csv", row.names=FALSE)
write.csv(sup_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$Table, file="LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich.csv", row.names=FALSE)
write.csv(sup_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$Table, file="WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich.csv", row.names=FALSE)
write.csv(sup_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$Table, file="RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich.csv", row.names=FALSE)
write.csv(sup_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$Table, file="LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich.csv", row.names=FALSE)
###

write.csv(sup_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$Table, file="WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich.csv", row.names=FALSE)
write.csv(sup_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$Table, file="RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich.csv", row.names=FALSE)
write.csv(sup_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$Table, file="LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich.csv", row.names=FALSE)
write.csv(sup_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$Table, file="WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich.csv", row.names=FALSE)
write.csv(sup_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$Table, file="RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich.csv", row.names=FALSE)
write.csv(sup_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$Table, file="LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich.csv", row.names=FALSE)
write.csv(sup_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$Table, file="WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich.csv", row.names=FALSE)
write.csv(sup_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$Table, file="RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich.csv", row.names=FALSE)
write.csv(sup_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$Table, file="LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich.csv", row.names=FALSE)
write.csv(sup_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$Table, file="WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich.csv", row.names=FALSE)
write.csv(sup_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$Table, file="RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich.csv", row.names=FALSE)
write.csv(sup_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$Table, file="LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich.csv", row.names=FALSE)

#### GO enrichment tables - export


write.table(RGL_Tbi_WB_FB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tbi_WB_FB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_WB_FB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tce_WB_FB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_WB_FB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tcm_WB_FB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_WB_FB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tpa_WB_FB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_WB_FB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tps_WB_FB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_RT_FB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tbi_RT_FB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_RT_FB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tce_RT_FB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_RT_FB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tcm_RT_FB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_RT_FB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tpa_RT_FB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_RT_FB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tps_RT_FB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_LG_FB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tbi_LG_FB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_LG_FB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tce_LG_FB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_LG_FB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tcm_LG_FB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_LG_FB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tpa_LG_FB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_LG_FB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tps_LG_FB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)

write.table(RGL_Tbi_WB_MB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tbi_WB_MB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_WB_MB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tce_WB_MB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_WB_MB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tcm_WB_MB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_WB_MB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tpa_WB_MB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_WB_MB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tps_WB_MB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_RT_MB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tbi_RT_MB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_RT_MB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tce_RT_MB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_RT_MB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tcm_RT_MB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_RT_MB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tpa_RT_MB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_RT_MB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tps_RT_MB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_LG_MB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tbi_LG_MB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_LG_MB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tce_LG_MB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_LG_MB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tcm_LG_MB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_LG_MB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tpa_LG_MB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_LG_MB_rankedlist_wFC_nr_enrich$allRes1_BP, "RGL_Tps_LG_MB_rankedlist_wFC_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)

write.table(RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)

write.table(RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich$allRes1_BP, "RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_first_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)

write.table(RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)

write.table(RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich$allRes1_BP, "RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_last_nr_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)



write.table(RGL_Tbi_WB_FB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tbi_WB_FB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_WB_FB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tce_WB_FB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_WB_FB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tcm_WB_FB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_WB_FB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tpa_WB_FB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_WB_FB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tps_WB_FB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_RT_FB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tbi_RT_FB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_RT_FB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tce_RT_FB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_RT_FB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tcm_RT_FB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_RT_FB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tpa_RT_FB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_RT_FB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tps_RT_FB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_LG_FB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tbi_LG_FB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_LG_FB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tce_LG_FB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_LG_FB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tcm_LG_FB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_LG_FB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tpa_LG_FB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_LG_FB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tps_LG_FB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)

write.table(RGL_Tbi_WB_MB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tbi_WB_MB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_WB_MB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tce_WB_MB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_WB_MB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tcm_WB_MB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_WB_MB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tpa_WB_MB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_WB_MB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tps_WB_MB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_RT_MB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tbi_RT_MB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_RT_MB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tce_RT_MB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_RT_MB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tcm_RT_MB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_RT_MB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tpa_RT_MB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_RT_MB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tps_RT_MB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_LG_MB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tbi_LG_MB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_LG_MB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tce_LG_MB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_LG_MB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tcm_LG_MB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_LG_MB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tpa_LG_MB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_LG_MB_rankedlist_wFC_DROSO_enrich$allRes1_BP, "RGL_Tps_LG_MB_rankedlist_wFC_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)

write.table(RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)

write.table(RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich$allRes1_BP, "RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_first_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)

write.table(RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tbi_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tce_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tcm_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tpa_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tps_WB_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tbi_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tce_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tcm_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tpa_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tps_RT_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tbi_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tce_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tcm_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tpa_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tps_LG_sig_wFCFB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)

write.table(RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tbi_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tce_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tcm_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tpa_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tps_WB_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tbi_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tce_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tcm_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tpa_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tps_RT_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tbi_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tce_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tcm_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tpa_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich$allRes1_BP, "RGL_Tps_LG_sig_wFCMB_rankedlist_positive_FC_last_DROSO_enrich.txt", sep = '\t', quote = FALSE, row.names = FALSE)




########################################################################################################################################################################
####### output session info

writeLines(capture.output(sessionInfo()), "TopGO.R_sessionInfo.txt")





