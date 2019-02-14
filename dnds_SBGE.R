###### dNdS pnps and SBGE

### libs

library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(lme4)
library(car)
library(MASS)
library(fitdistrplus)


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
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
 # [1] fitdistrplus_1.0-11 npsurv_0.4-0        lsei_1.2-0          survival_2.43-1     MASS_7.3-51.1       car_3.0-2           carData_3.0-2       lme4_1.1-19         Matrix_1.2-15      
# [10] RColorBrewer_1.1-2  cowplot_0.9.3       ggplot2_3.1.0      

# loaded via a namespace (and not attached):
 # [1] zip_1.0.0         Rcpp_1.0.0        cellranger_1.1.0  pillar_1.3.0      compiler_3.5.1    nloptr_1.2.1      plyr_1.8.4        bindr_0.1.1       tools_3.5.1       forcats_0.3.0    
# [11] digest_0.6.18     tibble_1.4.2      gtable_0.2.0      nlme_3.1-137      lattice_0.20-38   pkgconfig_2.0.2   rlang_0.3.0.1     openxlsx_4.1.0    curl_3.2          haven_1.1.2      
# [21] rio_0.5.10        bindrcpp_0.2.2    withr_2.1.2       dplyr_0.7.8       hms_0.4.2         grid_3.5.1        tidyselect_0.2.5  glue_1.3.0        data.table_1.11.8 R6_2.3.0         
# [31] readxl_1.1.0      foreign_0.8-71    minqa_1.2.4       purrr_0.2.5       magrittr_1.5      scales_1.0.0      splines_3.5.1     abind_1.4-5       assertthat_0.2.0  colorspace_1.3-2 
# [41] labeling_0.3      lazyeval_0.2.1    munsell_0.5.0     crayon_1.3.4     
# > 

#### remember that pN/pS are not really valid here as different species were calced from different N of indv

### data
setwd("Data/dNdS")
dat1_WB = read.csv("ALL_WB_10sp_disp_allsepar_wCPM_SB_asex_2_dnds_exp.csv")
dat1_RT = read.csv("ALL_RT_10sp_disp_allsepar_wCPM_SB_asex_2_dnds_exp.csv")
dat1_LG = read.csv("ALL_LG_10sp_disp_allsepar_wCPM_SB_asex_2_dnds_exp.csv")


### output 

setwd("../..")
dir.create("Output")
dir.create("Output/dNdS")
setwd("Output/dNdS")


dat1_WB_c <- na.omit(dat1_WB)
dat1_RT_c <- na.omit(dat1_RT)
dat1_LG_c <- na.omit(dat1_LG)

## make binomial dNdS
dat1_WB_c$dnds_bi = ifelse(dat1_WB_c$dnds > 0.0001 , 1, 0)
dat1_WB_c$dnds_bi_nonzero = ifelse(dat1_WB_c$dnds > 0.0001 , dat1_WB_c$dnds, NA)
dat1_RT_c$dnds_bi = ifelse(dat1_RT_c$dnds > 0.0001 , 1, 0)
dat1_RT_c$dnds_bi_nonzero = ifelse(dat1_RT_c$dnds > 0.0001 , dat1_RT_c$dnds, NA)
dat1_LG_c$dnds_bi = ifelse(dat1_LG_c$dnds > 0.0001 , 1, 0)
dat1_LG_c$dnds_bi_nonzero = ifelse(dat1_LG_c$dnds > 0.0001 , dat1_LG_c$dnds, NA)

dat1_WB_c$sexbias_ord  <- ordered(dat1_WB_c$sexbias, levels = c("female_biased", "male_biased", "Unbiased"))
dat1_RT_c$sexbias_ord  <- ordered(dat1_RT_c$sexbias, levels = c("female_biased", "male_biased", "Unbiased"))
dat1_LG_c$sexbias_ord  <- ordered(dat1_LG_c$sexbias, levels = c("female_biased", "male_biased", "Unbiased"))
dat1_WB_c$sexbias2_ord  <- ordered(dat1_WB_c$sexbias2, levels = c("female_biased", "male_biased", "Unbiased"))
dat1_RT_c$sexbias2_ord  <- ordered(dat1_RT_c$sexbias2, levels = c("female_biased", "male_biased", "Unbiased"))
dat1_LG_c$sexbias2_ord  <- ordered(dat1_LG_c$sexbias2, levels = c("female_biased", "male_biased", "Unbiased"))

dat1_WB_c$sp_ord = ordered(dat1_WB_c$sp, levels = c(
"Tbi", "Tte", "Tce", "Tms", "Tps", "Tdi", "Tcm", "Tsi", "Tpa", "Tge"))
dat1_RT_c$sp_ord = ordered(dat1_RT_c$sp, levels = c(
"Tbi", "Tte", "Tce", "Tms", "Tps", "Tdi", "Tcm", "Tsi", "Tpa", "Tge"))
dat1_LG_c$sp_ord = ordered(dat1_LG_c$sp, levels = c(
"Tbi", "Tte", "Tce", "Tms", "Tps", "Tdi", "Tcm", "Tsi", "Tpa", "Tge"))

dat1_WB_c$sp_type_ord  <- ordered(dat1_WB_c$sp_type, levels = c("SF", "AF"))
dat1_RT_c$sp_type_ord  <- ordered(dat1_RT_c$sp_type, levels = c("SF", "AF"))
dat1_LG_c$sp_type_ord  <- ordered(dat1_LG_c$sp_type, levels = c("SF", "AF"))

head(dat1_RT_c)

max_dnds <- max(max(dat1_WB$dnds, na.rm = TRUE),
max(dat1_RT$dnds, na.rm = TRUE),
max(dat1_LG$dnds, na.rm = TRUE))

#### Split

# male-biased/unbiased 

dat1_WB_MB2_c <- subset(dat1_WB_c, dat1_WB_c$sexbias2 != "female_biased")
dat1_WB_MB2_c$sexbias2 <- droplevels(dat1_WB_MB2_c$sexbias2)
dat1_RT_MB2_c <- subset(dat1_RT_c, dat1_RT_c$sexbias2 != "female_biased")
dat1_RT_MB2_c$sexbias2 <- droplevels(dat1_RT_MB2_c$sexbias2)
dat1_LG_MB2_c <- subset(dat1_LG_c, dat1_LG_c$sexbias2 != "female_biased")
dat1_LG_MB2_c$sexbias2 <- droplevels(dat1_LG_MB2_c$sexbias2)

# female-biased/unbiased

dat1_WB_FB2_c <- subset(dat1_WB_c, dat1_WB_c$sexbias2 != "male_biased")
dat1_WB_FB2_c$sexbias2 <- droplevels(dat1_WB_FB2_c$sexbias2)
dat1_RT_FB2_c <- subset(dat1_RT_c, dat1_RT_c$sexbias2 != "male_biased")
dat1_RT_FB2_c$sexbias2 <- droplevels(dat1_RT_FB2_c$sexbias2)
dat1_LG_FB2_c <- subset(dat1_LG_c, dat1_LG_c$sexbias2 != "male_biased")
dat1_LG_FB2_c$sexbias2 <- droplevels(dat1_LG_FB2_c$sexbias2)


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### boxplots by tiss_by_species (I.e. classification is on each sp independently)

plot_box_dnds_allsp = function(df,Tiss){
	exp_labels <- c("Convergent", "Lineage-\nspecific", "None" )

	P1 <- ggplot(df, aes(sp_ord, dnds)) + 
	theme_classic() +
	geom_boxplot(aes(fill = factor(sexbias)),position=position_dodge(0.6), width = 0.5, outlier.size = 0.5) +
	ylim(0,max(max_dnds, na.rm = TRUE) + 0.1) +
	ylab ("dN/dS") +
	xlab ("Gene-class by species") + 
	scale_fill_manual(values=c("firebrick2", "royalblue2", "grey")) + labs(fill='')  +
	ggtitle(paste(Tiss, "| SB class = FDR <0.05"))
	
	P2 <- ggplot(df, aes(sp_ord, dnds)) + 
	theme_classic() +
	geom_boxplot(aes(fill = factor(sexbias2)),position=position_dodge(0.6), width = 0.5, outlier.size = 0.5) +
	ylim(0,max(max_dnds, na.rm = TRUE) + 0.1) +
	ylab ("dN/dS") +
	xlab ("Gene-class by species") + 
	scale_fill_manual(values=c("firebrick2", "royalblue2", "grey")) + labs(fill='')  +
	ggtitle(paste(Tiss, "| SB class = FDR <0.05 FC > 2")	)
	
	out_plots = list("SB" = P1, "SBwFC" = P2)
	return(out_plots)
}


# plots

P_dnds_WB_allspwFC = plot_box_dnds_allsp(dat1_WB_c,"WB")$SBwFC
P_dnds_RT_allspwFC = plot_box_dnds_allsp(dat1_RT_c,"RT")$SBwFC
P_dnds_LG_allspwFC = plot_box_dnds_allsp(dat1_LG_c,"LG")$SBwFC


pdf("dnds_WB_allspwFC.pdf", width = 10, height = 5)
P_dnds_WB_allspwFC
dev.off()
getwd() ## where has my plot 

pdf("dnds_RT_allspwFC.pdf", width = 10, height = 5)
P_dnds_RT_allspwFC
dev.off()
getwd() ## where has my plot 

pdf("dnds_LG_allspwFC.pdf", width = 10, height = 5)
P_dnds_LG_allspwFC
dev.off()
getwd() ## where has my plot 


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### models


## horrible dists - so doing a binomial model and and non-zero model

### but what model for non-zero models?

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### SEX - biased

fit_n_WB_dnds  <- fitdist(as.numeric(na.omit(dat1_WB_c$dnds_bi_nonzero)), "norm")
fit_ln_WB_dnds  <- fitdist(as.numeric(na.omit(dat1_WB_c$dnds_bi_nonzero)), "lnorm")
fit_g_WB_dnds  <- fitdist(as.numeric(na.omit(dat1_WB_c$dnds_bi_nonzero)), "gamma")
fit_n_RT_dnds  <- fitdist(as.numeric(na.omit(dat1_RT_c$dnds_bi_nonzero)), "norm")
fit_ln_RT_dnds  <- fitdist(as.numeric(na.omit(dat1_RT_c$dnds_bi_nonzero)), "lnorm")
fit_g_RT_dnds  <- fitdist(as.numeric(na.omit(dat1_RT_c$dnds_bi_nonzero)), "gamma")
fit_n_LG_dnds  <- fitdist(as.numeric(na.omit(dat1_LG_c$dnds_bi_nonzero)), "norm")
fit_ln_LG_dnds  <- fitdist(as.numeric(na.omit(dat1_LG_c$dnds_bi_nonzero)), "lnorm")
fit_g_LG_dnds  <- fitdist(as.numeric(na.omit(dat1_LG_c$dnds_bi_nonzero)), "gamma")


png(filename = "dnds_dist_plots.png", width = 9, height = 14, units = "in", bg = "white", res = 300)
par(mfrow=c(3,2))
denscomp(list(fit_n_WB_dnds, fit_ln_WB_dnds, fit_g_WB_dnds), legendtext = c("norm", "lnorm", "gamma"))
qqcomp  (list(fit_n_WB_dnds, fit_ln_WB_dnds, fit_g_WB_dnds), legendtext = c("norm", "lnorm", "gamma"))
denscomp(list(fit_n_RT_dnds, fit_ln_RT_dnds, fit_g_RT_dnds), legendtext = c("norm", "lnorm", "gamma"))
qqcomp  (list(fit_n_RT_dnds, fit_ln_RT_dnds, fit_g_RT_dnds), legendtext = c("norm", "lnorm", "gamma"))
denscomp(list(fit_n_LG_dnds, fit_ln_LG_dnds, fit_g_LG_dnds), legendtext = c("norm", "lnorm", "gamma"))
qqcomp  (list(fit_n_LG_dnds, fit_ln_LG_dnds, fit_g_LG_dnds), legendtext = c("norm", "lnorm", "gamma"))
dev.off()
getwd() ## where has my plot gone....




### lood like I want a gamma for non-zero dnds 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### binomial

dat1_WB_c$sp_type <- relevel(dat1_WB_c$sp_type , "SF")
dat1_WB_c$sexbias2 <- relevel(dat1_WB_c$sexbias2 , "Unbiased")
dat1_RT_c$sp_type <- relevel(dat1_RT_c$sp_type , "SF")
dat1_RT_c$sexbias2 <- relevel(dat1_RT_c$sexbias2 , "Unbiased")
dat1_LG_c$sp_type <- relevel(dat1_LG_c$sp_type , "SF")
dat1_LG_c$sexbias2 <- relevel(dat1_LG_c$sexbias2 , "Unbiased")



### WB
## dnds
WB_dnds_m0  = glmer(dnds_bi ~ (1|OG_name), data = dat1_WB_c, family = "binomial") 
WB_dnds_m1  = glmer(dnds_bi ~ sp_type  +  (1|OG_name), data = dat1_WB_c, family = "binomial") 
WB_dnds_m2  = glmer(dnds_bi ~ sexbias2 +  (1|OG_name), data = dat1_WB_c, family = "binomial") 
WB_dnds_m3  = glmer(dnds_bi ~ sp_type +  sexbias2 +  (1|OG_name), data = dat1_WB_c, family = "binomial") 
WB_dnds_m4  = glmer(dnds_bi ~ sp_type *  sexbias2 +  (1|OG_name), data = dat1_WB_c, family = "binomial") 
WB_dnds_m4a = glm  (dnds_bi ~ sp_type *  sexbias2, data = dat1_WB_c, family = "binomial") 


drop1(WB_dnds_m4,test="Chisq")
drop1(WB_dnds_m3,.~.,test="Chisq")

WB_dnds_m4_m4a_LRT = 2 * (as.numeric(logLik(WB_dnds_m4)) - as.numeric(logLik(WB_dnds_m4a)))
WB_dnds_m4_m4a_LRT_p <- pchisq(WB_dnds_m4_m4a_LRT, df=1, lower.tail=F) 

writeLines(capture.output(drop1(WB_dnds_m3,.~.,test="Chisq")), "WB_dnds_m3_lrt_fixed.txt")
writeLines(capture.output(drop1(WB_dnds_m4,.~.,test="Chisq")), "WB_dnds_m4_lrt_fixed.txt")
writeLines(capture.output(c("lrt_OG_name = ", WB_dnds_m4_m4a_LRT,"p_lrt_OG_name", WB_dnds_m4_m4a_LRT_p )), "WB_dnds_m4_m4a_rand.txt")



### RT
## dnds
RT_dnds_m0  = glmer(dnds_bi ~ (1|OG_name), data = dat1_RT_c, family = "binomial") 
RT_dnds_m1  = glmer(dnds_bi ~ sp_type  +  (1|OG_name), data = dat1_RT_c, family = "binomial") 
RT_dnds_m2  = glmer(dnds_bi ~ sexbias2 +  (1|OG_name), data = dat1_RT_c, family = "binomial") 
RT_dnds_m3  = glmer(dnds_bi ~ sp_type +  sexbias2 +  (1|OG_name), data = dat1_RT_c, family = "binomial") 
RT_dnds_m4  = glmer(dnds_bi ~ sp_type *  sexbias2 +  (1|OG_name), data = dat1_RT_c, family = "binomial") 
RT_dnds_m4a = glm  (dnds_bi ~ sp_type *  sexbias2, data = dat1_RT_c, family = "binomial") 

drop1(RT_dnds_m4,test="Chisq")
drop1(RT_dnds_m3,.~.,test="Chisq")

RT_dnds_m4_m4a_LRT = 2 * (as.numeric(logLik(RT_dnds_m4)) - as.numeric(logLik(RT_dnds_m4a)))
RT_dnds_m4_m4a_LRT_p <- pchisq(RT_dnds_m4_m4a_LRT, df=1, lower.tail=F) 

writeLines(capture.output(drop1(RT_dnds_m3,.~.,test="Chisq")), "RT_dnds_m3_lrt_fixed.txt")
writeLines(capture.output(drop1(RT_dnds_m4,.~.,test="Chisq")), "RT_dnds_m4_lrt_fixed.txt")
writeLines(capture.output(c("lrt_OG_name = ", RT_dnds_m4_m4a_LRT,"p_lrt_OG_name", RT_dnds_m4_m4a_LRT_p )), "RT_dnds_m4_m4a_rand.txt")


### LG
## dnds
LG_dnds_m0  = glmer(dnds_bi ~ (1|OG_name), data = dat1_LG_c, family = "binomial") 
LG_dnds_m1  = glmer(dnds_bi ~ sp_type  +  (1|OG_name), data = dat1_LG_c, family = "binomial") 
LG_dnds_m2  = glmer(dnds_bi ~ sexbias2 +  (1|OG_name), data = dat1_LG_c, family = "binomial") 
LG_dnds_m3  = glmer(dnds_bi ~ sp_type +  sexbias2 +  (1|OG_name), data = dat1_LG_c, family = "binomial") 
LG_dnds_m4  = glmer(dnds_bi ~ sp_type *  sexbias2 +  (1|OG_name), data = dat1_LG_c, family = "binomial") 
LG_dnds_m4a = glm  (dnds_bi ~ sp_type *  sexbias2, data = dat1_LG_c, family = "binomial") 

drop1(LG_dnds_m4,test="Chisq")
drop1(LG_dnds_m3,.~.,test="Chisq")

LG_dnds_m4_m4a_LRT = 2 * (as.numeric(logLik(LG_dnds_m4)) - as.numeric(logLik(LG_dnds_m4a)))
LG_dnds_m4_m4a_LRT_p <- pchisq(LG_dnds_m4_m4a_LRT, df=1, lower.tail=F) 

writeLines(capture.output(drop1(LG_dnds_m3,.~.,test="Chisq")), "LG_dnds_m3_lrt_fixed.txt")
writeLines(capture.output(drop1(LG_dnds_m4,.~.,test="Chisq")), "LG_dnds_m4_lrt_fixed.txt")
writeLines(capture.output(c("lrt_OG_name = ", LG_dnds_m4_m4a_LRT,"p_lrt_OG_name", LG_dnds_m4_m4a_LRT_p )), "LG_dnds_m4_m4a_rand.txt")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### non-zero

### WB
## dnds
WB_dnds_m0_nz  = glmer(dnds_bi_nonzero ~ (1|OG_name), data = dat1_WB_c, family = "Gamma") 
WB_dnds_m1_nz   = glmer(dnds_bi_nonzero ~ sp_type  +  (1|OG_name), data = dat1_WB_c, family = "Gamma") 
WB_dnds_m2_nz   = glmer(dnds_bi_nonzero ~ sexbias2 +  (1|OG_name), data = dat1_WB_c, family = "Gamma") 
WB_dnds_m3_nz   = glmer(dnds_bi_nonzero ~ sp_type +  sexbias2 +  (1|OG_name), data = dat1_WB_c, family = "Gamma") 
WB_dnds_m4_nz   = glmer(dnds_bi_nonzero ~ sp_type *  sexbias2 +  (1|OG_name), data = dat1_WB_c, family = "Gamma") 
WB_dnds_m4a_nz  = glm  (dnds_bi_nonzero ~ sp_type *  sexbias2, data = dat1_WB_c, family = "Gamma") 

drop1(WB_dnds_m4_nz ,test="Chisq")
drop1(WB_dnds_m3_nz ,.~.,test="Chisq")

WB_dnds_m4_m4a_LRT_nz  = 2 * (as.numeric(logLik(WB_dnds_m4_nz )) - as.numeric(logLik(WB_dnds_m4a_nz )))
WB_dnds_m4_m4a_LRT_p_nz  <- pchisq(WB_dnds_m4_m4a_LRT_nz , df=1, lower.tail=F) 

writeLines(capture.output(drop1(WB_dnds_m3_nz ,.~.,test="Chisq")), "WB_dnds_m3_lrt_fixed_nz.txt")
writeLines(capture.output(drop1(WB_dnds_m4_nz ,.~.,test="Chisq")), "WB_dnds_m4_lrt_fixed_nz.txt")
writeLines(capture.output(c("lrt_OG_name = ", WB_dnds_m4_m4a_LRT_nz ,"p_lrt_OG_name", WB_dnds_m4_m4a_LRT_p_nz  )), "WB_dnds_m4_m4a_rand_nz.txt")

### RT
## dnds
RT_dnds_m0_nz  = glmer(dnds_bi_nonzero ~ (1|OG_name), data = dat1_RT_c, family = "Gamma") 
RT_dnds_m1_nz   = glmer(dnds_bi_nonzero ~ sp_type  +  (1|OG_name), data = dat1_RT_c, family = "Gamma") 
RT_dnds_m2_nz   = glmer(dnds_bi_nonzero ~ sexbias2 +  (1|OG_name), data = dat1_RT_c, family = "Gamma") 
RT_dnds_m3_nz   = glmer(dnds_bi_nonzero ~ sp_type +  sexbias2 +  (1|OG_name), data = dat1_RT_c, family = "Gamma") 
RT_dnds_m4_nz   = glmer(dnds_bi_nonzero ~ sp_type *  sexbias2 +  (1|OG_name), data = dat1_RT_c, family = "Gamma") 
RT_dnds_m4a_nz  = glm  (dnds_bi_nonzero ~ sp_type *  sexbias2, data = dat1_RT_c, family = "Gamma") 

drop1(RT_dnds_m4_nz ,test="Chisq")
drop1(RT_dnds_m3_nz ,.~.,test="Chisq")

RT_dnds_m4_m4a_LRT_nz  = 2 * (as.numeric(logLik(RT_dnds_m4_nz )) - as.numeric(logLik(RT_dnds_m4a_nz )))
RT_dnds_m4_m4a_LRT_p_nz  <- pchisq(RT_dnds_m4_m4a_LRT_nz , df=1, lower.tail=F) 

writeLines(capture.output(drop1(RT_dnds_m3_nz ,.~.,test="Chisq")), "RT_dnds_m3_lrt_fixed_nz.txt")
writeLines(capture.output(drop1(RT_dnds_m4_nz ,.~.,test="Chisq")), "RT_dnds_m4_lrt_fixed_nz.txt")
writeLines(capture.output(c("lrt_OG_name = ", RT_dnds_m4_m4a_LRT_nz ,"p_lrt_OG_name", RT_dnds_m4_m4a_LRT_p_nz  )), "RT_dnds_m4_m4a_rand_nz.txt")

### LG
## dnds
LG_dnds_m0_nz  = glmer(dnds_bi_nonzero ~ (1|OG_name), data = dat1_LG_c, family = "Gamma") 
LG_dnds_m1_nz   = glmer(dnds_bi_nonzero ~ sp_type  +  (1|OG_name), data = dat1_LG_c, family = "Gamma") 
LG_dnds_m2_nz   = glmer(dnds_bi_nonzero ~ sexbias2 +  (1|OG_name), data = dat1_LG_c, family = "Gamma") 
LG_dnds_m3_nz   = glmer(dnds_bi_nonzero ~ sp_type +  sexbias2 +  (1|OG_name), data = dat1_LG_c, family = "Gamma") 
LG_dnds_m4_nz   = glmer(dnds_bi_nonzero ~ sp_type *  sexbias2 +  (1|OG_name), data = dat1_LG_c, family = "Gamma") 
LG_dnds_m4a_nz  = glm  (dnds_bi_nonzero ~ sp_type *  sexbias2, data = dat1_LG_c, family = "Gamma") 

drop1(LG_dnds_m4_nz ,test="Chisq")
drop1(LG_dnds_m3_nz ,.~.,test="Chisq")

LG_dnds_m4_m4a_LRT_nz  = 2 * (as.numeric(logLik(LG_dnds_m4_nz )) - as.numeric(logLik(LG_dnds_m4a_nz )))
LG_dnds_m4_m4a_LRT_p_nz  <- pchisq(LG_dnds_m4_m4a_LRT_nz , df=1, lower.tail=F) 

writeLines(capture.output(drop1(LG_dnds_m3_nz ,.~.,test="Chisq")), "LG_dnds_m3_lrt_fixed_nz.txt")
writeLines(capture.output(drop1(LG_dnds_m4_nz ,.~.,test="Chisq")), "LG_dnds_m4_lrt_fixed_nz.txt")
writeLines(capture.output(c("lrt_OG_name = ", LG_dnds_m4_m4a_LRT_nz ,"p_lrt_OG_name", LG_dnds_m4_m4a_LRT_p_nz  )), "LG_dnds_m4_m4a_rand_nz.txt")




########################################################################################################################
#### using male-bias vs unbias

fit_n_WB_MB2_dnds  <- fitdist(as.numeric(na.omit(dat1_WB_MB2_c$dnds_bi_nonzero)), "norm")
fit_ln_WB_MB2_dnds  <- fitdist(as.numeric(na.omit(dat1_WB_MB2_c$dnds_bi_nonzero)), "lnorm")
fit_g_WB_MB2_dnds  <- fitdist(as.numeric(na.omit(dat1_WB_MB2_c$dnds_bi_nonzero)), "gamma")
fit_n_RT_MB2_dnds  <- fitdist(as.numeric(na.omit(dat1_RT_MB2_c$dnds_bi_nonzero)), "norm")
fit_ln_RT_MB2_dnds  <- fitdist(as.numeric(na.omit(dat1_RT_MB2_c$dnds_bi_nonzero)), "lnorm")
fit_g_RT_MB2_dnds  <- fitdist(as.numeric(na.omit(dat1_RT_MB2_c$dnds_bi_nonzero)), "gamma")
fit_n_LG_MB2_dnds  <- fitdist(as.numeric(na.omit(dat1_LG_MB2_c$dnds_bi_nonzero)), "norm")
fit_ln_LG_MB2_dnds  <- fitdist(as.numeric(na.omit(dat1_LG_MB2_c$dnds_bi_nonzero)), "lnorm")
fit_g_LG_MB2_dnds  <- fitdist(as.numeric(na.omit(dat1_LG_MB2_c$dnds_bi_nonzero)), "gamma")

png(filename = "dnds_dist_plots_MB2.png", width = 9, height = 14, units = "in", bg = "white", res = 300)
par(mfrow=c(3,2))
denscomp(list(fit_n_WB_MB2_dnds, fit_ln_WB_MB2_dnds, fit_g_WB_MB2_dnds), legendtext = c("norm", "lnorm", "gamma"))
qqcomp  (list(fit_n_WB_MB2_dnds, fit_ln_WB_MB2_dnds, fit_g_WB_MB2_dnds), legendtext = c("norm", "lnorm", "gamma"))
denscomp(list(fit_n_RT_MB2_dnds, fit_ln_RT_MB2_dnds, fit_g_RT_MB2_dnds), legendtext = c("norm", "lnorm", "gamma"))
qqcomp  (list(fit_n_RT_MB2_dnds, fit_ln_RT_MB2_dnds, fit_g_RT_MB2_dnds), legendtext = c("norm", "lnorm", "gamma"))
denscomp(list(fit_n_LG_MB2_dnds, fit_ln_LG_MB2_dnds, fit_g_LG_MB2_dnds), legendtext = c("norm", "lnorm", "gamma"))
qqcomp  (list(fit_n_LG_MB2_dnds, fit_ln_LG_MB2_dnds, fit_g_LG_MB2_dnds), legendtext = c("norm", "lnorm", "gamma"))
dev.off()
getwd() ## where has my plot gone....


### lood like I want a gamma for dnds 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### binomial

dat1_WB_MB2_c$sp_type <- relevel(dat1_WB_MB2_c$sp_type , "SF")
dat1_WB_MB2_c$sexbias2 <- relevel(dat1_WB_MB2_c$sexbias2 , "Unbiased")
dat1_RT_MB2_c$sp_type <- relevel(dat1_RT_MB2_c$sp_type , "SF")
dat1_RT_MB2_c$sexbias2 <- relevel(dat1_RT_MB2_c$sexbias2 , "Unbiased")
dat1_LG_MB2_c$sp_type <- relevel(dat1_LG_MB2_c$sp_type , "SF")
dat1_LG_MB2_c$sexbias2 <- relevel(dat1_LG_MB2_c$sexbias2 , "Unbiased")

str(dat1_WB_MB2_c)
### WB
## dnds

WB_MB2_dnds_m3  = glmer(dnds_bi ~ sp_type +  sexbias2 +  (1|OG_name), data = dat1_WB_MB2_c, family = "binomial") 
WB_MB2_dnds_m4  = glmer(dnds_bi ~ sp_type *  sexbias2 +  (1|OG_name), data = dat1_WB_MB2_c, family = "binomial") 
WB_MB2_dnds_m4a = glm  (dnds_bi ~ sp_type *  sexbias2, data = dat1_WB_MB2_c, family = "binomial") 

drop1(WB_MB2_dnds_m4,test="Chisq")
drop1(WB_MB2_dnds_m3,.~.,test="Chisq")

WB_MB2_dnds_m4_m4a_LRT = 2 * (as.numeric(logLik(WB_MB2_dnds_m4)) - as.numeric(logLik(WB_MB2_dnds_m4a)))
WB_MB2_dnds_m4_m4a_LRT_p <- pchisq(WB_MB2_dnds_m4_m4a_LRT, df=1, lower.tail=F) 

## writeLines(capture.output(summary(WB_MB2_dnds_m4)), "WB_MB2_dnds_m4_summary.txt")
writeLines(capture.output(drop1(WB_MB2_dnds_m3,.~.,test="Chisq")), "WB_MB2_dnds_m3_lrt_fixed.txt")
writeLines(capture.output(drop1(WB_MB2_dnds_m4,.~.,test="Chisq")), "WB_MB2_dnds_m4_lrt_fixed.txt")
writeLines(capture.output(c("lrt_OG_name = ", WB_MB2_dnds_m4_m4a_LRT,"p_lrt_OG_name", WB_MB2_dnds_m4_m4a_LRT_p )), "WB_MB2_dnds_m4_m4a_rand.txt")



### RT
## dnds

RT_MB2_dnds_m3  = glmer(dnds_bi ~ sp_type +  sexbias2 +  (1|OG_name), data = dat1_RT_MB2_c, family = "binomial") 
RT_MB2_dnds_m4  = glmer(dnds_bi ~ sp_type *  sexbias2 +  (1|OG_name), data = dat1_RT_MB2_c, family = "binomial") 
RT_MB2_dnds_m4a = glm  (dnds_bi ~ sp_type *  sexbias2, data = dat1_RT_MB2_c, family = "binomial") 

drop1(RT_MB2_dnds_m4,test="Chisq")
drop1(RT_MB2_dnds_m3,.~.,test="Chisq")

RT_MB2_dnds_m4_m4a_LRT = 2 * (as.numeric(logLik(RT_MB2_dnds_m4)) - as.numeric(logLik(RT_MB2_dnds_m4a)))
RT_MB2_dnds_m4_m4a_LRT_p <- pchisq(RT_MB2_dnds_m4_m4a_LRT, df=1, lower.tail=F) 

writeLines(capture.output(drop1(RT_MB2_dnds_m3,.~.,test="Chisq")), "RT_MB2_dnds_m3_lrt_fixed.txt")
writeLines(capture.output(drop1(RT_MB2_dnds_m4,.~.,test="Chisq")), "RT_MB2_dnds_m4_lrt_fixed.txt")
writeLines(capture.output(c("lrt_OG_name = ", RT_MB2_dnds_m4_m4a_LRT,"p_lrt_OG_name", RT_MB2_dnds_m4_m4a_LRT_p )), "RT_MB2_dnds_m4_m4a_rand.txt")


### LG
## dnds
LG_MB2_dnds_m3  = glmer(dnds_bi ~ sp_type +  sexbias2 +  (1|OG_name), data = dat1_LG_MB2_c, family = "binomial") 
LG_MB2_dnds_m4  = glmer(dnds_bi ~ sp_type *  sexbias2 +  (1|OG_name), data = dat1_LG_MB2_c, family = "binomial") 
LG_MB2_dnds_m4a = glm  (dnds_bi ~ sp_type *  sexbias2, data = dat1_LG_MB2_c, family = "binomial") 

drop1(LG_MB2_dnds_m4,test="Chisq")
drop1(LG_MB2_dnds_m3,.~.,test="Chisq")

LG_MB2_dnds_m4_m4a_LRT = 2 * (as.numeric(logLik(LG_MB2_dnds_m4)) - as.numeric(logLik(LG_MB2_dnds_m4a)))
LG_MB2_dnds_m4_m4a_LRT_p <- pchisq(LG_MB2_dnds_m4_m4a_LRT, df=1, lower.tail=F) 

writeLines(capture.output(drop1(LG_MB2_dnds_m3,.~.,test="Chisq")), "LG_MB2_dnds_m3_lrt_fixed.txt")
writeLines(capture.output(drop1(LG_MB2_dnds_m4,.~.,test="Chisq")), "LG_MB2_dnds_m4_lrt_fixed.txt")
writeLines(capture.output(c("lrt_OG_name = ", LG_MB2_dnds_m4_m4a_LRT,"p_lrt_OG_name", LG_MB2_dnds_m4_m4a_LRT_p )), "LG_MB2_dnds_m4_m4a_rand.txt")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### non-zero

### WB
## dnds

WB_MB2_dnds_m3_nz   = glmer(dnds_bi_nonzero ~ sp_type +  sexbias2 +  (1|OG_name), data = dat1_WB_MB2_c, family = "Gamma") 
WB_MB2_dnds_m4_nz   = glmer(dnds_bi_nonzero ~ sp_type *  sexbias2 +  (1|OG_name), data = dat1_WB_MB2_c, family = "Gamma") 
WB_MB2_dnds_m4a_nz  = glm  (dnds_bi_nonzero ~ sp_type *  sexbias2, data = dat1_WB_MB2_c, family = "Gamma") 

drop1(WB_MB2_dnds_m4_nz ,test="Chisq")
drop1(WB_MB2_dnds_m3_nz ,.~.,test="Chisq")

WB_MB2_dnds_m4_m4a_LRT_nz  = 2 * (as.numeric(logLik(WB_MB2_dnds_m4_nz )) - as.numeric(logLik(WB_MB2_dnds_m4a_nz )))
WB_MB2_dnds_m4_m4a_LRT_p_nz  <- pchisq(WB_MB2_dnds_m4_m4a_LRT_nz , df=1, lower.tail=F) 

writeLines(capture.output(drop1(WB_MB2_dnds_m3_nz ,.~.,test="Chisq")), "WB_MB2_dnds_m3_lrt_fixed_nz .txt")
writeLines(capture.output(drop1(WB_MB2_dnds_m4_nz ,.~.,test="Chisq")), "WB_MB2_dnds_m4_lrt_fixed_nz .txt")
writeLines(capture.output(c("lrt_OG_name = ", WB_MB2_dnds_m4_m4a_LRT_nz ,"p_lrt_OG_name", WB_MB2_dnds_m4_m4a_LRT_p_nz  )), "WB_MB2_dnds_m4_m4a_rand_nz.txt")

### RT
## dnds

RT_MB2_dnds_m3_nz   = glmer(dnds_bi_nonzero ~ sp_type +  sexbias2 +  (1|OG_name), data = dat1_RT_MB2_c, family = "Gamma") 
RT_MB2_dnds_m4_nz   = glmer(dnds_bi_nonzero ~ sp_type *  sexbias2 +  (1|OG_name), data = dat1_RT_MB2_c, family = "Gamma") 
RT_MB2_dnds_m4a_nz  = glm  (dnds_bi_nonzero ~ sp_type *  sexbias2, data = dat1_RT_MB2_c, family = "Gamma") 

drop1(RT_MB2_dnds_m4_nz ,test="Chisq")
drop1(RT_MB2_dnds_m3_nz ,.~.,test="Chisq")

RT_MB2_dnds_m4_m4a_LRT_nz  = 2 * (as.numeric(logLik(RT_MB2_dnds_m4_nz )) - as.numeric(logLik(RT_MB2_dnds_m4a_nz )))
RT_MB2_dnds_m4_m4a_LRT_p_nz  <- pchisq(RT_MB2_dnds_m4_m4a_LRT_nz , df=1, lower.tail=F) 

writeLines(capture.output(drop1(RT_MB2_dnds_m3_nz ,.~.,test="Chisq")), "RT_MB2_dnds_m3_lrt_fixed_nz .txt")
writeLines(capture.output(drop1(RT_MB2_dnds_m4_nz ,.~.,test="Chisq")), "RT_MB2_dnds_m4_lrt_fixed_nz .txt")
writeLines(capture.output(c("lrt_OG_name = ", RT_MB2_dnds_m4_m4a_LRT_nz ,"p_lrt_OG_name", RT_MB2_dnds_m4_m4a_LRT_p_nz  )), "RT_MB2_dnds_m4_m4a_rand_nz.txt")

### LG
## dnds

LG_MB2_dnds_m3_nz   = glmer(dnds_bi_nonzero ~ sp_type +  sexbias2 +  (1|OG_name), data = dat1_LG_MB2_c, family = "Gamma") 
LG_MB2_dnds_m4_nz   = glmer(dnds_bi_nonzero ~ sp_type *  sexbias2 +  (1|OG_name), data = dat1_LG_MB2_c, family = "Gamma") 
LG_MB2_dnds_m4a_nz  = glm  (dnds_bi_nonzero ~ sp_type *  sexbias2, data = dat1_LG_MB2_c, family = "Gamma") 

drop1(LG_MB2_dnds_m4_nz ,test="Chisq")
drop1(LG_MB2_dnds_m3_nz ,.~.,test="Chisq")

LG_MB2_dnds_m4_m4a_LRT_nz  = 2 * (as.numeric(logLik(LG_MB2_dnds_m4_nz )) - as.numeric(logLik(LG_MB2_dnds_m4a_nz )))
LG_MB2_dnds_m4_m4a_LRT_p_nz  <- pchisq(LG_MB2_dnds_m4_m4a_LRT_nz , df=1, lower.tail=F) 

writeLines(capture.output(drop1(LG_MB2_dnds_m3_nz ,.~.,test="Chisq")), "LG_MB2_dnds_m3_lrt_fixed_nz .txt")
writeLines(capture.output(drop1(LG_MB2_dnds_m4_nz ,.~.,test="Chisq")), "LG_MB2_dnds_m4_lrt_fixed_nz .txt")
writeLines(capture.output(c("lrt_OG_name = ", LG_MB2_dnds_m4_m4a_LRT_nz ,"p_lrt_OG_name", LG_MB2_dnds_m4_m4a_LRT_p_nz  )), "LG_MB2_dnds_m4_m4a_rand_nz.txt")



########################################################################################################################
#### using female-bias vs unbias

fit_n_WB_FB2_dnds  <- fitdist(as.numeric(na.omit(dat1_WB_FB2_c$dnds_bi_nonzero)), "norm")
fit_ln_WB_FB2_dnds  <- fitdist(as.numeric(na.omit(dat1_WB_FB2_c$dnds_bi_nonzero)), "lnorm")
fit_g_WB_FB2_dnds  <- fitdist(as.numeric(na.omit(dat1_WB_FB2_c$dnds_bi_nonzero)), "gamma")
fit_n_RT_FB2_dnds  <- fitdist(as.numeric(na.omit(dat1_RT_FB2_c$dnds_bi_nonzero)), "norm")
fit_ln_RT_FB2_dnds  <- fitdist(as.numeric(na.omit(dat1_RT_FB2_c$dnds_bi_nonzero)), "lnorm")
fit_g_RT_FB2_dnds  <- fitdist(as.numeric(na.omit(dat1_RT_FB2_c$dnds_bi_nonzero)), "gamma")
fit_n_LG_FB2_dnds  <- fitdist(as.numeric(na.omit(dat1_LG_FB2_c$dnds_bi_nonzero)), "norm")
fit_ln_LG_FB2_dnds  <- fitdist(as.numeric(na.omit(dat1_LG_FB2_c$dnds_bi_nonzero)), "lnorm")
fit_g_LG_FB2_dnds  <- fitdist(as.numeric(na.omit(dat1_LG_FB2_c$dnds_bi_nonzero)), "gamma")

png(filename = "dnds_dist_plots_FB2.png", width = 9, height = 14, units = "in", bg = "white", res = 300)
par(mfrow=c(3,2))
denscomp(list(fit_n_WB_FB2_dnds, fit_ln_WB_FB2_dnds, fit_g_WB_FB2_dnds), legendtext = c("norm", "lnorm", "gamma"))
qqcomp  (list(fit_n_WB_FB2_dnds, fit_ln_WB_FB2_dnds, fit_g_WB_FB2_dnds), legendtext = c("norm", "lnorm", "gamma"))
denscomp(list(fit_n_RT_FB2_dnds, fit_ln_RT_FB2_dnds, fit_g_RT_FB2_dnds), legendtext = c("norm", "lnorm", "gamma"))
qqcomp  (list(fit_n_RT_FB2_dnds, fit_ln_RT_FB2_dnds, fit_g_RT_FB2_dnds), legendtext = c("norm", "lnorm", "gamma"))
denscomp(list(fit_n_LG_FB2_dnds, fit_ln_LG_FB2_dnds, fit_g_LG_FB2_dnds), legendtext = c("norm", "lnorm", "gamma"))
qqcomp  (list(fit_n_LG_FB2_dnds, fit_ln_LG_FB2_dnds, fit_g_LG_FB2_dnds), legendtext = c("norm", "lnorm", "gamma"))
dev.off()
getwd() ## where has my plot gone....

### looks like I want a gamma for dnds 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### binomial

dat1_WB_FB2_c$sp_type <- relevel(dat1_WB_FB2_c$sp_type , "SF")
dat1_WB_FB2_c$sexbias2 <- relevel(dat1_WB_FB2_c$sexbias2 , "Unbiased")
dat1_RT_FB2_c$sp_type <- relevel(dat1_RT_FB2_c$sp_type , "SF")
dat1_RT_FB2_c$sexbias2 <- relevel(dat1_RT_FB2_c$sexbias2 , "Unbiased")
dat1_LG_FB2_c$sp_type <- relevel(dat1_LG_FB2_c$sp_type , "SF")
dat1_LG_FB2_c$sexbias2 <- relevel(dat1_LG_FB2_c$sexbias2 , "Unbiased")

str(dat1_WB_FB2_c)
### WB
## dnds

WB_FB2_dnds_m3  = glmer(dnds_bi ~ sp_type +  sexbias2 +  (1|OG_name), data = dat1_WB_FB2_c, family = "binomial") 
WB_FB2_dnds_m4  = glmer(dnds_bi ~ sp_type *  sexbias2 +  (1|OG_name), data = dat1_WB_FB2_c, family = "binomial") 
WB_FB2_dnds_m4a = glm  (dnds_bi ~ sp_type *  sexbias2, data = dat1_WB_FB2_c, family = "binomial") 

drop1(WB_FB2_dnds_m4,test="Chisq")
drop1(WB_FB2_dnds_m3,.~.,test="Chisq")

WB_FB2_dnds_m4_m4a_LRT = 2 * (as.numeric(logLik(WB_FB2_dnds_m4)) - as.numeric(logLik(WB_FB2_dnds_m4a)))
WB_FB2_dnds_m4_m4a_LRT_p <- pchisq(WB_FB2_dnds_m4_m4a_LRT, df=1, lower.tail=F) 

## writeLines(capture.output(summary(WB_FB2_dnds_m4)), "WB_FB2_dnds_m4_summary.txt")
writeLines(capture.output(drop1(WB_FB2_dnds_m3,.~.,test="Chisq")), "WB_FB2_dnds_m3_lrt_fixed.txt")
writeLines(capture.output(drop1(WB_FB2_dnds_m4,.~.,test="Chisq")), "WB_FB2_dnds_m4_lrt_fixed.txt")
writeLines(capture.output(c("lrt_OG_name = ", WB_FB2_dnds_m4_m4a_LRT,"p_lrt_OG_name", WB_FB2_dnds_m4_m4a_LRT_p )), "WB_FB2_dnds_m4_m4a_rand.txt")



### RT
## dnds

RT_FB2_dnds_m3  = glmer(dnds_bi ~ sp_type +  sexbias2 +  (1|OG_name), data = dat1_RT_FB2_c, family = "binomial") 
RT_FB2_dnds_m4  = glmer(dnds_bi ~ sp_type *  sexbias2 +  (1|OG_name), data = dat1_RT_FB2_c, family = "binomial") 
RT_FB2_dnds_m4a = glm  (dnds_bi ~ sp_type *  sexbias2, data = dat1_RT_FB2_c, family = "binomial") 

drop1(RT_FB2_dnds_m4,test="Chisq")
drop1(RT_FB2_dnds_m3,.~.,test="Chisq")

RT_FB2_dnds_m4_m4a_LRT = 2 * (as.numeric(logLik(RT_FB2_dnds_m4)) - as.numeric(logLik(RT_FB2_dnds_m4a)))
RT_FB2_dnds_m4_m4a_LRT_p <- pchisq(RT_FB2_dnds_m4_m4a_LRT, df=1, lower.tail=F) 

writeLines(capture.output(drop1(RT_FB2_dnds_m3,.~.,test="Chisq")), "RT_FB2_dnds_m3_lrt_fixed.txt")
writeLines(capture.output(drop1(RT_FB2_dnds_m4,.~.,test="Chisq")), "RT_FB2_dnds_m4_lrt_fixed.txt")
writeLines(capture.output(c("lrt_OG_name = ", RT_FB2_dnds_m4_m4a_LRT,"p_lrt_OG_name", RT_FB2_dnds_m4_m4a_LRT_p )), "RT_FB2_dnds_m4_m4a_rand.txt")


### LG
## dnds
LG_FB2_dnds_m3  = glmer(dnds_bi ~ sp_type +  sexbias2 +  (1|OG_name), data = dat1_LG_FB2_c, family = "binomial") 
LG_FB2_dnds_m4  = glmer(dnds_bi ~ sp_type *  sexbias2 +  (1|OG_name), data = dat1_LG_FB2_c, family = "binomial") 
LG_FB2_dnds_m4a = glm  (dnds_bi ~ sp_type *  sexbias2, data = dat1_LG_FB2_c, family = "binomial") 

drop1(LG_FB2_dnds_m4,test="Chisq")
drop1(LG_FB2_dnds_m3,.~.,test="Chisq")

LG_FB2_dnds_m4_m4a_LRT = 2 * (as.numeric(logLik(LG_FB2_dnds_m4)) - as.numeric(logLik(LG_FB2_dnds_m4a)))
LG_FB2_dnds_m4_m4a_LRT_p <- pchisq(LG_FB2_dnds_m4_m4a_LRT, df=1, lower.tail=F) 

writeLines(capture.output(drop1(LG_FB2_dnds_m3,.~.,test="Chisq")), "LG_FB2_dnds_m3_lrt_fixed.txt")
writeLines(capture.output(drop1(LG_FB2_dnds_m4,.~.,test="Chisq")), "LG_FB2_dnds_m4_lrt_fixed.txt")
writeLines(capture.output(c("lrt_OG_name = ", LG_FB2_dnds_m4_m4a_LRT,"p_lrt_OG_name", LG_FB2_dnds_m4_m4a_LRT_p )), "LG_FB2_dnds_m4_m4a_rand.txt")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### non-zero

### WB
## dnds

WB_FB2_dnds_m3_nz   = glmer(dnds_bi_nonzero ~ sp_type +  sexbias2 +  (1|OG_name), data = dat1_WB_FB2_c, family = "Gamma") 
WB_FB2_dnds_m4_nz   = glmer(dnds_bi_nonzero ~ sp_type *  sexbias2 +  (1|OG_name), data = dat1_WB_FB2_c, family = "Gamma") 
WB_FB2_dnds_m4a_nz  = glm  (dnds_bi_nonzero ~ sp_type *  sexbias2, data = dat1_WB_FB2_c, family = "Gamma") 

drop1(WB_FB2_dnds_m4_nz ,test="Chisq")
drop1(WB_FB2_dnds_m3_nz ,.~.,test="Chisq")

WB_FB2_dnds_m4_m4a_LRT_nz  = 2 * (as.numeric(logLik(WB_FB2_dnds_m4_nz )) - as.numeric(logLik(WB_FB2_dnds_m4a_nz )))
WB_FB2_dnds_m4_m4a_LRT_p_nz  <- pchisq(WB_FB2_dnds_m4_m4a_LRT_nz , df=1, lower.tail=F) 

writeLines(capture.output(drop1(WB_FB2_dnds_m3_nz ,.~.,test="Chisq")), "WB_FB2_dnds_m3_lrt_fixed_nz .txt")
writeLines(capture.output(drop1(WB_FB2_dnds_m4_nz ,.~.,test="Chisq")), "WB_FB2_dnds_m4_lrt_fixed_nz .txt")
writeLines(capture.output(c("lrt_OG_name = ", WB_FB2_dnds_m4_m4a_LRT_nz ,"p_lrt_OG_name", WB_FB2_dnds_m4_m4a_LRT_p_nz  )), "WB_FB2_dnds_m4_m4a_rand_nz.txt")

### RT
## dnds

RT_FB2_dnds_m3_nz   = glmer(dnds_bi_nonzero ~ sp_type +  sexbias2 +  (1|OG_name), data = dat1_RT_FB2_c, family = "Gamma") 
RT_FB2_dnds_m4_nz   = glmer(dnds_bi_nonzero ~ sp_type *  sexbias2 +  (1|OG_name), data = dat1_RT_FB2_c, family = "Gamma") 
RT_FB2_dnds_m4a_nz  = glm  (dnds_bi_nonzero ~ sp_type *  sexbias2, data = dat1_RT_FB2_c, family = "Gamma") 

drop1(RT_FB2_dnds_m4_nz ,test="Chisq")
drop1(RT_FB2_dnds_m3_nz ,.~.,test="Chisq")

RT_FB2_dnds_m4_m4a_LRT_nz  = 2 * (as.numeric(logLik(RT_FB2_dnds_m4_nz )) - as.numeric(logLik(RT_FB2_dnds_m4a_nz )))
RT_FB2_dnds_m4_m4a_LRT_p_nz  <- pchisq(RT_FB2_dnds_m4_m4a_LRT_nz , df=1, lower.tail=F) 

writeLines(capture.output(drop1(RT_FB2_dnds_m3_nz ,.~.,test="Chisq")), "RT_FB2_dnds_m3_lrt_fixed_nz .txt")
writeLines(capture.output(drop1(RT_FB2_dnds_m4_nz ,.~.,test="Chisq")), "RT_FB2_dnds_m4_lrt_fixed_nz .txt")
writeLines(capture.output(c("lrt_OG_name = ", RT_FB2_dnds_m4_m4a_LRT_nz ,"p_lrt_OG_name", RT_FB2_dnds_m4_m4a_LRT_p_nz  )), "RT_FB2_dnds_m4_m4a_rand_nz.txt")

### LG
## dnds

LG_FB2_dnds_m3_nz   = glmer(dnds_bi_nonzero ~ sp_type +  sexbias2 +  (1|OG_name), data = dat1_LG_FB2_c, family = "Gamma") 
LG_FB2_dnds_m4_nz   = glmer(dnds_bi_nonzero ~ sp_type *  sexbias2 +  (1|OG_name), data = dat1_LG_FB2_c, family = "Gamma") 
LG_FB2_dnds_m4a_nz  = glm  (dnds_bi_nonzero ~ sp_type *  sexbias2, data = dat1_LG_FB2_c, family = "Gamma") 

drop1(LG_FB2_dnds_m4_nz ,test="Chisq")
drop1(LG_FB2_dnds_m3_nz ,.~.,test="Chisq")

LG_FB2_dnds_m4_m4a_LRT_nz  = 2 * (as.numeric(logLik(LG_FB2_dnds_m4_nz )) - as.numeric(logLik(LG_FB2_dnds_m4a_nz )))
LG_FB2_dnds_m4_m4a_LRT_p_nz  <- pchisq(LG_FB2_dnds_m4_m4a_LRT_nz , df=1, lower.tail=F) 

writeLines(capture.output(drop1(LG_FB2_dnds_m3_nz ,.~.,test="Chisq")), "LG_FB2_dnds_m3_lrt_fixed_nz .txt")
writeLines(capture.output(drop1(LG_FB2_dnds_m4_nz ,.~.,test="Chisq")), "LG_FB2_dnds_m4_lrt_fixed_nz .txt")
writeLines(capture.output(c("lrt_OG_name = ", LG_FB2_dnds_m4_m4a_LRT_nz ,"p_lrt_OG_name", LG_FB2_dnds_m4_m4a_LRT_p_nz  )), "LG_FB2_dnds_m4_m4a_rand_nz.txt")



########################################################################################################################################################################
####### output session info

writeLines(capture.output(sessionInfo()), "dnds_SBGE.R_sessionInfo.txt")
