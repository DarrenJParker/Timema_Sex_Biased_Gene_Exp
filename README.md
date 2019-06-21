# Timema_Sex_Biased_Gene_Exp


This is the repository for the collected scripts used in the study:

>Parker, D. J., Bast, J., Jalvingh, K., Dumas, Z., Robinson-Rechavi, M., Schwander, T. (2019). Sex-biased gene expression is repeatedly masculinized in asexual females. bioRxiv: doi: https://doi.org/10.1101/553172

Currently under review.

# Components

## DATA

* Contains read counts, GO terms, and dN/dS estimates for input to the scripts below. 

## SCRIPTS

### Differential expression analyses

**_Main analyses using orthologs between sexual and asexual sister species_**

* Get Sex-biased genes:
    * Sex_bias_edgeR.R
* Get Sex-limited genes:
    * Sex_bias_edgeR_Sex_limited_expression.R

* bring expression data together:

```
mkdir Output/DE_joined	
for sp in Tbi Tce Tcm Tpa Tps; do
	for tiss in WB RT LG; do
	
	python3 sex_bias_asex_edgeR_tidier.py \
		-i Output/DE -s \
		"TTT_lrt_"$sp"_sex_bias_"$tiss".csv" -a sex_asex -o "Output/DE_joined/"$sp"_"$tiss"_RBBH_disp_allsepar" \
		 -x $sp"_"$tiss"_F_FPKM_RBBH_1.csv"
	done
done
```
* plot
    * Sex_bias_plotsetc.R

* Test for a relationship between the shift in sex-biased gene expression in asexual females and sexual-asexual divergence time
    * divergence_and_change_in_SB.R

**_Analyses using orthologs between all sexual and asexual species_**

* Get Sex-biased genes
    * Sex_bias_edgeR_10sp_orths.R

* bring expression data together

```
mkdir Output/DE_joined_10sp
python3 sex_bias_asex_edgeR_tidier.py -i Output/DE_10sp/ -s \
    TTT_lrt_Tbi_sex_bias_WB.csv,TTT_lrt_Tce_sex_bias_WB.csv,TTT_lrt_Tcm_sex_bias_WB.csv,TTT_lrt_Tpa_sex_bias_WB.csv,TTT_lrt_Tps_sex_bias_WB.csv \
	-a sex_asex -o Output/DE_joined_10sp/ALL_WB_10sp_disp_allsepar_wCPM -c cpm_df_WB_10sp_ALL.csv

python3 sex_bias_asex_edgeR_tidier.py -i Output/DE_10sp/ -s \
    TTT_lrt_Tbi_sex_bias_RT.csv,TTT_lrt_Tce_sex_bias_RT.csv,TTT_lrt_Tcm_sex_bias_RT.csv,TTT_lrt_Tpa_sex_bias_RT.csv,TTT_lrt_Tps_sex_bias_RT.csv \
	-a sex_asex -o Output/DE_joined_10sp/ALL_RT_10sp_disp_allsepar_wCPM -c cpm_df_RT_10sp_ALL.csv

python3 sex_bias_asex_edgeR_tidier.py -i Output/DE_10sp/ -s \
    TTT_lrt_Tbi_sex_bias_LG.csv,TTT_lrt_Tce_sex_bias_LG.csv,TTT_lrt_Tcm_sex_bias_LG.csv,TTT_lrt_Tpa_sex_bias_LG.csv,TTT_lrt_Tps_sex_bias_LG.csv \
	-a sex_asex -o Output/DE_joined_10sp/ALL_LG_10sp_disp_allsepar_wCPM -c cpm_df_LG_10sp_ALL.csv
```

* plot:
    * Sex_bias_plotsetc_10sp_orths.R	
    * N_SB_genes_10sp.R

*Analyses using virgin females*

* get SB genes, cluster virgin and mated samples
    * sex_bias_edgeR_withVIfemales.R


* bring expression data together:

```
mkdir Output/DE_joined_Virgin	
for sp in Tbi Tce Tcm Tpa Tps; do
	for tiss in WB; do
	
	python3 sex_bias_asex_edgeR_tidier.py \
		-i Output/DE_Virgin -s \
		"TTT_lrt_"$sp"_sex_bias_"$tiss".csv" -a sex_asex -o "Output/DE_joined_Virgin/"$sp"_"$tiss"_RBBH_disp_allsepar"
	done
done
```

* plot
    * Sex_bias_plotsetc_withVIfemales.R


**_Analyses without filtering genes with low expression in asexual females_**

* Get Sex-biased genes:
    * Sex_bias_edgeR_nocpmfilteronasex.R

* bring expression data together:

```
mkdir Output/DE_joined_nocpmfilteronasex
for sp in Tbi Tce Tcm Tpa Tps; do
	for tiss in WB RT LG; do
	
	python3 sex_bias_asex_edgeR_tidier.py \
		-i Output/DE_nocpmfiltonasex -s \
		"TTT_lrt_"$sp"_sex_bias_"$tiss".csv" -a sex_asex -o "Output/DE_joined_nocpmfilteronasex/"$sp"_"$tiss"_RBBH_disp_allsepar" 
	done
done
```

* plot
    * Sex_bias_plotsetc_nocpmfilteronasex.R

**_Analyses using full transcriptome references_**

* Using the sexual references
    * Sex_bias_edgeR_sexual_ref.R

* bring expression data together:

```
mkdir Output/DE_joined_sexual_ref	
for sp in Tbi Tce Tcm Tpa Tps; do
	for tiss in WB RT LG; do
	
	python3 sex_bias_asex_edgeR_tidier.py \
		-i Output/DE_sexual_ref -s \
		"TTT_lrt_"$sp"_sex_bias_"$tiss".csv" -a sex_asex -o "Output/DE_joined_sexual_ref/"$sp"_"$tiss"_longest_iso_disp_allsepar" \
		 -x $sp"_"$tiss"_N_FPKM_longest_iso_1.csv"
	done
done
```

* plot
    * Sex_bias_plotsetc_sexual_ref.R


* Using the asexual references
    * Sex_bias_edgeR_asexual_ref.R

* bring expression data together:

```
mkdir Output/DE_joined_asexual_ref	
for sp in Tbi Tce Tcm Tpa Tps; do
	for tiss in WB RT LG; do
	
	python3 sex_bias_asex_edgeR_tidier.py \
		-i Output/DE_asexual_ref -s \
		"TTT_lrt_"$sp"_sex_bias_"$tiss".csv" -a sex_asex -o "Output/DE_joined_asexual_ref/"$sp"_"$tiss"_longest_iso_disp_allsepar" \
		 -x $sp"_"$tiss"_N_FPKM_longest_iso_1.csv"
	done
done
```

* plot
    * Sex_bias_plotsetc_asexual_ref.R


### GO term analyses

* TopGO.R | for the main GO term analyses
* TopGO_10sp.R | for the GO term analyses using only the 10 species orthologs

### dN/dS analyses

* dnds_SBGE.R

### Additional scripts

* **B2G_to_topGO.py** | Script for converting Blast2GO output into a format usable by topGO
* **Get_GO_term_parent_and_child_overlap_adjuster.py** | Script for dealing with topographically close GO terms 
* **Get_GO_term_parent_and_child_overlap_adjuster_test_data_expl.R** | Script explaining the logic of Get_GO_term_parent_and_child_overlap_adjuster.py
* **super_exact_test_multitest_corrector.py** | Script to correct the p-values of SuperExactTest for multiple tests
* **super_exact_test_table_parser.py** | Script to tidy up the output of super_exact_test_multitest_corrector.py


# Infomation on running scripts

## General 

* All scripts should be run from the directory they are in. Output directories will be created to store Output as the code is run. 
* All python scripts were made using python 3.5. All contain help information which can be displayed by specifying no command line arguments.


# Abbreviations

## Species names:

Species name | Abbreviation | Reproductive mode 
--- | --- | --- 
*Timema bartmani* | Tbi | sexual 
*Timema tahoe* | Tte | asexual
*Timema cristinae* | Tce | sexual 
*Timema monikensis* | Tms | asexual
*Timema poppensis* | Tps | sexual 
*Timema douglasi* | Tdi | asexual
*Timema californicum* | Tcm | sexual 
*Timema shepardi* | Tsi | asexual
*Timema podura* | Tpa | sexual 
*Timema genevievae* | Tge | asexual

## Tissues:

Tissue | Abbreviation 
--- | --- 
Whole-body| WB
Reproductive tract | RT
Legs | LG




