##### sex_bias_asex_edgeR_tidier.py


import sys
import os
import getopt
import decimal

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:o:s:a:f:c:x:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

in_dir_name			 = "NOTHINGSET"
SB_file_names        = "NOTHINGSET"    ### sexual female vs sexual male
Asex_file_name_sub   = "NOTHINGSET"    ### asexual female vs sexual female file name sub for the "sex_bias" part of the sex_bias file names
output_filename_base = "NOTHINGSET"
FDR_cutoff           = "NOTHINGSET"
cpm_filename         = "NOTHINGSET"
FPKM_filename         = "NOTHINGSET"


#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** sex_bias_asex_edgeR_tidier_5.py ****\n")
		print("\ntakes 2 egdeR outputfiles from Timema samples and an expression file and joins them all together: \
			1: where sex bias has been compared per contig in sexuals, 2: of a different comparsion (e.g. sex/asex) 3 (opt) cpm file 4 (opt) FPKM file")
	
		print("NOTE it outputs genes alphabetical order. If gene is missing from a sample, values are set to NA")

		print("\n**** USAGE **** \n")
		print("\nSee README file \n")		
	

		sys.exit(2)
		
	elif opt in ('-i'):
		in_dir_name = arg
	elif opt in ('-o'):
		output_filename_base = arg
	elif opt in ('-s'):
		SB_file_names = arg
	elif opt in ('-a'):
		Asex_file_name_sub = arg
	elif opt in ('-f'):
		FDR_cutoff = arg
	elif opt in ('-c'):
		cpm_filename = arg
	elif opt in ('-x'):
		FPKM_filename = arg		
	else:
		print("i dont know")
		sys.exit(2)


if FDR_cutoff == "NOTHINGSET":
	print("\nFDR cutoff for calling sex-bias set to 0.05. Use -f to alter this\n")
	FDR_cutoff = 0.05
else:
	print("\nFDR cutoff set to " + str(FDR_cutoff) +"\n")
	
	
if cpm_filename == "NOTHINGSET":
	print("\nSkipping cpm in the output. Specify -c to add this\n")
else:
	print("Using CPM values from " + cpm_filename)

### GET ALL SB gene names

all_gene_names = set()
sample_name_list = []

sb_filename_list = SB_file_names.split(",")
for f in sb_filename_list:
	SB_file_path = os.path.join(in_dir_name, f)
	SB_file = open(SB_file_path)
	#print(f)
	sample_name = f.split("_")[2] + "_" + f.split("_")[5]
	sample_name = sample_name.split(".")[0]
	sample_name_list.append(sample_name)
	#print(sample_name)
	line_N = 0
	for line in SB_file:
		line_N = line_N + 1
		if line_N != 1:
			line = line.rstrip("\n").split(",")
			gene_name  = line[0]
			all_gene_names.add(gene_name)
	SB_file.close()

print("\nNumber of unique genes in input files: " + str(len(all_gene_names)) + "\n")

### get SB info to dict

SB_dict_2 = {}

sb_filename_list = SB_file_names.split(",")
for f in sb_filename_list:
	SB_file_path = os.path.join(in_dir_name, f)
	SB_file = open(SB_file_path)
	#print(f)
	sample_name = f.split("_")[2] + "_" + f.split("_")[5]
	sample_name = sample_name.split(".")[0]
	#print(sample_name)
	line_N = 0
	for line in SB_file:
		line_N = line_N + 1
		if line_N > 1:
			line = line.rstrip("\n").split(",")
			
			gene_name  = line[0] + "SPLTTTYNNNN____HH_" + sample_name
			log2_FC      = decimal.Decimal(line[1])
			FDR_val      = decimal.Decimal(line[5])
			bias = ""
			
			if FDR_val < FDR_cutoff and log2_FC > 0: #### female biased
				bias = "female_biased"
			elif FDR_val < FDR_cutoff and log2_FC < 0: #### male biased
				bias = "male_biased"
			elif FDR_val < FDR_cutoff and log2_FC == 0: #### just in case
				bias = "ERROR?"
			else:
				bias = "Unbiased"
			
			
			### add to dict
			
			new_rec = str(log2_FC) + "," + str(FDR_val) + "," + bias
			SB_dict_2[gene_name] = new_rec 
			# print(gene_name)
			# print(new_rec)
	SB_file.close()

#####################
## same for asex 

SA_dict_2 = {}

for f in sb_filename_list:
	asex_name = f.replace("sex_bias", Asex_file_name_sub)
	#print(asex_name)
	asex_file_path = os.path.join(in_dir_name, asex_name)
	asex_file = open(asex_file_path)

	sample_name = f.split("_")[2] + "_" + f.split("_")[5]
	sample_name = sample_name.split(".")[0]
	#print(sample_name)
	line_N = 0
	for line in asex_file:
		#print(line)
		line_N = line_N + 1
		if line_N > 1:
			line = line.rstrip("\n").split(",")
			
			gene_name  = line[0] + "SPLTTTYNNNN____HH_" + sample_name
			log2_FC      = decimal.Decimal(line[1])
			FDR_val      = decimal.Decimal(line[5])

			### add to dict

			new_rec = str(log2_FC) + "," + str(FDR_val)
			SA_dict_2[gene_name] = new_rec
	
	asex_file.close()



####################################################################################################################
#####  cpm 

cpm_dict = {}
cpm_header = ""

if cpm_filename != "NOTHINGSET":
	cpm_file_path = os.path.join(in_dir_name, cpm_filename)
	cpm_file = open(cpm_file_path)
	
	line_N = 0 	
	for line in cpm_file:
		line_N = line_N + 1

		line = line.rstrip("\n").replace('"', '').split(",")

		if line_N == 1:
			cpm_header = ""
			for c in range(1,len(line)):
				cpm_header = cpm_header + ",cpm_" + line[c]			
			

		else:
			
			gene_name = line[0]
			
			rest_of_line = ""
			for c in range(1,len(line)):
				rest_of_line = rest_of_line + "," + line[c]
		


			cpm_dict[gene_name] = rest_of_line
		
			# print(line)
			# print(rest_of_line)


####################################################################################################################
#####  FPKM 

FPKM_dict = {}
FPKM_header = ""

if FPKM_filename != "NOTHINGSET":
	FPKM_file_path = os.path.join(in_dir_name, FPKM_filename)
	FPKM_file = open(FPKM_file_path)
	
	line_N = 0 	
	for line in FPKM_file:
		line_N = line_N + 1

		line = line.rstrip("\n").replace('"', '').split(",")

		if line_N == 1:
			FPKM_header = ""
			for c in range(1,len(line)):
				FPKM_header = FPKM_header + ",FPKM_" + line[c]			
			

		else:
			
			gene_name = line[0]
			
			rest_of_line = ""
			for c in range(1,len(line)):
				rest_of_line = rest_of_line + "," + line[c]
		


			FPKM_dict[gene_name] = rest_of_line
		
			# print(line)
			# print(rest_of_line)


####################################################################################################################
##### export

## order genes:

all_gene_names_l = sorted(list(all_gene_names))

output_filename = output_filename_base + "_SB_asex.csv"
output_file = open(output_filename, "w")


### make header

out_header = ""
for s in sample_name_list:
	out_header = out_header + "," + s + "_log2FC_SB" + "," + s + "_FDR_SB" + "," + s + "_sexbias" + "," + s + "_log2FC_SA" + "," + s + "_FDR_SA"
out_header = "genename" + out_header

if cpm_filename != "NOTHINGSET":
	out_header = out_header + cpm_header

if FPKM_filename != "NOTHINGSET":
	out_header = out_header + FPKM_header

output_file.write(out_header + "\n")


for el in all_gene_names_l:
	
	gene_rec = ""
	
	for s in sample_name_list:
		gene_name = el + "SPLTTTYNNNN____HH_" + s
		SB_rec = SB_dict_2.get(gene_name)
		if SB_rec == None:
			SB_rec = "NA,NA,NA"
		#print(SB_rec)
		
		
		SA_rec = SA_dict_2.get(gene_name)
		if SA_rec == None:
			SA_rec = "NA,NA"
		#print(SA_rec)
				
		gene_rec = gene_rec + "," + SB_rec + "," + SA_rec
		
	if cpm_filename != "NOTHINGSET":
		cpm_rec = cpm_dict.get(el)
		gene_rec = gene_rec + cpm_rec

	if FPKM_filename != "NOTHINGSET":
		FPKM_rec = FPKM_dict.get(el)
		gene_rec = gene_rec + FPKM_rec

	output_file.write(el + gene_rec + "\n")


print("Output file: " + output_filename)

print("\n\n\nFinished, Philip Raven\n\n")
	



