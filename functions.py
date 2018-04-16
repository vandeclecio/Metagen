#!/usr/bin/envs python
#-*- encoding: utf-8 -*-
import os, re

def sampleMetaData(sampleFile, coluna):
	LIST = []
	for line in sampleFile:
		line = line.strip().split("\t")
		if bool(not re.search("^#", line[coluna])):
			LIST.append(line[coluna])

	return LIST	


def rename(CODE, DIR_INPUT, DIR_OUT, layout): # verificar saidas de arquivo paired  e single end sem merge 
	print("\n\n---------------------Renomeando arquivos---------------------\n\n")
	for ID in CODE:
		if(layout == 0):
			if(os.path.exists(DIR_OUT + ID + "_00_L001_R1_001.fastq.gz")):
				print("Arquivo existente")
			else:
				os.system("mv " + DIR_INPUT + ID + "_trimmed_no_hits.fastq.gz " + DIR_INPUT + ID + "_00_L001_R1_001.fastq.gz")
		else:
			if(layout == 2):
				if(os.path.exists(DIR_OUT + ID + "_00_L001_R1_001.fastq.gz")):
					print("Arquivo existente")
				else:
					os.system("mv " + DIR_INPUT + ID + "_merged_trimmed_no_hits.fastq.gz " + DIR_INPUT + ID + "_00_L001_R1_001.fastq.gz")
			elif(layout == 1):		
				if(os.path.exists(DIR_OUT + ID + "_00_L001_R1_001.fastq.gz") and os.path.exists(DIR_OUT + ID + "_00_L001_R2_001.fastq.gz")):
					print("Arquivos existentes")
				else:
					os.system("mv " + DIR_INPUT + ID + "_1_val_1.fq.gz " + DIR_INPUT + ID + "_00_L001_R1_001.fastq.gz")
					os.system("mv " + DIR_INPUT + ID + "_2_val_2.fq.gz " + DIR_INPUT + ID + "_00_L001_R2_001.fastq.gz")

def QCfastqDump(CODE, DIR_INPUT, DIR_OUT, layout):
	print("\n\n---------------------FASTQ DUMP---------------------\n\n")
	for ID in CODE:
		if(layout == 0):
			if(os.path.exists(DIR_OUT + ID + ".fastq.gz")):
				print("Arquivo existente")
			else:
				os.system("fastq-dump --gzip " + DIR_INPUT + ID + ".sra -O " + DIR_OUT)
		else:
			if(os.path.exists(DIR_OUT + ID + "_1.fastq.gz") and os.path.exists(DIR_OUT + ID + "_2.fastq.gz")):
				print("Arquivos existentes")
			else:
				os.system("fastq-dump --gzip --split-files " + DIR_INPUT + ID + ".sra -O " + DIR_OUT)


def QCfastqc(CODE, DIR_INPUT, DIR_OUT,layout):
	print("\n\n---------------------FASTQC---------------------\n\n")
	for ID in CODE:
		print("fastqc " + ID)
		if(not layout):
			if(not os.path.exists(DIR_OUT + ID + "_fastqc.html")):
				os.system("fastqc -f fastq " + DIR_INPUT + ID + ".fastq.gz -o " + DIR_OUT)
			else:
				print("Arquivo " + DIR_OUT + ID + " existente.")	
		else:
			if(not os.path.exists(DIR_OUT + ID + "_1_fastqc.html") and not os.path.exists(DIR_OUT + ID + "_2_fastqc.html")):
				os.system("fastqc -f fastq " + DIR_INPUT + ID + "_1.fastq.gz " + DIR_INPUT + ID + "_2.fastq.gz -o " + DIR_OUT)
			else:
				print("Arquivo " + DIR_OUT + ID + " existente.")

def QCSeqPrep(CODE, DIR_INPUT, DIR_OUT, layout):
	print("\n\n---------------------SEQPREP---------------------\n\n")
	for ID in CODE:
		print("SeqPrep " + ID)
		if(layout == 2):
			if(not os.path.exists(DIR_OUT + ID + "_merged.fastq.gz")):
				os.system("SeqPrep -f " + DIR_INPUT + ID + "_1.fastq.gz -r " + DIR_INPUT + ID + "_2.fastq.gz -1 " + DIR_INPUT + ID + "_1_.fastq.gz -2 "\
				 + DIR_INPUT + ID + "_2_.fastq.gz -3 " + DIR_OUT + ID + "_1.discarted.fastq.gz -4 " + DIR_OUT + ID + "_2.discarted.fastq.gz -s "\
				 + DIR_OUT + ID + "_merged.fastq.gz -E " + DIR_OUT + ID + ".aling.fastq.gz")
			else:
				print("Arquivo " + DIR_OUT + ID + " existente.")

			if(not os.path.exists(DIR_OUT + "FASTQC/" + ID + "_merged_fastqc.html")):	
				os.system("fastqc -f fastq " + DIR_INPUT + ID + "_merged.fastq.gz -o " + DIR_OUT + "FASTQC/")
			else:
				print("Arquivo " + DIR_OUT + ID + " existente.")	
		else:
			if(not layout):
				print("As sequencias são single-end")
			else:
			 	print("Não foi realizado merge nas sequencias")

def QCtrimGalore(CODE, DIR_INPUT, DIR_OUT, layout):
	print("\n\n---------------------TRIM_GALORE---------------------\n\n")
	for ID in CODE:
		print("Trim_Galore " + ID)
		if(not layout):
			if(not os.path.exists(DIR_OUT + ID + "_trimmed.fq.gz")):
				os.system("trim_galore --illumina --stringency 5 --fastqc --length 20 " + DIR_INPUT + ID + ".fastq.gz -o " + DIR_OUT)
			else:
				print("Arquivo " + DIR_OUT + ID + " existente.")	
		# Arquivo de entrada são os merged
		else:
			if(layout == 2):
				if(not os.path.exists(DIR_OUT + ID + "_merged_trimmed.fq.gz")):
					os.system("trim_galore --stringency 5 --illumina --fastqc --length 20 " + DIR_INPUT + ID + "_merged.fastq.gz -o " + DIR_OUT)
				else:
					print("Arquivo " + DIR_OUT + ID + " existente.")
					
			elif(layout == 1):
				if(not os.path.exists(DIR_OUT + ID + "_1_val_1.fq.gz") and not os.path.exists(DIR_OUT + ID + "_2_val_2.fq.gz")):
					print("Peired end")
					os.system("trim_galore --paired --stringency 5 --illumina --fastqc --length 20 " + DIR_INPUT + ID + "_1.fastq.gz " + DIR_INPUT + ID + "_2.fastq.gz -o " + DIR_OUT)
				else:
					print("Arquivo " + DIR_OUT + ID + " existente.")	

def QCfastqScreen(CODE, DIR_INPUT, DIR_OUT, layout, screenConf):
	print("\n\n---------------------FASTQ_SCREEN---------------------\n\n")
	if(screenConf == " "):
		print("Arquivo de configuração padrão")
	else:
		print("Arquivo de configuração inserido pelo usuario")
	for ID in CODE:
		print("Fastq_Screen " + ID)
		if(not layout):
			if(not os.path.exists(DIR_OUT + ID + "_trimmed_no_hits.fastq.gz")):
				os.system("fastq_screen --nohits --subset 0" + screenConf + " " + DIR_INPUT + ID + "_trimmed.fq.gz --outdir " + DIR_OUT )	
			else:
				print("Arquivo " + DIR_OUT + ID + " existente.")	
		else:
			if(layout == 2):
				if(not os.path.exists(DIR_OUT + ID + "_merged_trimmed_no_hits.fastq.gz")):
					os.system("fastq_screen --nohits --subset 0" + screenConf + " " + DIR_INPUT + ID + "_merged_trimmed.fq.gz --outdir " + DIR_OUT)
				else:
					print("Arquivo " + DIR_OUT + ID + " existente.")	
			else:
				if(os.path.exists(DIR_OUT + ID + "_1_val_1.fq.gz") and os.path.exists(DIR_OUT + ID + "_2_val_2.fq.gz")):
					os.system("fastq_screen --nohits --subset 0" + screenConf + " " + DIR_INPUT + ID + "_1_val_1.fq.gz " + DIR_INPUT + ID + "_2_val_2.fq.gz --outdir " + DIR_OUT)
				else:
					print("Arquivo " + DIR_OUT + ID + " existente.")
						
def QIIMEconfig(CODE, DIR_INPUT, DIR_OUT, layout):
	if(not os.path.exists(DIR_OUT + "config/")):
		os.system("mkdir " + DIR_OUT + "config/")
	if(not os.path.exists(DIR_OUT + "config/otu_settings.txt")):
		os.system("echo 'pick_otus:enable_rev_strand_match True' > " + DIR_OUT + "config/otu_settings.txt")
	if(not os.path.exists(DIR_OUT + "config/seqprep.txt")):
		os.system("echo 'join_paired_ends:pe_join_method SeqPrep' > " + DIR_OUT + "config/seqprep.txt")
	if(not os.path.exists(DIR_OUT + "alpha_params.txt")):
		os.system("echo alpha_diversity:metrics shannon,PD_whole_tree,chao1,observed_otus,shannon > " + DIR_OUT + "config/alpha_params.txt")

"""
ace, berger_parker_d, brillouin_d, chao1, chao1_ci, dominance, doubles, enspie, equitability, esty_ci, fisher_alpha, gini_index, goods_coverage, heip_e, kempton_taylor_q, margalef, mcintosh_d, mcintosh_e, menhinick, michaelis_menten_fit, observed_otus, observed_species, osd, simpson_reciprocal, robbins, shannon, simpson, simpson_e, singles, strong, PD_whole_tree
"""		