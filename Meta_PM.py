#!/usr/bin/env python
#-*- encoding:utf-8 -*-

import os
import sys
import multiprocessing

# Pega a quantidade de processadores total do computador menos 1
T = multiprocessing.cpu_count() - 1

# FUNCIONAL = "2.TAXONOMIA_FUNCIONAL/"
# PM_META = "2.TAXONOMIA_FUNCIONAL/PM-META/"

# os.system("mkdir " + FUNCIONAL)

# if(os.path.exists("PM-META/")):
# 	print("PM-META/")
# else:
# 	os.system("mkdir " + PM_META)

#metaFile = open(sys.argv[2], 'r')
# seqsFile = open(FUNCIONAL+ "seqs.list", 'w')
# metaNew = open(FUNCIONAL + "meta.txt", 'w')

def AF_parallelFiles(DIR_INPUT, DIR_OUT):
	for DIR in os.listdir(DIR_INPUT):
		if(sys.argv[1] == "-s"):
			ID = DIR.split("_")
			if(ID[1] == "merged" and ID[3] == "no"):
				os.system("gunzip " + DIR_INPUT + ID[0] + "_merged_trimmed_no_hits.fastq.gz")
				os.system("ln " + DIR_INPUT + ID[0] + "_merged_trimmed_no_hits.fastq " + PM_META + ID[0] + "_merged_trimmed_no_hits.fastq")
				seqsFile.writelines(ID[0] + "\t" + PM_META + ID[0] + "_merged_trimmed_no_hits.fastq" + "\n")
			elif(ID[1] == "trimmed" and ID[2] == "no"):
				os.system("gunzip " + DIR_INPUT + ID[0] + "_trimmed_no_hits.fastq.gz")
				os.system("ln " + DIR_INPUT + ID[0] + "_trimmed_no_hits.fastq " + PM_META + ID[0] + "_trimmed_no_hits.fastq")
				seqsFile.writelines(ID[0] + "\t" + PM_META + ID[0] + "_trimmed_no_hits.fastq" + "\n")	
			

		else:
			ID = DIR.split("_")
			FQ = DIR.split(".")
			if(FQ[1] == "fq"):
				print("gunzip " + DIR_INPUT + ID[0] + "_[12]_no_hits.fq.gz")
				
				os.system("gunzip " + DIR_INPUT + ID[0] + "_1_no_hits.fastq.gz")
				os.system("gunzip " + DIR_INPUT + ID[0] + "_2_no_hits.fastq.gz")

				os.system("ln " + DIR_INPUT + ID[0] + "_1_no_hits.fastq " + DIR_OUT)
				os.system("ln " + DIR_INPUT + ID[0] + "_2_no_hits.fastq " + DIR_OUT)

				seqsFile.writelines(ID[0] + "\t" + DIR_OUT + ID[0] + "_1_no_hits.fastq" + "\n")
				seqsFile.writelines(ID[0] + "\t" + DIR_OUT + ID[0] + "_2_no_hits.fastq" + "\n")

	for line in metaFile:
	        metaNew.writelines(line)

	metaFile.close()
	seqsFile.close()
	metaNew.close()


""" Função para executar o Parallel Meta
"""
def AF_parallel(DIR_INPUT, DIR_OUT, CODE, layout):
	AF_parallelFiles(DIR_INPUT, DIR_OUT)
	if(sys.argv[1] == "-s"):
		os.system("PM-pipeline -i " + FUNCIONAL + "seqs.list -m " + FUNCIONAL + "meta.txt -t " + str(T) + " -o " + FUNCIONAL + "PM-pipeline/")
	else:
		os.system("PM-pipeline -i " + FUNCIONAL + "seqs.list -m " + FUNCIONAL + "meta.txt -t " + str(T) + " -P 0 -o " + FUNCIONAL + "PM-pipeline/")
