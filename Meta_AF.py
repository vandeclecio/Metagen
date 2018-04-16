#!/usr/bin/envs python
# -*- encoding: utf-8 -*-

"""
	Analise Funcional

"""

import sys
import os
import multiprocessing
# Pega a quantidade de processadores total do computador menos 2
T = multiprocessing.cpu_count() - 1

""" Função para gerar a analise taxonomica utilizando o Kaiju
	@param 	DIR,INPUT, Diretorio contendo as amostras
	@param 	DIR_OUT, Diretorio para saida das amostras
	@param 	CODE, lista com os identificadores das amostras gerada pela função sampleMetaData
	@param 	layout, se as amostras são single-end ou paired-end
	@param 	DATABASE_DMP, base de dados do kaiju
	@param 	DATABASE_DMP, base de dados do kaiju 
""" 
def AT_kaiju(DIR_INPUT, DIR_OUT, CODE, layout, DATABASE_DMP, DATABASE_FMI):
	print("\n\n---------------------KAIJU---------------------\n\n")
	for ID in CODE:
		if(layout == "-s"):
			if(os.path.exists(DIR_INPUT + ID + "_00_L001_R1_001.fastq.gz")):
				os.system("kaiju -z " + str(T) + " -t " + DATABASE_DMP + "nodes.dmp -f " + DATABASE_FMI + "kaiju_db.fmi -i " + DIR_INPUT + ID + "_00_L001_R1_001.fastq.gz -o " + DIR_OUT + ID + ".kaiju.out")
		else: # Saber diferenciar amostras r1 e r2
			os.system("kaiju -z " + str(T) + " -t " + DATABASE_DMP + "nodes.dmp -f " + DATABASE_FMI + "kaiju_db.fmi -i " + DIR_INPUT + ID + " -j " + DIR_INPUT + ID + " -o " + ID + ".kaiju.out")	

	os.system("kaiju2krona -t " + DATABASE_DMP + "nodes.dmp -n " + DATABASE_DMP + "names.dmp -i " + ID + ".kaiju.out -o " + ID + ".kaiju.out.krona")
	os.system("ktImportText -o " + ID + ".kaiju.out.html " + ID + ".kaiju.out.krona")

	print("\n\n---------------------KRONA---------------------\n\n")
	KRONA = ""
	for ID in CODE:
		KRONA += ID + ".kaiju.out.krona "
		
	os.system("ktImportText -o all.kaiju.out.html " + KRONA)
