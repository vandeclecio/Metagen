#!/usr/bin/env python
#-*- encoding:utf-8 -*-
import os, re

## Função para renomear as amostras
# @param 	sampleFile	arquivo contendo o identificador e os detalhes de cada amostra
# @param 	coluna	coluna rsponsavel para nomear a lista
# @return 	Lista contendo os identificadores das amostras
def sampleMetaData(sampleFile,	coluna):
	LIST = []
	for line in sampleFile:
		line = line.strip().split("\t")
		if bool(not re.search("^#", line[coluna])):
			LIST.append(line[coluna])

	return LIST


## Função criar os diretorios das etapas do pipeline
# @param 	sampleFile	arquivo contendo o identificador e os detalhes de cada amostra
# @param 	coluna	coluna rsponsavel para nomear a lista
# @return 	Lista contendo os identificadores das amostras
def createDiretory(DIR, output):
	if(DIR == 1):
		DADO_BRUTO = output + "0.DADO_BRUTO/"
		DADO_BRUTO_FQ = output + "0.DADO_BRUTO/FASTQC/"
		CONTROLE_QUALIDADE =  output + "1.CONTROLE_QUALIDADE/"
		TRIMAGEM = output + "1.CONTROLE_QUALIDADE/1.TRIMAGEM/"
		CONTAMINANTES = output + "1.CONTROLE_QUALIDADE/2.CONTAMINANTES/"
		os.system("mkdir " + DADO_BRUTO + " " + CONTROLE_QUALIDADE)
		os.system("mkdir " + TRIMAGEM + " " + CONTAMINANTES + " " + DADO_BRUTO_FQ)
	elif(DIR == 2):
		TAXONOMI = output + "2.ANALISE_TAXONOMIA/"
		os.system("mkdir " + TAXONOMI)
	elif(DIR == 3):
		FUNCIONAL = output + "3.ANALISE_FUNCIONAL/"
		os.system("mkdir " + FUNCIONAL)
	else:
		print("Não foi possivel criar o diretorio!!")
		exit()	

