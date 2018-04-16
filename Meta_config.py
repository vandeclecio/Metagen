#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import os, re, sys, argparse
sys.path.append("Functions")
from Meta_QC import *
from Meta_AT import *
from Meta_AF import *
from Meta_PM import *
from Meta_functions import *

parser = argparse.ArgumentParser(description="Metagenomics Pipeline - version(1.0)")
parser.add_argument("-config", 	dest="config", 		nargs="?", 	required=True, default="config_meta",const="config_meta", help="Arquivo de configuração")
args = parser.parse_args()

if(args.config == "config_meta"):
	argumentos = []
	file = open(args.config)
	for line in file:
		line = line.strip().split("=")
		argumentos.append(line[1])
		
	mapping = argumentos[0]
	inputdir = argumentos[1]
	layout = argumentos[2]
	output = argumentos[3]
	screenConf = argumentos[4]
	seqprep = argumentos[5]
	instrument = argumentos[6]
	barcode = argumentos[7]
	otus = argumentos[8]
	refSeq = argumentos[9]
	taxonomy = argumentos[10]
	tree = argumentos[11]
	kaiju = argumentos[12]
	typesample = argumentos[13]

sampleFile = open(mapping)
sampleName = mapping
DATABASE_DMP = kaiju
DATABASE_FMI = kaiju
DIR_INPUT = inputdir


if(output == ''):
	output = "Metagenomics"
DIR_OUT = output + "/"

if(layout == "single"):
	layout = 0
else:
	layout = 1
	if(seqprep == "merge"):
		layout = 2		

if(screenConf == ''):
	screenConf = " "
else:
	screenConf = " --conf " + screenConf

if(instrument == "illumina"):
	instrument = "illumina"
else:
	instrument = "454"

if(barcode == "multiplex"):
	barcode = "M"
else:
	barcode = "D"

if(otus == "open"):
	otus = 0
elif(otus == "close"):
	otus = 1
elif(otus == "deNovo"):
	otus = 2

if(refSeq == ' '):
	reference = " "
else:
	reference = " -r " + refSeq

if(taxonomy == ''):
	taxonomy = " "
else:
	taxonomy = " -t " + taxonomy
	

if(tree == ''):
	tree = " --nonphylogenetic_diversity "
else:
	tree = " -t " + tree
	
if(not os.path.exists(DIR_OUT)):
	os.system("mkdir " + DIR_OUT)
elif(not os.path.exists(DIR_INPUT)):
	print("Diretorio das amostras não encontrado")
	exit()

print("Iniciando Processo")
print("\n--------------------------------------\n")	


CODE = sampleMetaData(sampleFile, 0)
sampleFile.seek(0,0)
NAME = sampleMetaData(sampleFile, 0)
sampleFile.close()
dic_name = dict(zip(CODE, NAME))

DADO_BRUTO = DIR_OUT + "0.DADO_BRUTO/"
DADO_BRUTO_FQ = DIR_OUT + "0.DADO_BRUTO/FASTQC/"
CONTROLE_QUALIDADE =  DIR_OUT + "1.CONTROLE_QUALIDADE/"
TRIMAGEM = DIR_OUT + "1.CONTROLE_QUALIDADE/1.TRIMAGEM/"
CONTAMINANTES = DIR_OUT + "1.CONTROLE_QUALIDADE/2.CONTAMINANTES/"

## Controle de qualidade das amoastras
createDiretory(1, DIR_OUT)
# lista, string, string, inteiro
QCfastqDump(CODE, DIR_INPUT, DADO_BRUTO, layout)

QCfastqc(CODE, DADO_BRUTO, DADO_BRUTO_FQ, layout)

QCSeqPrep(CODE, DADO_BRUTO, DADO_BRUTO, layout)

QCtrimGalore(CODE, DADO_BRUTO, TRIMAGEM, layout)

# lista, string, string, inteiro, string
QCfastqScreen(CODE, TRIMAGEM, CONTAMINANTES, layout, screenConf)

rename(CODE, CONTAMINANTES, CONTAMINANTES, layout)

print("\n\nConcluido processo de controle de qualidade\n\n")

# """
# 	Analise Taxonomica
# """
print("\n\nIniciando analise Funcional\n\n")

createDiretory(2, DIR_OUT)
DIR_INPUT = CONTAMINANTES
DIR_OUT = DIR_OUT + "2.ANALISE_TAXONOMIA/"

if(typesample == "16S"):
	QIIMEconfig(CODE, DIR_INPUT, DIR_OUT)
	# O comando não aceita diretorios com acentos
	if(layout == 1):
		os.system("multiple_join_paired_ends.py -i " + DIR_INPUT + " -o " + DIR_INPUT + " -p " + DIR_OUT + "config/seqprep.txt")

	#join_paired_ends.py -f $PWD/forward_reads.fastq -r $PWD/reverse_reads.fastq -o $PWD/fastq-join_joined

	validate_mapping(DIR_INPUT, sampleName)

	joinFastq, sample = sample_names_ids(DIR_INPUT, dic_name, layout)

	split_library(DIR_INPUT, DIR_OUT, joinFastq, sample, barcode)

	DIR_INPUT = DIR_OUT
	count_seqs(DIR_INPUT)

	pick_operational_taxonomic_units(DIR_INPUT, DIR_INPUT, reference, taxonomy, otus)

	biom_summarize(DIR_INPUT, otus)

	score = sampling_depth(DIR_INPUT)

	core_diversity_analyses(DIR_INPUT, DIR_INPUT, sampleName, tree, score, otus)

	make_emperor_plot(DIR_INPUT, score, sampleName)

	pick_otus(DIR_INPUT, DIR_OUT, method)

	pick_rep_set(DIR_INPUT, DIR_OUT, fasta, reference)

	assign_taxonomy(DIR_INPUT, DIR_OUT, taxonomy, reference)

	align_seqs()

	filter_ilignment()


	# fastq-join
	# Para instalar sudo apt-get install ea-utils
	# 

	if(instrument == "454"):
		#if(os.system("split_libraries.py -m Fasting_Map.txt -f Fasting_Example.fna -q Fasting_Example.qual -o split_library_output"))
		os.system("summarize_taxa_through_plots.py -i otus/otu_table.biom -o taxa_summary -m Fasting_Map.txt")
		os.system("make_otu_heatmap.py -i taxa_summary/otu_table_L3.biom -o taxa_summary/otu_table_L3_heatmap.pdf -c Treatment -m Fasting_Map.txt")
		os.system("alpha_rarefaction.py -i otus/otu_table.biom -m Fasting_Map.txt -o arare -p alpha_params.txt -t otus/rep_set.tre")
		os.system("beta_diversity_through_plots.py -i otus/otu_table.biom -m Fasting_Map.txt -o bdiv_even146 -t otus/rep_set.tre -e 146")
		os.system("jackknifed_beta_diversity.py -i otus/otu_table.biom -t otus/rep_set.tre -m Fasting_Map.txt -o jack -e 110")
		os.system("make_bootstrapped_tree.py -m jack/unweighted_unifrac/upgma_cmp/master_tree.tre -s jack/unweighted_unifrac/upgma_cmp/jackknife_support.txt -o jack/unweighted_unifrac/upgma_cmp/jackknife_named_nodes.pdf")
		os.system("make_emperor.py -i bdiv_even146/unweighted_unifrac_pc.txt -m Fasting_Map.txt -t taxa_summary/otu_table_L3.txt --n_taxa_to_keep 5 -o biplots")
	"""
split_libraries.py -m Fasting_Map.txt -f Fasting_Example.fna -q Fasting_Example.qual -o split_library_output
	ou 
split_libraries.py -m Fasting_Map_reverse_primers.txt -f Fasting_Example.fna -q Fasting_Example.qual -z truncate_only -o split_library_output_revprimers

pick_de_novo_otus.py -i split_lsample_names_ids(DIR_INPUT)ibrary_output/seqs.fna -o otus

## Caso seja necessario remover os ruido do sequenciamento podemos utilizar o seguinte comando .
denoiser.py -i 454Reads.sff.txt -f seqs.fna -v -o Outdir

biom summarize-table -i otus/otu_table.biom > otus/otu_table.txt 

summarize_taxa_through_plots.py -i otus/otu_table.biom -o taxa_summary -m Fasting_Map.txt

make_otu_heatmap.py -i taxa_summary/otu_table_L3.biom -o taxa_summary/otu_table_L3_heatmap.pdf -c Treatment -m Fasting_Map.txt

alpha_rarefaction.py -i otus/otu_table.biom -m Fasting_Map.txt -o arare -p alpha_params.txt -t otus/rep_set.tre

beta_diversity_through_plots.py -i otus/otu_table.biom -m Fasting_Map.txt -o bdiv_even146 -t otus/rep_set.tre -e 146

jackknifed_beta_diversity.py -i otus/otu_table.biom -t otus/rep_set.tre -m Fasting_Map.txt -o jack -e 110

make_bootstrapped_tree.py -m jack/unweighted_unifrac/upgma_cmp/master_tree.tre -s jack/unweighted_unifrac/upgma_cmp/jackknife_support.txt -o jack/unweighted_unifrac/upgma_cmp/jackknife_named_nodes.pdf

make_emperor.py -i bdiv_even146/unweighted_unifrac_pc.txt -m Fasting_Map.txt -t taxa_summary/otu_table_L3.txt --n_taxa_to_keep 5 -o biplots

# Caso as amostras sejam paired inicialmente utilizamos esse comando que faz a mesma coisa que o seqPrep
join_paired_ends.py -f $PWD/forward_reads.fastq -r $PWD/reverse_reads.fastq -o $PWD/fastq-join_joined
"""


# echo alpha_diversity:metrics shannon,PD_whole_tree,chao1,observed_otus > alpha_params.txt
else:
	AT_kaiju(DIR_INPUT, DIR_OUT, CODE, layout, DATABASE_DMP, DATABASE_FMI)

print("Finalizando processo")
print("\n--------------------------------------\n")


"""
	Analise Funcional
"""
print("\n\nIniciando analise taxonomica\n\n")
print(DIR_INPUT)
createDiretory(3, DIR_INPUT)

DIR_INPUT = DIR_OUT + "3.ANALISE_FUNCIONAL/"

