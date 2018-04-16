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
parser.add_argument("-config", dest="config", nargs="?", required=False, default="config_meta",const="config_meta", help="Arquivo de configuração")
parser.add_argument("-m", dest="mapping", 	nargs="?", 	required=True,	help="File containing the information of the samples to be processed.")
parser.add_argument("-d", dest="dir", 		nargs="?", 	required=True,	default="A", const="A", help="Diretory containing samples")
parser.add_argument("-l", dest="layout",	nargs="?", 	required=True,	default="single", const="single", choices=["single", "paired"], help="single-end sequences paired-end sequences")
parser.add_argument("-o", dest="output", 	nargs="?", 	required=False,	default="Metagenomics/", const="Metagenomics/", help="Diretory for output")
parser.add_argument("-c", dest="screenConf",nargs="?",	required=False,	default="F", const="F", help="Arquivo de configuração do fastq_screen")
parser.add_argument("-s", dest="seqprep",	nargs="?",	required=False,	default="merge", const="merge", choices=["merge", "normal"],help="Merge nas sequencias paired-end")
parser.add_argument("-i", dest="instrument",nargs="?", 	required=False,	default="illumina", const="illumina", choices=["illumina", "454"], help="Instrument usign for samples")
parser.add_argument("-b", dest="barcode", 	nargs="?", 	required=False,	default="demultiplex", const="demultiplex", choices=["multiplex", "demultiplex"], help="multiplex or demultiplex")
parser.add_argument("-p", dest="otus",		nargs="?",	required=False,	default="open", const="open", choices=["open", "close", "deNovo"], help="http://qiime.org/tutorials/otu_picking.html")
parser.add_argument("-r", dest="refSeq",	nargs="?", 	required=False, default="97_otus", const="97_otus", help="The reference sequences")
parser.add_argument("-t", dest="taxonomy",	nargs="?",	required=False, default="97_otu_taxonomy", const="97_otu_taxonomy", help="The taxonomy map")
parser.add_argument("-a", dest="tree",		nargs="?",	required=False, default="97_otus_tree", const="97_otus_tree", help="Tree for pick_closed_reference_otus")
parser.add_argument("-k", dest="kaiju",		nargs="?", 	required=True,	help="Informe o banco de dados da analise taxonomica")
parser.add_argument("-type", dest="type",	nargs="?", 	required=True,	default="16S", const="16S", choices=["16S","WGS"], help="16S ou WGS")
args = parser.parse_args()


if(len(sys.argv) < 5):
	print("Informe corretamente os paramentros")
	exit()

if(args.layout == "single"):
	layout = 0
else:
	layout = 1
	if(args.seqprep == "merge"):
		layout = 2		

if(args.screenConf == "F"):
	screenConf = " "
else:
	screenConf = " --conf " + args.screenConf

if(args.instrument == "illumina"):
	instrument = "illumina"
else:
	instrument = "454"

if(args.barcode == "multiplex"):
	barcode = "M"
else:
	barcode = "D"

if(args.otus == "open"):
	otus = 0
elif(args.otus == "close"):
	otus = 1
elif(args.otus == "deNovo"):
	otus = 2

if(args.refSeq != "97_otus"):
	reference = " -r " + args.refSeq
else:	
	reference = " "

if(args.taxonomy != "97_otu_taxonomy"):
	taxonomy = " -t " + args.taxonomy
else:
	taxonomy = " "

if(args.tree != "97_otus_tree"):
	tree = " -t " + args.tree
else:
	tree = " --nonphylogenetic_diversity "


DATABASE_DMP = args.kaiju
DATABASE_FMI = args.kaiju


DIR_INPUT = args.dir
DIR_OUT = args.output

if(not os.path.exists(DIR_OUT)):
	os.system("mkdir " + DIR_OUT)
elif(not os.path.exists(DIR_INPUT)):
	print("Diretorio das amostras não encontrado")
	exit()

print("Iniciando Processo")
print("\n--------------------------------------\n")	
sampleFile = open(args.mapping)
sampleName = args.mapping

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

"""
	Analise Taxonomica
"""
print("\n\nIniciando analise Funcional\n\n")

createDiretory(2, DIR_OUT)
DIR_INPUT = CONTAMINANTES
DIR_OUT = DIR_OUT + "2.ANALISE_TAXONOMIA/"

QIIMEconfig(CODE, DIR_INPUT, DIR_INPUT)

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
print("Finalizando processo")
print("\n--------------------------------------\n")


"""
	Analise Funcional
"""
print("\n\nIniciando analise taxonomica\n\n")
print(DIR_INPUT)
createDiretory(3, DIR_INPUT)

DIR_INPUT = DIR_OUT + "3.ANALISE_FUNCIONAL/"

AT_kaiju(DIR_INPUT, DIR_OUT, CODE, layout, DATABASE_DMP, DATABASE_FMI)