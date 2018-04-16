#!/usr/bin/envs python
#-*- encoding: utf-8 -*-
import os, re
import multiprocessing
T = multiprocessing.cpu_count() - 2

# Analise taxonomia

"""
"""
def summarize_table(DIR_INPUT, DIR_OUT, sampleName, category="-"):
	print("\n--------------------------------------\n")
	if(category != "-"):
		os.system("summarize_taxa_through_plots.py -i " + DIR_INPUT + " -o " + DIR_OUT + " -m " + sampleName + " -c " + category)
	else:
		os.system("summarize_taxa_through_plots.py -i " + DIR_INPUT + " -o " + DIR_OUT + " -m " + sampleName)	
"""
"""
def beta_diversity(DIR_INPUT, DIR_OUT, sampleName, tree, score):
	print("\n--------------------------------------\n")
	os.system("beta_diversity_through_plots.py -i " + DIR_INPUT + " -o " + DIR_OUT + "bdiv_even/" + score + " -t " + tree + " -m " + sampleName + " -e " + score + " -aO " + srt(T))

"""
"""
def alpha_diversity(DIR_INPUT, DIR_OUT, sampleName, tree, score):
	print("\n--------------------------------------\n")
	os.system("alpha_rarefaction.py -i " + DIR_INPUT + " -o " + DIR_OUT + "arare_max" + score + "/ -t " + tree + " -m " + sampleName + " -e " + score + " -aO " + str(T))

""" Função faz as analizes de diversidade alfa é beta de forna "automatica".
"""
def core_diversity_analyses(DIR_INPUT, DIR_OUT, sampleName, tree, score, otus):
	print("\n--------------------------------------\n")
	print("Analise de diversidade")
	if(otus == 0):
		print("core_diversity_analyses.py -i " + DIR_INPUT + "otus/otu_table_mc2.biom -o " + DIR_OUT + "cdout -m " + sampleName + " -t " + DIR_OUT + "otus/rep_set.tre -e " + score + " -aO " + str(T))
		os.system("core_diversity_analyses.py -i " + DIR_INPUT + "otus/otu_table_mc2.biom -o " + DIR_OUT + "cdout -m " + sampleName + " -t " + DIR_OUT + "otus/rep_set.tre -e " + score + " -aO " + str(T))

	elif(otus == 1):
		print("core_diversity_analyses.py -i " + DIR_INPUT + "otus/otu_table.biom -o " + DIR_OUT + "cdout -m " + sampleName + tree + " -e " + score + " -p " + DIR_INPUT + "config/alpha_params.txt -aO " + str(T))
		os.system("core_diversity_analyses.py -i " + DIR_INPUT + "otus/otu_table.biom -o " + DIR_OUT + "cdout -m " + sampleName + tree + " -e " + score + " -p " + DIR_INPUT + "config/alpha_params.txt -aO " + str(T))
		os.system("core_diversity_analyses.py -i " + DIR_INPUT + "otus/otu_table.biom -o " + DIR_OUT + "cdout --recover_from_failure -c \"Treatment\" -m " + sampleName + tree + " -e " + score + " -p " + DIR_INPUT + "config/alpha_params.txt -aO " + str(T))

	elif(otus == 2):
		print("core_diversity_analyses.py -i " + DIR_INPUT + "otus/otu_table.biom -o " + DIR_OUT + "cdout -m " + sampleName + " -t " + DIR_OUT + "otus/rep_set.tre -e " + score + " -p " + DIR_INPUT + "config/alpha_params.txt -aO " + str(T))
		os.system("core_diversity_analyses.py -i " + DIR_INPUT + "otus/otu_table.biom -o " + DIR_OUT + "cdout -m " + sampleName + " -t " + DIR_OUT + "otus/rep_set.tre -e " + score + " -p " + DIR_INPUT + "config/alpha_params.txt -aO " + str(T))
		# Verificar as caracteristicas que podem ser utilizadas (percorrer primeira linha do samplemeta e pegar as caracteristicas)
		os.system("core_diversity_analyses.py -i " + DIR_INPUT + "otus/otu_table.biom -o " + DIR_OUT + "cdout --recover_from_failure -c \"Treatment\" -m " + sampleName + " -t " + DIR_OUT + "otus/rep_set.tre -e " + score + " -p " + DIR_INPUT + "config/alpha_params.txt -aO " + str(T))

""" Para fazer download das OTU's referencia: ftp://greengenes.microbio.me/greengenes_release/gg_13_5/ ou ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/
"""
def pick_operational_taxonomic_units(DIR_INPUT, DIR_OUT, reference, taxonomy, otus):
	print("\n--------------------------------------\n")
	print("Pick Operational Taxonomic Units")
	if(otus == 0):
		print("pick_open_reference_otus")
		os.system("pick_open_reference_otus.py -i " + DIR_INPUT + "slout/seqs.fna -o " + DIR_OUT + "otus/ -aO " + str(T) + " -p " + DIR_OUT + "config/otu_settings.txt " + reference)

	elif(otus == 1):
		print("pick_closed_reference_otus")
		os.system("pick_closed_reference_otus.py -i " + DIR_INPUT + "slout/seqs.fna -o " + DIR_OUT + "otus/ -aO " + str(T) + " -p " + DIR_OUT + "config/otu_settings.txt " + reference + taxonomy)

	elif(otus == 2):
		print("pick_de_novo_otus")
		os.system("pick_de_novo_otus.py -i " + DIR_INPUT + "slout/seqs.fna -o " + DIR_OUT + "/otus/ -aO " + str(T) + " -p " + DIR_OUT + "config/otu_settings.txt")

""" Verifica qual amostra tem a menor sequencia apartir da saida da função biom_summarize()
	e a partir dela define um aprofundidade minima para analise das amostras
"""
def sampling_depth(DIR_OUT):
	print("\n--------------------------------------\n")
	print("Verificando sequencias")
	scoreEqual = open(DIR_OUT + "otus/otu_table_stats.txt", "r")
	for e in scoreEqual:
		e = e.split(":")
		if(e[0] == " Min"):
			score = e[1].strip().split(".")
			score = score[0]

	scoreEqual.close()
	return score

""" Nomeia as amostra de acordo com o sampleName do sample-meta.
"""
def sample_names_ids(DIR_INPUT, dic_name, layout):
	print("\n--------------------------------------\n")
	print("Nomeando amostras")
	joinFastq = ""
	sample = ""
	i = 1
	for ID in dic_name:
		if(layout == 0):
			joinFastq += DIR_INPUT + ID + "_00_L001_R1_001.fastq.gz"
		else:
			joinFastq += DIR_INPUT + ID + "_00_L001_R1_001/seqprep_assembled.fastq.gz"
		
		sample += dic_name[ID]
		if(i != len(dic_name)):
			i+=1
			joinFastq += ","
			sample += ","
	return joinFastq, sample

""" Verifica o sample-metadata.txt e gera um arquivo .html informando os problemas no arquivo.
	validate_mapping_file.py -m sample-meta.txt -o mapping_output/
"""
def validate_mapping(DIR_OUT, sampleName):
	print("\n--------------------------------------\n")
	print("Verificando mapping file")
	os.system("validate_mapping_file.py -m " + sampleName + " -o " + DIR_OUT + "mapping_output")

""" Gera a visualização do PCoA em 3D
	make_emperro.py -i cdout/bdiv_even000/*_unifrac_pc.txt -o cdout/bdiv_even000/*weighted_unifrac_pcoa_plot/ -m sample-meta.txt --custo_axes Month
"""
def make_emperor_plot(DIR_OUT, score, sampleName):
	print("\n--------------------------------------\n")
	# Verificar qual valor vode ser os da fleg --custom_axes
	os.system("make_emperor.py -i " + DIR_OUT + "cdout/bdiv_even" + score + "/weighted_unifrac_pc.txt -o " + DIR_OUT + "cdout/bdiv_even" + score + "/weighted_unifrac_emperor_pcoa_plot/ -m " + sampleName + " --custom_axes Month")
	os.system("make_emperor.py -i " + DIR_OUT + "cdout/bdiv_even" + score + "/unweighted_unifrac_pc.txt -o " + DIR_OUT + "cdout/bdiv_even" + score + "/unweighted_unifrac_emperor_pcoa_plot/ -m " + sampleName + " --custom_axes Month")

""" Criar um arquivo mostrando a quantidade de sequencias em cada amostra.
	biom summarize-table -i otus/otus_table.biom -o otus/otus_table_stats.txt
"""
def biom_summarize(DIR_OUT, otus):
	print("\n--------------------------------------\n")
	print("Summarize")
	if(otus == 0):
		os.system("biom summarize-table -i " + DIR_OUT + "otus/otu_table_mc2.biom  -o " + DIR_OUT + "otus/otu_table_stats.txt")
	elif(otus == 1):
		os.system("biom summarize-table -i " + DIR_OUT + "otus/otu_table.biom -o " + DIR_OUT + "otus/otu_table_stats.txt")
	elif(otus == 2):
		os.system("biom summarize-table -i " + DIR_OUT + "otus/otu_table.biom -o " + DIR_OUT + "otus/otu_table_stats.txt")	
""" Conta o número de sequencias no arquivo de saida da função split_library
	count_seqs.py -i slout/seqs.fna > slout/count_seqs.txt
"""
def count_seqs(DIR_INPUT):
	print("\n--------------------------------------\n")
	print("Contando quantidade de reads")
	os.system("count_seqs.py -i " + DIR_INPUT + "slout/seqs.fna > " + DIR_INPUT + "slout/count_seqs.txt")

""" Função que reuni todos as amostras em um unico arquivo e nomeia os ids 
	split_libraries_fastq.py -i seqA.fastq.gz,seqB.fastq.gz -o slout/ --sample_ids sampleA,sampleB --barcode_type 'not-barcoded' --phred_offset 33 -q 19 -r 1000
"""
def split_library(DIR_INPUT, DIR_OUT, joinFastq, sample, barcode):
	print("\n--------------------------------------\n")
	print("Reunindo arquivos")
	if(barcode == "D"):
		os.system("split_libraries_fastq.py -i " + joinFastq + " -o " + DIR_OUT + "slout/ --sample_ids " + sample + " --barcode_type \'not-barcoded\' --phred_offset 33 -q 19 -r 1000")
	else:
		print("split_libraries_fastq.py -i lane1_read1.fastq.gz -o " + DIR_OUT + "slout_q20/ -b lane1_barcode.fastq.gz --rev_comp_mapping_barcodes -m map.txt --store_qual_scores -q 19")
		print("Ainda não foi feito ...")
		exit()

""" Função para gerar o PCoA plot 2d usando o arquivo de coordenadas 
	principal gerado pela realização de medidas de diversidade beta de uma tabela OTU.
	parametros: arquivo de coordenadas, diretorio de saida, sample-metadata, categoria
	exemplo: make_2d_plots.py -i bdiv_even/unweighted_unifrac_pc.txt -m sample-metadata.txt -b 'tissue' -o output/ --scree
"""
def make_2d_plot(DIR_INPUT, DIR_OUT, sampleName, category):
	print("\n--------------------------------------\n")
	print("Criando PCoA plot 2D")
	os.system("make_2d_plots.py -i " + DIR_INPUT + " -m " + sampleName + " --scree -o " + DIR_OUT + "/2d_plot")
	os.system("make_2d_plots.py -i " + DIR_INPUT + " -m " + sampleName + " -b \'" + category + "\' --scree -o " + DIR_OUT + "/2d_plot_" + category)

""" Gera uma tabela heatmap correspondente a abundância relativa de uma OTU
"""
def make_heatmap(DIR_INPUT, DIR_OUT, sampleName, category="-", imageType="pdf"):
	print("\n--------------------------------------\n")
	print("Criando heatmap")
	if(category != "-"):
		os.system("make_otu_heatmap.py -i " + DIR_INPUT + " -o " + DIR_OUT + "heatmap -c " + category + " -m " + sampleName + " -g " + imageType)
	else:
		os.system("make_otu_heatmap.py -i " + DIR_INPUT + " -o " + DIR_OUT + "heatmap -m " + sampleName + " -g " + imageType)

def jackknifed_beta_diversity(DIR_INPUT, DIR_OUT, sampleName, referencia, score):
	print("\n--------------------------------------\n")
	print("Criando matriz de distancia")
	os.system("jackknifed_beta_diversity.py -i " + DIR_INPUT + " -t " + referencia + " -m " + sampleName +" -o " + DIR_OUT + "jackknife/ -e " + score)


def QIIMEconfig(CODE, DIR_INPUT, DIR_OUT):
	print("\n--------------------------------------\n")
	print("Criando arquivos de configuração")
	if(not os.path.exists(DIR_OUT + "config/")):
		os.system("mkdir " + DIR_OUT + "config/")
	if(not os.path.exists(DIR_OUT + "config/otu_settings.txt")):
		os.system("echo 'pick_otus:enable_rev_strand_match True' > " + DIR_OUT + "config/otu_settings.txt")
	if(not os.path.exists(DIR_OUT + "config/seqprep.txt")):
		os.system("echo 'join_paired_ends:pe_join_method SeqPrep' > " + DIR_OUT + "config/seqprep.txt")
	if(not os.path.exists(DIR_OUT + "alpha_params.txt")):
		os.system("echo alpha_diversity:metrics chao1,observed_otus,shannon,PD_whole_tree > " + DIR_OUT + "config/alpha_params.txt")
		#os.system("echo alpha_diversity:metrics PD_whole_tree,chao1,observed_otus,shannon > " + DIR_OUT + "config/alpha_params.txt")
"""
ace,berger_parker_d,brillouin_d,chao1,chao1_ci,dominance,doubles,enspie,equitability,esty_ci,fisher_alpha,gini_index,goods_coverage,heip_e,kempton_taylor_q,margalef,mcintosh_d,mcintosh_e,menhinick,michaelis_menten_fit,observed_otus,observed_species,osd,simpson_reciprocal,robbins,shannon,simpson,simpson_e,singles,strong,PD_whole_tree
"""

"""	Atribui sequências semelhantes a unidades taxonômicas operacionais, ou OTUs, agrupando sequências com base em um limite de similaridade definido pelo usuário
Exemplo: pick_otus.py -i slout/seqs.fna -o 1.picked_otus
"""
def pick_otus(DIR_INPUT, DIR_OUT, method="-"):
	if(method != "-"):
		os.system("pick_otus.py -i " + DIR_INPUT + " -o " + DIR_OUT + " -m " + method)
	else:
		os.system("pick_otus.py -i " + DIR_INPUT + " -o " + DIR_OUT)

"""	Depois de selecionar OTUs, você pode escolher um conjunto representativo de seqüências. Para cada OTU, você acabará 
com uma sequência que pode ser usada em análises subsequentes
Exemplo: pick_rep_set.py -i seqs_otus.txt -f seqs.fna -o rep_set1.fna
"""
def pick_rep_set(DIR_INPUT, DIR_OUT, fasta, reference='-'):
	if(reference != "-"):
		os.system("pick_rep_set.py -i " + DIR_INPUT + " -f " + fasta + " -r " + reference + " -o " + DIR_OUT)
	else:
		os.system("pick_rep_set.py -i " + DIR_INPUT + " -f " + fasta + " -o " + DIR_OUT)	

"""	Código para atribuir taxonomia, usando várias técnicas.
Exemplo: assign_taxonomy.py -i repr_set_seqs.fasta -r ref_seq_set.fna -t id_to_taxonomy.txt
"""
def assign_taxonomy(DIR_INPUT, DIR_OUT, taxonomy, reference):
	print("Olá mundo!")

"""	Alinha as sequências em um arquivo FASTA entre si ou para um alinhamento de seqüência de modelo, dependendo do método escolhido
Exemplo:
"""
def align_seqs():
	print("Olá mundo!")

"""	Este script irá remover posições que são lacunas em todas as seqüências 
(comuns para PyNAST, pois as seqüências típicas cobrem apenas 200-400 bases e estão sendo alinhadas contra o gene 16S completo)
Exemplo:
"""
def filter_ilignment():
	print("Olá mundo!")

"""	Restaurar um arquivo do fluxgram no formato .sff.txt, que é o resultado do sffinfo.
Exemplo: 
"""
def denoise_wrapper(DIR_INPUT, DIR_OUT, fasta, map_fname):
	print("Olá mundo!")

"""	produz esta árvore a partir de um alinhamento de seqüência múltiplas.
Exemplo:
"""
def make_phylogeny():
	print("Ola mundo!")

"""	Tabula o número de vezes que um OTU é encontrado em cada amostra e adiciona as previsões taxonômicas para
cada OTU na última coluna se um arquivo de taxonomia for fornecido.
Exemplo:
"""
def make_otu_table():
	print("Olá mundo!")

"""	Realiza o merge de todos os arquivos fasta das sequências criando um arquivo table_otu.biom
"""
def merge_otu_table():
	print("não implementada")