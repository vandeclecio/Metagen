# Validate the mapping file
 * validate_mapping_file.py -o vmf-map/ -m map.tsv

# Use split-libraries.py to split the input files according to barcode
 * split_libraries_fastq.py -o slout/ -i forward_reads.fastq.gz -b barcodes.fastq.gz -m map.tsv
 * count_seqs.py -i slout/seqs.fna

# Pick Operational Taxonomic Units

 * ## pick_de_novo_otus.py
    * pick_de_novo_otus.py -i $PWD/seqs.fna -o $PWD/uclust_otus/
 * ## pick_closed_reference_otus.py
    * pick_closed_reference_otus.py -i $PWD/seqs.fna -r $PWD/refseqs.fna -o $PWD/otus_w_tax/ -t $PWD/taxa.txt
 * ## pick_open_reference_otus.py
    * pick_open_reference_otus.py -o otus/ -i slout/seqs.fna -p ../uc_fast_params.txt

 * biom summarize-table -i otus/otu_table.biom

# Run diversity analyses
 * core_diversity_analyses.py -o cdout/ -i otus/otu_table.biom -m map.tsv -t otus/rep_set.tre -e 1114
 * core_diversity_analyses.py -o cdout/ --recover_from_failure -c "SampleType,DaysSinceExperimentStart" -i otus/otu_table.biom -m map.tsv -t otus/rep_set.tre -e 1114
 * make_emperor.py -i cdout/bdiv_even1114/weighted_unifrac_pc.txt -o cdout/bdiv_even1114/weighted_unifrac_emperor_pcoa_plot -m map.tsv --custom_axes DaysSinceExperimentStart
 * make_emperor.py -i cdout/bdiv_even1114/unweighted_unifrac_pc.txt -o cdout/bdiv_even1114/unweighted_unifrac_emperor_pcoa_plot -m map.tsv --custom_axes DaysSinceExperimentStart

# Calculate Alpha Diversity
 * alpha_rarefaction.py -i otu_table.biom -o alpha_output_folder -m mapping_file.txt -t rep_tree.tre
 * _This command will generate a new text file in a current working directory._ echo "alpha_diversty:metrics observed,shannon,simpson,pd_whole_tree,chao1" >> alpha_parameters.txt
 * compare_alpha_diversity.py -i alpha_output/alpha_div_collated/PD_whole_tree.txt -o alpha_pdwholetree_stats -m mapping_file.txt -t nonparametric -c SampleType
 * compare_alpha_diversity.py -i alpha_output/alpha_div_collated/chao1.txt -o alpha_chao1_stats -m mapping_file.txt -t nonparametric -c SampleType
 * compare_alpha_diversity.py -i alpha_output/alpha_div_collated/observed_otus.txt -o alpha_observed_otus_stats -m mapping_file.txt -t nonparametric -c SampleType
 * add_alpha_to_mapping_file.py -i alpha_div_collated/PD_whole_tree.txt -m mapping_file.txt -o mapping_file_with_alpha.txt
 * alpha_rarefaction.py -i otu_table.biom -o alpha_output_folder -m mapping_file.txt -t rep_tree.tre -p parameters.txt


# Calculate Beta Diversity
 * beta_diversity_through_plots.py -i otu_table.biom -o bdiv_plots/ -m mapping_file.txt -t rep_set.tre -p parameters.txt
 * compare_categories.py -i bdiv_plots/unweight_unifrac_dm.txt -o bdiv_stats_adonis_unweighted/ -m mapping_file.txt -c SampleType --method adonis
 * compare_categories.py -i bdiv_plots/weight_unifrac_dm.txt -o bdiv_stats_adonis_weighted/ -m mapping_file.txt -c SampleType --method adonis
 * make_2d_plots.py -i bdiv_plots/unweight_unifrac_pc.txt -o bdiv_2d_plot/ -m mapping_file.tx
 * make_distance_boxplots.py -d bdiv_plots/unweighted_unifrac_dm.txt -o bdiv_plots/unweighted_distance_boxplot -m mapping_file.txt -f "SampleType" --save_raw_data
 * make_emperor.py -i bdiv_plots/unweighted_unifrac_pc.txt -o bdiv_plots_axes -m mapping_file.txt --custom_axes Day
 * make_emperor.py -i bdiv_plots/unweighted_unifrac_pc.txt -o bdiv_plots_vector -m mapping_file.txt --add_vectors Subject
 * make_emperor.py -i bdiv_plots/unweighted_unifrac_pc.txt -o bdiv_plots_colorby -m mapping_file.txt --color_by "Subject&&Day"
 * summarize_taxa.py -i otu_table.biom -o summarize_taxa
 * make_emperor.py -i bdiv_plots/unweighted_unifrac_pc.txt -o bdiv_plots_biplot -m mapping_file.txt --taxa_fp summarize_taxa/otu_table_L6.txt -n 5 --biplot_fp biplot.txt

# Test whether any OTU is significantly associated with a particular experimental category

# Testing for differences in OTU abundance at different taxonomic levels

# Next step
 * ## Comparing Distance Matrices
    * compare_distance_matrices.py --method=mantel -i unweighted_unifrac_dm.txt,PH_dm.txt -o mantel_out -n 999
 * ## Running Supervised Learning
    * supervised_learning.py -i otu_table.biom -m Fasting_Map.txt -c Treatment -o ml -v
 * ##  Performing Procrustes Analysis
    * transform_coordinate_matrices.py -i illumina/precomputed-output/cdout/bdiv_even1114/unweighted_unifrac_pc.txt,illumina/precomputed-output/cdout/bdiv_even1114/weighted_unifrac_pc.txt -r 999 -o procrustes_results/
    * make_emperor.py -c -i procrustes_results/ -o procrustes_results/plots/ -m illumina/map.tsv --custom_axes DaysSinceExperimentStart