#!/usr/bin/env bash

main_dir='../../'
screen_data_path='./result/example-analyzed_data/2_base_PRAUC0.8_scRNA_seq_data_with_thv70_and_111_neuropal_reporter.csv'
neuron_data_path='./result/example-analyzed_data/confirm_Nmin_between_NeuroPALs/neuron_sets_and_recombination-union/neuron_sets_original-20-0.5-1-Both-False.xlsx'
gene_screened_data_path='./code/evaluation/example/union_gene_combinations_for_N20_combn11.log'

python evaluation_gene_combinations.py ${main_dir} ${gene_screened_data_path} ${neuron_data_path} ${screen_data_path}
echo "ALL FINISHED"
