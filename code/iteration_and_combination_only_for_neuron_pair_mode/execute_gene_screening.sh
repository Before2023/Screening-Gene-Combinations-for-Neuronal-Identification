#!/usr/bin/env bash

for Nmin in 6
do
    main_dir="../../"
    at_least_gene_n=2
    screen_data_path="./result/example-analyzed_data/2_base_PRAUC0.8_scRNA_seq_data_with_thv70_and_111_neuropal_reporter.csv"
    fneurons_data_path="./result/example-analyzed_data/confirm_Nmin_between_NeuroPALs/neuron_sets_and_recombination-union/example-Final_neuron_sets-${Nmin}-0.5-1-Both-False.xlsx"
    data_dir="./result/example-gene_screening/with_207_scrnaseq_genes/N${Nmin}-AtLeast_Gene${at_least_gene_n}"
    pool_n=2
    echo $Nmin
    echo " 1. To generate_gene_combinations"
    sub_pool_n=${pool_n}
    for pool_i in `eval echo {1..${sub_pool_n}}`
    do
        echo ${pool_i}
        nohup python generate_gene_combinations.py ${main_dir} ${screen_data_path} ${fneurons_data_path} ${data_dir} ${at_least_gene_n} ${pool_n} ${pool_i} > ./nohup_log/nohup_generate_gene_combinations_${pool_i}.log 2>&1 &
    done
    wait

    echo " 2. To summary unfold statistics files"
        nohup python sum_statistics_files.py ${main_dir} ${data_dir} ${sub_pool_n} > ./nohup_log/nohup_files_summary.log 2>&1 &
    wait

    echo " 3. To execute minimal union by logVersion"
    max_gene_combn=39
    file0_first_gene_combn=2
    file1_first_gene_combn=2
    file0_max_readrow_n=50000
    file1_max_readrow_n=10000
    nohup python minimal_union_logVersion.py ${main_dir} ${screen_data_path} ${data_dir} ${file0_first_gene_combn} ${file1_first_gene_combn} ${file0_max_readrow_n} ${file1_max_readrow_n} ${max_gene_combn} > ./nohup_log/nohup_minimal_union.log 2>&1 &
    wait
done




