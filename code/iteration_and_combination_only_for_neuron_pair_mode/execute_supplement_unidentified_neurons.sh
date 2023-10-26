#!/usr/bin/env bash

for Nmin in 6
do
    main_dir="../../"
    screen_data_path="./result/example-analyzed_data/2_base_PRAUC0.8_scRNA_seq_data_with_thv70_and_111_neuropal_reporter.csv"
    data_dir="./code/iteration_and_combination_only_for_neuron_pair_mode/example-supplement_unidentified_neurons"
    fneurons_data_path="${data_dir}/unidentified_neuron_sets.xlsx"
    at_least_gene_n=2
    pool_n=1
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

    echo " 3. To supplement minimal union by logVersion"
    incmp_gene_comb_data_dir="${data_dir}/incomplete_identified_gene_combinations"
    max_gene_combn=39
    file0_first_gene_combn=2
    file1_first_gene_combn=2
    file0_max_readrow_n=50000
    file1_max_readrow_n=10000
    start_idn=36
    nohup python minimal_union_logVersion_supplemented.py ${main_dir} ${incmp_gene_comb_data_dir} ${screen_data_path} ${data_dir} ${file0_first_gene_combn} ${file1_first_gene_combn} ${file0_max_readrow_n} ${file1_max_readrow_n} ${start_idn} ${max_gene_combn} > ./nohup_log/nohup_minimal_union_supplemented.log 2>&1 &
    wait
done




