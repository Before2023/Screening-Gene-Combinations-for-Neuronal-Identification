from minimal_union_logVersion import *
import warnings
warnings.filterwarnings("ignore")


def start_union_data_collect(incmp_gene_comb_data_dir, save_minimal_union_dir):
    data_path_list = walk_files(incmp_gene_comb_data_dir)
    for data_path in data_path_list:
        if "idn" in data_path.split("/")[-1]:
            to_path = save_minimal_union_dir + "/%s" % data_path.split("/")[-1]
            val = os.system("cp %s %s" % (data_path, to_path))


# union gene combinations of central neurons
def union_gene_combinations_central_neurons_complement(iterindex, gene_n, read_format_path,
                                                       save_union_path_format, save_stat_data_path,
                                                       file0_first_gene_combn=2, file0_max_readrow_n=np.inf,
                                                       file1_first_gene_combn=2, file1_max_readrow_n=100000,
                                                       single_csr_max_rown=10000, max_gene_combn=39, pool_n=1,
                                                       start_idn=3970):
    stat_dict = {'central neuron id': [],
                 'gene combination info. of file1': [],
                 'first gene combn. selected of file1': [],
                 'total combn. selected of file1': [],
                 'total combn. of file1': [],

                 'gene combination info. of file2': [],
                 'first gene combn. selected of file2': [],
                 'total combn. selected of file2': [],
                 'total combn. of file2': [],

                 'union file code': [],
                 'gene combination info. of union file': [],
                 'total combn. of union file': []}
    union_file_code = ''
    idn = start_idn
    for central_neuron_id in iterindex:
        central_neuron = "%s_central_neurons" % central_neuron_id
        # if central_neuron_id == 3167:
        idn += 1
        print("\n%s -- %s"%(idn, central_neuron))
        csr0_dict = get_csr_matrix_dict(central_neuron, read_format_path,
                                        max_gene_combn, gene_n, file0_first_gene_combn,
                                        single_csr_max_rown=single_csr_max_rown,
                                        max_readrow_n=file0_max_readrow_n)
        csr1_dict = get_csr_matrix_dict(idn-1, save_union_path_format,
                                        max_gene_combn, gene_n, file1_first_gene_combn,
                                        single_csr_max_rown=single_csr_max_rown,
                                        max_readrow_n=file1_max_readrow_n)
        gen_combn_info1 = statistics_gene_comb_info_from_file(idn-1, save_union_path_format,
                                                              max_gene_combn, file1_first_gene_combn,
                                                              save_union_path_format)
        print("len(csr0_dict):", len(csr0_dict), "len(csr1_dict):", len(csr1_dict))
        main_two_neuron_sets_union2(csr0_dict, csr1_dict, save_union_path_format, max_gene_combn,
                                    file1_first_gene_combn, pool_n, idn)
        # statistics
        union_file_code = union_file_code + "-%s" % (str(central_neuron_id))
        gen_combn_info0 = statistics_gene_comb_info_from_file(central_neuron, read_format_path,
                                                              max_gene_combn, file0_first_gene_combn,
                                                              save_union_path_format)
        gen_combn_info = statistics_gene_comb_info_from_file(idn, save_union_path_format,
                                                              max_gene_combn, file1_first_gene_combn,
                                                             save_union_path_format)
        stat_dict['central neuron id'].append(central_neuron_id)
        stat_dict['gene combination info. of file1'].append(gen_combn_info0)
        stat_dict['first gene combn. selected of file1'].append(file0_first_gene_combn)
        stat_dict['total combn. selected of file1'].append(sum([v.shape[0] for k, v in csr0_dict.items()]))
        stat_dict['total combn. of file1'].append(sum(list(gen_combn_info0.values())))

        stat_dict['gene combination info. of file2'].append(gen_combn_info1)
        stat_dict['first gene combn. selected of file2'].append(file1_first_gene_combn)
        stat_dict['total combn. selected of file2'].append(sum([v.shape[0] for k, v in csr1_dict.items()]))
        stat_dict['total combn. of file2'].append(sum(list(gen_combn_info1.values())))

        stat_dict['union file code'].append(union_file_code)
        stat_dict['gene combination info. of union file'].append(gen_combn_info)
        stat_dict['total combn. of union file'].append(sum(list(gen_combn_info.values())))
        # break
        # if idn == 3791:
        #     break
    # save statistics
    stat_data = pd.DataFrame(stat_dict)
    stat_data.to_excel(save_stat_data_path, index=False)


if __name__ == '__main__':
    import sys
    import os

    (main_dir, incmp_gene_comb_data_dir, screen_data_path, data_dir, file0_first_gene_combn, file1_first_gene_combn,
     file0_max_readrow_n, file1_max_readrow_n, start_idn, max_gene_combn) = sys.argv[1:]
    print('-----------------------------------------------------------------------')
    print('**************** execute: union of gene combinatins ********************')

    os.chdir(main_dir)
    max_gene_combn = int(max_gene_combn)
    file0_first_gene_combn = int(file0_first_gene_combn)
    file1_first_gene_combn = int(file1_first_gene_combn)
    start_idn = int(start_idn)
    file0_max_readrow_n = int(file0_max_readrow_n)
    file1_max_readrow_n = int(file1_max_readrow_n)
    single_csr_max_rown = 1000

    # 可选参数
    screen_data = pd.read_csv(screen_data_path)
    gene_n = screen_data.shape[0]
    # 可变参数
    save_minimal_union_dir = data_dir + "/minimal_union_%s_%s-%s" % (file0_first_gene_combn, file1_first_gene_combn, file1_max_readrow_n)
    mkdir(save_minimal_union_dir)
    start_union_data_collect(incmp_gene_comb_data_dir, save_minimal_union_dir)  # supplemented file1
    stat_data_path = data_dir + "/unfold_gene_combinations/unfold_statistics_files_summary.xlsx"  #
    stat_data = pd.read_excel(stat_data_path)
    iterindex = stat_data['central neuron id'].tolist()  # supplemented index
    read_unfold_gene_comb_format_path = data_dir + "/unfold_gene_combinations/gene_combination_for_%s_%s.log"  # supplemented file2
    save_union_path_format = save_minimal_union_dir + "/union_gene_combinations_for_idn%s_combn%s.log"
    save_stat_data_path = save_minimal_union_dir + "/statistics_union_gene_combinations.xlsx"
    union_gene_combinations_central_neurons_complement(iterindex, gene_n, read_unfold_gene_comb_format_path,
                                                       save_union_path_format, save_stat_data_path,
                                                       file0_first_gene_combn, file0_max_readrow_n,
                                                       file1_first_gene_combn, file1_max_readrow_n,
                                                       single_csr_max_rown,
                                                       max_gene_combn, 3, start_idn)


