import pandas as pd
import warnings
warnings.filterwarnings("ignore")


def is_Exist_file(path):
    import os
    if os.path.exists(path):
        os.remove(path)


def mkdir(path):
    import os
    path = path.strip()
    path = path.rstrip("\\")
    isExists = os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        print(path + ' 创建成功')


def walk_files(data_dir):
    import os
    g = os.walk(data_dir)
    data_path_list = []
    for path, dir_list, file_list in g:
        for file_name in file_list:
            data_path_list.append(os.path.join(path, file_name))
    return data_path_list


def gene_dict(seq_data_nbase):
    import copy
    a = copy.deepcopy(seq_data_nbase)
    a.sort_values(by='gene_id', ascending=True, inplace=True)
    a.reset_index(drop=True, inplace=True)
    a['id'] = a.index
    a.index = a['gene_id']
    gene_code_dict = dict(a['id'])
    return gene_code_dict


def reverse_gene_code_dict(gene_code_dict):
    rev_gene_code_dict = {}
    for k, v in gene_code_dict.items():
        rev_gene_code_dict[v] = k
    return rev_gene_code_dict


def dataframe_columns_to_list(data):
    gene_combn = data.shape[1]
    data_dict = {"gene combination": [],
                 "gene number": []}
    data = data.astype(int)
    for i in data.index:
        values = data.loc[i, :].tolist()
        values.sort(reverse=False)
        data_dict['gene combination'].append(str(values))
        data_dict['gene number'].append(gene_combn)
    data = pd.DataFrame(data_dict)
    return data


# check genes from these 39 genes in fneurons data
def get_gene_combination_screened(raw_screen_data, gene_screened_data,
                                  gene_combn_index="all", code_transform=True):
    gene_code_dict = gene_dict(raw_screen_data)
    gene_combs_screened = []
    if code_transform:
        rev_gene_code_dict = reverse_gene_code_dict(gene_code_dict)
        if gene_combn_index == 'all':
            for i, row in gene_screened_data.iterrows():
                gene_comb = [rev_gene_code_dict[x] for x in eval(row['gene combination'])]
                gene_combs_screened.append(gene_comb)
        else:
            gene_combs_screened.append([rev_gene_code_dict[x] for x in eval(gene_screened_data.loc[gene_combn_index, :]['gene combination'])])
    else:
        code_gene_dict = {}
        for gene, code in gene_code_dict.items():
            code_gene_dict[code] = gene
        if gene_combn_index == 'all':
            for i, row in gene_screened_data.iterrows():
                gene_comb = [code_gene_dict[x] for x in eval(row['gene combination'])]
                gene_combs_screened.append(gene_comb)
        else:
            gene_combs_screened.append([code_gene_dict[x] for x in eval(gene_screened_data.loc[gene_combn_index, :]['gene combination'])])
    return gene_combs_screened


def evaluate_gene_list_screened(fneurons, screen_data):
    stat_dict = {'T.B. neurons': [],
                 'Num. T.B. neurons': [],
                 'Num. neurons to be discernable': []
                }
    for i in fneurons.index:
        neurons = eval(fneurons.loc[i, 'T.B. neurons'])
        neruon_gene_mat = screen_data[neurons]
        stat_dict['T.B. neurons'].append(neurons)
        stat_dict['Num. T.B. neurons'].append(len(neurons))
        stat_dict['Num. neurons to be discernable'].append(neruon_gene_mat.T.drop_duplicates().shape[0])
    stat_data = pd.DataFrame(stat_dict)
    stat_data['abs_diff'] = stat_data.apply(lambda row: abs(row['Num. T.B. neurons'] - row['Num. neurons to be discernable']), axis=1)
    gene_combination = screen_data['gene_id'].tolist()
    stat_data['gene_combination_screened'] = str(gene_combination)
    stat_data['num. of gene combination'] = len(gene_combination)
    # print("before merging: stat_data.shape", stat_data.shape)
    stat_data['T.B. neurons'] = stat_data['T.B. neurons'].apply(lambda x: str(x))
    stat_data = pd.merge(fneurons, stat_data, how='left', on=['T.B. neurons', 'Num. T.B. neurons'])
    # print("after merging: stat_data.shape", stat_data.shape)
    #
    del_cols = ['Num. represent genes', 'Central neuron index',
            'etime (s)', 'etime (h)', 'random_index']
    for dcol in del_cols:
        if dcol in stat_data.columns:
            del stat_data[dcol]
        else:
            pass
    return stat_data


# 找出无法区分的神经元对
def find_undiscernable_neuron_pair(stat_data, screen_data):
    unable_dicern_neuron_pairs = []
    for i in stat_data.loc[stat_data['abs_diff']!=0, :].index:
        neurons = eval(stat_data.loc[i, 'T.B. neurons'])
        neruon_gene_mat = screen_data[neurons]
        a = neruon_gene_mat.T
        if (a[a.duplicated(keep=False)].shape[0] != 0) & (str(tuple(a[a.duplicated(keep=False)].index)) not in unable_dicern_neuron_pairs):
            undis_neuron_pair = list(a[a.duplicated(keep=False)].index)
            undis_neuron_pair.sort(reverse=False)
            unable_dicern_neuron_pairs.append(str(tuple(undis_neuron_pair)))
    return unable_dicern_neuron_pairs


# evaluation tool
def evaluate_one_gene_combinatin_discernable(fneurons, screen_data):
    from itertools import combinations
    total_neuron_pairs = []
    unable_dicern_neuron_pairs = []
    unable_dicern_neuron_pairs_dict = {}
    for i in fneurons.index:
        neurons = eval(fneurons.loc[i, 'T.B. neurons'])
        neruon_gene_matT = screen_data[neurons].T
        if neruon_gene_matT[neruon_gene_matT.duplicated(keep=False)].shape[0] != 0:
            for neuron_pair in combinations(neurons, 2):
                neuron_pair = list(neuron_pair)
                neuron_pair.sort(reverse=False)
                temp_matT = neruon_gene_matT.loc[neuron_pair, :]
                if temp_matT.drop_duplicates().shape[0] != temp_matT.shape[0]:
                    unable_dicern_neuron_pairs.append(str(neuron_pair))
                    if str(neurons) not in unable_dicern_neuron_pairs_dict:
                        unable_dicern_neuron_pairs_dict[str(neurons)] = [str(neuron_pair)]
                    else:
                        unable_dicern_neuron_pairs_dict[str(neurons)].append(str(neuron_pair))
                else:
                    pass
                total_neuron_pairs.append(str(neuron_pair))
        else:
            for neuron_pair in combinations(neurons, 2):
                neuron_pair = list(neuron_pair)
                neuron_pair.sort(reverse=False)
                total_neuron_pairs.append(str(neuron_pair))
    unable_dicern_neuron_pairs = list(set(unable_dicern_neuron_pairs))
    total_neuron_pairs = list(set(total_neuron_pairs))
    # statistics
    eval_dict = {'gene combination': [screen_data['gene_id'].tolist()],
                 'num. of gene combination': [len(screen_data['gene_id'].tolist())],
                 'num. of total neuron pairs': [len(total_neuron_pairs)],
                 'num. of undiscernable neuron pairs': [len(unable_dicern_neuron_pairs)],
                 'ratio of dicrenable neuron pairs': [1- round(len(unable_dicern_neuron_pairs)/len(total_neuron_pairs), 3)],
                 'undiscernable neuron pairs': [unable_dicern_neuron_pairs],
                 'original': [unable_dicern_neuron_pairs_dict]
                }
    eval_data = pd.DataFrame(eval_dict)
    return eval_data, unable_dicern_neuron_pairs


def get_gene_screened_39reporter_data_with_complement(gene_comb_screened_mode):
    gene_screened_data_dict = {}
    for mode, data_dir_list in gene_comb_screened_mode.items():
        data_list = []
        for one_data_dir in data_dir_list:
            data_path_list = walk_files(one_data_dir)
            for data_path in data_path_list:
                if ("idn3790" not in data_path) and ("idn" in data_path):
                    data = pd.read_csv(data_path, sep='\t', header=None)
                    data = dataframe_columns_to_list(data)
                    data_list.append(data)
                else:
                    pass
        data = pd.concat(data_list, axis=0)
        data.drop_duplicates(inplace=True)
        data.reset_index(drop=True, inplace=True)
        gene_screened_data_dict[mode] = data
    return gene_screened_data_dict


def count_total_neuron_pairs(neuron_data):
    from itertools import combinations
    total_neuron_pairs = []
    for i in neuron_data.index:
        neurons = eval(neuron_data.loc[i, 'T.B. neurons'])
        for neuron_pair in combinations(neurons, 2):
            neuron_pair = list(neuron_pair)
            neuron_pair.sort(reverse=False)
            total_neuron_pairs.append(str(neuron_pair))
    total_neuron_pairs = list(set(total_neuron_pairs))
    return len(total_neuron_pairs)


def count_neuron_pairs_undiscernable(stat_data, screen_data):
    from itertools import combinations
    unable_dicern_neuron_pairs = []
    for i in stat_data.loc[stat_data['abs_diff']!=0, :].index:
        neurons = eval(stat_data.loc[i, 'T.B. neurons'])
        neruon_gene_matT = screen_data[neurons].T
        for neuron_pair in combinations(list(neruon_gene_matT[neruon_gene_matT.duplicated(keep=False)].index), 2):
            neuron_pair = list(neuron_pair)
            neuron_pair.sort(reverse=False)
            temp_matT = neruon_gene_matT.loc[neuron_pair, :]
            if temp_matT.drop_duplicates().shape[0] != temp_matT.shape[0]:
                unable_dicern_neuron_pairs.append(str(neuron_pair))
            else:
                pass
    unable_dicern_neuron_pairs = list(set(unable_dicern_neuron_pairs))
    return unable_dicern_neuron_pairs


# #evaluation for gene list screened -- final neuron data
# find all duplicated neurons in one neuron sets
def split_index_for_multiprocess(all_index, pool_n, pool_i):
    if len(all_index) % pool_n == 0:
        splitn = int(len(all_index) / pool_n)
    else:
        splitn = int(len(all_index) / pool_n) + 1
    index_belong_i = all_index[(pool_i - 1) * splitn:pool_i * splitn]
    return index_belong_i


def get_gene_combinations_screened(gene_screened_data_dir_dict, gene_screened_mode, screen_data_path):
    raw_screen_data = pd.read_csv(screen_data_path)
    gene_screened_data_dict = get_gene_screened_39reporter_data_with_complement(gene_screened_data_dir_dict)
    gene_screened_data = gene_screened_data_dict[gene_screened_mode]
    gene_combs_screened = get_gene_combination_screened(raw_screen_data, gene_screened_data,
                                                        gene_combn_index="all", code_transform=False)
    return gene_combs_screened


def main_evaluation_gene_combinations(neuron_data_path, gene_screened_data_path, screen_data_path, save_dir):
    raw_screen_data = pd.read_csv(screen_data_path)
    gene_screened_data = pd.read_csv(gene_screened_data_path, sep='\t', header=None)
    gene_screened_data = dataframe_columns_to_list(gene_screened_data)
    gene_combs_screened = get_gene_combination_screened(raw_screen_data, gene_screened_data,
                                                        gene_combn_index="all", code_transform=False)
    # read
    eval_data_dict = {"neuron data": [],
                      "gene combination": [],
                      "gene number": [],
                      "gene screened mode": [],
                      "T.B. neurons from one neuron set undiscernable": [],
                      "number of neuron sets undiscernable": [],
                      "Neuron belong from": [],
                      'neuron pairs undiscernable': [],
                      "rato of neuron pairs discernable": [],
                      "number of neuron pairs undiscernable": [],
                      "number of total neuron pairs": []
                     }
    for one_gene_combination in gene_combs_screened:
        neuron_data = pd.read_excel(neuron_data_path)
        total_neuron_pn = count_total_neuron_pairs(neuron_data)
        screen_data = raw_screen_data.loc[raw_screen_data['gene_id'].isin(one_gene_combination), :]
        stat_data = evaluate_gene_list_screened(neuron_data, screen_data)
        unable_dicern_neuron_sets = find_undiscernable_neuron_pair(stat_data, screen_data)
        unable_dicern_neuron_pairs = count_neuron_pairs_undiscernable(stat_data, screen_data)
        eval_data_dict['neuron data'].append(neuron_data_path.split('/')[-1])
        eval_data_dict['gene combination'].append(one_gene_combination)
        eval_data_dict['gene number'].append(len(one_gene_combination))
        eval_data_dict['gene screened mode'].append(gene_screened_data_path.split('/')[-1])
        eval_data_dict['T.B. neurons from one neuron set undiscernable'].append(unable_dicern_neuron_sets)
        eval_data_dict['number of neuron sets undiscernable'].append(len(unable_dicern_neuron_sets))
        eval_data_dict['Neuron belong from'].append(stat_data.loc[stat_data['abs_diff'] != 0, :]['Neuron belong'].tolist())
        eval_data_dict['neuron pairs undiscernable'].append(unable_dicern_neuron_pairs)
        eval_data_dict['rato of neuron pairs discernable'].append(round(1-(len(unable_dicern_neuron_pairs)/total_neuron_pn), 3))
        eval_data_dict['number of neuron pairs undiscernable'].append(len(unable_dicern_neuron_pairs))
        eval_data_dict['number of total neuron pairs'].append(total_neuron_pn)
    eval_data = pd.DataFrame(eval_data_dict)
    mkdir(save_dir)
    save_eval_data_path = save_dir + "/evaluation_results_for_gene_combinations.xlsx"
    eval_data.to_excel(save_eval_data_path, index=False)
    return eval_data


if __name__ == '__main__':
    import sys
    import os

    print('-----------------------------------------------------------------------')
    print('**************** execute: gene combinatins ********************')
    main_dir, gene_screened_data_path, neuron_data_path, screen_data_path = sys.argv[1:]

    os.chdir(main_dir)
    gene_screened_data_dir = "/".join(gene_screened_data_path.split('/')[:-1])
    save_dir = gene_screened_data_dir + "/evaluation_gene_combinations"
    eval_data = main_evaluation_gene_combinations(neuron_data_path, gene_screened_data_path,
                                                  screen_data_path, save_dir)

