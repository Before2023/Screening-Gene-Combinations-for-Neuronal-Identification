import os
import time
import pandas as pd
import numpy as np
from itertools import combinations
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


# 归类：在 neurons 中相同基因表达模式的基因集
def classify_genes(data0, neurons):
    data = data0[['gene_id'] + neurons]
    data.sort_values(by=neurons, ascending=False, inplace=True)
    data.reset_index(drop=True, inplace=True)
    neuron_n = len(neurons)
    rep_gene_expr_mode = data[neurons].loc[0, :]
    rep_expr_gene = data.loc[0, 'gene_id']
    rep_gene_expr_info = {}
    for index, row in data.iterrows():
        if index == 0:
            pass
        else:
            gene_expr_mode = data[neurons].loc[index, :]
            expr_gene = data.loc[index, 'gene_id']
            if rep_expr_gene not in rep_gene_expr_info:
                if (sum(rep_gene_expr_mode) == neuron_n) | (sum(rep_gene_expr_mode) == 0):
                    rep_gene_expr_mode = gene_expr_mode
                    rep_expr_gene = expr_gene
                    continue
                else:
                    rep_gene_expr_info[rep_expr_gene] = [rep_expr_gene]
            else:
                pass
            if rep_gene_expr_mode.equals(gene_expr_mode):
                rep_gene_expr_info[rep_expr_gene].append(expr_gene)
            else:
                rep_gene_expr_mode = gene_expr_mode
                rep_expr_gene = expr_gene
                if (sum(rep_gene_expr_mode) == neuron_n) | (sum(rep_gene_expr_mode) == 0):
                    pass
                else:
                    rep_gene_expr_info[rep_expr_gene] = [rep_expr_gene]
    rep_data = data.loc[data['gene_id'].isin(list(rep_gene_expr_info.keys())), :]
    rep_data.sort_values(by='gene_id', inplace=True)
    rep_data.reset_index(drop=True, inplace=True)
    return rep_data, rep_gene_expr_info


# 归类：在 neurons 中相同基因表达模式的基因集
def classify_genes2(data0, neurons):
    data = data0[neurons]
    data['sum'] = data.sum(axis=1)
    neuron_n = len(neurons)
    keep_index = data.loc[(data['sum'] != 0) & (data['sum'] != neuron_n), :].index.tolist()
    rep_data = data0.loc[keep_index, ['gene_id'] + neurons]
    rep_data.reset_index(drop=True, inplace=True)
    return rep_data


# gene screening
def split_index_for_multiprocess(all_index, pool_n, pool_i):
    if len(all_index) % pool_n == 0:
        splitn = int(len(all_index) / pool_n)
    else:
        splitn = int(len(all_index) / pool_n) + 1
    index_belong_i = all_index[(pool_i - 1) * splitn:pool_i * splitn]
    return index_belong_i


# seq data --gene id
def gene_dict(seq_data_nbase):
    import copy
    a = copy.deepcopy(seq_data_nbase)
    a.sort_values(by='gene_id', ascending=True, inplace=True)
    a.reset_index(drop=True, inplace=True)
    a['id'] = a.index
    a.index = a['gene_id']
    gene_code_dict = dict(a['id'])
    return gene_code_dict


# Only for neuron pair mode
# to generate gene combinatins for each central neurons
def main_generate_gene_combinations2(fneurons_data, screen_data, at_least_gene_n, save_dir_format,
                                     pool_n, pool_i):
    # save gene_code
    gene_code_dict = gene_dict(screen_data)
    gene_code_data = screen_data[['gene_id']]
    gene_code_data['code'] = gene_code_data['gene_id'].apply(lambda x: gene_code_dict[x])
    gene_code_data.to_excel("/".join(save_dir_format.split("/")[:-1]) + "/gene_code_correspondence.xlsx", index=False)
    # check: number of each T.B. neurons
    fneurons_data['Num. T.B. neurons'] = fneurons_data['T.B. neurons'].apply(lambda x: len(eval(x)))
    if fneurons_data['Num. T.B. neurons'].value_counts()[2] == fneurons_data.shape[0]:  # neuron pair mode
        # print("GOOD: Neuron Pair Mode.")
        pass
    else:  # sort fneuron_data to neuron pair mode
        neuron_pair_list = []
        for i, row in fneurons_data.iterrows():
            neurons = eval(row['T.B. neurons'])
            for neuron_pair in combinations(neurons, 2):
                neuron_pair = list(neuron_pair)
                neuron_pair.sort(reverse=False)
                neuron_pair_list.append(neuron_pair)
        neuron_pair_list = list(set(neuron_pair_list))
        fneurons_dict = {'T.B. neurons': neuron_pair_list, 'Num. T.B. neurons': [2]*len(neuron_pair_list)}
        fneurons_data = pd.DataFrame(fneurons_dict)
        fneurons_data['order'] = fneurons_data.index
    # get part pool_i
    all_index = fneurons_data['order'].tolist()
    iterindex = split_index_for_multiprocess(all_index, pool_n, pool_i)
    stat_dict = {"data name": [],
                 "gene combination info.": [],
                 "Num. of total gene combination": []}
    for i in iterindex:
        neurons = eval(fneurons_data.loc[i, 'T.B. neurons'])
        central_neurons = "%s_central_neurons" % i
        # gene classification
        rep_screen_data_nbase = classify_genes2(screen_data, neurons)
        rep_screen_data_nbase['code'] = rep_screen_data_nbase['gene_id'].apply(lambda x: gene_code_dict[x])
        gene_code_list = rep_screen_data_nbase['code'].tolist()
        # write
        save_path = save_dir_format + "/gene_combination_for_%s_%s.log" % (central_neurons, at_least_gene_n)
        is_Exist_file(save_path)
        if len(gene_code_list) >= at_least_gene_n:
            gene_comb_info = {at_least_gene_n: 0}
            with open(save_path, 'a') as a:
                for gene_tuple in combinations(gene_code_list, at_least_gene_n):
                    temp_gene_list = list(gene_tuple)
                    temp_gene_list.sort(reverse=False)
                    a.write(('%s\t' * (len(temp_gene_list)-1) + '%s\n') % tuple(temp_gene_list))
                    gene_comb_info[at_least_gene_n] += 1
            data_name = "union_represent_gene_combinations%s_%s.log" % (central_neurons, at_least_gene_n)
            stat_dict['data name'].append(data_name)
            stat_dict['gene combination info.'].append(gene_comb_info)
            stat_dict['Num. of total gene combination'].append(gene_comb_info[at_least_gene_n])
        else:
            pass
    # save stat
    if len(stat_dict) != 0:
        save_stat_path = save_dir_format + "/0statistics_neuron_sets_gene_combinations_pool_%s.xlsx" % pool_i
        is_Exist_file(save_stat_path)
        stat_data = pd.DataFrame(stat_dict)
        stat_data.to_excel(save_stat_path, index=False)
    else:
        pass


if __name__ == '__main__':
    import sys
    print('-----------------------------------------------------------------------')
    print('**************** execute: gene combinatins ********************')
    main_dir, screen_data_path, fneurons_data_path, data_dir, at_least_gene_n, pool_n, pool_i = sys.argv[1:]

    # data path & parameters
    os.chdir(main_dir)
    at_least_gene_n = int(at_least_gene_n)
    pool_n = int(pool_n)
    pool_i = int(pool_i)
    mkdir(data_dir)
    time.sleep(pool_i)

    # read
    screen_data = pd.read_csv(screen_data_path)
    fneurons_data = pd.read_excel(fneurons_data_path)
    # 固定参数
    save_unfold_gene_comb_dir = data_dir + "/unfold_gene_combinations"
    mkdir(save_unfold_gene_comb_dir)
    main_generate_gene_combinations2(fneurons_data, screen_data, at_least_gene_n, save_unfold_gene_comb_dir,
                                    pool_n, pool_i)


