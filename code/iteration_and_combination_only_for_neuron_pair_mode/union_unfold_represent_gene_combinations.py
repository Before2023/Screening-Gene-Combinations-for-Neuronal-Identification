from scipy.sparse import vstack
from scipy.sparse import csr_matrix
import time
from generate_gene_combinations import *
import warnings
warnings.filterwarnings("ignore")


# step 1: deduplication of gene combinations including lower same gene combinations
# step 2: unfold genes with the same expression of the represent gene among the neuron set
# comb_data0: lower gene combination number
# comb_data1: higner gene combination number
# comb0_n: gene combination number of comb_data0
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


def combination_deduplicate(comb_data0, comb_data1):
    # matrix multiplication
    matrix_one, matrix_two = np.array(comb_data0.T), np.array(comb_data1)
    mdata0 = np.dot(matrix_two, matrix_one)
    mdata = pd.DataFrame(mdata0 - matrix_one.sum(axis=0))
    # print(mdata.shape)
    mdata.replace(0, np.nan, inplace=True)
    dhcomb_index = mdata[mdata.isnull().T.any()].index.tolist()
    khcomb_index = list(set(mdata.index.tolist()) - set(dhcomb_index))
    kcomb_data = comb_data1.loc[khcomb_index, :]
    # comb_data = pd.concat([comb_data0, comb_data1.loc[khcomb_index, :]], axis=0)
    kcomb_data.drop_duplicates(inplace=True)
    kcomb_data.reset_index(drop=True, inplace=True)  # shape: [combinatin num, gene list]
    return kcomb_data


def combination_deduplicate_for_sparse_matrix(csr0, csr1):
    # matrix multiplication
    csr_multiply = csr1.dot(csr0.T)
    mcsr = csr_multiply - vstack([csr_matrix(csr0.T.sum(axis=0))] * csr_multiply.shape[0])
    row_indexes, col_indexes = mcsr.nonzero()
    row_series = pd.Series(row_indexes).value_counts()
    keep_indexes = row_series[row_series == csr0.shape[0]]  # 某行非零个数等于 csr0 的行数（组合数）
    csr1_keep = []
    for i in keep_indexes.index:
        csr1_keep.append(csr1[i])
    if len(csr1_keep) != 0:
        kcsr1 = vstack(csr1_keep)
    else:
        kcsr1 = csr_matrix((0, 0), dtype=np.int8)
    return kcsr1


# colindex to csr_matrix
def one_colindex_tocsr(row, k, gene_n):
    colindex = []
    for t in range(k):
        colindex.append(row[t])
    indices = np.array(colindex)
    indptr = np.array([0, len(indices)])
    data = np.array([1]*len(indices))
    csr = csr_matrix((data, indices, indptr), shape=(1, gene_n))
    return csr


def colindex_to_csr_matrix(data_path, gene_n, start_row_n=0, max_readrow_n=1000000):
    array_dict = {}
    if max_readrow_n < np.inf:
        c = 0
        with open(data_path, 'r') as f:
            for line in f:
                if (c >= start_row_n) & (c < (start_row_n + max_readrow_n)):
                    line_s = line.strip('\n').split('\t')
                    if len(line_s) not in array_dict:
                        array_dict[len(line_s)] = [line_s]
                    else:
                        array_dict[len(line_s)].append(line_s)
                else:
                    pass
                c += 1
    else:
        with open(data_path, 'r') as f:
            for line in f:
                line_s = line.strip('\n').split('\t')
                if len(line_s) not in array_dict:
                    array_dict[len(line_s)] = [line_s]
                else:
                    array_dict[len(line_s)].append(line_s)
    # read
    all_csr_list = []
    for combn, array_list in array_dict.items():
        data = pd.DataFrame(array_list).astype(int)
        csr_list = data.apply(lambda row: one_colindex_tocsr(row, combn, gene_n), axis=1).tolist()
        del data
        all_csr_list += csr_list
    if len(all_csr_list) != 0:
        csr = vstack(all_csr_list)
    else:
        csr = csr_matrix((0, 0), dtype=np.int8)
    return csr


def count_file_rown(data_path):
    temp_path = data_path + '-file_count'
    is_Exist_file(temp_path)
    val = os.system("wc -l %s > %s" % (data_path, temp_path))
    with open(temp_path, 'r') as f:
        for line in f:
            file_count = int(line.strip(' ').split(' ')[0])
            break
    is_Exist_file(temp_path)
    return file_count


def gene_combinations_union(data_path1, csr0, save_path, gene_n, max_readrow_n=100000):
    csr1 = colindex_to_csr_matrix(data_path1, gene_n, start_row_n=0, max_readrow_n=np.inf)
    partn1 = int(csr1.shape[0]/max_readrow_n) + 1
    partn0 = int(csr0.shape[0]/max_readrow_n) + 1
    kacsr1_list = []
    for i in range(partn1):
        kcsr1 = csr1[i*max_readrow_n:(i+1)*max_readrow_n]
        for j in range(partn0):
            csr0_chunk = csr0[j*max_readrow_n:(j+1)*max_readrow_n]
            if (kcsr1.shape[0] != 0) & (csr0_chunk.shape[0] != 0):
                kcsr1 = combination_deduplicate_for_sparse_matrix(csr0_chunk, kcsr1)
            else:
                pass
                # break
        if kcsr1.shape[0] != 0:
            kacsr1_list.append(kcsr1)
    if len(kacsr1_list) != 0:
        kacsr1 = vstack(kacsr1_list)
        with open(save_path, 'a') as a:
            for index in range(kacsr1.shape[0]):
                col_index = kacsr1[index].nonzero()[1].tolist()
                col_index.sort(reverse=False)
                a.write(("%s\t" * (len(col_index) - 1) + '%s\n') % tuple(col_index))
        csr0 = vstack([csr0, kacsr1])
    else:
        pass
    return csr0


def gene_combinations_union3(csr1, csr0, save_path, max_readrow_n=100000, write_mode=True):
    partn1 = int(csr1.shape[0]/max_readrow_n) + 1
    partn0 = int(csr0.shape[0]/max_readrow_n) + 1
    kacsr1_list = []
    for i in range(partn1):
        kcsr1 = csr1[i*max_readrow_n:(i+1)*max_readrow_n]
        for j in range(partn0):
            csr0_chunk = csr0[j*max_readrow_n:(j+1)*max_readrow_n]
            if (kcsr1.shape[0] != 0) & (csr0_chunk.shape[0] != 0):
                kcsr1 = combination_deduplicate_for_sparse_matrix(csr0_chunk, kcsr1)
            else:
                pass
                # break
        if kcsr1.shape[0] != 0:
            kacsr1_list.append(kcsr1)
    if len(kacsr1_list) != 0:
        kacsr1 = vstack(kacsr1_list)
        csr0 = vstack([csr0, kacsr1])
        if write_mode:
            with open(save_path, 'a') as a:
                for index in range(kacsr1.shape[0]):
                    col_index = kacsr1[index].nonzero()[1].tolist()
                    col_index.sort(reverse=False)
                    a.write(("%s\t" * (len(col_index) - 1) + '%s\n') % tuple(col_index))
    else:
        kacsr1 = csr0[:0]
    return csr0, kacsr1


# for multiprocessing
def gene_combinations_union2(parameter_tuple):
    (csr, csr0, save_path, max_readrow_n) = parameter_tuple
    partn1 = int(csr.shape[0]/max_readrow_n) + 1
    partn0 = int(csr0.shape[0]/max_readrow_n) + 1
    kacsr1_list = []
    for i in range(partn1):
        kcsr1 = csr[i*max_readrow_n:(i+1)*max_readrow_n]
        for j in range(partn0):
            csr0_chunk = csr0[j*max_readrow_n:(j+1)*max_readrow_n]
            if (kcsr1.shape[0] != 0) & (csr0_chunk.shape[0] != 0):
                kcsr1 = combination_deduplicate_for_sparse_matrix(csr0_chunk, kcsr1)
            else:
                break
        if kcsr1.shape[0] != 0:
            kacsr1_list.append(kcsr1)
    if len(kacsr1_list) != 0:
        kacsr1 = vstack(kacsr1_list)
        with open(save_path, 'a') as a:
            for index in range(kacsr1.shape[0]):
                col_index = kacsr1[index].nonzero()[1].tolist()
                col_index.sort(reverse=False)
                a.write(("%s\t" * (len(col_index) - 1) + '%s\n') % tuple(col_index))


def union_gene_combinations_for_a_central_neuron(central_neurons, rep_comb_data_path_format,
                                                 save_rep_comb_data_path_format, gene_n,
                                                 max_comb_n=9, max_readrow_n=100000):
    # print("union_gene_combinations_for_a_central_neuron")
    s = 0
    for combn in range(1, (max_comb_n+1), 1):
        rep_comb_data_path = rep_comb_data_path_format % (central_neurons, combn)
        if not os.path.exists(rep_comb_data_path):
            pass
        else:
            save_rep_comb_data_path = save_rep_comb_data_path_format % (central_neurons, combn)
            is_Exist_file(save_rep_comb_data_path)
            if s == 0:
                csr0 = colindex_to_csr_matrix(rep_comb_data_path, gene_n, start_row_n=0, max_readrow_n=np.inf)
                val = os.system("cp %s %s" % (rep_comb_data_path, save_rep_comb_data_path))
                s = 1
            else:
                csr0 = gene_combinations_union(rep_comb_data_path, csr0, save_rep_comb_data_path,
                                               gene_n, max_readrow_n)


# step 2: unfold represent gene combination
# Cartesian prodcut
def cartesian_product_for_repgenes(rep_comb_genes):
    from itertools import product

    if len(rep_comb_genes) == 1:
        result = [(x,) for x in rep_comb_genes[0]]
    elif len(rep_comb_genes) == 2:
        result = product(rep_comb_genes[0], rep_comb_genes[1])
    elif len(rep_comb_genes) == 3:
        result = product(rep_comb_genes[0], rep_comb_genes[1], rep_comb_genes[2])
    elif len(rep_comb_genes) == 4:
        result = product(rep_comb_genes[0], rep_comb_genes[1], rep_comb_genes[2], rep_comb_genes[3])
    elif len(rep_comb_genes) == 5:
        result = product(rep_comb_genes[0], rep_comb_genes[1], rep_comb_genes[2], rep_comb_genes[3], rep_comb_genes[4])
    elif len(rep_comb_genes) == 6:
        result = product(rep_comb_genes[0], rep_comb_genes[1], rep_comb_genes[2], rep_comb_genes[3], rep_comb_genes[4],
                         rep_comb_genes[5])
    elif len(rep_comb_genes) == 7:
        result = product(rep_comb_genes[0], rep_comb_genes[1], rep_comb_genes[2], rep_comb_genes[3], rep_comb_genes[4],
                         rep_comb_genes[5], rep_comb_genes[6])
    elif len(rep_comb_genes) == 8:
        result = product(rep_comb_genes[0], rep_comb_genes[1], rep_comb_genes[2], rep_comb_genes[3], rep_comb_genes[4],
                         rep_comb_genes[5], rep_comb_genes[6], rep_comb_genes[7])
    elif len(rep_comb_genes) == 9:
        result = product(rep_comb_genes[0], rep_comb_genes[1], rep_comb_genes[2], rep_comb_genes[3], rep_comb_genes[4],
                         rep_comb_genes[5], rep_comb_genes[6], rep_comb_genes[7], rep_comb_genes[8])
    elif len(rep_comb_genes) == 10:
        result = product(rep_comb_genes[0], rep_comb_genes[1], rep_comb_genes[2], rep_comb_genes[3], rep_comb_genes[4],
                         rep_comb_genes[5], rep_comb_genes[6], rep_comb_genes[7], rep_comb_genes[8], rep_comb_genes[9])
    elif len(rep_comb_genes) == 11:
        result = product(rep_comb_genes[0], rep_comb_genes[1], rep_comb_genes[2], rep_comb_genes[3], rep_comb_genes[4],
                         rep_comb_genes[5], rep_comb_genes[6], rep_comb_genes[7], rep_comb_genes[8], rep_comb_genes[9],
                         rep_comb_genes[10])
    elif len(rep_comb_genes) == 12:
        result = product(rep_comb_genes[0], rep_comb_genes[1], rep_comb_genes[2], rep_comb_genes[3], rep_comb_genes[4],
                         rep_comb_genes[5], rep_comb_genes[6], rep_comb_genes[7], rep_comb_genes[8], rep_comb_genes[9],
                         rep_comb_genes[10], rep_comb_genes[11])
    elif len(rep_comb_genes) == 13:
        result = product(rep_comb_genes[0], rep_comb_genes[1], rep_comb_genes[2], rep_comb_genes[3], rep_comb_genes[4],
                         rep_comb_genes[5], rep_comb_genes[6], rep_comb_genes[7], rep_comb_genes[8], rep_comb_genes[9],
                         rep_comb_genes[10], rep_comb_genes[11], rep_comb_genes[12])
    elif len(rep_comb_genes) == 14:
        result = product(rep_comb_genes[0], rep_comb_genes[1], rep_comb_genes[2], rep_comb_genes[3], rep_comb_genes[4],
                         rep_comb_genes[5], rep_comb_genes[6], rep_comb_genes[7], rep_comb_genes[8], rep_comb_genes[9],
                         rep_comb_genes[10], rep_comb_genes[11], rep_comb_genes[12], rep_comb_genes[13])
    elif len(rep_comb_genes) == 15:
        result = product(rep_comb_genes[0], rep_comb_genes[1], rep_comb_genes[2], rep_comb_genes[3], rep_comb_genes[4],
                         rep_comb_genes[5], rep_comb_genes[6], rep_comb_genes[7], rep_comb_genes[8], rep_comb_genes[9],
                         rep_comb_genes[10], rep_comb_genes[11], rep_comb_genes[12], rep_comb_genes[13], rep_comb_genes[14])
    else:
        result = []
        print("Error: combination number of represent genes is more than 15, please chech again.")
    return result


# unfold represent genes
def unfold_represent_gene_combinations(central_neurons, rep_comb_data_path_format, rep_gene_digit_info,
                                       save_comb_unfold_path_format, save_stat_path, max_comb_n=9):
    for combn in range(1, (max_comb_n + 1), 1):
        rep_comb_data_path = rep_comb_data_path_format % (central_neurons, combn)
        save_comb_unfold_path = save_comb_unfold_path_format % (central_neurons, combn)
        is_Exist_file(save_comb_unfold_path)
        if not os.path.exists(rep_comb_data_path):
            pass
        else:
            data_name = rep_comb_data_path.split("/")[-1]
            stat_dict = {"data name": [data_name],
                         "gene combination info.": [{}]}
            total_combn = 0
            with open(save_comb_unfold_path, 'a') as a:
                with open(rep_comb_data_path, 'r') as f:
                    for line in f:
                        rep_gene_codes = ''.join(line.split("\n")).split("\t")
                        rep_comb_genes = [rep_gene_digit_info[int(x)] for x in rep_gene_codes]
                        result = cartesian_product_for_repgenes(rep_comb_genes)
                        for gene_comb_tuple in result:
                            a.write(("%s\t"*(len(gene_comb_tuple) - 1) + "%s\n") % gene_comb_tuple)
                            if len(gene_comb_tuple) not in stat_dict['gene combination info.'][0]:
                                stat_dict['gene combination info.'][0][len(gene_comb_tuple)] = 1
                            else:
                                stat_dict['gene combination info.'][0][len(gene_comb_tuple)] += 1
                            total_combn += 1
            # save statistics
            stat_dict['Num. of total gene combination'] = [total_combn]
            stat_data = pd.DataFrame(stat_dict)
            if os.path.exists(save_stat_path):
                stat_data0 = pd.read_excel(save_stat_path)
                stat_data = pd.concat([stat_data0, stat_data], axis=0)
            stat_data.to_excel(save_stat_path, index=False)


def main_union_and_unfold_for_represent_gene_combinations(fneurons_data, screen_data, max_comb_n,
                                                          save_central_neurons_comb_dir, save_union_rep_gene_dir,
                                                          save_unfold_gene_comb_dir, pool_n, pool_i):
    # iteration
    all_index = fneurons_data['order'].tolist()
    iterindex = split_index_for_multiprocess(all_index, pool_n, pool_i)
    # save
    rep_comb_data_path_format = save_central_neurons_comb_dir + "/%s_combination%s.log"
    save_rep_comb_data_path_format = save_union_rep_gene_dir + "/union_represent_gene_combinations%s_%s.log"
    save_comb_unfold_path_format = save_unfold_gene_comb_dir + "/gene_combination_for_%s_%s.log"
    save_stat_path = save_unfold_gene_comb_dir + "/0statistics_neuron_sets_gene_combinations_pool_%s.xlsx" % pool_i
    is_Exist_file(save_stat_path)
    # execute
    gene_n = screen_data.shape[0]
    gene_code_dict = gene_dict(screen_data)
    for i in iterindex:
        neurons = eval(fneurons_data.loc[i, 'T.B. neurons'])
        central_neurons = "%s_central_neurons" % i
        print(central_neurons)
        start = time.time()
        # gene classification & combination
        # gene classification
        rep_screen_data_nbase, rep_gene_expr_info = classify_genes(screen_data, neurons)
        rep_gene_digit_info = {}
        for rep_gene, genes in rep_gene_expr_info.items():
            rep_gene_digit_info[gene_code_dict[rep_gene]] = [gene_code_dict[x] for x in genes]
        print("union_gene_combinations_for_a_central_neuron ...")
        union_gene_combinations_for_a_central_neuron(central_neurons, rep_comb_data_path_format,
                                                     save_rep_comb_data_path_format,
                                                     gene_n, max_comb_n, max_readrow_n=10000)
        end = time.time()
        print("Using time:", end - start)
        print("unfold_represent_gene_combinations ... ")
        unfold_represent_gene_combinations(central_neurons, save_rep_comb_data_path_format, rep_gene_digit_info,
                                           save_comb_unfold_path_format, save_stat_path, max_comb_n)
        print("Using time:", time.time() - end)


# if __name__ == '__main__':
    # import sys
    # main_dir, screen_data_path, fneurons_data_path, data_dir, max_comb_n, pool_n, pool_i = sys.argv[1:]
    # print('-----------------------------------------------------------------------')
    # print('**************** execute: union of gene combinatins ********************')
    # # data path & parameters
    # os.chdir(main_dir)
    # max_comb_n = int(max_comb_n)
    # pool_n = int(pool_n)
    # pool_i = int(pool_i)
    # # read
    # screen_data = pd.read_csv(screen_data_path)
    # fneurons_data = pd.read_excel(fneurons_data_path)
    # # 固定参数
    # save_central_neurons_comb_dir = data_dir + "/central_neuron_combinations"
    # save_union_rep_gene_dir = data_dir + "/union_represent_gene_combinations"
    # save_unfold_gene_comb_dir = data_dir + "/unfold_gene_combinations"
    # mkdir(save_central_neurons_comb_dir)
    # mkdir(save_union_rep_gene_dir)
    # mkdir(save_unfold_gene_comb_dir)
    # print("combinations & unfolding...")
    # main_union_and_unfold_for_represent_gene_combinations(fneurons_data, screen_data, max_comb_n,
    #                                                       save_central_neurons_comb_dir, save_union_rep_gene_dir,
    #                                                       save_unfold_gene_comb_dir, pool_n, pool_i)

