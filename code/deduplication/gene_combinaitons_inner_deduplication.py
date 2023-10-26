import time
import os
import pandas as pd
import numpy as np
from scipy.sparse import vstack
from scipy.sparse import csr_matrix
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


def gene_combinations_union(csr1, csr0, save_path, max_readrow_n=10000):
    is_Exist_file(save_path)
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


if __name__ == '__main__':
    main_dir = "../../"
    os.chdir(main_dir)
    # 可变参数
    start_combn = 11
    other_combn_order = range((start_combn+1), 14, 1)
    max_readrow_n = 10000
    screen_data_path = "./result/example-analyzed_data/" \
                       "2_base_PRAUC0.8_scRNA_seq_data_with_thv70_and_111_neuropal_reporter.csv"
    data_path_format = "./code/deduplication/example/union_gene_combinations_for_N20_combn%s.log"
    save_dir = "./code/deduplication/example/deduplication"
    mkdir(save_dir)
    # 固定参数
    save_path_format = save_dir + '/' + data_path_format.split("/")[-1]
    # initial
    data0_path = data_path_format % start_combn
    save_data0_path = save_path_format % start_combn
    is_Exist_file(save_data0_path)
    data0 = pd.read_csv(data0_path, sep='\t', header=None)
    data0.drop_duplicates(inplace=True)
    data0.to_csv(save_data0_path, index=False, header=False)
    # deduplication
    screen_data = pd.read_csv(screen_data_path)
    gene_n = screen_data.shape[0]
    csr0 = colindex_to_csr_matrix(data0_path, gene_n, start_row_n=0, max_readrow_n=np.inf)
    for combn in other_combn_order:
        data1_path = data_path_format % combn
        save_path = save_path_format % combn
        csr1 = colindex_to_csr_matrix(data1_path, gene_n, start_row_n=0, max_readrow_n=np.inf)
        csr0 = gene_combinations_union(csr1, csr0, save_path, max_readrow_n=max_readrow_n)



