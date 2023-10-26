from multiprocessing.dummy import Pool
from union_unfold_represent_gene_combinations import *
import numpy as np
import warnings
warnings.filterwarnings("ignore")


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


def get_csr_matrix_dict(central_neuron, data_path_format, max_comb_n, gene_n, file_first_gene_combn,
                        single_csr_max_rown=10000, max_readrow_n=1000000):
    csr_dict = {}
    read_file_list = []
    csr_count = 0
    s = 0
    for combn in range(1, (max_comb_n + 1), 1):
        data_path = data_path_format % (central_neuron, combn)
        if os.path.exists(data_path):
            if len(read_file_list) < file_first_gene_combn:
                if s == 0:
                    csr = colindex_to_csr_matrix(data_path, gene_n, start_row_n=0, max_readrow_n=np.inf)
                    csr_count += csr.shape[0]
                    partn = int(csr.shape[0]/single_csr_max_rown) + 1
                    for i in range(partn):
                        if csr[i*single_csr_max_rown:(i+1)*single_csr_max_rown].shape[0] != 0:
                            csr_dict[i] = csr[i*single_csr_max_rown:(i+1)*single_csr_max_rown]
                        else:
                            break
                    del csr
                    read_file_list.append(data_path)
                    s = 1
                else:
                    if csr_count < max_readrow_n:
                        left_readrow_n = max_readrow_n - csr_count
                        csr = colindex_to_csr_matrix(data_path, gene_n, start_row_n=0, max_readrow_n=left_readrow_n)
                        csr_count += csr.shape[0]
                        partn = int(csr.shape[0] / single_csr_max_rown) + 1
                        start_i = max(list(csr_dict.keys())) + 1
                        for i in range(partn):
                            if csr[i * single_csr_max_rown:(i + 1) * single_csr_max_rown].shape[0] != 0:
                                csr_dict[start_i+i] = csr[i * single_csr_max_rown:(i + 1) * single_csr_max_rown]
                            else:
                                break
                        del csr
                        read_file_list.append(data_path)
                    else:
                        break
            else:
                break
    return csr_dict


def get_union_combn2(parameters):
    csr_pair_codes, csr0_dict, csr1_dict, file1_first_gene_combn, save_h5_path = parameters
    print("part csr_pair_codes:", len(csr_pair_codes))
    min_combn_list = []
    csr_pair_info = {}
    is_Exist_file(save_h5_path)
    store = pd.HDFStore(save_h5_path)
    for pair_code in csr_pair_codes:
        csr0_chunk = csr0_dict[pair_code[0]]
        csr1_chunk = csr1_dict[pair_code[1]]
        if (csr0_chunk.shape[0] != 0) & (csr1_chunk.shape[0] != 0):
            csr_combn = (csr0_chunk.dot(csr1_chunk.T)) + \
                        (csr0_chunk.dot((csr_matrix(np.ones(csr1_chunk.shape)) - csr1_chunk).T)) + \
                        ((csr_matrix(np.ones(csr0_chunk.shape)) - csr0_chunk).dot(csr1_chunk.T))
            a = list(set(csr_combn.toarray().reshape(1, -1)[0].astype(int).tolist()))
            a.sort(reverse=False)
            one_min_combn_tuple = tuple(a[:file1_first_gene_combn])
            min_combn_list.extend(one_min_combn_tuple)
            csr_pair_info[pair_code] = one_min_combn_tuple
            store['%s-%s' % (str(pair_code), str(one_min_combn_tuple))] = pd.DataFrame(csr_combn.toarray()).astype(int)
            # print("pair_code:", pair_code, csr0_chunk.shape, csr1_chunk.shape)
            # print(store['%s-%s' % (str(pair_code), str(one_min_combn_tuple))].shape)
            del csr_combn
    store.close()
    min_combn_list = list(set(min_combn_list))
    min_combn_list.sort(reverse=False)
    min_combn_list = min_combn_list[:file1_first_gene_combn]
    store = pd.HDFStore(save_h5_path)
    for pair_code, one_min_combn_tuple in csr_pair_info.items():
        if len(set(one_min_combn_tuple) & set(min_combn_list)) != 0:
            pass
        else:
            del store['%s-%s' % (str(pair_code), str(one_min_combn_tuple))]
    store.close()


# crs add
def apply_csr_addition(row, csr0, csr1):
    csr0_array = csr0[row['csr0_index']].toarray()
    csr1_array = np.array(pd.DataFrame(csr1.toarray()).loc[row['csr1_index_list'], :])
    csr_add = csr_matrix(pd.DataFrame(csr1_array + csr0_array).replace(2, 1).drop_duplicates())
    return csr_add


# first, find index of minimal union AND csr add
def find_minimal_union_index_and_csr_add(parameters):
    print("Function: find_minimal_union_index_and_csr_add ...")
    (pair_code_min_combn_dict, csr0_dict, csr1_dict,
     read_pooli_h5_path_format, save_union_h5_path_format, idn, whole_min_combn_list) = parameters
    for combn in whole_min_combn_list:
        save_union_h5_path = save_union_h5_path_format % (idn, combn)
        is_Exist_file(save_union_h5_path)

    print("part pair_code_min_combn_dict:", len(pair_code_min_combn_dict))
    pooli_info = {}
    for str_one_pair_code, min_combn_info in pair_code_min_combn_dict.items():
        pool_i, str_one_min_combn_tuple, share_min_combn_tuple = min_combn_info
        if pool_i not in pooli_info:
            pooli_info[pool_i] = [(str_one_pair_code, str_one_min_combn_tuple, share_min_combn_tuple)]
        else:
            pooli_info[pool_i].append((str_one_pair_code, str_one_min_combn_tuple, share_min_combn_tuple))
    #
    for pool_i, pair_csr_info_list in pooli_info.items():
        read_pooli_h5_path = read_pooli_h5_path_format % (idn, pool_i)
        store = pd.HDFStore(read_pooli_h5_path)
        for pair_csr_info in pair_csr_info_list:
            str_one_pair_code, str_one_min_combn_tuple, share_min_combn_tuple = pair_csr_info
            df_combn = store['%s-%s' % (str_one_pair_code, str_one_min_combn_tuple)]
            csr0 = csr0_dict[eval(str_one_pair_code)[0]]
            csr1 = csr1_dict[eval(str_one_pair_code)[1]]
            print("pool_i - str_one_pair_code:", pool_i, str_one_pair_code)
            print("csr0.shape:", csr0.shape, "csr1.shape:", csr1.shape)
            print(df_combn.shape)

            # csr0 中包含 csr1 的数据: 意味着这些数据在两个neuron sets 中都可识别，降低计算量考虑可单独提出计算
            print("csr0 中包含 csr1 的数据 ...")
            df_temp0 = pd.DataFrame(np.array(df_combn) - np.array(csr0.sum(axis=1)).reshape(csr0.shape[0], 1))
            both_index0 = df_temp0.replace(0, np.nan)[df_temp0.replace(0, np.nan).isnull().T.any()].index.tolist()
            keep_index0 = [i for i in df_temp0.index if i not in both_index0]
            csr0_both = csr_matrix(pd.DataFrame(csr0.toarray()).loc[both_index0, :])
            csr0 = csr_matrix(pd.DataFrame(csr0.toarray()).loc[keep_index0, :])
            df_combn = df_combn.loc[keep_index0, :]
            # print("keep_index0:", len(keep_index0))
            # print("both_index0:", len(both_index0))
            # print("df_combn.shape:", df_combn.shape, "csr0.shape:", csr0.shape)
            del df_temp0

            # csr1 中包含 csr0 的数据: 意味着这些数据在两个neuron sets 中都可识别，降低计算量考虑可单独提出计算
            print("csr1 中包含 csr0 的数据 ...")
            df_temp1 = pd.DataFrame(np.array(df_combn) - np.array(csr1.T.sum(axis=0)))
            s = (df_temp1 == 0)
            result1 = s.agg(lambda s: s.index[s].values, axis=1).to_frame(name='csr1_both_index_list')
            result1['csr1_both_index_list'] = result1['csr1_both_index_list'].apply(lambda x: list(x))
            both_index1 = list(set(result1['csr1_both_index_list'].sum()))
            keep_index1 = [i for i in range(csr1.shape[0]) if i not in both_index1]
            csr1_both = csr_matrix(pd.DataFrame(csr1.toarray()).loc[both_index1, :])
            csr1 = csr_matrix(pd.DataFrame(csr1.toarray()).loc[keep_index1, :])
            df_combn = df_combn.loc[:, keep_index1].reset_index(drop=True)
            df_combn.columns = [x for x in range(df_combn.shape[1])]
            # print("keep_index1:", len(keep_index1))
            # print("both_index1:", len(both_index1))
            # print("df_combn.shape:", df_combn.shape, "csr1.shape:", csr1.shape)
            del df_temp1
            del result1
            # crs add
            print("crs add")
            print("df_combn.shape:", df_combn.shape, "; csr0_both.shape:", csr0_both.shape,
                  "; csr1_both.shape:", csr1_both.shape)
            # print("csr0[0].nonzero()[1]:", len(csr0[0].nonzero()[1]))
            # print("csr1[0].nonzero()[1]:", len(csr1[0].nonzero()[1]))
            if csr0.shape[0] < csr1.shape[0]:
                csr_add_list = [vstack([csr0[c]] * csr1.shape[0]) + csr1 for c in
                                range(csr0.shape[0])]
            else:
                csr_add_list = [vstack([csr1[c]] * csr0.shape[0]) + csr0 for c in
                                range(csr1.shape[0])]
            if len(csr_add_list) != 0:
                csr_add = vstack(csr_add_list)
                del csr_add_list
            else:
                csr_add = csr0[:0]
            df = pd.DataFrame(csr_add.toarray()).replace(2, 1).drop_duplicates()
            print("csr_add.shape:", csr_add.shape, "after dropping duplicates, df.shape:", df.shape)
            del csr_add
            csr0_both_df = pd.DataFrame(csr0_both.toarray())
            csr1_both_df = pd.DataFrame(csr1_both.toarray())
            a = df.sum(axis=1)
            a0 = csr0_both_df.sum(axis=1)
            a1 = csr1_both_df.sum(axis=1)
            share_min_combn_list = list(share_min_combn_tuple)
            share_min_combn_list.sort(reverse=False)
            for combn in share_min_combn_list:
                csr_add = vstack([csr_matrix(df.loc[a == combn, :]),
                                  csr_matrix(csr0_both_df.loc[a0 == combn, :]),
                                  csr_matrix(csr1_both_df.loc[a1 == combn, :])])
                csr_add = csr_matrix(pd.DataFrame(csr_add.toarray()).drop_duplicates())
                if csr_add.shape[0] != 0:
                    prev_combn_list = [x for x in whole_min_combn_list if x <= combn]
                    pcsr0_list = []
                    for pcombn in prev_combn_list:
                        prev_union_h5_path = save_union_h5_path_format % (idn, pcombn)
                        if os.path.exists(prev_union_h5_path):
                            pstore = pd.HDFStore(prev_union_h5_path)
                            if '/df' in pstore.keys():
                                pcsr0_list.append(csr_matrix(pstore['/df']))
                            pstore.close()
                    save_union_h5_path = save_union_h5_path_format % (idn, combn)
                    temp_store = pd.HDFStore(save_union_h5_path)
                    if len(pcsr0_list) != 0:
                        pcsr0 = vstack(pcsr0_list)
                        del pcsr0_list
                        pcsr0, kcsr_add = gene_combinations_union3(csr_add, pcsr0, save_path='',
                                                                   max_readrow_n=10000, write_mode=False)
                        if (kcsr_add.shape[0] != 0) & ('/df' in temp_store.keys()):
                            temp_store['/df'] = pd.concat([temp_store['/df'], pd.DataFrame(kcsr_add.toarray())],
                                                          axis=0).drop_duplicates().astype(int)
                        elif kcsr_add.shape[0] != 0:
                            temp_store['/df'] = pd.DataFrame(kcsr_add.toarray()).astype(int)
                        else:
                            pass
                    else:
                        temp_store['/df'] = pd.DataFrame(csr_add.toarray()).astype(int)
                    temp_store.close()
            del df
            del csr0_both_df
            del csr1_both_df
        store.close()


# 执行多进程
def do_multiprocessing(argument_list, function, pool_n=10):
    import time
    s = time.time()
    print('Parent process %s.' % os.getpid())
    p = Pool(pool_n)
    for parameter_tuple in argument_list:
        p.apply_async(function, args=(parameter_tuple,))
    print('Waiting for all subprocesses done...')
    p.close()
    p.join()
    e = time.time()
    print('All using time: %.3f' % (e - s))
    print('Make_corpus_by_pool: All subprocesses done.')


def main_two_neuron_sets_union2(csr0_dict, csr1_dict, save_union_path_format, max_gene_combn,
                                file1_first_gene_combn, pool_n, idn):
    print("main_two_neuron_sets_union2 ...")
    pair_codes = []
    for i in csr0_dict.keys():
        for j in csr1_dict.keys():
            pair_codes.append((i, j))
    # 筛选出 union 后满足 file1_first_gene_combn 文件的 pair_codes
    print("\n* 1. To get combintion nnumber of minimal union ...")
    print("sum of pair_codes:", len(pair_codes))
    print("csr0_dict:", len(csr0_dict), [(k, v.shape) for k, v in csr0_dict.items()])
    print("csr1_dict:", len(csr1_dict), [(k, v.shape) for k, v in csr1_dict.items()])
    partn = int(len(pair_codes) / pool_n) + 1
    argument_list = []
    save_pooli_dir = "/".join(save_union_path_format.split("/")[:-1])
    for pool_i in range(pool_n):
        chunk_pair_codes = pair_codes[pool_i*partn:(pool_i+1)*partn]
        if len(chunk_pair_codes) != 0:
            chunk_csr0_dict = {}
            chunk_csr1_dict = {}
            for one_pair_code in chunk_pair_codes:
                i, j = one_pair_code
                if i not in chunk_csr0_dict:
                    chunk_csr0_dict[i] = csr0_dict[i]
                else:
                    pass
                if j not in chunk_csr1_dict:
                    chunk_csr1_dict[j] = csr1_dict[j]
                else:
                    pass
            save_pooli_h5_path = save_pooli_dir + "/get_union_combination_matrix_idn%s_pool%s.h5" % (idn, pool_i + 1)
            parameter_tuple = (chunk_pair_codes, chunk_csr0_dict, chunk_csr1_dict, file1_first_gene_combn, save_pooli_h5_path)
            argument_list.append(parameter_tuple)
    pool_n1 = len(argument_list)
    print("pool n:", pool_n1)
    if pool_n1 > 1:
        do_multiprocessing(argument_list, get_union_combn2, pool_n1)
    elif pool_n1 == 1:
        get_union_combn2(argument_list[0])
    else:
        pass
    del argument_list
    del parameter_tuple
    del chunk_csr0_dict
    del chunk_csr1_dict
    print("over: get_union_combn2.")

    # find whole minimal file1_first_gene_combn gene combinations
    print("\n* 2. to obtain combination number of whole minimal union ...")
    whole_min_combn_dict = {}
    whole_min_combn_list = []
    for pool_i in range(pool_n):
        save_pooli_h5_path = save_pooli_dir + "/get_union_combination_matrix_idn%s_pool%s.h5" % (idn, pool_i + 1)
        if os.path.exists(save_pooli_h5_path):
            store = pd.HDFStore(save_pooli_h5_path)
            for key in store.keys():
                str_one_pair_code, str_one_min_combn_tuple = tuple(key.strip('/').split('-'))
                whole_min_combn_list.extend(list(eval(str_one_min_combn_tuple)))
                if str_one_min_combn_tuple not in whole_min_combn_dict:
                    whole_min_combn_dict[str_one_min_combn_tuple] = [(pool_i + 1, str_one_pair_code)]
                else:
                    whole_min_combn_dict[str_one_min_combn_tuple].append((pool_i + 1, str_one_pair_code))
            store.close()
        else:
            pass
    whole_min_combn_list = list(set(whole_min_combn_list))
    whole_min_combn_list.sort(reverse=False)
    min_combn_list = whole_min_combn_list[:file1_first_gene_combn]
    print("min_combn_list:", min_combn_list)
    print("union of whole min. combn list with each pair csr's")
    min_combn_dict = {}
    for str_one_min_combn_tuple, str_one_pair_code_list in whole_min_combn_dict.items():
        share_min_combn_tuple = tuple(set(eval(str_one_min_combn_tuple)) & set(min_combn_list))
        if len(share_min_combn_tuple) != 0:
            for (pool_i, str_one_pair_code) in str_one_pair_code_list:
                min_combn_dict[str_one_pair_code] = (pool_i, str_one_min_combn_tuple, share_min_combn_tuple)
        else:
            del_pooli_dict = {}
            for (pool_i, str_one_pair_code) in str_one_pair_code_list:
                if pool_i not in del_pooli_dict:
                    del_pooli_dict[pool_i] = ['%s-%s' % (str_one_pair_code, str_one_min_combn_tuple)]
                else:
                    del_pooli_dict[pool_i].append('%s-%s' % (str_one_pair_code, str_one_min_combn_tuple))
            for pool_i in del_pooli_dict.keys():
                save_pooli_h5_path = save_pooli_dir + "/get_union_combination_matrix_idn%s_pool%s.h5" % (idn, pool_i)
                store = pd.HDFStore(save_pooli_h5_path)
                for key in del_pooli_dict[pool_i]:
                    del store[key]
                store.close()

    # get index of minimal combinations for csr pair AND add
    print("\n* 3. To get index of minimal combinations for csr pair AND add ...")
    print("sum of min_combn_dict:", len(min_combn_dict))
    list_str_pair_codes = list(min_combn_dict.keys())
    partn1 = int(len(list_str_pair_codes) / pool_n) + 1
    argument_list = []
    for pool_j in range(pool_n):
        chunk_str_pair_codes = list_str_pair_codes[pool_j*partn1:(pool_j+1)*partn1]
        chunk_pair_code_min_combn_dict = {}
        chunk_csr0_dict = {}
        chunk_csr1_dict = {}
        for str_one_pair_code in chunk_str_pair_codes:
            chunk_pair_code_min_combn_dict[str_one_pair_code] = min_combn_dict[str_one_pair_code]
            i, j = eval(str_one_pair_code)
            if i not in chunk_csr0_dict:
                chunk_csr0_dict[i] = csr0_dict[i]
            else:
                pass
            if j not in chunk_csr1_dict:
                chunk_csr1_dict[j] = csr1_dict[j]
            else:
                pass
        read_pooli_h5_path_format = save_pooli_dir + "/get_union_combination_matrix_idn%s_pool%s.h5"
        save_poolj_union_h5_path_format = save_union_path_format[:-4] + ("_pool-%s" % (pool_j + 1)) + '.h5'
        parameter_tuple = (chunk_pair_code_min_combn_dict, chunk_csr0_dict, chunk_csr1_dict,
                           read_pooli_h5_path_format, save_poolj_union_h5_path_format, idn, min_combn_list)
        argument_list.append(parameter_tuple)
    pool_n2 = len(argument_list)
    print("pool_n2:", pool_n2)
    if pool_n2 > 1:
        do_multiprocessing(argument_list, find_minimal_union_index_and_csr_add, pool_n2)
    elif pool_n2 == 1:
        find_minimal_union_index_and_csr_add(argument_list[0])
    else:
        pass
    del argument_list
    del parameter_tuple
    del chunk_csr1_dict
    del chunk_csr0_dict
    print("it's over: find_minimal_union_index_and_csr_add")

    # 合并去重之合并
    print("\n* 4. To union & drop duplicates")
    print("(1) union by combn")
    csr_add_dict = {}  # union by combn
    for combn in min_combn_list:
        df_list = []
        for i in range(pool_n):
            read_pooj_union_h5_path = (save_union_path_format[:-4] + "_pool-%s") % (idn, combn, (i + 1))+'.h5'
            if os.path.exists(read_pooj_union_h5_path):
                store = pd.HDFStore(read_pooj_union_h5_path)
                if '/df' in store.keys():
                    df_list.append(store['df'])
                else:
                    pass
                store.close()
        if len(df_list) != 0:
            csr_add_dict[combn] = csr_matrix(pd.concat(df_list, axis=0).drop_duplicates())
            print(combn, csr_add_dict[combn].shape)
    # 合并去重之去重
    print("(2) drop duplicates by function: gene_combinations_union3")
    if len(csr_add_dict) != 0:
        combn_sort = list(csr_add_dict.keys())
        combn_sort.sort(reverse=False)
        csr0 = csr_add_dict[combn_sort[0]]
        print(combn_sort[0], csr0.shape)
        save_union_path = save_union_path_format % (idn, combn_sort[0])
        is_Exist_file(save_union_path)
        # write
        with open(save_union_path, 'a') as a:
            for index in range(csr0.shape[0]):
                col_index = csr0[index].nonzero()[1].tolist()
                col_index.sort(reverse=False)
                a.write(("%s\t" * (len(col_index) - 1) + '%s\n') % tuple(col_index))
        for combn in combn_sort[1:]:
            csr1 = csr_add_dict[combn]
            save_union_path = save_union_path_format % (idn, combn)
            is_Exist_file(save_union_path)
            csr0, kcsr1 = gene_combinations_union3(csr1, csr0, save_union_path, max_readrow_n=10000)
            print(combn, kcsr1.shape)

    # 中间文件删除
    print("\n* 5. To delete temporary files.")
    for i in range(pool_n):
        read_pooli_h5_path = save_pooli_dir + "/get_union_combination_matrix_idn%s_pool%s.h5" % (idn, i + 1)
        is_Exist_file(read_pooli_h5_path)
        for combn in min_combn_list:
            read_pooj_union_h5_path = (save_union_path_format[:-4] + "_pool-%s") % (idn, combn, (i + 1)) + '.h5'
            is_Exist_file(read_pooj_union_h5_path)
    # to delete previous files
    if (idn != 3791) and (idn != 14745):
        for combn in range(1, max_gene_combn+1, 1):
            prev_union_combn_path = save_union_path_format % (idn-1, combn)
            is_Exist_file(prev_union_combn_path)


def intialize_save_union_path(save_union_path_format, intial_central_neuron, read_format_path,
                              max_comb_n, file1_first_gene_combn, idn):
    for combn in range(1, (max_comb_n + 1), 1):
        save_union_combn_path = save_union_path_format % (idn, combn)
        is_Exist_file(save_union_combn_path)
    data_path_list = []
    for combn in range(1, (max_comb_n + 1), 1):
        save_union_combn_path = save_union_path_format % (idn, combn)
        data_path = read_format_path % (intial_central_neuron, combn)
        if os.path.exists(data_path):
            if len(data_path_list) < file1_first_gene_combn:
                data_path_list.append(data_path)
                val = os.system("cp %s %s" % (data_path, save_union_combn_path))
            else:
                break


def statistics_gene_comb_info_from_file(central_neuron, read_format_path,
                                        max_gene_combn, first_gene_combn,
                                        save_union_path_format):
    save_dir = "/".join(save_union_path_format.split("/")[:-1])
    data_path_list = []
    gen_combn_info = {}
    for combn in range(1, (max_gene_combn + 1), 1):
        data_path = read_format_path % (central_neuron, combn)
        if os.path.exists(data_path):
            if len(data_path_list) < first_gene_combn:
                data_path_list.append(data_path)
                stat_rown_data_path = save_dir + '/' + data_path.split('/')[-1] + '_file_count'
                val = os.system("wc -l %s > %s" % (data_path, stat_rown_data_path))
                with open(stat_rown_data_path, 'r') as f:
                    for line in f:
                        file_count = int(line.strip(" ").split(" ")[0])
                        break
                gen_combn_info[combn] = file_count
                is_Exist_file(stat_rown_data_path)
    return gen_combn_info


# union gene combinations of central neurons
def union_gene_combinations_central_neurons(iterindex, gene_n, read_format_path,
                                            save_union_path_format, save_stat_data_path,
                                            file0_first_gene_combn=2, file0_max_readrow_n=np.inf,
                                            file1_first_gene_combn=2, file1_max_readrow_n=100000,
                                            single_csr_max_rown=10000, max_gene_combn=39, pool_n=1):
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
    s = 0
    union_file_code = ''
    idn = 0
    for central_neuron_id in iterindex:
        # if central_neuron_id == 255:
        central_neuron = "%s_central_neurons" % central_neuron_id
        idn += 1
        print("\n%s -- %s"%(idn, central_neuron))
        if s == 0:
            csr0_dict = get_csr_matrix_dict(central_neuron, read_format_path,
                                            max_gene_combn, gene_n, file0_first_gene_combn,
                                            single_csr_max_rown=single_csr_max_rown,
                                            max_readrow_n=100)
            if len(csr0_dict) != 0:
                union_file_code = str(central_neuron_id)
                print("Start central_neuron ID:", central_neuron_id)
                intialize_save_union_path(save_union_path_format, central_neuron, read_format_path,
                                          max_gene_combn, file1_first_gene_combn, idn)
                # statistics
                gen_combn_info0 = statistics_gene_comb_info_from_file(central_neuron, read_format_path,
                                                                      max_gene_combn, file0_first_gene_combn,
                                                                      save_union_path_format)
                gen_combn_info1 = statistics_gene_comb_info_from_file(idn, save_union_path_format,
                                                                      max_gene_combn, file1_first_gene_combn,
                                                                      save_union_path_format)
                stat_dict['central neuron id'].append(central_neuron_id)
                stat_dict['gene combination info. of file1'].append(gen_combn_info0)
                stat_dict['first gene combn. selected of file1'].append(file0_first_gene_combn)
                stat_dict['total combn. selected of file1'].append(sum(list(gen_combn_info0.values())))
                stat_dict['total combn. of file1'].append(sum(list(gen_combn_info0.values())))

                stat_dict['gene combination info. of file2'].append(gen_combn_info1)
                stat_dict['first gene combn. selected of file2'].append(file1_first_gene_combn)
                stat_dict['total combn. selected of file2'].append(sum(list(gen_combn_info1.values())))
                stat_dict['total combn. of file2'].append(sum(list(gen_combn_info1.values())))

                stat_dict['union file code'].append(union_file_code)
                stat_dict['gene combination info. of union file'].append(gen_combn_info1)
                stat_dict['total combn. of union file'].append(sum(list(gen_combn_info1.values())))
                s = 1
            else:
                pass
        else:
            # print("Union central_neuron ID:", central_neuron_id)
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
        # if idn == 100:
        #     break
    # save statistics
    stat_data = pd.DataFrame(stat_dict)
    stat_data.to_excel(save_stat_data_path, index=False)


if __name__ == '__main__':
    import sys
    import os

    (main_dir, screen_data_path, data_dir, file0_first_gene_combn, file1_first_gene_combn,
     file0_max_readrow_n, file1_max_readrow_n, max_gene_combn) = sys.argv[1:]
    print('-----------------------------------------------------------------------')
    print('**************** execute: union of gene combinations ********************')
    # data path & parameters
    os.chdir(main_dir)
    max_gene_combn = int(max_gene_combn)
    file0_first_gene_combn = int(file0_first_gene_combn)  # 2
    file1_first_gene_combn = int(file1_first_gene_combn)  # 3
    file0_max_readrow_n = int(file0_max_readrow_n)  # 50000
    file1_max_readrow_n = int(file1_max_readrow_n)  # 10000
    single_csr_max_rown = 1000  # 1000
    screen_data = pd.read_csv(screen_data_path)
    gene_n = screen_data.shape[0]

    # 固定参数
    unfold_gene_comb_dir = data_dir + "/unfold_gene_combinations"
    save_minimal_union_dir = data_dir + "/minimal_union_%s_%s" % (file0_first_gene_combn,
                                                                  file1_first_gene_combn)
    mkdir(save_minimal_union_dir)
    # read & save
    read_unfold_gene_comb_format_path = unfold_gene_comb_dir + "/gene_combination_for_%s_%s.log"
    save_union_path_format = save_minimal_union_dir + "/union_gene_combinations_for_idn%s_combn%s.log"
    save_stat_data_path = save_minimal_union_dir + "/statistics_union_gene_combinations.xlsx"
    # 确定 central neuron id 的 union order
    unfold_data_path = unfold_gene_comb_dir + "/unfold_statistics_files_summary.xlsx"
    unfold_data = pd.read_excel(unfold_data_path)
    iterindex = unfold_data['central neuron id'].tolist()
    # union gene combinations of central neurons
    print("==============================================")
    print("union gene combinations of central neurons")
    union_gene_combinations_central_neurons(iterindex, gene_n, read_unfold_gene_comb_format_path,
                                            save_union_path_format, save_stat_data_path,
                                            file0_first_gene_combn, file0_max_readrow_n,
                                            file1_first_gene_combn, file1_max_readrow_n,
                                            single_csr_max_rown, max_gene_combn, 3)

