from divide_neuron_sets_by_distance import *


def compute_and_plot_PR_AUC(data, neuron_or_gene, save_dir):
    # 导入数据模拟模块
    from sklearn.metrics import precision_recall_curve  # 计算精确召回曲线
    from sklearn.metrics import auc

    # 创建ROC曲线
    y_test = np.array(data['reporter']).reshape(-1, 1)
    y_pred = np.array(data['seq']).reshape(-1, 1)
    # 创建PR曲线
    precision, recall, thresholds = precision_recall_curve(y_test, y_pred)
    auc_score = auc(recall, precision)

    # 计算F分数
    fscore = (2 * precision * recall) / (precision + recall)

    # 找到最佳阈值
    index = np.argmax(fscore)
    thresholdOpt = round(thresholds[index], ndigits=4)
    fscoreOpt = round(fscore[index], ndigits=4)
    recallOpt = round(recall[index], ndigits=4)
    precisionOpt = round(precision[index], ndigits=4)
    # print('Best Threshold: {} with F-Score: {}'.format(thresholdOpt, fscoreOpt))
    # print('Recall: {}, Precision: {}'.format(recallOpt, precisionOpt))
    best_threshlod_tuple = (auc_score, thresholdOpt, fscoreOpt, recallOpt, precisionOpt)
    # # 创建数据图
    # from plotnine import *
    # import plotnine
    # 绘制ROC曲线
    # df_recall_precision = pd.DataFrame({'Precision':precision[:-1],
    #                                     'Recall':recall[:-1],
    #                                     'Threshold':thresholds})
    # plotnine.options.figure_size = (8, 4.8)
    # p = (
    #     ggplot(data = df_recall_precision)+
    #     geom_point(aes(x = 'Recall',
    #                    y = 'Precision'),
    #                size = 0.4)+
    #     # 最佳阈值
    #     geom_point(aes(x = recallOpt,
    #                    y = precisionOpt),
    #                color = '#981220',
    #                size = 4)+
    #     geom_line(aes(x = 'Recall',
    #                   y = 'Precision'))+
    #     # 注释
    #     geom_text(aes(x = recallOpt,
    #                   y = precisionOpt),
    #               label = 'Optimal threshold \n for class: {}'.format(thresholdOpt),
    #               nudge_x = 0.18,
    #               nudge_y = 0,
    #               size = 10,
    #               fontstyle = 'italic')+
    #     labs(title = 'Recall Precision Curve {}'.format(round(auc_score, 3)))+
    #     xlab('Recall')+
    #     ylab('Precision')+
    #     theme_minimal()
    # )
    # p.save(filename = '{}/{}_select_reporter_threshold.png'.format(save_dir, neuron_or_gene), units = 'in', dpi=300)
    # p.show()
    return best_threshlod_tuple


def reporter_threshold_selection_based_on_same_neuron(seq_data, prom_data, rep_neurons_info, save_dir):
    best_thv_data_list = []
    for neuron_abb, neuron in rep_neurons_info.items():
        temp_prom = merge_compare_seq_reporter_for_same_neuron(neuron_abb, neuron, seq_data, prom_data)
        count_1 = temp_prom['reporter'].sum()
        count_0 = temp_prom.shape[0] - count_1
        # best threshold selection of reporter expression
        best_threshlod_tuple = compute_and_plot_PR_AUC(temp_prom, "neuron_%s" % neuron_abb, save_dir)
        auc_score, thresholdOpt, fscoreOpt, recallOpt, precisionOpt = best_threshlod_tuple
        temp_best_thv_data = pd.DataFrame({"Neuron": [neuron_abb],
                                           "0_count": [count_0],
                                           "1_count": [count_1],
                                           "PR_AUC": [auc_score],
                                           "Best_threshold": [thresholdOpt],
                                           "F1_score": [fscoreOpt],
                                           "Recall": [recallOpt],
                                           "Precision": [precisionOpt]
                                           })
        best_thv_data_list.append(temp_best_thv_data)
    #     break
    # concat
    best_thv_data = pd.concat(best_thv_data_list, axis=0)
    best_thv_data.reset_index(drop=True)
    return best_thv_data


# reporter threshold seletion based on the same gene
def reporter_threshold_selection_based_on_same_gene(seq_data, prom_data, rep_neurons_info, save_dir):
    seqgene_ids = seq_data['gene_id'].tolist()
    prgene_ids = prom_data['gene_id'].tolist()
    common_gene_ids = list(set(seqgene_ids) & set(prgene_ids))
    best_thv_data_list = []
    for gene_id in common_gene_ids:
        temp_prom = merge_compare_seq_reporter_for_same_gene(gene_id, seq_data, prom_data, rep_neurons_info)
        count_1 = temp_prom['reporter'].sum()
        count_0 = temp_prom.shape[0] - count_1
        # best threshold selection of reporter expression
        best_threshlod_tuple = compute_and_plot_PR_AUC(temp_prom, "gene_%s" % gene_id, save_dir)
        auc_score, thresholdOpt, fscoreOpt, recallOpt, precisionOpt = best_threshlod_tuple
        temp_best_thv_data = pd.DataFrame({"Gene_id": [gene_id],
                                           "0_count": [count_0],
                                           "1_count": [count_1],
                                           "PR_AUC": [auc_score],
                                           "Best_threshold": [thresholdOpt],
                                           "F1_score": [fscoreOpt],
                                           "Recall": [recallOpt],
                                           "Precision": [precisionOpt]
                                           })
        best_thv_data_list.append(temp_best_thv_data)
    # concat
    best_thv_data = pd.concat(best_thv_data_list, axis=0)
    best_thv_data.reset_index(drop=True)
    return best_thv_data


# obtain rna-seq screened data & promoter screened data
# sorting promoter data
def sort_promoter_data(prom_data, gn_data):
    prom_data.index = prom_data.Neuron
    del prom_data['Neuron class']
    del prom_data['Neurotransmitter identity']
    del prom_data['Neuron']
    prom_data = prom_data.T
    neurons = prom_data.columns.tolist()
    prom_data['gene_name'] = prom_data.index
    prom_data = pd.merge(prom_data, gn_data[['gene_id', 'gene_name']], how='inner', on='gene_name')
    prom_data = prom_data[['gene_id', 'gene_name'] + neurons]
    prom_data.reset_index(drop=True, inplace=True)
    del prom_data['gene_name']
    print(prom_data.shape)
    return prom_data


# 同类型神经元 promoter expression 的差异
# check difference of promoter expression between Left and right of a neuron
def get_neuron_LR_info(nn_data):
    neuron_abbr_dict = {}
    for index, row in nn_data.iterrows():
        seq_neuron_name = row['Neuron abbr. from Seq']
        pal_neuron_name = row['Neurons']
        if seq_neuron_name in neuron_abbr_dict:
            neuron_abbr_dict[seq_neuron_name].append(pal_neuron_name)
        else:
            neuron_abbr_dict[seq_neuron_name] = [pal_neuron_name]
    return neuron_abbr_dict


# number of difference genes
def get_number_diff_genes(array_1, array_2):
    diff_array = array_1 - array_2
    return sum(abs(diff_array))


# check the difference
def check_diff_promoter_expression(promoter_data, neuron_abbr_dict):
    from itertools import combinations

    prom_gene_num = promoter_data.shape[0]
    _DVLR_diff_neurons = []
    rep_neurons_info = {}
    for neuron_abbr, neuron_list in neuron_abbr_dict.items():
        if len(neuron_list) != 1:
            for neuron_1, neuron_2 in combinations(neuron_list, 2):
                if (neuron_1 in promoter_data.columns) & (neuron_2 in promoter_data.columns):
                    temp_data = promoter_data[neuron_1].equals(promoter_data[neuron_2])
                    if temp_data:
                        pass
                    else:
                        diff_num = get_number_diff_genes(promoter_data[neuron_1], promoter_data[neuron_2])
                        _DVLR_diff_neurons.append(
                            {neuron_abbr: (neuron_1, neuron_2, "%s/%s" % (diff_num, prom_gene_num))})
                else:
                    pass
        elif len(neuron_list) == 1:
            pass
        else:
            raise ("Error: %s No corresponding Neurons." % neuron_abbr)
        for neuron in neuron_list:
            if neuron in promoter_data.columns:
                rep_neurons_info[neuron_abbr] = neuron
                break
            else:
                pass
    # statistics
    # differences of represent neurons
    rep_neurons_diff = []
    rep_neurons_diff_info = {}
    for n1, n2 in combinations(list(rep_neurons_info.values()), 2):
        _diff = promoter_data[n1].equals(promoter_data[n2])
        if _diff:
            rep_neurons_diff.append(0)
            rep_neurons_diff_info[(n1, n2)] = "0/%s" % prom_gene_num
            pass
        else:
            diff_num = get_number_diff_genes(promoter_data[n1], promoter_data[n2])
            rep_neurons_diff.append(diff_num)
            rep_neurons_diff_info[(n1, n2)] = "%s/%s" % (diff_num, prom_gene_num)
    # statistics for similar neurons
    _L_R_diff = []
    _D_V_diff = []
    _D_LR_diff = []
    _V_LR_diff = []
    _L_R_diff_info = {}
    _D_V_diff_info = {}
    _D_LR_diff_info = {}
    _V_LR_diff_info = {}
    for diff_neuron_dict in _DVLR_diff_neurons:
        neuron_abbr, comb_neurons = list(diff_neuron_dict.keys())[0], list(diff_neuron_dict.values())[0]
        n1, n2, diff_num = comb_neurons
        _n1 = n1.replace(neuron_abbr, '')
        _n2 = n2.replace(neuron_abbr, '')
        if ("D" in _n1) | ("D" in _n2):
            if ("V" in _n1) | ("V" in _n2):
                _D_V_diff.append(int(diff_num.split("/")[0]))
                _D_V_diff_info[(n1, n2)] = diff_num
            else:
                _D_LR_diff.append(int(diff_num.split("/")[0]))
                _D_LR_diff_info[(n1, n2)] = diff_num
        elif ("V" in _n1) | ("V" in _n2):
            _V_LR_diff.append(int(diff_num.split("/")[0]))
            _V_LR_diff_info[(n1, n2)] = diff_num
        else:
            _L_R_diff.append(int(diff_num.split("/")[0]))
            _L_R_diff_info[(n1, n2)] = diff_num
    # DataFrame
    #     print(rep_neurons_diff)
    diff_neurons_stat = pd.DataFrame({"Num. of represent neurons": [len(rep_neurons_diff)],
                                      "Mean of diff. genes among represent neurons": [
                                          "%s/%s" % (int(np.mean(rep_neurons_diff)), prom_gene_num)],
                                      "Median of diff. genes among represent neurons": [
                                          "%s/%s" % (np.median(rep_neurons_diff), prom_gene_num)],
                                      "Num. of _L_R_diff neurons": [len(_L_R_diff)],
                                      "Mean of diff. genes among _L_R_diff neurons": [
                                          "%s/%s" % (int(np.mean(_L_R_diff)), prom_gene_num)],
                                      "Median of diff. genes among _L_R_diff neurons": [
                                          "%s/%s" % (np.median(_L_R_diff), prom_gene_num)],
                                      "Num. of _D_V_diff neurons": [len(_D_V_diff)],
                                      "Mean of diff. genes among _D_V_diff neurons": [
                                          "%s/%s" % (int(np.mean(_D_V_diff)), prom_gene_num)],
                                      "Median of diff. genes among _D_V_diff neurons": [
                                          "%s/%s" % (np.median(_D_V_diff), prom_gene_num)],
                                      "Num. of _D_LR_diff neurons": [len(_D_LR_diff)],
                                      "Mean of diff. genes among _D_LR_diff neurons": [
                                          "%s/%s" % (int(np.mean(_D_LR_diff)), prom_gene_num)],
                                      "Median of diff. genes among _D_LR_diff neurons": [
                                          "%s/%s" % (np.median(_D_LR_diff), prom_gene_num)],
                                      "Num. of _V_LR_diff neurons": [len(_V_LR_diff)],
                                      "Mean of diff. genes among _V_LR_diff neurons": [
                                          "%s/%s" % (int(np.mean(_V_LR_diff)), prom_gene_num)],
                                      "Median of diff. genes among _V_LR_diff neurons": [
                                          "%s/%s" % (np.median(_V_LR_diff), prom_gene_num)],
                                      })
    diff_neuron_genes_info = {"represent neurons": rep_neurons_diff_info,
                              "_L_R diff neurons": _L_R_diff_info,
                              "_D_V_diff neurons": _D_V_diff_info,
                              "_D_LR_diff neurons": _D_LR_diff_info,
                              "_V_LR_diff neurons": _V_LR_diff_info
                              }
    return _DVLR_diff_neurons, rep_neurons_info, diff_neuron_genes_info, diff_neurons_stat


# 选取 RNA-Seq 中与 promoter reporter expression 高度一致的基因做下游分析
# （1）基于同个基因不同神经元之间的 reporter 阈值选择
# （2）RNA-Seq 与 reporter expression 高度一致性
# merge and compare the corresponding relation between seq and reporter for same neuron
def merge_compare_seq_reporter_for_same_neuron(neuron_abb, neuron, seq_data, prom_data):
    seq_temp = seq_data[['gene_id', neuron_abb]]
    pr_temp = prom_data[['gene_id', neuron]]
    seq_temp.rename(columns={neuron_abb: 'seq'}, inplace=True)
    pr_temp.rename(columns={neuron: 'reporter'}, inplace=True)
    temp = pd.merge(seq_temp, pr_temp, how='left', on='gene_id')
    nan_rep = temp.iloc[:, 2].isnull()
    nan_index = nan_rep[nan_rep == True].index
    prom_index = nan_rep[nan_rep == False].index
    temp_nan = temp.loc[nan_index, :]  # only in scRNA-Seq
    temp_prom = temp.loc[prom_index, :]  # reporter 0/1 data
    temp_prom.sort_values(by='reporter', ascending=False, inplace=True)
    temp_prom.reset_index(drop=True, inplace=True)
    return temp_prom


# compare and merge data to check the corresponding relation between seq and reporter for a gene
def merge_compare_seq_reporter_for_same_gene(gene_id, seq_data, prom_data, rep_neurons_info):
    neurons_abb = list(rep_neurons_info.keys())
    reverse_rep_neurons_info = {}
    for neuron_abb, neuron in rep_neurons_info.items():
        reverse_rep_neurons_info[neuron] = neuron_abb
    pr_neurons = list(reverse_rep_neurons_info.keys())
    seq_temp = seq_data.loc[seq_data['gene_id'] == gene_id, neurons_abb]
    pr_temp = prom_data.loc[prom_data['gene_id'] == gene_id, pr_neurons]
    seq_temp.index = ['seq'] * seq_temp.shape[0]
    pr_temp.index = ['reporter'] * pr_temp.shape[0]
    pr_temp.rename(columns=reverse_rep_neurons_info, inplace=True)
    seq_tempT, pr_tempT = seq_temp.T, pr_temp.T
    seq_tempT['Neuron'] = seq_tempT.index
    pr_tempT['Neuron'] = pr_tempT.index
    try:
        seq_tempT = seq_tempT[['Neuron', 'seq']]
    except KeyError as e:
        print("seq_tempT.columns:", seq_tempT.columns)
        print("gene_id:", gene_id)
        raise (e)
    temp_prom = pd.merge(seq_tempT, pr_tempT, how='inner', on='Neuron')
    temp_prom.sort_values(by='reporter', ascending=False, inplace=True)
    temp_prom.reset_index(drop=True, inplace=True)
    return temp_prom


# Adjust best threshold of reporter expression: best threshlod must be non-0
def adjust_0threshlod(thv_data, seq_data, prom_path, gene_name_path, neuron_name_path, pr_auc_thv=0.6):
    adthv_data = thv_data.loc[thv_data['PR_AUC'] >= pr_auc_thv, :]
    print("Num. of genes used for downstream analysis:", adthv_data.shape[0])
    # Difference anlysis of gene expresion between
    gn_data = pd.read_excel(gene_name_path)
    nn_data = pd.read_excel(neuron_name_path)
    prom_data = pd.read_csv(prom_path)
    prom_data = sort_promoter_data(prom_data, gn_data)
    neuron_abbr_dict = get_neuron_LR_info(nn_data)  # neuron_abb: all neurons in promoter dataset
    _DVLR_diff_neurons, rep_neurons_info, diff_neuron_genes_info, diff_neurons_stat = check_diff_promoter_expression(
        prom_data, neuron_abbr_dict)
    #
    adbest_thvs = []
    adthv_data.reset_index(drop=True, inplace=True)
    for index, row in adthv_data.iterrows():
        gene_id = row['Gene_id']
        count_0 = row['0_count']
        count_1 = row['1_count']
        best_thv = row['Best_threshold']
        if best_thv > 0:
            adbest_thvs.append(best_thv)
            pass
        else:
            if 0 in [count_0, count_1]:
                adbest_thvs.append(np.nan)
            else:
                # compare and merge data to check the corresponding relation between seq expression and reporter for a gene
                temp_data = merge_compare_seq_reporter_for_same_gene(gene_id, seq_data, prom_data, rep_neurons_info)
                adbest_thv = temp_data.loc[temp_data['seq'] != 0, :]['seq'].min()
                adbest_thvs.append(adbest_thv)
    # add adjusted best threshold column
    adthv_data['Adjusted_best_threshold'] = adbest_thvs
    adthv_data.index = adthv_data['Gene_id']
    gene_id_info = dict(adthv_data['Adjusted_best_threshold'])
    return adthv_data, gene_id_info


# Get gene screening data based on the seq data
def get_gene_screening_data_from_seq(seq_data, thv_data,
                                     prom_path, gene_name_path, ganglia_name_path,
                                     pr_auc_thv, n_base, neurons):
    try:
        if len(neurons) == len(set(list(seq_data.columns)) & set(neurons)):
            pass
        else:
            raise ("KeyError: there are neurons that not in seq_data. Please check again.")
    except TypeError as e:
        print("KeyError: there are neurons that not in seq_data. Please check again.")
        raise e
    # seq data & thv data
    seq_data1 = seq_data[['gene_id'] + neurons]
    adthv_data, gene_id_info = adjust_0threshlod(thv_data, seq_data,
                                                 prom_path, gene_name_path, ganglia_name_path,
                                                 pr_auc_thv=pr_auc_thv)
    # n base
    seq_data1.index = seq_data1.gene_id
    del seq_data1['gene_id']
    seq_data1T = seq_data1.T
    gene_ids = []  # for next analysis
    for gene_id, adbest_thv in gene_id_info.items():
        if adbest_thv == np.nan:
            pass
        else:
            seq_data1T[gene_id] = seq_data1T[gene_id].apply(
                lambda x: n_base - 1 if (int(x // adbest_thv) >= n_base - 1) else int(x // adbest_thv))
            # seq_data1T[gene_id] = seq_data1T[gene_id].apply(lambda x: int(x//adbest_thv))
            gene_ids.append(gene_id)
    seq_data_nbase = seq_data1T[gene_ids].T
    seq_data_nbase.reset_index(drop=False, inplace=True)
    return seq_data_nbase


# Get gene screening data based on the promoter
def get_gene_screening_data_from_promter(prom_path, gene_name_path, ganglia_name_path, neurons):
    # Difference anlysis of gene expresion between
    gn_data = pd.read_excel(gene_name_path)
    nn_data = pd.read_excel(ganglia_name_path)
    prom_data = pd.read_csv(prom_path)
    prom_data = sort_promoter_data(prom_data, gn_data)
    neuron_abbr_dict = get_neuron_LR_info(nn_data)  # neuron_abb: all neurons in promoter dataset
    _DVLR_diff_neurons, rep_neurons_info, diff_neuron_genes_info, diff_neurons_stat = check_diff_promoter_expression(
        prom_data, neuron_abbr_dict)
    #
    reverse_rep_neurons_info = {}
    for neuron_abb, neuron in rep_neurons_info.items():
        reverse_rep_neurons_info[neuron] = neuron_abb
    pr_neurons = list(reverse_rep_neurons_info.keys())
    prom_data = prom_data[['gene_id'] + pr_neurons]
    prom_data.rename(columns=reverse_rep_neurons_info, inplace=True)
    comm_neurons = list(set(list(prom_data.columns)) & set(neurons))
    prom_data = prom_data[['gene_id'] + comm_neurons]
    return prom_data


# 通过分析得：
# 1. _L & _R: 在 promoter 数据中一致, 除 AWCL & AWCR 之外（存在 6/796 genes 差异）；
# 2. _D & _V: 在 promoter 数据中存在个位数基因表达差异；
# 3. _D & _L/R:在 promoter 数据中存在个位数基因表达差异；
# 4. _V & _L/R:在 promoter 数据中存在个位数基因表达差异；
if __name__ == '__main__':
    main_dir = "../../"
    os.chdir(main_dir)

    # For scRNA + Promoter data
    # read data
    scrna_path = "./data/RNA-Seq/L4_10x_R_data_files/021821_medium_threshold2.csv"
    prom_path = "./data/promoter/cegenelist_head_expression_filter.csv"
    gene_name_path = "./data/sorted_worm_gene_names.xlsx"
    neuron_name_path = "./data/ganglia_neuron_name.xlsx"
    save_dir = "./result/example-analyzed_data"
    mkdir(save_dir)
    scrna_data = pd.read_csv(scrna_path)
    prom_data = pd.read_csv(prom_path)
    gn_data = pd.read_excel(gene_name_path)
    nn_data = pd.read_excel(neuron_name_path)

    # rename scRNA-Seq data
    del scrna_data['Unnamed: 0']
    del scrna_data['gene_name']
    scrna_data.rename(columns={'Wormbase_ID': 'gene_id'}, inplace=True)
    scrna_data.sort_values(by='gene_id', ascending=True, inplace=True)
    scrna_data.reset_index(drop=True, inplace=True)

    # Difference anlysis of gene expresion between
    prom_data = sort_promoter_data(prom_data, gn_data)
    neuron_abbr_dict = get_neuron_LR_info(nn_data)  # neuron_abb: all neurons in promoter dataset
    _DVLR_diff_neurons, rep_neurons_info, diff_neuron_genes_info, diff_neurons_stat = check_diff_promoter_expression(
                                                                                          prom_data, neuron_abbr_dict)
    # compute the threshold
    best_thv_data_bog = reporter_threshold_selection_based_on_same_gene(scrna_data, prom_data, rep_neurons_info,
                                                                        save_dir)
    # save
    best_thv_data_bog.to_excel(save_dir + "/based_on_same_gene_reporter_best_threshold_selection.xlsx", index=False)


