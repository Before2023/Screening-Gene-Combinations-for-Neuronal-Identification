import os
import numpy as np
import pandas as pd
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


# to get neurons and ganglia info. from position based on Ganglia
def get_neurons_from_position(data_neup, data_ganglia, same_data=False):
    '''
    neuropal_path = "./data from ltk/position_7.csv"
    ganglia_name_path="./ganglia_neuron_name.xlsx"
    data_neup = pd.read_csv(neuropal_path)
    data_ganglia = pd.read_excel(ganglia_name_path)
    data_neurons, ganglions = get_neurons_from_position(data_neup, data_ganglia)
    '''
    if not same_data:
        data_neup.rename(columns={'Neuron': 'Neurons'}, inplace=True)
        data_neurons = pd.merge(data_ganglia, data_neup, how='inner', on='Neurons')
        data_neurons.dropna(subset=['Ganglia', 'Neurons'], axis=0, inplace=True)
        ganglions = data_neurons[['Ganglia', 'Neurons']].groupby(by='Ganglia').count()
        ganglions.sort_values(by='Neurons', ascending=False, inplace=True)
        # print(data_ganglia.shape, data_neurons.shape)
        return data_neurons, ganglions
    else:
        data_neup.dropna(subset=['Ganglia', 'Neurons'], axis=0, inplace=True)
        ganglions = data_neup[['Ganglia', 'Neurons']].groupby(by='Ganglia').count()
        ganglions.sort_values(by='Neurons', ascending=False, inplace=True)
        return data_neup, ganglions


# only for asymmetrical side
def get_asymmetrical_neurons(data_neurons, side='left', only_in_ganglia=True):
    import copy
    data_copy = copy.deepcopy(data_neurons)
    if only_in_ganglia:
        sort_columns = ['Ganglia', 'Neuron abbr. from Seq', 'Neurons']
    else:
        sort_columns = ['Neuron abbr. from Seq', 'Neurons']
    data_copy.sort_values(by=sort_columns, ascending=True, inplace=True)
    if side == 'left':
        data_copy.drop_duplicates(subset=sort_columns[:-1], keep='first', inplace=True)
    elif side == 'right':
        data_copy.drop_duplicates(subset=sort_columns[:-1], keep='last', inplace=True)
    else:
        pass
    data_copy.reset_index(drop=True, inplace=True)
    return data_copy


# compute spatial distance between two dataset points
def array_from_columns(data, x='x', y='y', z='z'):
    data_array = np.array(data[[x, y, z]])
    return data_array


def compute_spatial_distance(data1_array, data2_array):
    from scipy.spatial.distance import cdist
    dist_array = cdist(data1_array, data2_array, metric='euclidean')
    return dist_array


def distance_between_ganglions(data, ganglia_name_pair, same_ganglia=False):
    '''
    # ganglion_name1 = "Anterior Pharyngeal Bulb"
    # ganglion_name2 = "Posterior Pharyngeal Bulb"
    ganglion_name1 = "Lateral Ganglion"
    ganglion_name2 = "Lateral Ganglion"
    # the same ganglia
    ganglia_name_pair = (ganglion_name1, ganglion_name2)
    dist_data, min_dist_dict = distance_between_ganglions(data, ganglia_name_pair, same_ganglia=True)
    '''
    ganglion_name1, ganglion_name2 = ganglia_name_pair
    datag1 = data.loc[data['Ganglia'] == ganglion_name1, :]
    datag2 = data.loc[data['Ganglia'] == ganglion_name2, :]
    datag1.reset_index(drop=True, inplace=True)
    datag2.reset_index(drop=True, inplace=True)
    datag1_array = array_from_columns(datag1)
    datag2_array = array_from_columns(datag2)
    # compute dist
    dist_array = compute_spatial_distance(datag1_array, datag2_array)
    # add row and columns
    dist_data = pd.DataFrame(dist_array)
    dist_data.index = datag1['Neurons']
    dist_data.columns = datag2['Neurons']
    if not same_ganglia:
        pass
    else:
        dist_array = np.array(dist_data.replace(0, np.inf))
    # the minimum of distances between ganglions
    row_min_index = np.argmin(dist_array, axis=0)
    min_points = {}
    for index in range(len(row_min_index)):
        min_p = dist_array[row_min_index[index], index]
        neuron1 = datag1['Neurons'][row_min_index[index]]
        neuron2 = datag2['Neurons'][index]
        if min_p not in min_points:
            min_points[min_p] = [(neuron1, neuron2)]
        else:
            min_points[min_p].append((neuron1, neuron2))
    #
    min_dist = np.min(list(min_points.keys()))
    min_dist_neuron_comp = min_points[min_dist]
    min_dist_dict = {min_dist: min_dist_neuron_comp}
    return dist_data, min_dist_dict


# get distance set between ganglions
def get_distance_between_ganglia(data):
    data.dropna(subset=['Ganglia', 'Neurons'], axis=0, inplace=True)
    ganglia_cdist_data = pd.DataFrame()
    ganglions = list(data['Ganglia'].unique())
    for ganglia_pair in combinations(ganglions, 2):
        #         print(ganglia_pair)
        dist_data, min_dist_dict = distance_between_ganglions(data, ganglia_pair, same_ganglia=False)
        temp_data = pd.DataFrame({'Ganglia_1': [ganglia_pair[0]],
                                  'Gangali_2': [ganglia_pair[1]],
                                  'Distance': [list(min_dist_dict.keys())[0]],
                                  'From neuron pair': [str(list(min_dist_dict.values())[0])]
                                  })
        # concat
        ganglia_cdist_data = pd.concat([ganglia_cdist_data, temp_data], axis=0)
    # sort & reset index
    ganglia_cdist_data.sort_values(by='Distance', ascending=True, inplace=True)
    ganglia_cdist_data.reset_index(drop=True, inplace=True)
    return ganglia_cdist_data


# list sort return string
def list_sort(a_list):
    a_list.sort(reverse=False)
    return " vs ".join(a_list)


# get distance set between neurons in ganglia
# only_in_ganglia == True, meaning that only computing distances of neurons in the same ganglia
# among multiplen ganglions
def get_distance_between_neurons(data, only_in_ganglia=True):
    data.dropna(subset=['Ganglia', 'Neurons'], axis=0, inplace=True)
    dist_data_info = {}
    if only_in_ganglia:
        for ganglia in data['Ganglia'].unique():
            ganglia_pair = (ganglia, ganglia)
            dist_data, min_dist_dict = distance_between_ganglions(data, ganglia_pair, same_ganglia=False)
            dist_data_info[ganglia_pair] = dist_data
    else:
        import copy
        temp_data = copy.deepcopy(data)
        temp_data['Ganglia'] = 'All ganglions'
        ganglia_pair = ('All ganglions', 'All ganglions')
        dist_data, min_dist_dict = distance_between_ganglions(temp_data, ganglia_pair, same_ganglia=True)
        dist_data_info[ganglia_pair] = dist_data
    # parse
    dist_data_list = []
    for ganglia_pair, dist_data in dist_data_info.items():
        # print(ganglia_pair, dist_data.shape)
        for col in dist_data.columns:
            temp_data = dist_data[[col]]
            temp_data['Ganglia_1'] = ganglia_pair[0]
            temp_data['Ganglia_2'] = ganglia_pair[1]
            temp_data['Neuron_1'] = temp_data.index
            temp_data['Neuron_2'] = col
            temp_data.rename(columns={col: 'Distance'}, inplace=True)
            temp_data = temp_data[['Ganglia_1', 'Ganglia_2',
                                   'Neuron_1', 'Neuron_2',
                                   'Distance'
                                   ]]
            dist_data_list.append(temp_data)
    # concat
    dist_data = pd.concat(dist_data_list, axis=0)
    # drop same neurons
    dist_data.reset_index(drop=True, inplace=True)
    dist_data = dist_data.loc[dist_data['Neuron_1'] != dist_data['Neuron_2'], :]
    # drop neuron_pair duplicates
    if dist_data.shape[0] != 0:
        dist_data['neuron_pair'] = dist_data.apply(lambda row: list_sort([row['Neuron_1'], row['Neuron_2']]),
                                                   axis=1)
        dist_data.drop_duplicates(subset=['Distance', 'neuron_pair'], keep='first', inplace=True)
        del dist_data['neuron_pair']
        dist_data.sort_values(by=['Distance', 'Ganglia_1', 'Neuron_1'], ascending=True, inplace=True)
        dist_data.reset_index(drop=True, inplace=True)
    else:  # dist_data 无值
        pass
    return dist_data


# compute neuron's diameter
def get_neuron_diameter_table(data):
    all_unique_neurons = list(set(data['Neuron_1'].tolist() + data['Neuron_2'].tolist()))
    index = []
    for neuron in all_unique_neurons:
        min_index = data.loc[(data['Neuron_1'] == neuron) | (data['Neuron_2'] == neuron), :]['Distance'].idxmin()
        index.append(min_index)
    data_table = data.loc[index, :]
    data_table.drop_duplicates(inplace=True)
    data_table.sort_values(by='Distance', ascending=True, inplace=True)
    data_table.reset_index(drop=True, inplace=True)
    return data_table


# to compute angle between two vectors
def angle2(v1, v2):
    x = np.array(v1)
    y = np.array(v2)
    # 分别计算两个向量的模：
    module_x = np.sqrt(x.dot(x))
    module_y = np.sqrt(y.dot(y))
    # 计算两个向量的点积
    dot_value = x.dot(y)
    # 计算夹角的cos值：
    cos_theta = dot_value / (module_x * module_y)
    # 求得夹角（弧度制）：
    angle_radian = np.arccos(cos_theta)
    # 转换为角度值：
    angle_value = angle_radian * 180 / np.pi
    return angle_value


# compute spatial vector
def compute_spatial_vector(s_point, e_point):
    s_point = np.array(s_point)
    e_point = np.array(e_point)
    return e_point - s_point


def find_tb_neurons(central_data, candid_table, Nmin=25, sigma=0.5, sita=60):
    # for central data
    central_neuron = central_data.loc[0, 'Central_neuron']
    coord_central = np.array(central_data.loc[0, ['x', 'y', 'z']])
    # for other candidate neurons
    # To delete neurons with more distance than Nmin
    drop_neuron_table1 = candid_table.loc[candid_table['Distance'] > Nmin, :]
    drop_neuron_table1['Description'] = 'More than Nmin'
    candid_table1 = candid_table.loc[candid_table['Distance'] <= Nmin, :]
    candid_table1.sort_values(by=['Distance'], ascending=True, inplace=True)
    candid_table1.reset_index(drop=True, inplace=True)
    # to iterate candidate neurons
    if candid_table1.shape[0] > 0:
        # The first boundary neuron
        bound_neuron = candid_table1.loc[0, 'Neurons']
        bound_dist = candid_table1.loc[0, 'Distance']
        # equal neurons with bound_neuron using the parameter sigma
        bound_neurons_index = candid_table1.loc[candid_table1['Distance'] <= bound_dist + sigma, :].index.tolist()
        # print("bound_neurons_index:", bound_neurons_index)
        # to filter other candidate neurons using the angle occupied by topological boundary neurons
        left_index = list(set(candid_table1.index.tolist()) - set(bound_neurons_index))
        drop_index = []
        for index in left_index:
            # to judge candidate neurons not in the regions of topological neurons
            for bound_index in bound_neurons_index:
                coord_bound = np.array(candid_table1.loc[bound_index, ['x', 'y', 'z']])
                coord_candid = np.array(candid_table1.loc[index, ['x', 'y', 'z']])
                v_central2bound = compute_spatial_vector(coord_central, coord_bound)
                v_central2candid = compute_spatial_vector(coord_central, coord_candid)
                angle = angle2(v_central2bound, v_central2candid)
                if angle > sita / 2:  # keep candidate neurons outside the regions of topological neurons
                    pass
                else:
                    drop_index.append(index)
                    break
        next_candid_neuron_index = list(set(left_index) - set(drop_index))
        next_candid_table = candid_table1.loc[next_candid_neuron_index, :]
        bound_neuron_table = candid_table1.loc[bound_neurons_index, :]
        bound_neuron_table['Description'] = "Topological Boundary Neuron"
        drop_neuron_table2 = candid_table1.loc[drop_index, :]
        drop_neuron_table2['Description'] = 'Inside regions of the topological neurons'
        drop_neuron_table = pd.concat([drop_neuron_table1, drop_neuron_table2], axis=0)
        if next_candid_table.shape[0] > 0:
            status = 0
        else:
            status = 1
    else:
        status = 1
        next_candid_table = candid_table1
        bound_neuron_table = candid_table1
        bound_neuron_table['Description'] = []
        drop_neuron_table = drop_neuron_table1
    done_neuron_table = pd.concat([bound_neuron_table, drop_neuron_table], axis=0)
    # reset index
    next_candid_table.reset_index(drop=True, inplace=True)
    done_neuron_table.reset_index(drop=True, inplace=True)
    return status, next_candid_table, done_neuron_table


# main function: to get topological boundary of a neuron
def main_get_topological_boundary(central_neuron, data, neuron_cdist_data,
                                  Nmin=25, sigma=0.5, sita=60):
    # Initail candidate neuron table for the central Neuron
    central_data = data.loc[data['Neurons'] == central_neuron, :]
    central_data.reset_index(drop=True, inplace=True)
    candid_table = neuron_cdist_data.loc[(neuron_cdist_data['Neuron_1'] == central_neuron) |
                                         (neuron_cdist_data['Neuron_2'] == central_neuron), :]
    candid_table['Neuron_2'] = candid_table.apply(
        lambda row: row['Neuron_1'] if row['Neuron_1'] != central_neuron else row['Neuron_2'],
        axis=1)
    candid_table['Neuron_1'] = central_neuron
    candid_table.rename(columns={'Neuron_1': 'Central_neuron',
                                 'Neuron_2': 'Neurons'}, inplace=True)
    central_data.rename(columns={'Neurons': 'Central_neuron'}, inplace=True)
    candid_table = pd.merge(candid_table, data[['Neurons', 'x', 'y', 'z']], how='left', on=['Neurons'])
    # iterate candidate neurons
    done_neuron_table_list = []
    r = 0
    status = 0
    while status == 0:
        r += 1
        status, candid_table, done_neuron_table = find_tb_neurons(central_data, candid_table,
                                                                  Nmin=Nmin, sigma=sigma, sita=sita)
        done_neuron_table_list.append(done_neuron_table)
    #
    done_neuron_table = pd.concat(done_neuron_table_list, axis=0)
    done_neuron_table.sort_values(by='Description', ascending=False, inplace=True)
    done_neuron_table.reset_index(drop=True, inplace=True)
    return done_neuron_table


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


# Firstly, to identify topological boundary of a Neuron
# update functin: 先计算距离后根据 L/R，D/V 神经元类型去重
def get_ganglia_neurons_info(data_neup, data_ganglia, neuron_asymmetry='left', only_in_ganglia=True,
                             Nmin=25, sigma=0.5, sita=60):
    data, ganglions = get_neurons_from_position(data_neup, data_ganglia)
    # update
    data = get_asymmetrical_neurons(data, side=neuron_asymmetry,
                                    only_in_ganglia=only_in_ganglia)  # only one side neurons
    # select one ganglia
    neuron_cdist_data = get_distance_between_neurons(data, only_in_ganglia=only_in_ganglia)  # compute dist bet neurons
    diam_table = get_neuron_diameter_table(neuron_cdist_data)  # compute neuron diameter
    # TO GET TOPOLOGICAL BOUNDARY OF A NEURON
    # For a central Neuron
    bound_neuron_table_list = []
    for ganglia in diam_table['Ganglia_1'].unique():
        temp_diam_table = diam_table.loc[diam_table['Ganglia_1'] == ganglia, :]
        temp_diam_table.reset_index(drop=True, inplace=True)
        temp_all_unique_neuron_list = []
        for index, row in temp_diam_table.iterrows():
            neuron1, neuron2 = row['Neuron_1'], row['Neuron_2']
            if neuron1 not in temp_all_unique_neuron_list:
                temp_all_unique_neuron_list.append(neuron1)
            else:
                pass
            if neuron2 not in temp_all_unique_neuron_list:
                temp_all_unique_neuron_list.append(neuron2)
            else:
                pass
        for central_neuron in temp_all_unique_neuron_list:
            done_neuron_table = main_get_topological_boundary(central_neuron, data, neuron_cdist_data,
                                                              Nmin=Nmin, sigma=sigma, sita=sita)
            bound_neuron_table_list.append(done_neuron_table)
    # concat
    if len(bound_neuron_table_list) != 0:
        bound_neuron_table = pd.concat(bound_neuron_table_list, axis=0)
        bound_neuron_table.reset_index(drop=True, inplace=True)
    else:
        bound_neuron_table = pd.DataFrame({'Ganglia_1': [],
                                           'Ganglia_2': [],
                                           'Central_neuron': [],
                                           'Neurons': [],
                                           'Distance': [],
                                           'x': [],
                                           'y': [],
                                           'z': [],
                                           'Description': []})
    return data, bound_neuron_table, diam_table


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


# to sum up for obtaining neuron sets and statistics of common neurons between NeuronPALs based on Nmin
# get neuron set based on Nmin for NeuroPAL different datasets
def get_neuron_sets_based_on_Nmin(data_ganglia, neuropal_path_format, neuropal_codes,
                                  neuron_asymmetry="both", only_in_ganglia=False,
                                  Nmin=15, sigma=0.5, sita=5):
    backup_data_dict = {}
    for nc in neuropal_codes:
        backup_data_dict['central_neuron(%s)' % nc] = []
        backup_data_dict['TB_neurons(%s)' % nc] = []
    NeuroPAL_neuronset_dict = {}
    max_cneruons_n = 0
    for nc in neuropal_codes:
        neuropal_path = neuropal_path_format % (nc)
        data_neup = pd.read_csv(neuropal_path)
        # 所有神经元搜索，不区分神经节来源
        data, bound_neuron_table, diam_table = get_ganglia_neurons_info(data_neup, data_ganglia,
                                                                        neuron_asymmetry=neuron_asymmetry,
                                                                        only_in_ganglia=only_in_ganglia,
                                                                        Nmin=Nmin, sigma=sigma, sita=sita)
        # Get gene screening data based on the seq data
        neuron_name_dict = {}
        for i, row in data.iterrows():
            neuron_name_dict[row['Neurons']] = row['Neuron abbr. from Seq']
        ganglia_iter_info = {}  # {ganglia: [TB_neuron_set with dropping duplicates]}
        neurons_ited_info = {}  # {central neuron: [all TB_neurons]}
        for ganglia in bound_neuron_table['Ganglia_1'].unique():
            gang_data = bound_neuron_table.loc[bound_neuron_table['Ganglia_1'] == ganglia, :]
            ganglia_iter_info[ganglia] = []
            for central_neuron in gang_data['Central_neuron'].unique():
                neurons = gang_data.loc[(gang_data['Central_neuron'] == central_neuron) &
                                        (gang_data['Description'] == "Topological Boundary Neuron"), :
                          ]['Neurons'].tolist()
                neurons = [neuron_name_dict[x] for x in neurons] + [neuron_name_dict[central_neuron]]
                neurons = list(set(neurons))
                neurons.sort(reverse=False)
                neurons_str = str(neurons)
                ganglia_iter_info[ganglia].append(neurons_str)
                if neurons_str not in neurons_ited_info:
                    neurons_ited_info[neurons_str] = [central_neuron]
                else:
                    neurons_ited_info[neurons_str].append(central_neuron)
        # drop duplicates in ganglia_iter_info (to reduce computation)
        nganglia_iter_info = {}
        for ganglia, neurons_str_list in ganglia_iter_info.items():
            ddp_v = list(set(neurons_str_list))
            for neurons_str in ddp_v:
                central_neurons_tuple = tuple(neurons_ited_info[neurons_str])
                nganglia_iter_info[(ganglia, central_neurons_tuple)] = eval(neurons_str)
        ganglia_iter_info = nganglia_iter_info
        NeuroPAL_neuronset_dict[nc] = ganglia_iter_info
        # save backup
        for k, v in ganglia_iter_info.items():
            backup_data_dict['central_neuron(%s)' % nc].append(str(k))
            backup_data_dict['TB_neurons(%s)' % nc].append(str(v))
        max_cneruons_n = max([max_cneruons_n, len(ganglia_iter_info)])
    #
    for k, v in backup_data_dict.items():
        backup_data_dict[k] = v + [np.nan] * (max_cneruons_n - len(v))
    backup_data = pd.DataFrame(backup_data_dict)
    # get commone neurons between NeuroPAL position()*10
    data_neurp0 = pd.read_csv(neuropal_path_format % neuropal_codes[0])
    data_neurp0.dropna(subset=['Neuron'], axis=0, inplace=True)
    CNeurons = list(data_neurp0['Neuron'].unique())
    print(neuropal_codes[0], len(CNeurons))
    for nc in neuropal_codes[1:]:
        temp = pd.read_csv(neuropal_path_format % nc)
        temp.dropna(subset=['Neuron'], axis=0, inplace=True)
        Neurons = list(temp['Neuron'].unique())
        CNeurons = list(set(CNeurons) & set(Neurons))
        print(nc, len(Neurons))
    print("Number of common neurons: ", len(CNeurons))
    data_ganglia.dropna(subset=['Ganglia', 'Neurons'], axis=0, inplace=True)
    neuron_name_dict = {}
    for i, row in data_ganglia.iterrows():
        neuron_name_dict[row['Neurons']] = row['Neuron abbr. from Seq']
    CNeurons = list(set([neuron_name_dict[x] for x in CNeurons]))
    # print("using Neuron Seq name (num. of CNeurons):", len(CNeurons))
    return NeuroPAL_neuronset_dict, CNeurons, backup_data


