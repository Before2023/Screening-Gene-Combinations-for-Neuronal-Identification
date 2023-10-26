from divide_neuron_sets_by_distance import *


# 取 central neuron 在各数据集中 neuron set 的并
def get_one_ganglia_iter_info(read_dir, fix_code, Nmin, sigma, sita, neuron_asymmetry, only_in_ganglia):
    sheet_name = "%s_%s_%s_%s_%s" % (Nmin, sigma, sita, neuron_asymmetry, str(only_in_ganglia))
    neuron_sets_backup_path = read_dir + "/neuron_sets_based_on_Nmin_for_NeuroPALs-%s.xlsx" % sheet_name
    backup_data = pd.read_excel(neuron_sets_backup_path)
    col1 = 'central_neuron(%s)' % fix_code
    col2 = 'TB_neurons(%s)' % fix_code
    temp = backup_data[[col1, col2]]
    temp.dropna(inplace=True)
    temp[col1] = temp[col1].apply(lambda x: eval(x))
    temp[col2] = temp[col2].apply(lambda x: eval(x))
    temp.index = temp[col1]
    ganglia_iter_info = dict(temp[col2])
    return ganglia_iter_info


def get_union_neuron_sets_under_Nmin(data_ganglia, save_excel_dir, neuropal_codes,
                                     Nmin_list, sigma, sita, neuron_asymmetry, only_in_ganglia):
    backup_data_dict = {}
    for Nmin in Nmin_list:
        backup_data_dict['central_neuron(Nmin%s)' % Nmin] = []
        backup_data_dict['TB_neurons(Nmin%s)' % Nmin] = []
    Nall_cneuron_union_info = {}
    max_cneruons_n = 0
    for Nmin in Nmin_list:
        NeuroPAL_neuronset_dict = {}
        for nc in neuropal_codes:
            ganglia_iter_info = get_one_ganglia_iter_info(save_excel_dir, nc, Nmin, sigma, sita, neuron_asymmetry,
                                                          only_in_ganglia)
            NeuroPAL_neuronset_dict[nc] = ganglia_iter_info
        # union
        all_cneuron_union_info = {}
        all_neurons = data_ganglia['Neurons'].tolist()
        for aneuron in all_neurons:
            all_cneuron_union_info[aneuron] = []
        ganglia = ''  # only for creating a variable
        for aneuron in all_neurons:
            for nc in neuropal_codes:
                for neuron_belong, neurons in NeuroPAL_neuronset_dict[nc].items():
                    ganglia, central_neurons_tuple = neuron_belong
                    if aneuron not in central_neurons_tuple:
                        pass
                    else:
                        all_cneuron_union_info[aneuron] += neurons
            all_cneuron_union_info[aneuron] = list(set(all_cneuron_union_info[aneuron]))
        Nall_cneuron_union_info[Nmin] = all_cneuron_union_info
        # drop duplicates
        ganglia_iter_info = {ganglia: []}  # {ganglia: [TB_neuron_set with dropping duplicates]}
        neurons_ited_info = {}  # {central neuron: [all TB_neurons]}
        for central_neuron, neurons in all_cneuron_union_info.items():
            if len(neurons) == 0:
                pass
            else:
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
        # save backup
        for k, v in ganglia_iter_info.items():
            backup_data_dict['central_neuron(Nmin%s)' % Nmin].append(str(k))
            backup_data_dict['TB_neurons(Nmin%s)' % Nmin].append(str(v))
        max_cneruons_n = max([max_cneruons_n, len(ganglia_iter_info)])
    for k, v in backup_data_dict.items():
        backup_data_dict[k] = v + [np.nan] * (max_cneruons_n - len(v))
    backup_data = pd.DataFrame(backup_data_dict)
    return backup_data, Nall_cneuron_union_info


# 3. To cross Validation between NeuroPALs for confirming parameter Nmin-sigma-sita
def cross_validation_between_neuronpals(data_ganglia, read_dir, neuropal_codes, Nmin_list, sigma, sita,
                                        neuron_asymmetry, only_in_ganglia):
    stat_data_list = []
    for fix_code in neuropal_codes:
        print(fix_code, '... ...')
        left_codes = [x for x in neuropal_codes if x != fix_code]
        backup_data0, Nunion_neuron_sets_info = get_union_neuron_sets_under_Nmin(data_ganglia, read_dir,
                                                                                 left_codes, Nmin_list, sigma, sita,
                                                                                 neuron_asymmetry, only_in_ganglia)
        backup_data1, Na_neuron_sets_info = get_union_neuron_sets_under_Nmin(data_ganglia, read_dir, [fix_code],
                                                                             Nmin_list, sigma, sita, neuron_asymmetry,
                                                                             only_in_ganglia)
        for Nmin in Nmin_list:
            union_neuron_sets_info, a_neuron_sets_info = Nunion_neuron_sets_info[Nmin], Na_neuron_sets_info[Nmin]
            stat_data = compute_similarity_between_NeuroPAL(union_neuron_sets_info, a_neuron_sets_info)
            stat_data['Nmin'] = Nmin
            stat_data['fix_code'] = fix_code
            stat_data_list.append(stat_data)
    stat_data = pd.concat(stat_data_list, axis=0)
    return stat_data


# computing similarity for confirming Nmin
def compute_similarity_between_NeuroPAL(union_neuron_sets_info, a_neuron_sets_info):
    stat_dict = {'central_neuron': [],
                 'similarity': [],
                 'neurons_n_test': [],
                 'neurons_n_union': [],
                 'common_neurons_n': []
                 }
    for cneuron, neurons in a_neuron_sets_info.items():
        if cneuron not in union_neuron_sets_info:
            pass
        else:
            neurons1 = union_neuron_sets_info[cneuron]
            if (len(neurons) < 2) | (len(neurons1) < 2):
                pass
            else:
                comm_neurons = list(set(neurons) & set(neurons1))
                sim = round(len(comm_neurons) / len(neurons), 3)
                neurons_n = len(neurons)
                neurons1_n = len(neurons1)
                comm_n = len(comm_neurons)
                stat_dict['central_neuron'].append(cneuron)
                stat_dict['similarity'].append(sim)
                stat_dict['neurons_n_test'].append(neurons_n)
                stat_dict['neurons_n_union'].append(neurons1_n)
                stat_dict['common_neurons_n'].append(comm_n)
    #
    stat_data = pd.DataFrame(stat_dict)
    return stat_data


# orignal
# 得到 final Neuron sets
def get_neuron_sets(neurons_source, read_dir, Nmin=6, sigma=0.5, sita=1,
                    neuron_asymmetry="Both", only_in_ganglia=False):
    if neurons_source == "union":
        neuron_set_path = read_dir + "/neuron_sets_union_between_NeuroPALs.xlsx"
        if os.path.exists(neuron_set_path):
            col1 = 'central_neuron(Nmin%s)' % Nmin
            col2 = 'TB_neurons(Nmin%s)' % Nmin
        else:
            col1 = ''
            col2 = ''
            print("Neuron sets of Nmin=%s need to be generated by running \nFunciton (get_neuron_sets_based_on_Nmin)")
            print("AND Function (get_union_neuron_sets_under_Nmin)")
    else:
        neuron_set_path = read_dir + "/neuron_sets_based_on_Nmin_for_NeuroPALs-%s_%s_%s_%s_%s.xlsx" % (Nmin, sigma, sita,
                                                                                                        neuron_asymmetry,
                                                                                                        str(only_in_ganglia))
        if os.path.exists(neuron_set_path):
            pos = neurons_source.split("_")[-1]
            col1 = "central_neuron(%s)" % pos
            col2 = "TB_neurons(%s)" % pos
        else:
            col1 = ''
            col2 = ''
            print("Neuron sets of Nmin=%s need to be generated by running \nFunciton (get_neuron_sets_based_on_Nmin)")
    neuron_set_data = pd.read_excel(neuron_set_path)
    return neuron_set_data, col1, col2


## neuron set 的 split & combinations
def get_ganglia_iter_info2(neuron_set_data, col1, col2):
    temp = neuron_set_data[[col1, col2]]
    temp.dropna(inplace=True)
    temp[col1] = temp[col1].apply(lambda x: eval(x))
    temp[col2] = temp[col2].apply(lambda x: eval(x))
    temp.index = temp[col1]
    ganglia_iter_info = dict(temp[col2])
    return ganglia_iter_info


# drop duplicates of neuron sets
def neuron_belong_sum(cols):
    neuron_belongs = list(cols)
    return str(neuron_belongs)


# get recombinated neuron sets
def get_recombination_neuron_sets(neurons_source, Nmin, sigma, sita, neuron_asymmetry, only_in_ganglia,
                                  read_dir, save_dir):
    # print("\nNmin: %s ... ..." % Nmin)
    # after generating Neuron sets
    neuron_set_data, col1, col2 = get_neuron_sets(neurons_source, read_dir, Nmin, sigma, sita, neuron_asymmetry,
                                                  only_in_ganglia)
    ganglia_iter_info = get_ganglia_iter_info2(neuron_set_data, col1, col2)
    #
    params = (Nmin, sigma, sita, neuron_asymmetry, str(only_in_ganglia))
    str_params = ("%s" + "-%s" * 4) % params
    save_stat_data_path = save_dir + "/neuron_sets_original-%s.xlsx" % str_params
    save_final_neurons_path = save_dir + '/Final_neuron_sets-%s.xlsx' % str_params
    if os.path.exists(save_stat_data_path):
        pass
    else:
        print("================================\nNmin:", Nmin)
        stat_dict = {"Neuron belong": [],
                     "Num. T.B. neurons": [],
                     "T.B. neurons": []
                     }
        all_neurons_n = 0  # checkpoint
        rep_neurons = []  # checkpoint
        for neuron_belong, neurons in ganglia_iter_info.items():
            ganglia, central_neurons_tuple = neuron_belong
            all_neurons_n += len(central_neurons_tuple)
            rep_neurons += neurons
            if len(neurons) < 2:
                pass
            else:
                stat_dict['Neuron belong'].append(neuron_belong)
                stat_dict['Num. T.B. neurons'].append(len(neurons))
                stat_dict['T.B. neurons'].append(str(neurons))
        print("all neuron number:", all_neurons_n)
        stat_data = pd.DataFrame(stat_dict)
        # save
        groupby_data = stat_data.groupby(["T.B. neurons", "Num. T.B. neurons"]).agg(neuron_belong_sum)
        groupby_data.reset_index(drop=False, inplace=True)
        cols = ["Neuron belong", "Num. T.B. neurons", "T.B. neurons"]
        stat_data = groupby_data[cols]
        stat_data.to_excel(save_stat_data_path, index=False)
    if os.path.exists(save_final_neurons_path):
        fneurons_data = pd.read_excel(save_final_neurons_path)
    else:
        stat_data = pd.read_excel(save_stat_data_path)
        # groupby_data = stat_data.groupby(["T.B. neurons", "Num. T.B. neurons"]).agg(neuron_belong_sum)
        # groupby_data.reset_index(drop=False, inplace=True)
        # cols = ["Neuron belong", "Num. T.B. neurons", "T.B. neurons"]
        # groupby_data = groupby_data[cols]
        fneuron_dict = {'Neuron belong': [], 'Num. T.B. neurons': [], 'T.B. neurons': []}
        for i, row in stat_data.iterrows():
            neuron_belong = row['Neuron belong']
            neurons = eval(row['T.B. neurons'])
            for neuron_pair in combinations(neurons, 2):
                neuron_pair = list(neuron_pair)
                neuron_pair.sort(reverse=False)
                fneuron_dict['Neuron belong'].append(neuron_belong)
                fneuron_dict['Num. T.B. neurons'].append(len(neuron_pair))
                fneuron_dict['T.B. neurons'].append(str(list(neuron_pair)))
        fneurons_data = pd.DataFrame(fneuron_dict)
        fneurons_data.drop_duplicates(subset=['T.B. neurons'], keep='first', inplace=True)
        fneurons_data.reset_index(drop=True, inplace=True)
        fneurons_data['order'] = fneurons_data.index
        fneurons_data.to_excel(save_final_neurons_path, index=False)
    return fneurons_data

