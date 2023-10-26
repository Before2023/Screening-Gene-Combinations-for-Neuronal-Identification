import os
import pandas as pd
import warnings
warnings.filterwarnings("ignore")


def sum_dict(df):
    adict = {}
    for i in df.index:
        temp_dict = eval(df.loc[i, 'gene combination info.'])
        for k, v in temp_dict.items():
            adict[k] = v
    # del df['central neuron id']
    df['gene combination info.'] = str(adict)
    df.drop_duplicates(inplace=True)
    return df


def sum_unfold_statistics_files(data_dir, pool_n):
    result_list = []
    for pool_i in range(1, (pool_n + 1), 1):
        data_path = data_dir + "/0statistics_neuron_sets_gene_combinations_pool_%s.xlsx" % pool_i
        data = pd.read_excel(data_path)
        data['central neuron id'] = data['data name'].apply(lambda x: int(x.split("gene_combinations")[-1].split("_central_neurons")[0]))
        result = data[['central neuron id', 'gene combination info.']].groupby('central neuron id').apply(lambda df:
                                                                                                          sum_dict(df))
        if "central neuron id" in result.columns:
            result.reset_index(drop=True, inplace=True)
        else:
            result.reset_index(drop=False, inplace=True)
        result = result[['central neuron id', 'gene combination info.']]
        result_list.append(result)
    result = pd.concat(result_list, axis=0)
    result.reset_index(drop=True, inplace=True)
    result['minimal gene number'] = result['gene combination info.'].apply(lambda x: min(list(eval(x).keys())))
    result['combn of minimal gene comb.'] = result.apply(lambda row: eval(row['gene combination info.'])[row['minimal gene number']], axis=1)
    result['total combn'] = result['gene combination info.'].apply(lambda x: sum(list(eval(x).values())))
    result.sort_values(by='total combn', ascending=True, inplace=True)
    result.reset_index(drop=True, inplace=True)
    # save
    save_path = data_dir + '/unfold_statistics_files_summary.xlsx'
    result.to_excel(save_path, index=False)
    return result


if __name__ == '__main__':
    import sys
    main_dir, data_dir, pool_n = sys.argv[1:]

    os.chdir(main_dir)
    pool_n = int(pool_n)
    # 固定参数
    unfold_gene_comb_dir = data_dir + "/unfold_gene_combinations"
    sum_unfold_statistics_files(unfold_gene_comb_dir, pool_n)
