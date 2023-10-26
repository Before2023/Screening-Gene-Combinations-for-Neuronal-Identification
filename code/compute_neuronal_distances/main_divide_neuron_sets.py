import os
from union_neuron_sets_among_neuropals import *


class Parameter:
    neuropal_codes = [7, 11, 20, 38, 55, 56, 62, 64, 70, 76]
    neuropal_path_format = "./data/NeuroPAL/position_%s.csv"
    ganglia_name_path = "./data/ganglia_neuron_name.xlsx"
    save_excel_dir = "./result/example-analyzed_data/confirm_Nmin_between_NeuroPALs"
    save_final_dir_format = save_excel_dir + "/neuron_sets_and_recombination-%s/"

    def __init__(self, main_dir, Nmin_list, sigma, sita, neuron_asymmetry, only_in_ganglia, neurons_source):
        self.main_dir = main_dir
        self.neuron_asymmetry = neuron_asymmetry  # or 'left', 'right' for distinguishing L/R or D/V ...
        self.only_in_ganglia = only_in_ganglia  # Whether to divide the neuron set according to ganglia
        self.Nmin_list = Nmin_list  # The maximum of distance between neurons
        self.sigma = sigma  # the minimum of distance between equal status neurons
        self.sita = sita  # The angle occupied by a topological boundary neuron
        self.neurons_source = neurons_source  # to select final neuron sets from NeuroPAL union or only one position
        self.save_final_dir = self.save_final_dir_format % neurons_source  # to save final neuron sets

    def build_folder(self):
        os.chdir(self.main_dir)
        mkdir(self.save_excel_dir)
        mkdir(self.save_final_dir)


if __name__ == '__main__':
    main_dir = "../../"
    print(os.getcwd())
    Nmin_list = [6, 20]  # range(6, 30, 1) # The maximum of distance between neurons
    sigma = 0.5   # the minimum of distance between equal status neurons
    sita = 1  # The angle occupied by a topological boundary neuron
    neuron_asymmetry = 'Both'  # or 'left', 'right' for distinguishing L/R or D/V ...
    only_in_ganglia = False  # Whether to divide the neuron set according to ganglia
    neurons_source = 'union'  # to select final neuron sets from NeuroPAL union or only one position

    # GET PARAMETERS
    parameters = Parameter(main_dir, Nmin_list, sigma, sita, neuron_asymmetry, only_in_ganglia, neurons_source)
    parameters.build_folder()
    neuropal_codes = parameters.neuropal_codes
    neuropal_path_format = parameters.neuropal_path_format
    ganglia_name_path = parameters.ganglia_name_path
    save_excel_dir = parameters.save_excel_dir
    save_final_dir = parameters.save_final_dir  # to save final neuron sets

    # execute
    # 1. To get neuron set based on Nmin for each NeuroPAL dataset
    print("\n==================================")
    print("1. To get neuron set based on Nmin for each NeuroPAL dataset")
    data_ganglia = pd.read_excel(ganglia_name_path)
    for Nmin in Nmin_list:
        print("==================\nNmin: %s" % Nmin)
        print("To get neuron set based on Nmin for NeuroPAL different datasets ...")
        NeuroPAL_neuronset_dict, CNeurons, backup_data = get_neuron_sets_based_on_Nmin(data_ganglia,
                                                                                       neuropal_path_format,
                                                                                       neuropal_codes,
                                                                                       neuron_asymmetry=neuron_asymmetry,
                                                                                       only_in_ganglia=only_in_ganglia,
                                                                                       Nmin=Nmin, sigma=sigma, sita=sita)
        sheet_name = "%s_%s_%s_%s_%s" % (Nmin, sigma, sita, neuron_asymmetry, str(only_in_ganglia))
        save_neuron_sets_backup_path = save_excel_dir + "/neuron_sets_based_on_Nmin_for_NeuroPALs-%s.xlsx"%sheet_name
        backup_data.to_excel(save_neuron_sets_backup_path, index=False)

    # 2. To union neuron sets Among NeuroPALs
    print("\n==================================")
    print("2. To union neuron sets Among NeuroPALs")
    backup_data, Nall_cneuron_union_info = get_union_neuron_sets_under_Nmin(data_ganglia, save_excel_dir,
                                                                            neuropal_codes, Nmin_list, sigma, sita,
                                                                            neuron_asymmetry, only_in_ganglia)
    save_neuron_sets_backup_path = save_excel_dir + "/neuron_sets_union_between_NeuroPALs.xlsx"
    backup_data.to_excel(save_neuron_sets_backup_path, index=False)

    # 3. To cross Validation between NeuroPALs for confirming parameter Nmin-sigma-sita
    print("\n==================================")
    print("3. To cross Validation between NeuroPALs for confirming parameter Nmin-sigma-sita")
    stat_data = cross_validation_between_neuronpals(data_ganglia, save_excel_dir, neuropal_codes,
                                                    Nmin_list, sigma, sita,
                                                    neuron_asymmetry, only_in_ganglia)
    save_path = save_excel_dir + "/confirm_and_plot_similarity_between_union_and_fix_under_Nmin.xlsx"
    stat_data.to_excel(save_path, index=False)
    # groupby
    stat_groupby_data = stat_data[['Nmin', 'fix_code', 'similarity']].groupby(['Nmin', 'fix_code']).mean()
    stat_groupby_data.to_excel(save_excel_dir + "/confirm_and_plot_similarity_groupby.xlsx")

    # 4. To recombine the union neuron sets for obtaining final neuron sets
    print("\n==================================")
    print("4. To recombine the union neuron sets for obtaining final neuron sets")
    # 得到 final Neuron sets
    mkdir(save_final_dir)
    for Nmin in Nmin_list:
        fneurons_data = get_recombination_neuron_sets(neurons_source, Nmin, sigma, sita, neuron_asymmetry,
                                                      only_in_ganglia, save_excel_dir, save_final_dir)
    print("Finish.")

