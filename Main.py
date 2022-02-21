from FileReader import *

filepath_positions = "C:/Users/mabed/Documents/Travail/Etudes_techniques/fichiers_travail/positions_correspondance.txt"
filepath_reads = "C:/Users/mabed/Documents/Travail/Etudes_techniques/fichiers_travail/reads_statistics.txt"
filepath_nucleotypes = "C:/Users/mabed/Documents/Travail/Etudes_techniques/fichiers_travail/nucleotypes.txt"
filepath_mixtures = "C:/Users/mabed/Documents/Travail/Etudes_techniques/fichiers_travail/simulated_mixtures_composition.txt"

data_utils = build_data_utils(filepath_positions, filepath_reads, filepath_nucleotypes, filepath_mixtures)
mixtures_dict = data_utils[0]
nucleotypes = data_utils[1]
column_nucleotypes_dict = data_utils[2]
positions_errors = data_utils[3]
harp_dict = data_utils[4]

mixtures_list = get_mixtures_list(mixtures_dict)

# for mixture_name in mixtures_list:
#     G = generate_G_from_mix(mixture_name, mixtures_dict, nucleotypes, column_nucleotypes_dict)
#     print(G[0:10])
#
#     reads = reads_of_mix(mixture_name, harp_dict, positions_errors)
#     print(reads)

G_test = generate_G_from_mix('Tm1001', mixtures_dict, nucleotypes, column_nucleotypes_dict)
G_T_test = G_test.transpose()
reads_test = reads_of_mix('Tm1001', harp_dict, positions_errors)
