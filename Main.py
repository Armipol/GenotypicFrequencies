from FileReader import *

filepath_positions = "C:/Users/mabed/Documents/Travail/Etudes_techniques/fichiers_travail/positions_correspondance.txt"
filepath_reads = "C:/Users/mabed/Documents/Travail/Etudes_techniques/fichiers_travail/reads_statistics.txt"
filepath_nucleotypes = "C:/Users/mabed/Documents/Travail/Etudes_techniques/fichiers_travail/nucleotypes.txt"
filepath_mixtures = "C:/Users/mabed/Documents/Travail/Etudes_techniques/fichiers_travail/simulated_mixtures_composition.txt"

mixtures_dict = build_mixtures_dictionnary(filepath_mixtures)
mixtures_list = get_mixtures_list(mixtures_dict)

for mixture_name in mixtures_list:
    data = generate_data_for_mix(mixture_name, filepath_positions, filepath_reads, filepath_nucleotypes, filepath_mixtures)
    G = data[0]
    reads = data[1]