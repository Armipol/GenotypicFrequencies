import numpy as np
import csv


#traitement du fichier simulated_mixtures
def build_mixtures_dictionnary(filepath):
    f = open(filepath, 'r')
    mixtures_dict = {}
    #format mixture_dict = {‘mixture_name_1’: mixture_array_1, ‘mixture_name_2’: mixture_array_2, …}
    i = 0
    j = 1
    for line in csv.reader(f):
        if(i == 0):
            i += 1
        else:
            if j <= 12:
                lines = np.array(line)
                elements = lines[0].split(';')
                if j == 1:
                    mixture = []
                    previous_mix = elements[0]
                else:
                    if elements[0] != previous_mix:
                        print("problem : not a 12 genotypes mix")
                        # On vérifie qu'on a bien des mélanges de 12 génotypes à chaque fois
                    previous_mix = elements[0]
                mixture.append(elements[1])
                if j == 12:
                    mixtures_dict[elements[0]] = mixture
                j += 1
            if j == 13:
                j = 1
    print("\nmixtures dict built")
    f.close()
    return mixtures_dict


def encode_nucleotypes_0_1(filepath, nb_genotypes, nb_snips):
    nucleotypes = np.zeros((nb_snips, nb_genotypes))
    f = open(filepath, 'r')
    i = 0
    lines_nucleotypes_dict = {} # lines_dict est le dictionnaire des positions harp de la matrice G selon le numéro de ligne
    pairs_dict = {} # dictionnaire des paires parentes selon la position harp
    column_nucleotypes_dict = {} # dictionnaire des colonnes correspondant à chaque génotype
    count_errors = 0
    positions_errors = [] # liste des positions contenant des erreurs de lecture
    k = 0 # décalage d'indice lié au nombre de lignes (positions) supprimées (car elles contiennent des erreurs de lecture)
    for line in csv.reader(f):
        if i == 0: # traitement de l'en-tête
            nb_genotypes_in_file = len(line) - 2 # pour vérifier qu'il y a bien 96 génotypes
            print("Nbr de génotypes dans le fichier :", nb_genotypes_in_file)
            for j in range(2, len(line)):
                column_nucleotypes_dict[line[j]] = j - 2
            print(column_nucleotypes_dict)
        else:
            lines_nucleotypes_dict[i - 1] = int(line[0])
            ref = line[1]
            nucleotids_already_seen = []
            for j in range(2, len(line)):
                # A = 0, C = 1, G = 2, T = 3, autre : 4
                if line[j] == ref:
                    nucleotypes[i - 1][j - 2] = 1
                    nucleotids_already_seen.append(ref)
                else:
                    if line[j] != 'A' and line[j] != 'C' and line[j] != 'G' and line[j] != 'T':
                        # Il peut y avoir des M, S, Y, W...
                        # on met une valeur à 2 mais on supprimera les lignes (positions) correspondantes ensuite
                        nucleotypes[i - 1][j - 2] = 2
                        count_errors += 1
                        positions_errors.append(int(line[0]))
                    else:
                        nucleotypes[i - 1][j - 2] = 0
                        nucleotids_already_seen.append(line[j])
            nucleotids_already_seen = list(set(nucleotids_already_seen))
            if len(nucleotids_already_seen) >= 3:
                print("more than 2 nucleotids types seen", len(nucleotids_already_seen)) # (sans compter les autres que ACGT)
            pairs_dict[int(line[0])] = nucleotids_already_seen
        i += 1
    positions_errors = list(set(positions_errors))
    tuple_to_delete = ()
    for key in lines_nucleotypes_dict:
        if lines_nucleotypes_dict[key] in positions_errors:
            tuple_to_delete += (key,)
    print("tuple to delete",tuple_to_delete)
    nucleotypes = np.delete(nucleotypes, tuple_to_delete, axis=0)
    print("nb of reading errors in parents nucleotypes", count_errors)
    print("nb of positions that have reading errors", len(positions_errors))
    print("positions errors ", positions_errors)
    print("\nnucleotypes encoded")
    f.close()
    return [nucleotypes, lines_nucleotypes_dict, column_nucleotypes_dict, pairs_dict, count_errors, positions_errors]


def generate_G_from_mix(mixture_name, mixtures_dict, nucleotypes, column_nucleotypes_dict):
    genotypes = mixtures_dict[mixture_name]
    G = np.zeros((len(nucleotypes), len(genotypes)))
    print("len nucleotypes", len(nucleotypes))
    for i in range(len(nucleotypes)):
        for j in range(len(genotypes)):
            #utiliser np.matrix
            G[i][j] = int(nucleotypes[i][column_nucleotypes_dict[genotypes[j]]])
    return G



# CI-DESSOUS FONCTIONS LIEES AUX READS


def reads_statistics_reader(filepath):
    f = open(filepath, 'r')
    dict = {}
    #format dict = {'line_id_1' : line_dict_1, 'line_id_2' : line_dict_2, ...}
    #format line_dict = {'fragment_name': fragment_name, 'position_number': position_number, 'ref_nucleotid': ref_nucleotid, 'p_or_m': p_or_m, 'reads_dict': reads_dict}
    #format reads_dict : {'reads_array_1', 'reads_array_2', ...}
    #format reads_array : [genotype_id, nb_reads, nb_a, nb_c, nb_g, nb_t]
    # tester que somme vaut nbreads
    for line in zip(*[line for line in csv.reader(f)]):
        #print(line)
        lines = np.array(line)
        m = 0
        for i in range(len(lines)):
            line_dict = {'fragment_name': None, 'position_number': None, 'ref_nucleotid': None, 'p_or_m': None, 'reads_dict': None}
            line_array = lines[i].split('\t')
            fragment_name = line_array[0]
            position_number = int(line_array[1])
            ref_nucleotid = line_array[2]
            p_or_m = line_array[3]
            reads_dict = {}
            if p_or_m == 'M':
                k = 0
                for j in range(4, len(line_array)):
                    genotype_id = j - 3
                    nb_reads = int(line_array[j])
                    nb_a = 0
                    nb_c = 0
                    nb_g = 0
                    nb_t = 0
                    if ref_nucleotid == 'A':
                        nb_a = int(line_array[j])
                    if ref_nucleotid == 'C':
                        nb_c = int(line_array[j])
                    if ref_nucleotid == 'G':
                        nb_g = int(line_array[j])
                    if ref_nucleotid == 'T':
                        nb_t = int(line_array[j])
                    reads_array = []
                    reads_array.append(genotype_id)
                    reads_array.append(nb_reads)
                    reads_array.append(nb_a)
                    reads_array.append(nb_c)
                    reads_array.append(nb_g)
                    reads_array.append(nb_t)
                    reads_array_name = 'reads_array_' + str(k)
                    reads_dict[reads_array_name] = reads_array
                    k += 1

            else:
                k = 0
                for j in range(4, len(line_array)):
                    genotype_id = j - 3
                    left_decomposed = line_array[j].split('[')
                    nb_reads = int(left_decomposed[0])
                    right_decomposed = left_decomposed[1].split(']')
                    nucleotids_decomposed = right_decomposed[0].split('/')
                    nb_a = int(nucleotids_decomposed[0])
                    nb_c = int(nucleotids_decomposed[1])
                    nb_g = int(nucleotids_decomposed[2])
                    nb_t = int(nucleotids_decomposed[3])
                    reads_array = []
                    reads_array.append(genotype_id)
                    reads_array.append(nb_reads)
                    reads_array.append(nb_a)
                    reads_array.append(nb_c)
                    reads_array.append(nb_g)
                    reads_array.append(nb_t)
                    reads_array_name = 'reads_array_' + str(k)
                    reads_dict[reads_array_name] = reads_array
                    k += 1

            line_dict['fragment_name'] = fragment_name
            line_dict['position_number'] = position_number
            line_dict['ref_nucleotid'] = ref_nucleotid
            line_dict['p_or_m'] = p_or_m
            line_dict['reads_dict'] = reads_dict

            line_id = 'line_' + str(i)
            dict[line_id] = line_dict

            if k != 96:
                print("problem with number of genotypes for position", position_number)
    f.close()
    return dict

    # il faut utiliser une clé distincte du nom du fragment. Plusieurs snps sont en fait sur le même fragment (2031 lignes font
    # intervenir des fragments déjà apparus avant)

    # nb de lignes total : 5242
    # il y a le même nombre de reads_array pour toutes les lignes : 96

    # vérifier s'il y a des valeurs vides/nulles sur certaines lignes


def build_positions_dict(filepath):
    f = open(filepath, 'r')
    dict = {}
    for line in zip(*[line for line in csv.reader(f)]):
        lines = np.array(line)
        for i in range(1, len(lines)):
            line_array = lines[i].split(' ')
            line_array[0] = line_array[0][1:] # pour enlever le chevron
            if i == 1:
                print(line_array[0])
            line_array[1] = int(line_array[1])
            line_array[2] = int(line_array[2])
            dict[i] = {'fragment_name': line_array[0], 'basic_position': line_array[1], 'harp_position': line_array[2]}

    print("\npositions dictionnary built")
    f.close()
    return dict


def add_harp_positions(filepath_reads, positions_dict):
    n = 0
    dict = reads_statistics_reader(filepath_reads)
    harp_dict = {}
    # format : harp_dict = {'position_harp_1' : line_dict_1, 'position_harp_2' : line_dict_2, ...}
    j = 0
    for key in dict:
        line_dict = dict[key]
        fragment_name = line_dict.get('fragment_name')
        position_number = line_dict.get('position_number') # renommer position_ref
        for pos_key in positions_dict:
            if positions_dict[pos_key]['basic_position'] == position_number and positions_dict[pos_key]['fragment_name'] == fragment_name:
                corresp_harp_position = positions_dict[pos_key]['harp_position']
                harp_dict[corresp_harp_position] = line_dict

        # la position n'est pas une clé car on retrouve plusieurs fois une même position.
        # on doit mettre les numéros de positions à la référence HARP pour comparer avec les nucléotypes
        # il y a 5242 numéros harp
    print("\nharp positions added")
    return harp_dict


def check_sum_is_nbReads(harp_dict):
    i = 0
    for key in harp_dict:
        line_dict = harp_dict[key]
        reads_dict = line_dict.get('reads_dict')
        for reads_key in reads_dict:
            sum = reads_dict[reads_key][2] + reads_dict[reads_key][3] + reads_dict[reads_key][4] + reads_dict[reads_key][5]
            if sum != reads_dict[reads_key][1]:
                print("la somme des reads de nucléotides ne vaut pas le nbr de reads total")
                print("somme:", sum)
                print("nbReads:", reads_dict[reads_key][1])
                i += 1
    print("nombre de sommes différentes:", i)
    return i


def delete_reads_errors(harp_dict, pairs_dict):
    j = 0 # indice de la position harp
    for key in harp_dict:
        line_dict = harp_dict[key]
        reads_dict = line_dict.get('reads_dict')
        i = 0 # indice du mélange de génotypes, càd correspondant à un tableau de reads
        for reads_key in reads_dict:
            reads_array_name = 'reads_array_' + str(i)
            if ('A' not in pairs_dict[key]) and reads_dict[reads_array_name][2] != 0:
                # print("harp", key)
                # print("reads", i)
                # print(pairs_dict[key])
                # print(reads_dict[reads_array_name])
                reads_dict[reads_array_name][1] = reads_dict[reads_array_name][1] - reads_dict[reads_array_name][2]
                reads_dict[reads_array_name][2] = 0
                # print(reads_dict[reads_array_name])
            if ('C' not in pairs_dict[key]) and reads_dict[reads_array_name][3] != 0:
                reads_dict[reads_array_name][1] = reads_dict[reads_array_name][1] - reads_dict[reads_array_name][3]
                reads_dict[reads_array_name][3] = 0
            if ('G' not in pairs_dict[key]) and reads_dict[reads_array_name][4] != 0:
                reads_dict[reads_array_name][1] = reads_dict[reads_array_name][1] - reads_dict[reads_array_name][4]
                reads_dict[reads_array_name][4] = 0
            if ('T' not in pairs_dict[key]) and reads_dict[reads_array_name][5] != 0:
                reads_dict[reads_array_name][1] = reads_dict[reads_array_name][1] - reads_dict[reads_array_name][5]
                reads_dict[reads_array_name][5] = 0
            i += 1
        j += 1
    print("\nreads errors deleted")


def extract_reads_nb(harp_dict):
    reads_array = []
    for key in harp_dict:
        line_dict = harp_dict[key]
        reads_dict = line_dict['reads_dict']
        for reads_key in reads_dict:
            reads_array.append(reads_dict[reads_key][1])
    return reads_array


def reads_of_mix(mixture_name, harp_dict, positions_errors):
    reads = np.zeros((2,len(harp_dict.keys()) - len(positions_errors)))
    reads_array_name = 'reads_array_' + str(int(mixture_name[3:]) - 1)
    i = 0
    for position in harp_dict:
        if position not in positions_errors:
            ref = harp_dict[position]['ref_nucleotid']
            reads_array = harp_dict[position]['reads_dict'][reads_array_name]
            reads[0][i] = reads_array[1]
            if(ref == 'A'):
                reads[1][i] = reads_array[2]
            if (ref == 'C'):
                reads[1][i] = reads_array[3]
            if (ref == 'G'):
                reads[1][i] = reads_array[4]
            if (ref == 'T'):
                reads[1][i] = reads_array[5]
            i += 1

    return reads



#FONCTION OBSOLETE
def obtain_reads_proba(harp_dict, nb_genotypes):
    for key in harp_dict:
        line_dict = harp_dict[key]
        reads_dict = line_dict.get('reads_dict')
        proba_dict = {}
        # format proba_array : [genotype_id, nb_reads, proba_a, proba_c, proba_g, proba_t]
        for i in range(nb_genotypes):
            reads_array_name = 'reads_array_' + str(i)
            proba_array_i = reads_dict.get(reads_array_name).copy()
            if line_dict['p_or_m'] == 'P':
                if proba_array_i[1] != 0:
                    proba_array_i[2] = proba_array_i[2]/proba_array_i[1]
                    proba_array_i[3] = proba_array_i[3]/proba_array_i[1]
                    proba_array_i[4] = proba_array_i[4]/proba_array_i[1]
                    proba_array_i[5] = proba_array_i[5]/proba_array_i[1]
                else:
                    proba_array_i[2] = 0
                    proba_array_i[3] = 0
                    proba_array_i[4] = 0
                    proba_array_i[5] = 0
            else:
                proba_array_i[2] = 0
                proba_array_i[3] = 0
                proba_array_i[4] = 0
                proba_array_i[5] = 0
                if line_dict['ref_nucleotid'] == 'A':
                    proba_array_i[2] = 1
                if line_dict['ref_nucleotid'] == 'C':
                    proba_array_i[3] = 1
                if line_dict['ref_nucleotid'] == 'G':
                    proba_array_i[4] = 1
                if line_dict['ref_nucleotid'] == 'T':
                    proba_array_i[5] = 1
            proba_array_name = 'proba_array_' + str(i)
            proba_dict[proba_array_name] = proba_array_i
        line_dict['proba_dict'] = proba_dict
    return harp_dict


#FONCTION OBSOLETE. Encodage multiple. Pourrait s'avérer utile dans le cas multinomial ?
def encode_nucleotypes(filepath, nb_genotypes, nb_snips):
    G = np.zeros((nb_snips, nb_genotypes))
    f = open(filepath, 'r')
    i = 0
    lines_dict = {} # lines_dict est le dictionnaire des positions harp de la matrice G selon le numéro de ligne
    pairs_dict = {} # dictionnaire des paires parentes selon la position harp
    for line in csv.reader(f):
        if i == 0:
            nb_genotypes_in_file = len(line) - 2 # pour vérifier qu'il y a bien 96 génotypes
            print("Nbr de génotypes dans le fichier :", nb_genotypes_in_file)
        else:
            line_nb = i - 1
            lines_dict[line_nb] = int(line[0])
            # note  : le nucléotide de référence n'est pas exploité pour le moment
            nucleotids_already_seen = []
            for j in range(2, len(line)):
                # A = 0, C = 1, G = 2, T = 3, autre : 4
                if line[j] != 'A' and line[j] != 'C' and line[j] != 'G' and line[j] != 'T':
                    # Il peut y avoir des M, S, Y, W...
                    G[i - 1][j - 2] = 4
                if line[j] == 'A':
                    G[i - 1][j - 2] = 0
                    nucleotids_already_seen.append('A')
                if line[j] == 'C':
                    G[i - 1][j - 2] = 1
                    nucleotids_already_seen.append('C')
                if line[j] == 'G':
                    G[i - 1][j - 2] = 2
                    nucleotids_already_seen.append('G')
                if line[j] == 'T':
                    G[i - 1][j - 2] = 3
                    nucleotids_already_seen.append('T')
            nucleotids_already_seen = list(set(nucleotids_already_seen))
            if len(nucleotids_already_seen) >= 3:
                print("more than 2 nucleotids types seen", len(nucleotids_already_seen)) # (sans compter les autres que ACGT)
            pairs_dict[int(line[0])] = nucleotids_already_seen
        i += 1
    f.close()
    return [G, lines_dict, pairs_dict]


def get_mixtures_list(mixtures_dict):
    list = []
    for key in mixtures_dict:
        list.append(key)
    return list


def build_data_utils(filepath_positions, filepath_reads, filepath_nucleotypes, filepath_mixtures):
    positions_dict = build_positions_dict(filepath_positions)
    harp_dict = add_harp_positions(filepath_reads, positions_dict)

    nucleotypes_data = encode_nucleotypes_0_1(filepath_nucleotypes, 96, 5242)
    pairs_dict = nucleotypes_data[3]

    delete_reads_errors(harp_dict, pairs_dict)

    mixtures_dict = build_mixtures_dictionnary(filepath_mixtures)
    nucleotypes = nucleotypes_data[0]
    column_nucleotypes_dict = nucleotypes_data[2]
    positions_errors = nucleotypes_data[5]

    return [mixtures_dict, nucleotypes, column_nucleotypes_dict, positions_errors, harp_dict]


def check_if_position_is_error(fragment_name, basic_position, positions_dict, positions_errors):
    for key in positions_dict:
        if positions_dict[key]['fragment_name'] == fragment_name and positions_dict[key]['basic_position'] == basic_position:
            harp_position = positions_dict[key]['harp_position']
            print("corresp harp position", harp_position)
            if harp_position in positions_errors:
                print("position is error")
                return True
            else:
                print("position is not error")
                return False

positions_dict = build_positions_dict("D:/Rémi/Documents/IMT/3A/S10/EtudeTech/melange_simul_renom/positions_correspondance.txt")
harp_dict = add_harp_positions("D:/Rémi/Documents/IMT/3A/S10/EtudeTech/melange_simul_renom/SIMULS_MIXTURES_fauxBAM_vraiREADS/reads_statistics.txt", positions_dict)

# check_sum_is_nbReads(harp_dict)
# ATTENTION, on remarque que la somme des tirages pour les nucléotides ACGT n'est pas toujours égale au nbr de reads.
# Retirer le chiffre de nbReads en trop ?

nucleotypes_data = encode_nucleotypes_0_1("D:/Rémi/Documents/IMT/3A/S10/EtudeTech/melange_simul_renom/nucleotypes.txt", 96, 5242)
pairs_dict = nucleotypes_data[3]
positions_errors = nucleotypes_data[5]
