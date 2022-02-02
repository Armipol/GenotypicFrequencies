import numpy as np
import csv


def reads_statistics_reader(filepath):
    f = open(filepath, 'r')
    dict = {}
    #format dict = {'line_id_1' : line_dict_1, 'line_id_2' : line_dict_2, ...}
    #format line_dict = {'fragment_name': fragment_name, 'position_number': position_number, 'ref_nucleotid': ref_nucleotid, 'p_or_m': p_or_m, 'reads_dict': reads_dict}
    #format reads_dict : {'reads_array_1', 'reads_array_2', ...}
    #format reads_array : [genotype_id, nb_reads, nb_a, nb_c, nb_g, nb_t]
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
            if(p_or_m == 'M'):
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

            if(k != 96):
                print("problem with number of genotypes for position", position_number)
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
            if(i == 1):
                print(line_array[0])
            line_array[1] = int(line_array[1])
            line_array[2] = int(line_array[2])
            dict[i] = {'fragment_name': line_array[0], 'basic_position': line_array[1], 'harp_position': line_array[2]}

    print("\npositions dictionnary built")
    return dict


positions_dict = build_positions_dict("C:/Users/mabed/Documents/Travail/Etudes_techniques/fichiers_travail/positions_correspondance.txt")

def add_harp_positions(i, filepath_reads, positions_dict):
    n = 0
    dict = reads_statistics_reader(filepath_reads)
    new_dict = {}
    j = 0
    for key in dict:
        line_dict = dict[key]
        fragment_name = line_dict.get('fragment_name')
        position_number = line_dict.get('position_number')
        for pos_key in positions_dict:
            if (positions_dict[pos_key]['basic_position'] == position_number and positions_dict[pos_key]['fragment_name'] == fragment_name):
                corresp_harp_position = positions_dict[pos_key]['harp_position']
                new_dict[corresp_harp_position] = line_dict

        # la position n'est pas une clé car on retrouve plusieurs fois une même position.
        # on doit mettre les numéros de positions à la référence HARP pour comparer avec les nucléotypes
        # à priori il semble y avoir 5242 numéros harp, donc la matrice G devrait avoir 5242 lignes.
    print("\nharp positions added")
    return new_dict

harp_dict = add_harp_positions(1, "C:/Users/mabed/Documents/Travail/Etudes_techniques/fichiers_travail/reads_statistics.txt", positions_dict)
#format : dict = {'numero_harp_1' : line_dict_1, 'numero_harp_2' : line_dict_2, ...}

def obtain_reads_proba(harp_dict, nb_genotypes):
    for key in harp_dict:
        line_dict = harp_dict[key]
        reads_dict = line_dict.get('reads_dict')
        proba_dict = {}
        # format proba_array : [genotype_id, nb_reads, proba_a, proba_c, proba_g, proba_t]
        for i in range(nb_genotypes):
            reads_array_name = 'reads_array_' + str(i)
            proba_array_i = reads_dict.get(reads_array_name).copy()
            if(line_dict['p_or_m'] == 'P'):
                if(proba_array_i[1] != 0) :
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
                if(line_dict['ref_nucleotid'] == 'A'):
                    proba_array_i[2] = 1
                if (line_dict['ref_nucleotid'] == 'C'):
                    proba_array_i[3] = 1
                if (line_dict['ref_nucleotid'] == 'G'):
                    proba_array_i[4] = 1
                if (line_dict['ref_nucleotid'] == 'T'):
                    proba_array_i[5] = 1
            proba_array_name = 'proba_array_' + str(i)
            proba_dict[proba_array_name] = proba_array_i
        line_dict['proba_dict'] = proba_dict
    return harp_dict

dict_with_proba = obtain_reads_proba(harp_dict, 96)

def encode_nucleotype(filepath):
    G = []
    # A = 0, C = 1, G = 2, T = 3



#treat_one_mixture(0, "C:/Users/mabed/Documents/Travail/Etudes_techniques/fichiers_travail/reads_statistics.txt")

# les 4 probas dans des tableaux pour chaque position

#pratique : le num harp de la ligne 1 avec position 97 est 7783
