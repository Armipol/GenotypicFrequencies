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

def treat_one_mixture(i, filepath):
    mixture_dict = {}
    n = 0
    dict = reads_statistics_reader(filepath)
    for key in dict:
        line_dict = dict[key]
        position_number = line_dict.get('position_number')
        reads_dict = line_dict.get('reads_dict')
        reads_array_name = 'reads_array_' + str(i)
        reads_array_i = reads_dict.get(reads_array_name)
        if(mixture_dict.get(position_number) != None):
            print("already exists")
            n+=1
        print("n", n)
        # la position n'est pas une clé car on retrouve plusieurs fois une même position.

        mixture_dict[position_number] = reads_array_i


treat_one_mixture(0, "C:/Users/mabed/Documents/Travail/Etudes_techniques/fichiers_travail/reads_statistics.txt")

# les 4 probas dans des tableaux pour chaque position
