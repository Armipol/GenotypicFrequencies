import random
import pandas as pd

from FileReader import *

filepath_positions = "D:/Rémi/Documents/IMT/3A/S10/EtudeTech/melange_simul_renom/positions_correspondance.txt"
filepath_reads = "D:/Rémi/Documents/IMT/3A/S10/EtudeTech/melange_simul_renom/SIMULS_READS_MIXTURES_fauxBAM_fauxREADS/reads_statistics.txt"
filepath_nucleotypes = "D:/Rémi/Documents/IMT/3A/S10/EtudeTech/melange_simul_renom/nucleotypes.txt"
filepath_mixtures = "D:/Rémi/Documents/IMT/3A/S10/EtudeTech/melange_simul_renom/simulated_mixtures_composition.txt"

data_utils = build_data_utils(filepath_positions, filepath_reads, filepath_nucleotypes, filepath_mixtures)
mixtures_dict = data_utils[0]
nucleotypes = data_utils[1]
column_nucleotypes_dict = data_utils[2]
positions_errors = data_utils[3]
harp_dict = data_utils[4]
# mixtures_list = get_mixtures_list(mixtures_dict)

def generateLambda(nb_genotypes):
    freq_geno = np.random.random(size=nb_genotypes)
    somme = np.sum(freq_geno)
    freq_geno = np.round(freq_geno/somme,2)
    if sum(freq_geno) > 1:
        freq_geno[1] -= 0.01
    return freq_geno

indiceLigne = []
indiceCol = []
for i in range(12):
    indL = 'G' + str(i + 1)
    indiceLigne.append(indL)
for j in range(4404):
    indC = 'SNP' + str(j + 1)
    indiceCol.append(indC)
indiceLigne.append('fq1 attendue')
indiceLigne.append('nbReads')
indiceLigne.append('nb1')
indiceCol.append('freq_init')

def generationSemiReal(fileExcel,listeMelanges):
    for melange in listeMelanges:
        print("Mélange n° : ",melange)

        # Génération de G à partir du mélange demandé
        G = generate_G_from_mix("Tm10"+str(melange).zfill(2), mixtures_dict, nucleotypes, column_nucleotypes_dict)

        # Récupération des reads de ce mélange
        reads = reads_of_mix("Tm10"+str(melange).zfill(2), harp_dict, positions_errors)
        reads = reads[0]
        reads = np.array(reads)
        reads = reads.astype(int)

        # Génération aléatoire des fréquences
        freq = np.array(generateLambda(12))

        # Génération binomiale des 1 obtenus
        Pf1 = G @ freq
        nb1_observe = np.random.binomial(reads, Pf1)

        # On met en forme les données pour les stocker dans Excel
        reads_excel = np.append(reads,0) # On ajoute des 0 pour avoir le bon nombre de termes dans la matrice
        nb1_observe = np.append(nb1_observe,0)
        Matrice = np.c_[G.T, freq]
        Pf1_attendu_excel = np.append(Pf1, 0)
        Matrice = np.vstack([Matrice, Pf1_attendu_excel])
        Matrice = np.vstack([Matrice, reads_excel])
        Matrice = np.vstack([Matrice,nb1_observe])

        # Export vers Excel
        df = pd.DataFrame(Matrice, index=indiceLigne, columns=indiceCol)
        with pd.ExcelWriter(fileExcel, mode='a') as writer:
            df.to_excel(writer, sheet_name="Tm10"+str(melange).zfill(2))

generationSemiReal("outputSemiReal.xlsx",[5,8])
