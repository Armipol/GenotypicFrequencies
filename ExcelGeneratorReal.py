import random

from scipy import stats
import random as rd
import pandas as pd
import numpy as np
from math import sqrt

from FileReader import *

def generation_reelle(G, reads, nb_iter):
    for iter in range(nb_iter):
        #freq = np.array(generateLambda(nb_genomes))
        #G = generateGauto(nb_genomes,nb_snp)
        #reads = generateReads_observ(nb_snp,freq,G)
        # Pf1_attendu = G@freq
        matrice = G.T
        # Pf1_attendu_excel = np.append(Pf1_attendu,0)
        # Matrice = np.vstack([Matrice,Pf1_attendu_excel])
        reads_excel = np.c_[reads,[[0],[0]]]
        matrice = np.vstack([matrice,reads])
        print("matrice :", matrice[-1])
        indiceLigne = []
        indiceCol = []
        for i in range(len(G[0])):
            indL = 'G'+str(i+1)
            indiceLigne.append(indL)
        for j in range(len(G)):
            indC = 'SNP'+str(j+1)
            indiceCol.append(indC)
        print("indice ligne", indiceLigne)
        print("len indice colonne", len(indiceCol))
        print("ind col -2", indiceCol[-2])

        # indiceLigne.append('fq1 attendue')
        indiceLigne.append('nbReads')
        indiceLigne.append('nb1')
        df = pd.DataFrame(matrice,index=indiceLigne,columns=indiceCol)
        with pd.ExcelWriter('C:/Users/mabed/Documents/Travail/Etudes_techniques/resultats/output.xlsx', engine="openpyxl", mode="a") as writer:
            df.to_excel(writer,sheet_name='N°'+str(iter+1))


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

G_test = generate_G_from_mix('Tm1001', mixtures_dict, nucleotypes, column_nucleotypes_dict)
G_T_test = G_test.transpose()
reads_test = reads_of_mix('Tm1001', harp_dict, positions_errors)


generation_reelle(G_test, reads_test, 1)

#ci-dessous, à inclure dans TraitementExcelStats
#normalement j'ai pris les bonnes notations. À tester, notamment pas sûre du bon fonctionnement pour l'histogramme avec plt

# from sklearn.metrics import mean_squared_error, mean_absolute_error, max_error
# import matplotlib.pyplot as plt
#
# #Nelder-Mead
# print("RMSE Nelder-Mead:", sqrt(mean_squared_error(freq, min_NM.x)))
# print("MSE Nelder-Mead:", mean_squared_error(freq, min_NM.x))
# print("MAE Nelder-Mead:", mean_absolute_error(freq, min_NM.x))
# print("Max error Nelder-Mead:", max_error(freq, min_NM.x))
#
# plt.hist(sqrt(mean_squared_error(freq, min_NM.x)), bins=50)
# plt.show()
#
# #Log Nelder-Mead
# print("RMSE Log Nelder-Mead:", sqrt(mean_squared_error(freq, min_NMLog.x)))
# print("MSE Log Nelder-Mead:", mean_squared_error(freq, min_NMLog.x))
# print("MAE Log Nelder-Mead:", mean_absolute_error(freq, min_NMLog.x))
# print("Max error Log Nelder-Mead:", max_error(freq, min_NMLog.x))
#
# plt.hist(sqrt(mean_squared_error(freq, min_NMLog.x)), bins=50)
# plt.show()
#
# #TNC Log
# print("RMSE TNC Log:", sqrt(mean_squared_error(freq, min_TNCLog.x.x)))
# print("MSE TNC Log:", mean_squared_error(freq, min_TNCLog.x.x))
# print("MAE TNC Log:", mean_absolute_error(freq, min_TNCLog.x.x))
# print("Max error TNC Log:", max_error(freq, min_TNCLog.x.x))
#
# plt.hist(sqrt(mean_squared_error(freq, min_TNCLog.x.x)), bins=50)
# plt.show()




