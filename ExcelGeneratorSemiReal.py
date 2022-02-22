import random

from scipy import stats
import random as rd
import pandas as pd
import numpy as np

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
mixtures_list = get_mixtures_list(mixtures_dict)

# for mixture_name in mixtures_list:
#     G = generate_G_from_mix(mixture_name, mixtures_dict, nucleotypes, column_nucleotypes_dict)
#     print(G[0:10])
#
#     reads = reads_of_mix(mixture_name, harp_dict, positions_errors)
#     print(reads)

# positions_dict = build_positions_dict(filepath_positions)
# harp_dict = add_harp_positions(filepath_reads, positions_dict)
reads_possibles = extract_reads_nb(harp_dict)

def generateLambda(nb_genotypes):
    freq_geno = np.random.random(size=nb_genotypes)
    somme = np.sum(freq_geno)
    freq_geno = np.round(freq_geno/somme,2)
    if sum(freq_geno) < 1:
        freq_geno[1]+=0.01
    if sum(freq_geno) > 1:
        freq_geno[1] -= 0.01
    return freq_geno

# def generateGauto(nb_genotypes,nb_snips):
#     return np.random.randint(2, size=(nb_snips,nb_genotypes))
G_test = generate_G_from_mix('Tm1002', mixtures_dict, nucleotypes, column_nucleotypes_dict)
G_T_test = G_test.transpose()


def generateReads_observ(freq,G):
    reads_observ = np.array(random.choices(reads_possibles,k=len(G)))
    Pf1 = G @ freq
    nb1_observe = np.random.binomial(reads_observ, Pf1)
    return np.vstack((reads_observ,nb1_observe))

# def binomInverse(n,p):
#     x=0
#     rand = rd.random()
#     while(stats.binom.cdf(x,n,p)<rand):
#         x+=1
#     return x

def generationSemiReal(G,nb_iter):
    for iter in range(nb_iter):
        freq = np.array(generateLambda(12))
        # G = generateGauto(nb_genomes,nb_snp)
        reads = generateReads_observ(freq,G)
        Pf1_attendu = G@freq
        Matrice = np.c_[G.T,freq]
        Pf1_attendu_excel = np.append(Pf1_attendu,0)
        Matrice = np.vstack([Matrice,Pf1_attendu_excel])
        reads_excel = np.c_[reads,[[0],[0]]]
        Matrice = np.vstack([Matrice,reads_excel])
        indiceLigne = []
        indiceCol = []
        for i in range(12):
            indL = 'G'+str(i+1)
            indiceLigne.append(indL)
        for j in range(len(G)):
            indC = 'SNP'+str(j+1)
            indiceCol.append(indC)
        indiceLigne.append('fq1 attendue')
        indiceLigne.append('nbReads')
        indiceLigne.append('nb1')
        indiceCol.append('freq_init')
        df = pd.DataFrame(Matrice,index=indiceLigne,columns=indiceCol)
        with pd.ExcelWriter('outputSR.xlsx',mode='a') as writer:
            df.to_excel(writer,sheet_name='N°'+str(iter+1))

generationSemiReal(G_test,30)
