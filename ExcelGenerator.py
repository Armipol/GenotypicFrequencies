import random

from scipy import stats
import random as rd
import pandas as pd
import numpy as np

import FileReader

# positions_dict = FileReader.build_positions_dict("D:/Rémi/Documents/IMT/3A/S10/EtudeTech/melange_simul_renom/positions_correspondance.txt")
# harp_dict = FileReader.add_harp_positions("D:/Rémi/Documents/IMT/3A/S10/EtudeTech/melange_simul_renom/SIMULS_READS_MIXTURES_fauxBAM_fauxREADS/reads_statistics.txt", positions_dict)
# reads_possibles = FileReader.extract_reads_nb(harp_dict)

def generateLambda(nb_genotypes):
    freq_geno = np.zeros(nb_genotypes)
    for i in range(nb_genotypes):
        freq_geno[i] = round(rd.random(),2)
        print(freq_geno[i])
    print(freq_geno)
    somme = np.sum(freq_geno)
    freq_geno = np.round(freq_geno/somme,2)
    if sum(freq_geno) < 1:
        freq_geno[1]+=0.01
    if sum(freq_geno) > 1:
        freq_geno[1] -= 0.01
    return freq_geno

f = generateLambda(12)
print(sum(f))
def generateG(nb_genotypes,nb_snips):
    G = np.zeros((nb_snips,nb_genotypes))
    for i in range(nb_snips):
        for j in range(nb_genotypes):
            G[i][j] = rd.randint(0,1)
    return G

def generateReads_observ(nb_snips,freq,G):
    reads_observ = np.zeros((2,nb_snips))
    nb_reads = np.zeros(nb_snips)
    nb1_observe = np.zeros(nb_snips)
    Pf = G @ freq
    for i in range(nb_snips):
        read = rd.randint(0,10000)
        nb_reads[i] = read
        nb1_observe[i] = np.random.binomial(read,Pf[i])
    reads_observ[0] = nb_reads
    reads_observ[1] = nb1_observe
    return reads_observ

# def binomInverse(n,p):
#     x=0
#     rand = rd.random()
#     while(stats.binom.cdf(x,n,p)<rand):
#         x+=1
#     return x

def generationBasique(nb_genomes,nb_snp,nb_iter):
    for iter in range(nb_iter):
        freq = np.array(generateLambda(nb_genomes))
        G = generateG(nb_genomes,nb_snp)
        reads = generateReads_observ(nb_snp,freq,G)
        Pf1_attendu = G@freq
        Matrice = np.c_[G.T,freq]
        Pf1_attendu_excel = np.append(Pf1_attendu,0)
        Matrice = np.vstack([Matrice,Pf1_attendu_excel])
        reads_excel = np.c_[reads,[[0],[0]]]
        Matrice = np.vstack([Matrice,reads_excel])
        indiceLigne = []
        indiceCol = []
        for i in range(nb_genomes):
            indL = 'G'+str(i+1)
            indiceLigne.append(indL)
        for j in range(nb_snp):
            indC = 'SNP'+str(j+1)
            indiceCol.append(indC)
        indiceLigne.append('fq1 attendue')
        indiceLigne.append('nbReads')
        indiceLigne.append('nb1')
        indiceCol.append('freq_init')
        df = pd.DataFrame(Matrice,index=indiceLigne,columns=indiceCol)
        with pd.ExcelWriter('output.xlsx',mode='a') as writer:
            df.to_excel(writer,sheet_name='N°'+str(iter+1))

generationBasique(12,10,10)
# print(random.choices(reads_possibles,k=10))
