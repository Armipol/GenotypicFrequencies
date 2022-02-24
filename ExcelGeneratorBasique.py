import random as rd
import pandas as pd
import numpy as np
from openpyxl import Workbook


def generateLambda(nb_genotypes):
    freq_geno = np.zeros(nb_genotypes)
    for i in range(nb_genotypes):
        freq_geno[i] = round(rd.random(),2)
    somme = np.sum(freq_geno)
    freq_geno = np.round(freq_geno/somme,2)
    if sum(freq_geno) < 1:
        freq_geno[1]+=0.01
    if sum(freq_geno) > 1:
        freq_geno[1] -= 0.01
    return freq_geno

def generateG(nb_genotypes,nb_snips):
    return np.random.randint(2,size=(nb_genotypes,nb_snips))

def generateReads_observ(nb_snips,freq,G):
    reads_observ = np.zeros((2,nb_snips))
    Pf1 = G @ freq
    reads1 = np.random.randint(0,15,size=int(nb_snips/2))
    reads2 = np.random.randint(0,200,size=int(nb_snips/2))
    Pf1_split = np.split(Pf1,2)
    nb1_observe1 = np.random.binomial(reads1,Pf1_split[0])
    nb1_observe2 = np.random.binomial(reads2,Pf1_split[1])
    reads = np.concatenate((reads1,reads2))
    nb1_observe = np.concatenate((nb1_observe1,nb1_observe2))
    reads_observ[0] = reads
    reads_observ[1] = nb1_observe
    return reads_observ

def generationBasique(nb_genomes,nb_snp,nb_iter,fileExcel):
    # Création du fichier excel
    wb = Workbook()
    wb.save(fileExcel)
    for iter in range(nb_iter):
        # Génération des fréquences à retrouver
        freq = np.array(generateLambda(nb_genomes))

        # Génération de la matrice G
        G = generateG(nb_genomes,nb_snp)

        # Création des reads
        reads = generateReads_observ(nb_snp,freq,G.T)
        Pf1_attendu = G.T@freq

        # Mise en forme des données pour export vers Excel
        Matrice = np.c_[G,freq]
        Pf1_attendu_excel = np.append(Pf1_attendu,0)
        Matrice = np.vstack([Matrice,Pf1_attendu_excel])
        reads_excel = np.c_[reads,[[0],[0]]]
        Matrice = np.vstack([Matrice,reads_excel])

        # Création des titres des colonnes/lignes
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

        # Export vers Excel
        df = pd.DataFrame(Matrice,index=indiceLigne,columns=indiceCol)
        with pd.ExcelWriter(fileExcel,mode='a') as writer:
            df.to_excel(writer,sheet_name='N°'+str(iter+1))

generationBasique(12,1000,5,"outputBasique3.xlsx")
