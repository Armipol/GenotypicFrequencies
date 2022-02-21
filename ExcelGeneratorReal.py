import random

from scipy import stats
import random as rd
import pandas as pd
import numpy as np

import FileReader

positions_dict = FileReader.build_positions_dict("D:/Rémi/Documents/IMT/3A/S10/EtudeTech/melange_simul_renom/positions_correspondance.txt")
harp_dict = FileReader.add_harp_positions("D:/Rémi/Documents/IMT/3A/S10/EtudeTech/melange_simul_renom/SIMULS_READS_MIXTURES_fauxBAM_fauxREADS/reads_statistics.txt", positions_dict)
reads_possibles = FileReader.extract_reads_nb(harp_dict)

def generation_reelle(G, reads,nb_iter):
    for iter in range(nb_iter):
        #freq = np.array(generateLambda(nb_genomes))
        #G = generateGauto(nb_genomes,nb_snp)
        #reads = generateReads_observ(nb_snp,freq,G)
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
        with pd.ExcelWriter('output3.xlsx',mode='a') as writer:
            df.to_excel(writer,sheet_name='N°'+str(iter+1))

generation_reelle(12,5242,1)
