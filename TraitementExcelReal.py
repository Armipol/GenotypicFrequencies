from math import *
import numpy as np
import pandas as pd
from scipy import stats
import scipy.optimize
import sys
from openpyxl import load_workbook

########## Fichier de traitement pour les données réelles

## Méthodes de calcul

# Méthode pseudo-inversion directe de G : lambda = pseudoinv(G) * Pf1
# Cette fonction est équivalente à la closed-form d'Essr
def pinv_only(G_t, Pf1) :
    G = np.array(G_t).transpose()
    return np.linalg.pinv(G) @ Pf1

# Fonction de vraisemblance classique
# def likelihood(lam,reads,G):
#     lam = lam/np.sum(lam)
#     Pflambda = np.dot(G,lam)
#     vraisemblance = 1
#     for i in range(len(reads[1])):
#         proba = stats.binom.pmf(reads[1][i],reads[0][i],Pflambda[i])
#         vraisemblance = vraisemblance * proba
#     return -vraisemblance

# Fonction de vraisemblance logarithmique
def likelihoodLog(lam,reads,G):
    lam = lam/np.sum(lam)
    Pflambda = np.dot(G,lam)
    vraisemblance = 0
    for i in range(len(reads[1])):
        proba = stats.binom.pmf(reads[1][i],reads[0][i],Pflambda[i])
        if proba == 0:
            vraisemblance += sys.float_info.epsilon # Plus petit flottant possible en Python
        else :
            vraisemblance += np.log(proba)
    return -vraisemblance

## Traitement des données initiales :

def traitementDonneesReelles(fileExcel,listeMelanges):
    for melange in listeMelanges: # Nom des feuilles du excel = N° du mélange

        print("mélange N° :", melange)
        data = pd.read_excel(fileExcel,"Tm10"+str(melange).zfill(2))

        # Récupération de G
        G = pd.DataFrame(data.iloc[0:len(data.index)-2,1:len(data.columns)-1]).to_numpy()
        G_T = G.T
        nb_geno = np.shape(G_T)[1]
        nb_snps = np.shape(G_T)[0]

        # Récupération des reads
        reads = pd.DataFrame(data.iloc[[len(data.index)-2,len(data.index)-1],1:len(data.columns)-1]).to_numpy()

        # Pf1 observé qui servira à construire lam_init
        Pf1_observe = np.zeros(nb_snps)
        for i in range(nb_snps): # Obligé de tout parcourir car certaines positions ont 0 reads
            if reads[0][i] != 0:
                Pf1_observe[i] = reads[1][i] / reads[0][i]

        # Définition du lambda initial qui sera utilisé pour les optimisations
        lam_init = pinv_only(G,Pf1_observe)
        for i in range(lam_init.size):
            if lam_init[i] < 0 :
                lam_init[i]=0

        # Création des bounds pour l'utilisation dans les algorithmes
        bounds = ()
        for i in range(nb_geno):
            bounds = bounds + ((0,1),)

        ## Maximisation de la fonction de vraisemblance (minimisation de son opposé) :
        print("début calcul algorithmes")
        min_NMLog = scipy.optimize.minimize(likelihoodLog,lam_init,args=(reads,G_T),method= 'Nelder-Mead',bounds=bounds,options={'maxiter':50,'maxfev':50})
        print("algos ok")

        ## Traitement des résultats :

        # On crée une matrice à afficher sur excel
        lambda_differents = np.array([min_NMLog.x,lam_init])
        likelihoodsLog = np.array([likelihoodLog(min_NMLog.x,reads,G_T),likelihoodLog(lam_init,reads,G_T)])
        lambda_differents = lambda_differents.T
        matrice = np.vstack((lambda_differents,likelihoodsLog))
        colonnes = ["NMLog","LambdaMatrix"]

        # Export vers Excel
        book = load_workbook(fileExcel)
        writer = pd.ExcelWriter(fileExcel,engine='openpyxl')
        writer.book = book
        writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
        df = pd.DataFrame(matrice,columns=colonnes)
        df.to_excel(writer,"Tm10"+str(melange).zfill(2), startcol=3, startrow=nb_geno + 5, index=False)
        writer.save()

traitementDonneesReelles("outputReal.xlsx",[4,6])
