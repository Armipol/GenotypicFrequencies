from math import *
import numpy as np
import pandas as pd
from scipy import stats
import scipy.optimize
import sys
from openpyxl import load_workbook
from sklearn.metrics import mean_squared_error, mean_absolute_error, max_error

## Méthodes mathématiques

# Méthode pseudo-inversion directe de G : lambda = pseudoinv(G) * Pf1
# Cette fonction est équivalente à la closed-form d'Essr
def pinv_only(G_t, Pf1) :
    G = np.array(G_t).transpose()
    return np.linalg.pinv(G) @ Pf1

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
            # print(vraisemblance)
    return -vraisemblance

# Fonction de vraisemblance classique
# def likelihood(lam,reads,G):
#     lam = lam/np.sum(lam)
#     Pflambda = np.dot(G,lam)
#     # print(Pflambda)
#     vraisemblance = 1
#     for i in range(len(reads[1])):
#         proba = stats.binom.pmf(reads[1][i],reads[0][i],Pflambda[i])
#         vraisemblance = vraisemblance * proba
#     return -vraisemblance

### Traitement des données :

def traitementDonneesSemiReal(fileExcel,listeMelanges):
    print("départ")
    for melange in listeMelanges:
        data = pd.read_excel(fileExcel,"Tm10"+str(melange).zfill(2))

        # Récupération des fréquences initiales
        freq = pd.DataFrame(data['freq_init']).to_numpy()
        freq = np.reshape(freq.T, np.shape(freq)[0])
        freq = freq[0:len(freq) - 3]

        # Récupération de G
        G = pd.DataFrame(data.iloc[0:len(data.index)-3,1:len(data.columns)-1]).to_numpy()
        G_T = G.T

        # Récupération des reads
        reads = pd.DataFrame(data.iloc[[len(data.index) - 2, len(data.index) - 1], 1:len(data.columns) - 1]).to_numpy()
        nb_geno = np.shape(G_T)[1]
        nb_snps = np.shape(G_T)[0]

        # Pf1 observé qui servira à construire lam_init
        Pf1_observe = np.zeros(nb_snps)
        for i in range(nb_snps): # Obligé de tout parcourir pour les cas à nb_reads = 0
            if reads[0][i] !=0:
                Pf1_observe[i] = reads[1][i]/reads[0][i]

        # Définition du lambda initial qui sera utilisé pour les optimisations
        lam_init = pinv_only(G, Pf1_observe)
        for i in range(lam_init.size):
            if lam_init[i] < 0:
                lam_init[i] = 0

        # Création des bounds nécessaires aux algorithmes
        bounds = ()
        for i in range(nb_geno):
            bounds = bounds + ((0,1),)

        ### Maximisation de la fonction de vraisemblance (minimisation de son opposé) : tentative avec toutes les méthodes
        print("lancement algo")
        min_NMLog = scipy.optimize.minimize(likelihoodLog,lam_init,args=(reads,G_T),method= 'Nelder-Mead',bounds=bounds,options={'maxiter':50,'maxfev':50})
        print("algos ok")

        # On calcule les erreurs de chaque algorithme par rapport au lambda recherché
        lambda_differents = np.array([lam_init,min_NMLog.x])
        likelihoodsLog = np.zeros(len(lambda_differents))
        erreur_abs = np.zeros(len(lambda_differents))
        for i in range(len(lambda_differents)):
            lambda_algo = lambda_differents[i]
            likelihoodsLog[i] = likelihoodLog(lambda_algo,reads,G_T)
            erreurAbs = abs(freq - lambda_algo)
            erreur_abs[i] = np.mean(erreurAbs)
        erreur_abs_min = min(erreur_abs)
        indiceMinAbs = erreur_abs.tolist().index(erreur_abs_min)
        memes_algo_concat = ["LambdaMatrix", "NMLog"]
        algoMinAbs = memes_algo_concat[indiceMinAbs]

        # Log Nelder-Mead
        print("RMSE Log Nelder-Mead:", sqrt(mean_squared_error(freq, min_NMLog.x)))
        print("MSE Log Nelder-Mead:", mean_squared_error(freq, min_NMLog.x))
        print("MAE Log Nelder-Mead:", mean_absolute_error(freq, min_NMLog.x))
        print("Max error Log Nelder-Mead:", max_error(freq, min_NMLog.x))

        # Erreurs lamda matriciel
        print("RMSE Matriciel:", sqrt(mean_squared_error(freq, lam_init)))
        print("MSE Matriciel:", mean_squared_error(freq, lam_init))
        print("MAE Matriciel:", mean_absolute_error(freq, lam_init))
        print("Max error matriciel:", max_error(freq, lam_init))

        # On crée une matrice à afficher sur excel
        lambda_differents = lambda_differents.T
        matrice = np.vstack((lambda_differents, likelihoodsLog))
        matrice = np.vstack((matrice, erreur_abs))
        matrix = [erreur_abs_min, algoMinAbs]

        # On va écrire sur excel
        book = load_workbook(fileExcel)
        writer = pd.ExcelWriter(fileExcel, engine='openpyxl')
        writer.book = book
        writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
        df = pd.DataFrame(matrice, columns=memes_algo_concat)
        df.to_excel(writer,"Tm10"+str(melange).zfill(2), startcol=5, startrow=nb_geno + 5, index=False)
        df2 = pd.DataFrame(matrix,index = ['valeur','algorithme'],columns=['erreur absolue min'])
        df2.to_excel(writer,"Tm10"+str(melange).zfill(2), startrow=nb_geno+6,startcol=1)
        writer.save()

traitementDonneesSemiReal("outputSemiReal.xlsx",[5,8])
