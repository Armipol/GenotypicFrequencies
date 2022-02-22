### IMPORTS

from math import *
import numpy as np
import pandas as pd
from scipy import stats
import scipy.optimize
import sys
from openpyxl import load_workbook


### Implémentation des méthodes de traitement :

## Méthodes matricielles :

# Méthode pseudo-inversion directe de G : lambda = pseudoinv(G) * Pf1
# Attention ici à l'inversion de G qui n'est pas toujours possible ?

#Cette fonction est équivalente à la closed-form d'Essr
def pinv_only(G_t, Pf1) :
    G = np.array(G_t).transpose()
    return np.linalg.pinv(G) @ Pf1


## Méthodes de maximisation de la vraisemblance

# Fonction de vraisemblance classique
def likelihood(lam,reads,G):
    lam = lam/np.sum(lam)
    Pflambda = np.dot(G,lam)
    # print(Pflambda)
    vraisemblance = 1
    for i in range(len(reads[1])):
        proba = stats.binom.pmf(reads[1][i],reads[0][i],Pflambda[i])
        vraisemblance = vraisemblance * proba
    return -vraisemblance

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


### Traitement des données :
print("départ")

# nb_feuille = 10

# for f in range(1,nb_feuille+1):
# data = pd.read_excel("outputSR.xlsx","N°"+str(f))
data = pd.read_excel("outputReel.xlsx","N°1")
# print(data)
print("excel reader ok")

# Récupération des fréquences initiales
# freq = pd.DataFrame(data['freq_init']).to_numpy()
# freq = np.reshape(freq.T,np.shape(freq)[0])
# freq = freq[0:len(freq)-3]

# Récupération du Pf1
# Pf1 = pd.DataFrame(data.iloc[len(data.index)-3]).to_numpy()
# Pf1 = np.reshape(Pf1.T,np.shape(Pf1)[0])
# Pf1 = Pf1[1:len(Pf1)-1]

# Récupération de G
G = pd.DataFrame(data.iloc[0:len(data.index)-3,1:len(data.columns)-1]).to_numpy()
G = G.T

# Récupération des reads
reads = pd.DataFrame(data.iloc[[len(data.index)-2,len(data.index)-1],1:len(data.columns)-1]).to_numpy()
nb_geno = np.shape(G)[1]
nb_snp = np.shape(G)[0]

# Pf1 observé qui servira à construire lam_init
# Pf1_observe = np.zeros(nb_snp)
# for i in range(nb_snp):
#     Pf1_observe[i] = reads[1][i]/reads[0][i] # Attention au cas où 1obs = 0
# Pf1_observe = reads[1]/reads[0]
Pf1_observe = np.zeros(nb_snp)
for i in range(nb_snp):
    if reads[0][i] != 0:
        Pf1_observe[i] = reads[1][i]/reads[0][i]

# Définition du lambda initial qui sera utilisé pour les optimisations
lam_init = pinv_only(G.T,Pf1_observe)
for i in range(lam_init.size):
    if lam_init[i] < 0 :
        lam_init[i]=0
# print("lam_init",lam_init)
# lam_init = [0.2,0.2,0.2,0.2,0.2]

bounds = ()
for i in range(nb_geno):
    bounds = bounds + ((0,1),)

print("variables ok")
### Maximisation de la fonction de vraisemblance (minimisation de son opposé) : tentative avec toutes les méthodes
min_NM = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'Nelder-Mead',bounds=bounds,options={'maxiter':50,'maxfev':50})
# print("NM")
# print(min_NM.x)
min_NMLog = scipy.optimize.minimize(likelihoodLog,lam_init,args=(reads,G),method= 'Nelder-Mead',bounds=bounds,options={'maxiter':50,'maxfev':50})
# print("NMLog")
# print(min_NMLog.x)
# min_Powell = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'Powell',bounds=bounds)
# print("Powell")
# print(min_Powell.x)
# min_PowellLog = scipy.optimize.minimize(likelihoodLog,lam_init,args=(reads,G),method= 'Powell',bounds=bounds)
# print("PowellLog")
# print(min_PowellLog)
# print(min_PowellLog.x)
# min_CG = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'CG')
# print("CG")
# print(min_CG.x)
# min_CGLog = scipy.optimize.minimize(likelihoodLog,lam_init,args=(reads,G),method= 'CG')
# print("CGLog")
# print(min_CGLog.x)
# min_BFGS = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'BFGS')
# print("BFGS")
# print(min_BFGS.x)
# min_BFGSLog = scipy.optimize.minimize(likelihoodLog,lam_init,args=(reads,G),method= 'BFGS')
# print("BFGSLog")
# print(min_BFGSLog.x)
# min_LBFGSB = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'L-BFGS-B',bounds=bounds)
# print("LBFGSB")
# print(min_LBFGSB.x)
# min_LBFGSBLog = scipy.optimize.minimize(likelihoodLog,lam_init,args=(reads,G),method= 'L-BFGS-B',bounds=bounds)
# print("LBFGSBLog")
# print(min_LBFGSBLog.x)
# min_TNC = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'TNC')
# print("TNC")
# print(min_TNC.x)
min_TNCLog = scipy.optimize.minimize(likelihoodLog,lam_init,args=(reads,G),method= 'TNC')
# print("TNCLog")
# print(min_TNCLog.x)
# min_COBYLA = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'COBYLA')
# print("COBYLA")
# print(min_COBYLA.x)
# min_COBYLALog = scipy.optimize.minimize(likelihoodLog,lam_init,args=(reads,G),method= 'COBYLA')
# print("COBYLALog")
# print(min_COBYLALog.x)
# min_trust_constr = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'trust-constr',bounds=bounds)
# min_trust_constrLog = scipy.optimize.minimize(likelihoodLog,lam_init,args=(reads,G),method= 'trust-constr',bounds=bounds)
print("algos ok")
# Jacobian is required
# min_NewtonCG = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'Newton-CG',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
# min_dogleg = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'dogleg',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
# min_trust_ncg = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'trust-ncg',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)))
# min_trust_exact = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'trust-exact',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
# min_trust_krylov = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'trust-krylov',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)))
# min_SLSQP = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'SLSQP',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})



### Analyses

# def affichageErreur(freqOrigin,freqObs):
#     erreurAbs = abs(freqOrigin-freqObs)
#     erreurRelative = abs(freqOrigin-freqObs)
#     for i in range(len(erreurRelative)):
#         erreurRelative[i] = erreurRelative[i]/freq[i]
#     print("Erreur absolue : ", erreurAbs)
#     print("Erreur absolue moyenne % : ",round(np.mean(erreurAbs)*100,2))
#     print("Erreur relative : ", erreurRelative)
#     print("Erreur relative moyenne % : ", round(np.mean(erreurRelative)*100,2))
#
# print("\n")
# print("lambda à trouver : ",freq)
# print("\n")
#
# print("lambda initial donné pour minimisation = lambda prédit inversion simple : ", lam_init)
# affichageErreur(freq,lam_init)
# print("\n")
#
# # print("lambda prédit minimisation Essr : ", minimizeEssr(G, Pf1_observe))
# # affichageErreur(freq,minimizeEssr(G.T,Pf1_observe))
# # print("\n")
#
# print("lambda à trouver : ",freq)
# print("lambda de l'algo Nelder-Mead : ", min_NM.x/np.sum(min_NM.x))
# affichageErreur(freq,min_NM.x)
# print("lambda de l'algo Nelder-Mead Log: ", min_NMLog.x/np.sum(min_NMLog.x))
# affichageErreur(freq,min_NMLog.x)
# print("\n")
#
# print("lambda à trouver : ",freq)
# print("lambda de l'algo Powell : ",min_Powell.x/np.sum(min_Powell.x))
# affichageErreur(freq,min_Powell.x)
# print("lambda de l'algo Powell Log : ",min_PowellLog.x/np.sum(min_PowellLog.x))
# affichageErreur(freq,min_PowellLog.x)
# print("\n")
#
# print("lambda à trouver : ",freq)
# print("lambda de l'algo COBYLA : ", min_COBYLA.x/np.sum(min_COBYLA.x))
# affichageErreur(freq,min_COBYLA.x)
# print("lambda de l'algo COBYLA Log : ", min_COBYLALog.x/np.sum(min_COBYLALog.x))
# affichageErreur(freq,min_COBYLALog.x)
# print("\n")
#
# print("lambda à trouver : ",freq)
# print("lambda de l'algo CG : ", min_CG.x/np.sum(min_CG.x)) #Renvoie equiproportion
# affichageErreur(freq,min_CG.x)
# print("lambda de l'algo CG Log : ", min_CGLog.x/np.sum(min_CGLog.x)) #Renvoie equiproportion
# affichageErreur(freq,min_CGLog.x)
# print("\n")
#
# print("lambda à trouver : ",freq)
# print("lambda de l'algo BFGS : ", min_BFGS.x/np.sum(min_BFGS.x))  #Renvoie equiproportion
# affichageErreur(freq,min_BFGS.x)
# print("lambda de l'algo BFGS Log : ", min_BFGSLog.x/np.sum(min_BFGSLog.x))  #Renvoie equiproportion
# affichageErreur(freq,min_BFGSLog.x)
# print("\n")
#
# print("lambda à trouver : ",freq)
# print("lambda de l'algo L-BFGS-B : ", min_LBFGSB.x/np.sum(min_LBFGSB.x))  #Renvoie equiproportion
# affichageErreur(freq,min_LBFGSB.x)
# print("lambda de l'algo L-BFGS-B Log : ", min_LBFGSBLog.x/np.sum(min_LBFGSBLog.x))  #Renvoie equiproportion
# affichageErreur(freq,min_LBFGSBLog.x)
# print("\n")
#
# print("lambda à trouver : ",freq)
# print("lambda de l'algo TNC : ", min_TNC.x/np.sum(min_TNC.x))  #Renvoie equiproportion
# affichageErreur(freq,min_TNC.x)
# print("lambda de l'algo TNC Log : ", min_TNCLog.x/np.sum(min_TNCLog.x))  #Renvoie equiproportion
# affichageErreur(freq,min_TNCLog.x)
# print("\n")
#
# print("lambda à trouver : ",freq)
# print("lambda de l'algo trust-constr ", min_trust_constr.x/np.sum(min_trust_constr.x))  #Renvoie equiproportion
# affichageErreur(freq,min_trust_constr.x)
# print("lambda de l'algo trust-constr Log ", min_trust_constrLog.x/np.sum(min_trust_constrLog.x))  #Renvoie equiproportion
# affichageErreur(freq,min_trust_constrLog.x)
# print("\n")

### Traitement des données :

# Création d'un dictionnaire
# lambdas = {"NM":min_NM.x,"NMLog":min_NMLog.x,"Powell":min_Powell.x,"PowellLog":min_PowellLog.x,"COBYLA":min_COBYLA.x,"COBYLALog":min_COBYLALog.x,"CG":min_CG.x,"CGLog":min_CGLog.x,"BFGS":min_BFGS.x,"BFGSLog":min_BFGSLog.x,"L-BFGS-B":min_LBFGSB.x,"L-BFGS-BLog":min_LBFGSBLog.x,"TNC":min_TNC.x,"TNCLog":min_TNCLog.x} #,"trust-constr":min_trust_constr.x,"trust-constrLog":min_trust_constrLog.x}
lambdas = {"NM":min_NM.x,"NMLog":min_NMLog.x,"TNCLog":min_TNCLog.x} #On a supprimé Powell, Cobyla, CGLog, PowellLog et BFGSLog
# NM tourne pour NM,CobylaLog,CG,BFGS,LBFGSB,TNC

# lambdas2 = np.array([min_NM.x,min_NMLog.x,min_Powell.x,min_PowellLog.x,min_COBYLA.x,min_COBYLALog.x,min_CG.x,min_CGLog.x,min_BFGS.x,min_BFGSLog.x,min_LBFGSB.x,min_LBFGSBLog.x,min_TNC.x,min_TNCLog.x,min_trust_constr.x,min_trust_constrLog.x])

# Voyons quels algorithmes convergent vers le même lambda
memes_algo = [["NM"]]
lambda_differents = np.array([min_NM.x])
for cle,valeur in lambdas.items():
    if valeur.tolist() in lambda_differents.tolist():
        for i in range(len(memes_algo)):
            val = lambdas.get(memes_algo[i][0])
            if np.all(val == valeur) :
                memes_algo[i].append(cle)
    else:
        lambda_differents = np.vstack([lambda_differents,valeur])
        memes_algo.append([cle])
del memes_algo[0][0] # On supprime le premier "NM" qu'il était nécessaire de rajouter pour passer sans la boucle
lambda_differents = np.vstack((lam_init,lambda_differents)) # On rajoute le lambda initial
# print(memes_algo,len(memes_algo))
# print(lambda_differents,np.shape(lambda_differents))
print(memes_algo)
# On concatène les mêmes algos
memes_algo_concat = []
for i in range(len(memes_algo)):
    algos = memes_algo[i][0]
    for j in range(1,len(memes_algo[i])):
        algos += ', '+memes_algo[i][j]
    memes_algo_concat.append(algos)
memes_algo_concat = np.array(memes_algo_concat)
memes_algo_concat = np.concatenate((np.array(["LambdaMatrix"]),memes_algo_concat))
# memes_algo_concat = np.concatenate(["LambdaMatrix"],memes_algo_concat)
# print(memes_algo_concat)
# print(memes_algo_concat,len(memes_algo_concat))
print("memes algos ok")
# On calcule les erreurs de chaque algorithme par rapport au lambda recherché
likelihoods = np.zeros(len(lambda_differents))
likelihoodsLog = np.zeros(len(lambda_differents))
erreur_abs = np.zeros(len(lambda_differents))
erreur_rel = np.zeros(len(lambda_differents))
for i in range(len(lambda_differents)):
    lambda_algo = lambda_differents[i]
    likelihoods[i] = likelihood(lambda_algo,reads,G)*10**40
    likelihoodsLog[i] = likelihoodLog(lambda_algo,reads,G)
    # print(likelihoods[i])
    # erreurAbs = abs(freq - lambda_algo)
    # erreurRelative = abs(freq-lambda_algo)
    # for j in range(len(erreurRelative)):
    #     if freq[j] != 0:
    #         erreurRelative[j] = erreurRelative[j]/freq[j]
    #     else :
    #         erreurRelative[j]=inf
    # erreur_abs[i] = np.mean(erreurAbs)
    # erreur_rel[i] = np.mean(erreurRelative)
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

# On crée une matrice à afficher sur excel
lambda_differents = lambda_differents.T
matrice = np.vstack((lambda_differents,likelihoods))
matrice = np.vstack((matrice,likelihoodsLog))
# matrice = np.vstack((matrice,erreur_abs))
# matrice = np.vstack((matrice,erreur_rel))
#
# print("erreurs ok")
#
book = load_workbook('outputReel.xlsx')
writer = pd.ExcelWriter('outputReel.xlsx',engine='openpyxl')
writer.book = book
writer.sheets = dict((ws.title, ws) for ws in book.worksheets)

df = pd.DataFrame(matrice,columns=memes_algo_concat)
# print(df)
# df.to_excel(writer,"N°"+str(f),startcol=5,startrow=nb_geno+5,index=False)
df.to_excel(writer,"N°1",startcol=5,startrow=nb_geno+5,index=False)

# print(matrice)

# erreur_abs_min = min(erreur_abs)
# erreur_rel_min = min(erreur_rel)
# indiceMinAbs = erreur_abs.tolist().index(erreur_abs_min)
# indiceMinRel = erreur_rel.tolist().index(erreur_rel_min)
# algosMinAbs = memes_algo_concat[indiceMinAbs]
# algosMinRel = memes_algo_concat[indiceMinRel]

# matrix = [[erreur_abs_min,algosMinAbs],[erreur_rel_min,algosMinRel]]

# df2 = pd.DataFrame(matrix,index = ['erreur absolue min','erreur relative min'], columns = ['valeur','algorithmes'])
# df2.to_excel(writer,"N°"+str(f), startrow=nb_geno+6,startcol=1)
# df2.to_excel(writer,"N°3", startrow=nb_geno+6,startcol=1)

writer.save()
