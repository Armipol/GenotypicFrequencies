### IMPORTS

from math import *
import numpy as np
import random as rd
from scipy import stats
import scipy.optimize
import sys



### Définition des fonctions de génération :

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
    G = np.zeros((nb_genotypes,nb_snips))
    for i in range(nb_genotypes):
        for j in range(nb_snips):
            G[i][j] = rd.randint(0,1)
    return G

def generateReads_observ(nb_snips,freq,G):
    reads_observ = np.zeros((2,nb_snips))
    nb_reads = np.zeros(nb_snips)
    nb1_observe = np.zeros(nb_snips)
    for i in range(nb_snips):
        read = rd.randint(0,10000)
        nb_reads[i] = read
        Pf = G@freq
        nb1_observe[i] = np.random.binomial(read,Pf[i])
    reads_observ[0] = nb_reads
    reads_observ[1] = nb1_observe
    return reads_observ


### Implémentation des méthodes de traitement :

## Méthodes matricielles :


# Méthode pseudo-inversion directe de G : lambda = pseudoinv(G) * Pf1
# Attention ici à l'inversion de G qui n'est pas toujours possible ?

def pinv_only(G_t, Pf1) :
    G = np.array(G_t).transpose()
    return np.linalg.pinv(G) @ Pf1

#Méthode de minimisation de Essr pour trouver une valeur de lambda : utilisation de la closed-form
def minimizeEssr(G_t, Pf1) :
    G = np.array(G_t).transpose()
    G_nbraws = np.size(G, 0)
    Gd = np.c_[np.ones(G_nbraws), G]  # design matrix
    Gd_t = Gd.transpose()
    lam = np.linalg.inv(Gd_t @ Gd) @ Gd_t @ Pf1 # closed-form solution
    return lam[1:6]


## Méthodes de maximisation de la vraisemblance

# Fonction de vraisemblance classique
def likelihood(lam,reads,G):
    lam = lam/np.sum(lam)
    Pflambda = np.dot(G,lam)
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


### Génération des données

# Vecteur lambda à retrouver :
freq = np.array(generateLambda(5))

# Matrice des génotypes :
G = generateG(5,9)
print("Matrice G : ","\n", G)
G = G.T

# Matrice des lectures :
reads = generateReads_observ(9,freq,G)
print("Reads & Nb1 : ", reads)

# Pf1 observé qui servira à construire lam_init
Pf1_observe = np.zeros(9)
for i in range(9):
    Pf1_observe[i] = reads[1][i]/reads[0][i] # Attention au cas où 1obs = 0

# Définition du lambda initial qui sera utilisé pour les optimisations
lam_init = pinv_only(G.T,Pf1_observe)
# lam_init = [0.2,0.2,0.2,0.2,0.2]


### Traitement des données :

### Maximisation de la fonction de vraisemblance (minimisation de son opposé) : tentative avec toutes les méthodes
min_NM = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'Nelder-Mead',bounds=((0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
min_NMLog = scipy.optimize.minimize(likelihoodLog,lam_init,args=(reads,G),method= 'Nelder-Mead',bounds=((0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
min_Powell = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'Powell',bounds=((0,1),(0,1),(0,1),(0,1),(0,1)))
min_PowellLog = scipy.optimize.minimize(likelihoodLog,lam_init,args=(reads,G),method= 'Powell',bounds=((0,1),(0,1),(0,1),(0,1),(0,1)))
min_CG = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'CG')
min_CGLog = scipy.optimize.minimize(likelihoodLog,lam_init,args=(reads,G),method= 'CG')
min_BFGS = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'BFGS')
min_BFGSLog = scipy.optimize.minimize(likelihoodLog,lam_init,args=(reads,G),method= 'BFGS')
min_LBFGSB = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'L-BFGS-B',bounds=((0,1),(0,1),(0,1),(0,1),(0,1)))
min_LBFGSBLog = scipy.optimize.minimize(likelihoodLog,lam_init,args=(reads,G),method= 'L-BFGS-B',bounds=((0,1),(0,1),(0,1),(0,1),(0,1)))
min_TNC = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'TNC')
min_TNCLog = scipy.optimize.minimize(likelihoodLog,lam_init,args=(reads,G),method= 'TNC')
min_COBYLA = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'COBYLA')
min_COBYLALog = scipy.optimize.minimize(likelihoodLog,lam_init,args=(reads,G),method= 'COBYLA')
min_trust_constr = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'trust-constr',bounds=((0,1),(0,1),(0,1),(0,1),(0,1)))
min_trust_constrLog = scipy.optimize.minimize(likelihoodLog,lam_init,args=(reads,G),method= 'trust-constr',bounds=((0,1),(0,1),(0,1),(0,1),(0,1)))

# Jacobian is required
# min_NewtonCG = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'Newton-CG',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
# min_dogleg = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'dogleg',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
# min_trust_ncg = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'trust-ncg',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)))
# min_trust_exact = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'trust-exact',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
# min_trust_krylov = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'trust-krylov',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)))
# min_SLSQP = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'SLSQP',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})



### Analyses

def affichageErreur(freqOrigin,freqObs):
    erreurAbs = abs(freqOrigin-freqObs)
    erreurRelative = abs(freqOrigin-freqObs)
    for i in range(len(erreurRelative)):
        erreurRelative[i] = erreurRelative[i]/freq[i]
    print("Erreur absolue : ", erreurAbs)
    print("Erreur absolue moyenne % : ",round(np.mean(erreurAbs)*100,2))
    print("Erreur relative : ", erreurRelative)
    print("Erreur relative moyenne % : ", round(np.mean(erreurRelative)*100,2))

print("\n")
print("lambda à trouver : ",freq)
print("\n")

print("lambda initial donné pour minimisation = lambda prédit inversion simple : ", lam_init)
affichageErreur(freq,lam_init)
print("\n")

# print("lambda prédit minimisation Essr : ", minimizeEssr(G, Pf1_observe))
# affichageErreur(freq,minimizeEssr(G.T,Pf1_observe))
# print("\n")

print("lambda à trouver : ",freq)
print("lambda de l'algo Nelder-Mead : ", min_NM.x/np.sum(min_NM.x))
affichageErreur(freq,min_NM.x)
print("lambda de l'algo Nelder-Mead Log: ", min_NMLog.x/np.sum(min_NMLog.x))
affichageErreur(freq,min_NMLog.x)
print("\n")

print("lambda à trouver : ",freq)
print("lambda de l'algo Powell : ",min_Powell.x/np.sum(min_Powell.x))
affichageErreur(freq,min_Powell.x)
print("lambda de l'algo Powell Log : ",min_PowellLog.x/np.sum(min_PowellLog.x))
affichageErreur(freq,min_PowellLog.x)
print("\n")

print("lambda à trouver : ",freq)
print("lambda de l'algo COBYLA : ", min_COBYLA.x/np.sum(min_COBYLA.x))
affichageErreur(freq,min_COBYLA.x)
print("lambda de l'algo COBYLA Log : ", min_COBYLALog.x/np.sum(min_COBYLALog.x))
affichageErreur(freq,min_COBYLALog.x)
print("\n")

print("lambda à trouver : ",freq)
print("lambda de l'algo CG : ", min_CG.x/np.sum(min_CG.x)) #Renvoie equiproportion
affichageErreur(freq,min_CG.x)
print("lambda de l'algo CG Log : ", min_CGLog.x/np.sum(min_CGLog.x)) #Renvoie equiproportion
affichageErreur(freq,min_CGLog.x)
print("\n")

print("lambda à trouver : ",freq)
print("lambda de l'algo BFGS : ", min_BFGS.x/np.sum(min_BFGS.x))  #Renvoie equiproportion
affichageErreur(freq,min_BFGS.x)
print("lambda de l'algo BFGS Log : ", min_BFGSLog.x/np.sum(min_BFGSLog.x))  #Renvoie equiproportion
affichageErreur(freq,min_BFGSLog.x)
print("\n")

print("lambda à trouver : ",freq)
print("lambda de l'algo L-BFGS-B : ", min_LBFGSB.x/np.sum(min_LBFGSB.x))  #Renvoie equiproportion
affichageErreur(freq,min_LBFGSB.x)
print("lambda de l'algo L-BFGS-B Log : ", min_LBFGSBLog.x/np.sum(min_LBFGSBLog.x))  #Renvoie equiproportion
affichageErreur(freq,min_LBFGSBLog.x)
print("\n")

print("lambda à trouver : ",freq)
print("lambda de l'algo TNC : ", min_TNC.x/np.sum(min_TNC.x))  #Renvoie equiproportion
affichageErreur(freq,min_TNC.x)
print("lambda de l'algo TNC Log : ", min_TNCLog.x/np.sum(min_TNCLog.x))  #Renvoie equiproportion
affichageErreur(freq,min_TNCLog.x)
print("\n")

print("lambda à trouver : ",freq)
print("lambda de l'algo trust-constr ", min_trust_constr.x/np.sum(min_trust_constr.x))  #Renvoie equiproportion
affichageErreur(freq,min_trust_constr.x)
print("lambda de l'algo trust-constr Log ", min_trust_constrLog.x/np.sum(min_trust_constrLog.x))  #Renvoie equiproportion
affichageErreur(freq,min_trust_constrLog.x)
print("\n")
