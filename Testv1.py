### IMPORTS
from math import *
import numpy as np
import random as rd
from scipy import stats
import scipy.optimize

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
        nb1_observe[i] = binomInverse(read,Pf[i])
    reads_observ[0] = nb_reads
    reads_observ[1] = nb1_observe
    return reads_observ

def binomInverse(n,p):
    x=0
    rand = rd.random()
    while(stats.binom.cdf(x,n,p)<rand):
        x+=1
    return x

### Définition des matrices

print("Valeurs générées aléatoirement :")
freq = np.array(generateLambda(5))
G = generateG(5,9)
print("Matrice G : ","\n", G)
G = G.T
reads = generateReads_observ(9,freq,G)
print("Reads & Nb1 : ", reads)
Pf1_attendu = G@freq
lam_init = [0.2,0.2,0.2,0.2,0.2]
Pf_init =G@lam_init
print("fréquences initiales à retrouver : ", freq, sum(freq))
print("fréquences attendues par position (Pf1) : ", Pf1_attendu)
print("lambda initial donné pour minimisation : ", lam_init)


## Génotypes
# G1 = np.array([[1]*10,[0,0,0,0,0,1,1,1,1,1],[1,1,0,0,0,0,0,0,1,1]])
# G1 = G1.T
# G2 = np.array([[1]*17,[1,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0],[0,1,1,0,1,0,0,0,1,1,1,0,1,0,0,0,1],[0,1,1,1,1,0,0,1,0,1,1,1,1,1,1,0,1],[1,0,0,1,0,0,1,1,1,1,1,1,0,1,0,0,1]])
# G2 = G2.T

## Fréquences réelles
# lamb1 = np.array([0.6,0.3,0.1])
# lamb2 = np.array([0.17,0.21,0.47,0.06,0.09])
## Fréquences de 1 attendues
# Pf1 = np.dot(G1,lamb1)
# # print(Pf1)
# Pf2 = np.dot(G2,lamb2)
# # print(Pf2)
## Observations [[reads],[nb_1]
# reads1 = [[15,5,8,64,31,93,63,1,66,49],[9,4,4,31,22,86,59,1,66,49]]
# reads2 = [[129,574,36,816,685,425,133,752,526,456,809,706,482,744,903,716,305],[60,385,33,263,615,66,36,217,492,356,645,210,324,220,209,128,249]]
## Fréquences de départ pour la détermination du maximum de vraisemblance
# lam1 = [0.6,0.3,0.1]
# lam2 = [0.2]*5


### Implémentation des fonctions :

## Probabilité de tirer nb_1 parmi nb_tirages avec une probabilité de lamb_i
# def probaBinomiale(nb_1,nb_tirages,pflambda_i):
#     nb_1 = nb_1.astype(int)
#     nb_tirages = nb_tirages.astype(int)
#     proba = comb(nb_tirages,nb_1)*pow(pflambda_i,nb_1)*pow(1-pflambda_i,nb_tirages-nb_1)
#     return proba

## Fonction de vraisemblance
def likelihood(lam,reads,G):
    lam = lam/np.sum(lam)
    Pflambda = np.dot(G,lam)
    vraisemblance = 1
    for i in range(len(reads[1])):
        proba = stats.binom.pmf(reads[1][i],reads[0][i],Pflambda[i])
        vraisemblance = vraisemblance * proba
    # print(lam,vraisemblance)
    return -vraisemblance

### Maximisation de la fonction de vraisemblance (minimisation de son opposé)
min_NM = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'Nelder-Mead',bounds=((0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
min_Powell = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'Powell',bounds=((0,1),(0,1),(0,1),(0,1),(0,1)))
min_CG = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'CG')
min_BFGS = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'BFGS')
min_LBFGSB = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'L-BFGS-B',bounds=((0,1),(0,1),(0,1),(0,1),(0,1)))
min_TNC = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'TNC')
min_COBYLA = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'COBYLA')
min_trust_constr = scipy.optimize.minimize(likelihood,lam_init,args=(reads,G),method= 'trust-constr',bounds=((0,1),(0,1),(0,1),(0,1),(0,1)))

erreurC = abs(freq-(min_COBYLA.x/np.sum(min_COBYLA.x)))
for i in range(len(erreurC)):
    erreurC[i] = erreurC[i]/freq[i]

erreurNM = abs(freq-(min_NM.x/np.sum(min_NM.x)))
for i in range(len(erreurNM)):
    erreurNM[i] = erreurNM[i]/freq[i]

erreurP = abs(freq-(min_Powell.x/np.sum(min_Powell.x)))
for i in range(len(erreurP)):
    erreurP[i] = erreurP[i]/freq[i]

print("lambda à trouver : ",freq)
print("lambda de l'algo Nelder-Mead : ", min_NM.x/np.sum(min_NM.x))
print("Erreur moyenne avec l'algo Nelder-Mead (en %) : ",round(np.mean(erreurNM)*100,2))
print("lambda de l'algo Powell : ",min_Powell.x/np.sum(min_Powell.x))
print("Erreur moyenne avec l'algo Powell (en %) : ",round(np.mean(erreurP)*100,2))
print("lambda de l'algo COBYLA : ", min_COBYLA.x/np.sum(min_COBYLA.x))
print("Erreur moyenne avec l'algo COBYLA (en %) : ",round(np.mean(erreurC)*100,2))
print("lambda de l'algo CG : ", min_CG.x/np.sum(min_CG.x)) #Renvoie equiproportion
print("lambda de l'algo BFGS : ", min_BFGS.x/np.sum(min_BFGS.x))  #Renvoie equiproportion
print("lambda de l'algo L-BFGS-B : ", min_LBFGSB.x/np.sum(min_LBFGSB.x))  #Renvoie equiproportion
print("lambda de l'algo TNC : ", min_TNC.x/np.sum(min_TNC.x))  #Renvoie equiproportion
print("lambda de l'algo trust-constr ", min_trust_constr.x/np.sum(min_trust_constr.x))  #Renvoie equiproportion


# min1_NM = scipy.optimize.minimize(likelihood,lam1,args=(reads1, G1),method= 'Nelder-Mead',bounds=((0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
# print(min1_NM)
# min1_Powell = scipy.optimize.minimize(likelihood,lam1,args=(reads1,G1),method= 'Powell',bounds=((0,1),(0,1),(0,1)))
# print(min1_powell)

# min2_NM = scipy.optimize.minimize(likelihood,lam2,args=(reads2,G2),method= 'Nelder-Mead',bounds=((0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
# print(min2_NM)
# min2_Powell = scipy.optimize.minimize(likelihood,lam2,args=(reads2,G2),method= 'Powell',bounds=((0,1),(0,1),(0,1),(0,1),(0,1)))
# print(min2_Powell)

### Analyses

## Premières valeurs

# G_1 = np.linalg.pinv(G1)
# lam1_final_NelderMead = np.dot(G_1,min1_NM.x)
# lam1_final_Powell = np.dot(G_1,min1_powell.x)
# print("fréquences attendues à trouver : ", Pf1)
# print("fréquences attendues Powell : ", min1_powell.x)
# print("fréquences attendues Nelder-Mead : ", min1_NM.x)
# print("lambda à trouver : ",lamb1)
# print("lambda de l'algo Powell : ", min1_powell.x/np.sum(min1_powell.x)) #,sum(lam2_final_Powell))
# print("lambda de l'algo Nelder-Mead : ", min1_NM.x/np.sum(min1_NM.x)) #,sum(lam2_final_NelderMead))
# print(likelihood(lam1,reads1,G1))
# print(likelihood(min1_powell.x,reads1,G1))
# print(likelihood(min1_NM.x,reads1,G1))

# erreur_lam1_Powell = abs(lam1_final_Powell - lamb1)/lamb1
# print("erreur moyenne Powell = ", np.mean(erreur_lam1_Powell)) # En moyenne 9,98% d'erreur
# erreur_lam1_NM = abs(lam1_final_NelderMead - lamb1)/lamb1
# print("erreur moyenne Nelder-Mead = ", np.mean(erreur_lam1_NM)) # En moyenne erreur de 10%

# Deuxièmes valeurs

# G2_1 = np.linalg.pinv(G2)
# lam2_final_NelderMead = np.dot(G2_1,min2_NM.x)
# lam2_final_Powell = np.dot(G2_1,min2_Powell.x)
# print("fréquences attendues à trouver : ",Pf2)
# print("fréquences attendues Powell : ", min2_Powell.x)
# print("fréquences attendues Nelder-Mead : ", min2_NM.x)

# erreurNM = abs(lamb1-(min1_NM.x/np.sum(min1_NM.x)))
# for i in range(len(erreurNM)):
#     erreurNM[i] = erreurNM[i]/lamb1[i]
#
# erreurP = abs(lamb1-(min1_Powell.x/np.sum(min1_Powell.x)))
# for i in range(len(erreurP)):
#     erreurP[i] = erreurP[i]/lamb1[i]
#
# print("reads : ",reads1)
# print("lambda à trouver : ",lamb1)
# print("lambda de l'algo Nelder-Mead : ", min1_NM.x/np.sum(min1_NM.x))
# print("Erreur moyenne avec l'algo Nelder-Mead (en %) : ",round(np.mean(erreurNM)*100,2))
# print("lambda de l'algo Powell : ",min1_Powell.x/np.sum(min1_Powell.x))
# print("Erreur moyenne avec l'algo Powell (en %) : ",round(np.mean(erreurP)*100,2))

# print(likelihood(lamb2,reads2,G2))
# print(likelihood(min2_Powell.x,reads2,G2))
# print(likelihood(min2_NM.x,reads2,G2))
#
# erreur_lam2_Powell = abs(lam2_final_Powell - lamb2)/lamb2
# print("erreur moyenne Powell = ", np.mean(erreur_lam2_Powell)) # En moyenne erreur de 178%
# erreur_lam2_NM = abs(lam2_final_NelderMead - lamb2)/lamb2
# print("erreur moyenne Nelder-Mead = ", np.mean(erreur_lam2_NM)) # En moyenne erreur de 123%

## Essayons d'autres méthodes pour G2 :
# min2_CG = scipy.optimize.minimize(likelihood,lam2,args=(reads2,G2),method= 'CG')
# min2_BFGS = scipy.optimize.minimize(likelihood,lam2,args=(reads2,G2),method= 'BFGS')
# min2_LBFGSB = scipy.optimize.minimize(likelihood,lam2,args=(reads2,G2),method= 'L-BFGS-B',bounds=((0,1),(0,1),(0,1),(0,1),(0,1)))
# min2_TNC = scipy.optimize.minimize(likelihood,lam2,args=(reads2,G2),method= 'TNC')
# min2_COBYLA = scipy.optimize.minimize(likelihood,lam2,args=(reads2,G2),method= 'COBYLA')
# min2_trust_constr = scipy.optimize.minimize(likelihood,lam2,args=(reads2,G2),method= 'trust-constr',bounds=((0,1),(0,1),(0,1),(0,1),(0,1)))
#
# erreurC = abs(lamb2-(min2_COBYLA.x/np.sum(min2_COBYLA.x)))
# for i in range(len(erreurC)):
#     erreurC[i] = erreurC[i]/lamb2[i]
#
# erreurNM = abs(lamb2-(min2_NM.x/np.sum(min2_NM.x)))
# for i in range(len(erreurNM)):
#     erreurNM[i] = erreurNM[i]/lamb2[i]
#
# erreurP = abs(lamb2-(min2_Powell.x/np.sum(min2_Powell.x)))
# for i in range(len(erreurP)):
#     erreurP[i] = erreurP[i]/lamb2[i]
#
# print("reads : ",reads2)
# print("lambda à trouver : ",lamb2)
# print("lambda de l'algo Nelder-Mead : ", min2_NM.x/np.sum(min2_NM.x))
# print("Erreur moyenne avec l'algo Nelder-Mead (en %) : ",round(np.mean(erreurNM)*100,2))
# print("lambda de l'algo Powell : ",min2_Powell.x/np.sum(min2_Powell.x))
# print("Erreur moyenne avec l'algo Powell (en %) : ",round(np.mean(erreurP)*100,2))
# print("lambda à trouver : ",lamb2)
# print("lambda de l'algo COBYLA : ", min2_COBYLA.x/np.sum(min2_COBYLA.x)) ##
# print("Erreur moyenne avec l'algo COBYLA (en %) : ",round(np.mean(erreurC)*100,2))
# print("lambda de l'algo CG : ", min2_CG.x/np.sum(min2_CG.x)) #Renvoie equiproportion
# print("lambda de l'algo BFGS : ", min2_BFGS.x/np.sum(min2_BFGS.x))  #Renvoie equiproportion
# print("lambda de l'algo L-BFGS-B : ", min2_LBFGSB.x/np.sum(min2_LBFGSB.x))  #Renvoie equiproportion
# print("lambda de l'algo TNC : ", min2_TNC.x/np.sum(min2_TNC.x))  #Renvoie equiproportion
# print("lambda de l'algo trust-constr ", min2_trust_constr.x/np.sum(min2_trust_constr.x))  #Renvoie equiproportion

# Jacobian is required
# min2_NewtonCG = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'Newton-CG',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
# min2_dogleg = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'dogleg',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
# min2_trust_ncg = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'trust-ncg',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)))
# min2_trust_exact = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'trust-exact',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
# min2_trust_krylov = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'trust-krylov',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)))
# # min2_SLSQP = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'SLSQP',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
