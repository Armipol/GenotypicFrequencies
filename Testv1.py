### IMPORTS
from math import *
import numpy as np
import scipy.optimize

### Définition des matrices

## Génotypes
G1 = np.array([[1]*10,[0,0,0,0,0,1,1,1,1,1],[1,1,0,0,0,0,0,0,1,1]])
G1 = G1.T

G2 = np.array([[1]*17,[1,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0],[0,1,1,0,1,0,0,0,1,1,1,0,1,0,0,0,1],[0,1,1,1,1,0,0,1,0,1,1,1,1,1,1,0,1],[1,0,0,1,0,0,1,1,1,1,1,1,0,1,0,0,1]])
G2 = G2.T

## Fréquences réelles
lamb1 = np.array([0.6,0.3,0.1])
lamb2 = np.array([0.17,0.21,0.47,0.06,0.09])

## Fréquences de 1 attendues
Pf1 = np.dot(G1,lamb1)
# print(Pf1)
Pf2 = np.dot(G2,lamb2)
# print(Pf2)

## Observations [[reads],[nb_1]
reads1 = [[15,5,8,64,31,93,63,1,66,49],[9,4,4,31,22,86,59,1,66,49]]
reads2 = [[129,574,36,816,685,425,133,752,526,456,809,706,482,744,903,716,305],[60,385,33,263,615,66,36,217,492,356,645,210,324,220,209,128,249]]

## Fréquences de départ pour la détermination du maximum de vraisemblance
lam1 = [0.3]*10
lam2 = [0.8]*17


### Implémentation des fonctions :

## Probabilité de tirer nb_1 parmi nb_tirages avec une probabilité de lamb_i
def probaBinomiale(nb_1,nb_tirages,lamb_i):
    proba = comb(nb_tirages,nb_1)*pow(lamb_i,nb_1)*pow(1-lamb_i,nb_tirages-nb_1)
    return proba

## Fonction de vraisemblance
def likelihood(lam,reads):
    vraisemblance = 1
    for i in range(len(reads[1])):
        proba = probaBinomiale(reads[1][i],reads[0][i],lam[i])
        #print(reads[1][i],reads[0][i],lam[i],proba)
        vraisemblance = vraisemblance * proba
    return -vraisemblance

# print(likelihood(lam,reads))

### Maximisation de la fonction de vraisemblance (minimisation de son opposé)

min1_NM = scipy.optimize.minimize(likelihood,lam1,args=(reads1),method= 'Nelder-Mead',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
# print(min1_NM)
min1_powell = scipy.optimize.minimize(likelihood,lam1,args=(reads1),method= 'Powell',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)))
# print(min1_powell)

min2_NM = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'Nelder-Mead',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
# print(min2_NM)
min2_Powell = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'Powell',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)))
# print(min2_powell)

### Analyses

## Premières valeurs

# G_1 = np.linalg.pinv(G1)
# lam1_final_NelderMead = np.dot(G_1,min1_NM.x)
# lam1_final_Powell = np.dot(G_1,min1_powell.x)
# print("fréquences attendues à trouver : ", Pf1)
# print("fréquences attendues Powell : ", min1_powell.x)
# print("fréquences attendues Nelder-Mead : ", min1_NM.x)
# print("lambda à trouver : ",lamb1)
# print("lambda de l'algo Powell : ", lam1_final_Powell) #,sum(lam2_final_Powell))
# print("lambda de l'algo Nelder-Mead : ", lam1_final_NelderMead) #,sum(lam2_final_NelderMead))
#
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
# print("lambda à trouver : ",lamb2)
# print("lambda de l'algo Powell : ", lam2_final_Powell,sum(lam2_final_Powell))
# print("lambda de l'algo Nelder-Mead : ", lam2_final_NelderMead,sum(lam2_final_NelderMead))
#
# erreur_lam2_Powell = abs(lam2_final_Powell - lamb2)/lamb2
# print("erreur moyenne Powell = ", np.mean(erreur_lam2_Powell)) # En moyenne erreur de 178%
# erreur_lam2_NM = abs(lam2_final_NelderMead - lamb2)/lamb2
# print("erreur moyenne Nelder-Mead = ", np.mean(erreur_lam2_NM)) # En moyenne erreur de 123%

## Essayons d'autres méthodes pour G2 :
# min2_CG = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'CG')
# min2_BFGS = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'BFGS')
# min2_LBFGSB = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'L-BFGS-B',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)))
# min2_TNC = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'TNC')
# min2_COBYLA = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'COBYLA')
# min2_trust_constr = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'trust-constr',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)))

# Jacobian is required
# min2_NewtonCG = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'Newton-CG',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
# min2_dogleg = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'dogleg',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
# min2_trust_ncg = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'trust-ncg',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)))
# min2_trust_exact = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'trust-exact',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})
# min2_trust_krylov = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'trust-krylov',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)))
# # min2_SLSQP = scipy.optimize.minimize(likelihood,lam2,args=(reads2),method= 'SLSQP',bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),options={'maxiter':5000,'maxfev':5000})

# print("fréquences attendues CG : ", min2_CG.x)
# print("fréquences attendues BFGS : ", min2_BFGS.x)
# print("fréquences attendues L-BFGS-B : ", min2_LBFGSB.x)
# print("fréquences attendues TNC : ", min2_TNC.x)
# print("fréquences attendues COBYLA : ", min2_COBYLA.x)
# print("fréquences attendues trust-constr ", min2_trust_constr.x)

# import autograd.numpy as np
# from autograd import grad, jacobian
#
# x = np.array([5,3], dtype=float)
#
# def cost(X): # X = (x,y)
#     return X[0]**2 / X[1] - np.log(X[1])
# #f1(x,y) = x²
# #f2(x,y) = 1/y
# #f3(x,y) = -log(y)
# jacobian_cost = jacobian(cost)
# x = np.array([1,1])
# print(jacobian_cost(np.array([x,x,x])))
#
# help(jacobian)

from jax import jacfwd, jacrev

W = [1]*len(lam1)
J = jacfwd(likelihood)(W)
print("jacfwd result, with shape", J.shape)
print(J)

J = jacrev(likelihood)(W)
print("jacrev result, with shape", J.shape)
print(J)
