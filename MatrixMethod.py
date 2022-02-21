import numpy as np
import scipy

import Generator

#G = np.array([[1, 0, 1], [1,0,1], [1,0,0], [1,0,0], [1,0,0], [1,1,0], [1,1,0], [1,1,0], [1,1,1], [1,1,1]])
#Pf1 = np.array([0.7,0.7,0.6,0.6,0.6,0.9,0.9,0.9,1,1])
#lambda à trouver : [0.6, 0.3, 0.1]

#lambda_reel = [0.2,0.3,0.1,0.3,0.1]
#G = np.array([[1, 0, 1, 0, 1], [1, 0, 1, 0, 0], [1, 0, 0, 0, 1], [1, 0, 0, 0, 0], [1,0,0,1,0], [1,1,0,0,0], [1,1,0,1,1], [1,1,0,0,0], [1,1,1,1,1], [1,1,1,0,0]])

#avec entre 1 et 100 reads
#Pf1 = np.array([0.373134, 0.291667, 0.285714, 0.173913, 0.52, 0.506849, 0.9, 0.428571, 1, 0.65333])

#avec un plus grand nombre de reads (entre 1 et 1000)
#Pf2 = np.array([0.416666667,0.313207547,0.313765182,0.197916667,0.515151515,0.493761141,0.928571429,0.516129032,1,0.651639344])

#valeurs à faire varier
nb_genotypes = 5
nb_snips = 8

G_t = Generator.generateG(nb_genotypes, nb_snips)
print(G_t)
lam_reel = Generator.generateLambda(nb_genotypes)
freq = Generator.generateTrue_frequencies(G_t, lam_reel)
print("freq ", freq)
obs = Generator.generateReads_observ(nb_snips, freq)
print("obs ", obs)
Pf1_list = []
for i in range(len(obs[0])) :
    Pf1_list.append(obs[1][i]/obs[0][i])

print("pf1", Pf1_list)
Pf1 = np.array(Pf1_list)
print("lambda réel ",lam_reel)



#Méthode de minimisation de Essr pour trouver une valeur de lambda
def minimizeEssr(G_t, Pf1) :
    G = np.array(G_t).transpose()
    G_nbraws = np.size(G, 0)
    Gd = np.c_[np.ones(G_nbraws), G]  # design matrix
    Gd_t = Gd.transpose()
    lam = np.linalg.inv(Gd_t @ Gd) @ Gd_t @ Pf1 # closed-form solution
    return lam[1:6]

#Méthode pseudo-inversion directe de G
# Lambda = pseudoinv(G) * Pf1
def pinv_only(G_t, Pf1) :
    G = np.array(G_t).transpose()
    return np.linalg.pinv(G) @ Pf1

def erreur(lam_reel, lam_predicted) :
    erreur = abs(lam_reel - lam_predicted)
    for i in range(len(erreur)):
        if(lam_reel[i] != 0) :
            erreur[i] = erreur[i]/lam_reel[i]
    print("erreur : ",erreur)
    print("erreur moyenne :", round(np.mean(erreur)*100,2))
    return round(np.mean(erreur)*100,2)

lam_predicted1 = pinv_only(G_t, Pf1)
lam_predicted2 = minimizeEssr(G_t, Pf1)
# print("lambda réel", lam_reel)
print("lambda prédit inversion simple", lam_predicted1)
print("lambda prédit minimisation essr", lam_predicted2)
erreur(lam_reel, lam_predicted1)
erreur(lam_reel, lam_predicted2)
