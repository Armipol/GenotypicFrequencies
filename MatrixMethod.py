import numpy as np

#G = np.array([[1, 0, 1], [1,0,1], [1,0,0], [1,0,0], [1,0,0], [1,1,0], [1,1,0], [1,1,0], [1,1,1], [1,1,1]])
#Pf1 = np.array([0.7,0.7,0.6,0.6,0.6,0.9,0.9,0.9,1,1])
#lambda à trouver : [0.6, 0.3, 0.1]

lambda_reel = [0.2,0.3,0.1,0.3,0.1]
G = np.array([[1, 0, 1, 0, 1], [1, 0, 1, 0, 0], [1, 0, 0, 0, 1], [1, 0, 0, 0, 0], [1,0,0,1,0], [1,1,0,0,0], [1,1,0,1,1], [1,1,0,0,0], [1,1,1,1,1], [1,1,1,0,0]])

#avec entre 1 et 100 reads
Pf1 = np.array([0.373134, 0.291667, 0.285714, 0.173913, 0.52, 0.506849, 0.9, 0.428571, 1, 0.65333])

#avec un plus grand nombre de reads (entre 1 et 1000)
Pf2 = np.array([0.416666667,0.313207547,0.313765182,0.197916667,0.515151515,0.493761141,0.928571429,0.516129032,1,0.651639344])


print(G.shape)
print(Pf1.shape)

#Méthode de minimisation de Essr pour trouver une valeur de lambda

G_nbraws = np.size(G, 0)
Gd = np.c_[np.ones(G_nbraws), G]  # design matrix

Gd_t = Gd.transpose()

print("shape1 ", (np.linalg.pinv(Gd_t @ Gd) @ Gd_t).shape)
print("shape2 ", (Gd_t @ Gd).shape)

lam = np.linalg.pinv(Gd_t @ Gd) @ Gd_t @ Pf2
print("lambda: {}".format(lam))
print("sol: {}".format(Gd @ lam))
#prediction for any new vector
x = [0,0,0,0,0,0]
pf_predicted = []
for i in range(10) :
    x[0] = 1
    x[1] = G[i][0]
    x[2] = G[i][1]
    x[3] = G[i][2]
    x[4] = G[i][3]
    x[5] = G[i][4]
    print("x : ", x)
    print("prediction: {}".format(x @ lam))
    pf_predicted.append(x @ lam)

print(pf_predicted)

lam_predicted = np.linalg.pinv(G) @ pf_predicted

print(lam_predicted)

#Méthode inversion de G
# Lambda = pseudoinv(G) * Pf1

print("prediction lambda simple : ", np.linalg.pinv(G) @ Pf2)

erreurP = abs(lambda_reel-np.linalg.pinv(G) @ Pf2)
for i in range(len(erreurP)):
    erreurP[i] = erreurP[i]/lambda_reel[i]

print("erreur : ",erreurP)
print("moyenne erreur ", round(np.mean(erreurP)*100,2))