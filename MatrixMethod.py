import numpy as np

#X = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [0, 0, 0, 0, 0, 1, 1, 1, 1, 1], [1, 1, 0, 0, 0, 0, 0, 0, 1, 1]])
G = np.array([[1, 0, 1], [1,0,1], [1,0,0], [1,0,0], [1,0,0], [1,1,0], [1,1,0], [1,1,0], [1,1,1], [1,1,1]])
Pf1 = np.array([0.7,0.7,0.6,0.6,0.6,0.9,0.9,0.9,1,1])

print(G.shape)
print(Pf1.shape)

G_nbraws = np.size(G, 0)
Gd = np.c_[np.ones(G_nbraws), G]  # design matrix

Gd_t = Gd.transpose()
lam = np.linalg.pinv(Gd_t @ Gd) @ Gd_t @ Pf1
print("lambda: {}".format(lam))
#print("sol: {}".format(Gd @ theta))
# prediction for any new vector
#x = np.array([1.0, 8.0, 15.0])
#print("prediction: {}".format(x @ theta))