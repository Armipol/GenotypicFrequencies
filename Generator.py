import numpy as np
import random as rd
from scipy import stats

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
    G = []
    for i in range(nb_genotypes):
        geno = []
        for j in range(nb_snips):
            geno.append(rd.randint(0,1))
        G.append(geno)
    return G

def generateTrue_frequencies(G, lam) :
    freq = []
    for j in range(len(G[0])):
        ones_freq = 0
        for i in range(len(G)):
            if (G[i][j] == 1):
                ones_freq += lam[i]
        freq.append(ones_freq)
    return freq

def generateReads_observ(nb_snips,freq):
    reads_observ = []
    nb_reads = []
    nb1_observe = []
    for j in range(nb_snips):
        read = rd.randint(0,10000)
        nb_reads.append(read)
        nb1_observe.append(binomInverse(read,freq[j]))
    reads_observ.append(nb_reads)
    reads_observ.append(nb1_observe)
    return reads_observ

def binomInverse(n,p):
    x=0
    rand = rd.random()
    while(stats.binom.cdf(x,n,p)<rand):
        x+=1
    return x

# freq = np.array(generateLambda(7))
# G = np.array(generateG(7,5))
# print(G)
# G = G.T
# reads = generateReads_observ(5,freq)
# print(freq, sum(freq))
# print(G@freq)
# print(reads)

# def generator(nb_genotypes,nb_snips):
#     frequences = generateLambda(nb_genotypes)
#     G = generateG(nb_genotypes,nb_snips)
#     Pf1_attendues = G & frequences
