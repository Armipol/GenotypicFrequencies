import pandas as pd

from FileReader import *

## Filepaths
filepath_positions = "D:/Rémi/Documents/IMT/3A/S10/EtudeTech/melange_simul_renom/positions_correspondance.txt"
filepath_reads = "D:/Rémi/Documents/IMT/3A/S10/EtudeTech/melange_simul_renom/SIMULS_READS_MIXTURES_fauxBAM_fauxREADS/reads_statistics.txt"
filepath_nucleotypes = "D:/Rémi/Documents/IMT/3A/S10/EtudeTech/melange_simul_renom/nucleotypes.txt"
filepath_mixtures = "D:/Rémi/Documents/IMT/3A/S10/EtudeTech/melange_simul_renom/simulated_mixtures_composition.txt"

## Data
data_utils = build_data_utils(filepath_positions, filepath_reads, filepath_nucleotypes, filepath_mixtures)
mixtures_dict = data_utils[0]
nucleotypes = data_utils[1]
column_nucleotypes_dict = data_utils[2]
positions_errors = data_utils[3]
harp_dict = data_utils[4]
mixtures_list = get_mixtures_list(mixtures_dict)

## Creation titres
indiceLigne = []
indiceCol = []
for i in range(12):
    indL = 'G'+str(i+1)
    indiceLigne.append(indL)
indiceLigne.append('nbReads')
indiceLigne.append('nb1')
for j in range(4404):
    indC = 'SNP'+str(j+1)
    indiceCol.append(indC)

## Lancement génération du Excel

def generationReal(fileExcel,listeMelanges): # Ici, listeMelanges doit être de la forme [1,4,11,...]
    for melange in listeMelanges:
        print("Mélange N° : ",melange)

        # On récupère G
        G_melange = generate_G_from_mix("Tm10"+str(melange).zfill(2), mixtures_dict, nucleotypes, column_nucleotypes_dict)

        # On récupère les reads
        reads_melange = reads_of_mix("Tm10"+str(melange).zfill(2), harp_dict, positions_errors)

        # On met en forme pour Excel
        matrice = G_melange.T
        matrice = np.vstack([matrice,reads_melange])

        # On exporte
        df = pd.DataFrame(matrice,index=indiceLigne,columns=indiceCol)
        with pd.ExcelWriter(fileExcel, engine="openpyxl", mode="a") as writer:
            df.to_excel(writer,sheet_name="Tm10"+str(melange).zfill(2))

generationReal("outputReal.xlsx",[4,6])
