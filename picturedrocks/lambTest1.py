import sys

from singlecell import SingleCell, genEpsList


import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt


# Import the data.  
# 
# * ```df``` is all the Gene Expression Data, cells x genes
# * ```clustsdf``` is cluster labels
# * ```jdf``` is 2112 x 17156 cells x {genes, cluster}


bigData = False
normalized = True


if bigData:
    
    if (normalized):
        X = np.load("/home/ahsvargo/6mergedNormX.npy")
        y = np.load("/home/ahsvargo/6mergedy.npy")

        newrsums = np.sum(X, axis=1)
        ref = newrsums[0]
        if not np.all( np.abs(np.sum(X, axis=1) - ref) < 0.000000001):
            print("Warning: normalization failed.")

    else:
        X = np.load("/home/ahsvargo/6merged.npy")
        y = np.array([X[:,-1]]).T  # clusters
        X = X[:,:-1]   # expression levels


else:

# April dataset
#    df = pd.read_csv("/home/ahsvargo/markerData/03-sertoli/GeneExpressions_python_read.txt", 
#            delim_whitespace=True).set_index("x").T
#    clustsdf = pd.read_csv("/home/ahsvargo/markerData/03-sertoli/ClusterIDs_python_read.txt", 
#            sep=" ", header=None)

# July dataset
    df = pd.read_csv("/mnt/c/Users/Alexander Vargo/Downloads/GeneExpressions_python_read.txt",
            delim_whitespace=True).set_index("x").T
    clustsdf = pd.read_csv("/mnt/c/Users/Alexander Vargo/Downloads/ClusterIDs_python_read.txt", 
            delim_whitespace=True, header=None)

    clustsdf.columns = ["Gene", "Cluster"]
    clustsdf = clustsdf.set_index("Gene")
    df = df.join(clustsdf).sort_values("Cluster") # joint data frame

    df = df.as_matrix()
    X = df[:,:-1]   # expression levels
    y = np.array([df[:,-1]]).T  # clusters
    y = y -1

    if (normalized):
        rowsums = np.reshape(np.sum(X, axis=1), (X.shape[0], 1))
        norm = np.median(rowsums)
        X = (1.0*X)*norm/(1.0*rowsums)

myCell = SingleCell(X,y)
myCell.simpleCS(1.2)



# # some parameters
# num_cells, num_genes = X.shape
# K = y.max()  # number of clusters
# 
# print("Data uploaded sucessfully!", flush=True)
# 
# # The following tests a bunch of values of all of the parameters for a fixed value of $\lambda$.
# # 
# # Data are saved to named files in the current working directory: check the script for more specifics
# 
# # more parameters
# currLamb=float(sys.argv[1])
# numClusts=K
# 
# epsList = gz.genEpsList(maxEps=10.0**(-3))
# #lambMesh = list(np.linspace(0.01,0.1,num=10,endpoint=True))
# #alphaMesh = (np.linspace(0.01, 0.1, num=10, endpoint=True))
# 
# lambMesh = list(np.linspace(0.1,1.0,num=10,endpoint=True))
# alphaMesh = (np.linspace(0.05, 0.5, num=10, endpoint=True))
# 
# print(lambMesh, flush=True)
# print(alphaMesh, flush=True)
# 
# # store lambFac, alpha, #correct, and #support_genes
# ovrResults = np.zeros([len(alphaMesh)*len(lambMesh),4])
# fName = 'ovr-lamb{}.dat'.format(currLamb)
# rFile = open(fName, 'wb')
# 
# 
# for eye in range(len(lambMesh)):
#     
#     lambFac = lambMesh[eye]
#     
#     for jay in range(len(alphaMesh)):
#         
#         alpha = alphaMesh[jay]
#         index = len(lambMesh)*eye + jay
#         # index in the final list
#         
#         print("Working on alpha = {}, lambFac = {}".format(alpha, lambFac), flush=True)
#         
#         
#         # This finds the data for our current parameters
#         dubs,margins,epsInd = gz.runOvr1(y,X, currLamb=currLamb, lambFac=lambFac, 
# 			alpha=alpha, epsList=epsList, numClusts=numClusts)
#         
#         print("Found markers", flush=True)
# 
#         # the actual classification
#         classes = np.argmax(margins, axis=0) + 1
#         currCorrect = y[ classes[:,np.newaxis]==y ].shape[0]
#         
#         # find the support genes by truncating according to the best value of epsilion
#         bestEps = epsList[epsInd]
#         for w in dubs:
#             w[abs(w) < bestEps] = 0
#   
# 
#         # write the support genes to another data file
#         geneFile = "ovrGenes-lamb{}-lFac{}-alpha{}.dat".format(currLamb,lambFac,alpha)
#         gFile = open(geneFile, 'w')
#         support_genes = [np.nonzero(w)[0] for w in dubs]
#         for genes in support_genes:
#             for gene in genes:
#                 gFile.write("{} ".format(gene))
#             gFile.write("\n")
#         
#         gFile.close()
#         
#         ovrResults[index] = np.array([lambFac, alpha, currCorrect, 
# 		np.concatenate(support_genes).shape[0]])
#     
#     
# # write the results for the current value of lambFac
# np.savetxt(rFile, ovrResults, fmt=['%1.2f', '%1.2f', '%5i', '%5i'])
# rFile.close()
