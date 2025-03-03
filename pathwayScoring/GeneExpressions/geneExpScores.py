import scanpy as sc
import numpy as np

from GeneSets.geneSetScores import GeneSetScore 
from GeneSets.geneSetObjects import GeneSet

from jsonParser import pjson

def score(adata,geneSetIndex,geneSetScore):

    expMatrix = np.array(adata.X)

    expScore = 0
    arrList = []
    arrBuffer = []

    for var in adata.var_names:
        if var in geneSetScore:
            arrList.append(geneSetScore[var])
            
        else:
            arrList.append(0.0)
    for rowIndex,row in enumerate(expMatrix):

        expScore = np.dot(arrList,row)
        arrBuffer.append(expScore)

    return arrBuffer

if __name__ == "__main__" :

    newjson = pjson("/home/sadigungor/pathwayScoring/big_genesets_relations.json")
    

    for i in newjson:

        newGeneSet = GeneSet(f"{i}",newjson[i])

        print(newGeneSet.getID)

        _geneNames = newGeneSet.getGeneNames

        newGeneSetScore = GeneSetScore(newGeneSet.getMatrix,_geneNames)

        adata = sc.read_h5ad("/home/sadigungor/Desktop/pathWayScoring/test_data/pbmc3k.h5ad")

        #print(score(adata,"GeneSet1",newGeneSetScore))
