import json
import numpy as np
from dataclasses import dataclass
import scanpy as sc
from icecream import ic
import matplotlib.pyplot as plt

"""

THIS CODE HERE IS A PROOF OF CONCEPT. IT DOES NOT REPRESENT ANY FINISHED OR OPTIMIZED ALGORITHM OR THE LIBRARY ITSELF.

Otherwise this would be a very sloppy work :P

"""
geneSetSayac = -1
geneIndexes = []
length = 0

@dataclass
class MatrixItem:
    value : float
    row : int
    column : int

jsonFile = "../test_data/big_genesets_relations.json"
GeneSetList = []
unique_genes = {}
with open(jsonFile) as f:
    js = json.load(f)

for i in js:
    print(i)
    geneSetSayac += 1
    ic(i)
    unique_identifiers = {}
    unique_genes = {}
    GeneNamesList = []
    GeneList = []
    sayac = 0
    geneSayac = 0
    n1 = 0
    score = 0

    for j in js[i]:
        n1 += 1
        GeneNamesList.append(j)
        for k in js[i][j]:
                
            unique_hash = hash(j) * hash(k) * hash(js[i][j][k])
            
            if unique_hash not in unique_identifiers:
                unique_identifiers[unique_hash] = sayac
                sayac += 1
            GeneList.append(MatrixItem(float(js[i][j][k]),n1-1,unique_identifiers[unique_hash]-1))


    n = len(unique_identifiers)
    matrix = np.zeros((n1,n))
    for item in GeneList:
        w,r,c = item.value,item.row,item.column
        matrix[r,c] = w
    for row_index,row in enumerate(matrix):
        cols = np.nonzero(row)[0]
        vals = row[cols]
        for col,val in zip(cols,vals):
            rows = np.nonzero(matrix[:, col])[0]
            
            for i in rows:
                if i != row_index:
                    edges = np.count_nonzero(matrix[i])
            score += val * edges
        
ic(matrix)
        if GeneNamesList[row_index] not in unique_genes:
            unique_genes[GeneNamesList[row_index]] = score
            geneSayac += 1
        #GeneSetList.append(MatrixItem(score,geneSetSayac,unique_genes[GeneNamesList[row_index]]))
        score = 0
        #ic(unique_genes[GeneNamesList[row_index]]-1)    
    
    #if len(unique_genes) > length:
        #length = len(unique_genes)
    geneIndexes.append(unique_genes)

adata = sc.read_h5ad("pbmc3k.h5ad")
expressionMatrix = np.array(adata.X)

#minval = np.min(expressionMatrix)
#maxval = np.max(expressionMatrix)
#expressionMatrix = (expressionMatrix - minval) / (maxval - minval)


def execute(genesetIndex,geneset):

    expScore = 0
    arrList = []
    arrBuffer = []
    ic(f"Geneset{genesetIndex}")
    with open("genesetSonuc33.txt","a") as f:
        arrList = []
        for var in adata.var_names:
            if var in geneset:
               arrList.append(geneset[var])
            else:
                arrList.append(0.0)# Some genes in the annotated data does not exist in some genesets(nor all of them). Is that a problem?

        for rowIndex,row in enumerate(expressionMatrix):
            expScore = np.dot(arrList,row) 
            f.write(f"Geneset{genesetIndex}\tCell{rowIndex}\t{expScore}\n")
            arrBuffer.append(expScore)
    
    return arrBuffer

arr = []  
for i,j in enumerate(geneIndexes):
    score = execute(i+1,j)
    arr.append(score)

arr = np.array(arr)

ic(arr)
