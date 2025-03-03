import numpy as np

from GeneSets.geneSetObjects import GeneSet
from jsonParser import pjson


class GeneSetScore(dict): # Holds the scores of a given geneset in dictionary format for easy accession.

    def __init__(self,matrix,geneNamesList):
        
        self.matrix = matrix

        self.geneNamesList = geneNamesList # Get the gene names list. To get the gene names, use pjson.getGeneNames() function. 

        self.unique_genes = {} # Dictionary to hold the genes and scores. 

        for row_index, row in enumerate(self.matrix): # The scoring algorithm.
            score = 0

            cols = np.nonzero(row)[0]

            vals = row[cols]

            for col,val in zip(cols,vals):

                rows = np.nonzero(self.matrix[:, col])[0]

                for i in rows:
                    if i != row_index:
                        edges = np.count_nonzero(self.matrix[i])
                score += val * edges

            if self.geneNamesList[row_index] not in self.unique_genes:
                self.unique_genes[self.geneNamesList[row_index]] = score

        self.update(self.unique_genes) # Update self to the unique_genes dict

if __name__ == "__main__" : 

    newjson = pjson("/home/sadigungor/pathwayScoring/deneme.json")
    
    y = 0
    for i in newjson:
        
        print(i)
        y += 1
        newGeneSet = GeneSet(f"GeneSet{y}",newjson[i])
        
        _geneNames = newGeneSet.getGeneNames
        print(newGeneSet.getAsJson)

        
        newGeneSetScore = GeneSetScore(newGeneSet.getMatrix,_geneNames)
        
        print(newGeneSetScore) 
        






