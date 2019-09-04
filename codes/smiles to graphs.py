# SMILES to graphs
# You can get an adjacency matrix with the function GetAdjacencyMatrix:

from rdkit import Chem

mol = Chem.GetMolFromSmiles(‘CCCCC’)
am = Chem.rdmolops.GetAdjacencyMatrix(mol) 

# After the adjacency matrix from rdkit you can do a lot with the networkx python library.