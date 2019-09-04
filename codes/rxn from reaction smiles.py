# Generating .rxn and/or .rdf fiiles from reaction smiles

from rdkit import Chem
from rdkit.Chem import AllChem

rxn = AllChem.ReactionFromSmarts('CC(=O)O.OCC>>CC(=O)OCC')
AllChem.ReactionToRxnBlock(rxn)