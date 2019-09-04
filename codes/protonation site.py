'''
>> I am trying to use rdkit descriptors to investigate protonation of 
    molecules in the gas phase while testing my procedure on Ritalin, 
    I try to calculate the protonation energy on different functional 
    groups in the moleculealthough my understanding tells me protonation 
    should occur on the secondary amine, energy calculations point at 
    the keto group as the favorable site.
    
>> you should try to modify the molecule conformation as little as possible 
   between the two calculations; see the example below.
'''

from rdkit import Chem
from rdkit.Chem import AllChem


test_ion1  =  Chem.MolFromSmiles('COC(=O)C(c1ccccc1)C1CCCC[NH2+]1')
test_ion1  =  Chem.AddHs(test_ion1)
AllChem.EmbedMolecule(test_ion1,  randomSeed=2)
prop1  =  AllChem.MMFFGetMoleculeProperties(test_ion1,
                                           mmffVariant="MMFF94")
ff1  =  AllChem.MMFFGetMoleculeForceField(test_ion1,  prop1)
ff1.Minimize(maxIts=1000)
print(ff1.CalcEnergy())


mn  =  test_ion1.GetSubstructMatch(Chem.MolFromSmarts("[NH2+]([H])[H]"))
n  =  mn[0]
na  =  test_ion1.GetAtomWithIdx(n)
mo  =  test_ion1.GetSubstructMatch(Chem.MolFromSmarts("O=CO"))
o  =  mo[0]
for  h  in  mn[1:]:
     test_ion2  =  Chem.RWMol(test_ion1)
     na.SetFormalCharge(0)
     test_ion2.RemoveBond(n,  h)
     test_ion2.RemoveAtom(h)
     oa  =  test_ion2.GetAtomWithIdx(o)
     oa.SetFormalCharge(1)
     oa.SetNoImplicit(True)
     oa.SetNumExplicitHs(1)
     Chem.SanitizeMol(test_ion2)
     test_ion2  =  Chem.AddHs(test_ion2,  addCoords=True)
     prop2  =  AllChem.MMFFGetMoleculeProperties(test_ion2,  mmffVariant="MMFF94")
     ff2  =  AllChem.MMFFGetMoleculeForceField(test_ion2,  prop2)
     print(ff2.CalcEnergy())