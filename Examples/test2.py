# File name: Example/test2.py
from GenTopo.Coord import PDBobj
from GenTopo.Graph import MolGraph

mol = PDBobj("test.pdb") 
graph = MolGraph(mol, guessImpropers=True)

# accessing bonds 
for bond in graph.bonds: 
    pass   #process as needed by user  

# accessing angles  
for angle in graph.angles:
    pass 

# accessing dihedrals   
for dihedral in graph.dihedrals:
    pass 

