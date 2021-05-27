from GenTopo.Coord import PDBobj
from GenTopo.Graph import MolGraph
from GenTopo.GMXTopo import Topo


def main():

    mol = PDBobj("test.pdb")
    graph = MolGraph(mol, guessImpropers=True)
    gmx = Topo(mol, graph)

    gmx.write()


main()
