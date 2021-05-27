from GenTopo.Graph import MolGraph

# Define bonds as python list 
bonds = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]

# Create graph from bond list 
graph = MolGraph(bonds)

# Write/Process as needed 
graph.write("graph.dat")

