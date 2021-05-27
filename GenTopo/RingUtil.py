from collections import defaultdict


class RingUtil:
    """
    A utility class for handling aromatic molecule.
    It uses depth first search (DFS) to  determine if an atom
    is a part of aromatic ring. It also
    uses breadth first search (BFS) to determine if a bond is  part of aromatic
    ring.
    """

    def __init__(self, bonds):

        self.adjList = defaultdict(list)
        self.visited = defaultdict()
        for inode, jnode in bonds:
            self.adjList[inode].append(jnode)
            self.adjList[jnode].append(inode)

            self.visited[inode] = False
            self.visited[jnode] = False

        self.nVerts = len(self.adjList)

    def reset(self):
        for v in self.visited:
            self.visited[v] = False

    def isRingMember(self, current):
        self.start = current
        self.reset()

        return self.DFS(current, parent=-1)

    def DFS(self, current, parent):

        self.visited[current] = True
        for v in self.adjList[current]:
            if not self.visited[v]:

                if self.DFS(current=v, parent=current):
                    return True
            else:
                if v != parent and v == self.start:
                    return True

        return False

    def isFormRing(self, inode, jnode):

        self.reset()

        self.adjList[inode].remove(jnode)
        self.adjList[jnode].remove(inode)

        self.bfsChain = []

        self.BFS(inode)

        self.adjList[inode].append(jnode)
        self.adjList[jnode].append(inode)

        if jnode in self.bfsChain:
            return True
        else:
            return False

    def BFS(self, s):

        queue = []

        queue.append(s)
        self.visited[s] = True

        while queue:

            s = queue.pop(0)
            self.bfsChain.append(s)

            for i in self.adjList[s]:
                if not self.visited[i]:
                    queue.append(i)
                    self.visited[i] = True
