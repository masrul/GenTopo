from GenTopo.RingUtil import RingUtil
from collections import defaultdict
import numpy as np


class DihedralEstimator:
    """
    This class calculates dihedral angle of a molecule.
    The class is initialized with a python-object which
    contains coordinate information.

    Example:
        mol1 = DihedralEstimator(mol_obj)
        dihedralAngle = mol1.get([1,2,3,4])
    """

    def __init__(self, mol):
        self.mol = mol

    def get(self, dihedral):
        i, j, k, l = dihedral
        x = (self.mol.x[i - 1], self.mol.x[j - 1], self.mol.x[k - 1], self.mol.x[l - 1])
        y = (self.mol.y[i - 1], self.mol.y[j - 1], self.mol.y[k - 1], self.mol.y[l - 1])
        z = (self.mol.z[i - 1], self.mol.z[j - 1], self.mol.z[k - 1], self.mol.z[l - 1])

        x = np.array(x)
        y = np.array(y)
        z = np.array(z)

        b1 = (x[1] - x[0], y[1] - y[0], z[1] - z[0])
        b2 = (x[2] - x[1], y[2] - y[1], z[2] - z[1])
        b3 = (x[3] - x[2], y[3] - y[2], z[3] - z[2])

        b12 = np.cross(b1, b2)
        b23 = np.cross(b2, b3)
        b123 = np.cross(b12, b23)

        lb2 = np.linalg.norm(b2)
        ub2 = b2 / lb2

        angle = np.arctan2(np.dot(b123, ub2), np.dot(b12, b23))
        angle = angle * 57.2958

        return angle


class ImproperDihedralGenerator:
    """
    This class is used to generate improper dihedral
    in aromatic molecule. It uses following heuristics,
        + Center atom of improper is a part of aromatic ring
        + Dihedral angle < 5 degree (or user defined)
    """

    def __init__(self, mol):
        self.adjList = defaultdict(list)
        for inode, jnode in mol.bonds:
            self.adjList[inode].append(jnode)
            self.adjList[jnode].append(inode)

        self.RingUtil = RingUtil(mol.bonds)
        self.dihedralEstimator = DihedralEstimator(mol)

    def gen(self, cutoff=5.0):
        self.impDihedrals = []

        for atom in self.adjList:
            if len(self.adjList[atom]) == 3 and self.RingUtil.isRingMember(atom):
                tempList = tuple([atom] + self.adjList[atom])
                dihedralAngle = self.dihedralEstimator.get(tempList)
                if dihedralAngle < cutoff:
                    self.impDihedrals.append(tempList)

        return self.impDihedrals
