from GenTopo.RingUtil import RingUtil
from GenTopo.ImproperDihedral import ImproperDihedralGenerator
from GenTopo.Coord import PDBobj
import copy


class MolGraph:

    """
    It represent molecule as graph,
    which takes PDBobj/plain list as input. It generates
    angles, dihedrals, 1-4 and improper dihedrals using
    connectivity information.
    """

    def __init__(self, inp, guessImpropers=False, onlyCyclic14s=False):

        if isinstance(inp, PDBobj):
            self.coordObj = inp
            self.bonds = None
        else:
            self.coordObj = None
            self.bonds = inp

        self.gen(guessImpropers, onlyCyclic14s)

    def gen(self, guessImpropers, onlyCyclic14s):
        # generates necessary internal coordinates

        if self.bonds:
            self.nBonds = len(self.bonds)
            self.atoms = set()
            for i in range(self.nBonds):
                iatom = self.bonds[i][0]
                jatom = self.bonds[i][1]
                self.atoms.add(iatom)
                self.atoms.add(jatom)
            self.atoms = list(self.atoms)
            self.nAtoms = len(self.atoms)
            print("Number of Atoms: %-5d" % self.nAtoms)
        else:
            self.nAtoms = self.coordObj.nAtoms
            print("Number of Atoms: %-5d" % self.nAtoms)
            self.genBonds()

        self.nAngles = 0
        self.nDihedrals = 0
        self.nImDihedrals = 0
        self.nOneFours = 0

        self.genAngles()
        self.genDihedrals()
        self.genOneFours(onlyCyclic14s)

        if guessImpropers and self.coordObj:
            self.genImDihedrals()

    def genBonds(self):

        if self.coordObj.nBonds == 0:
            raise RuntimeError("Input file does not contains connectivity")
        else:
            self.bonds = copy.deepcopy(self.coordObj.bonds)
            self.nBonds = len(self.bonds)
        print("Number of Bonds: %-5d" % self.nBonds)

    def genAngles(self):
        self.angles = self.getNext(self.bonds)
        self.angles.sort()
        self.nAngles = len(self.angles)

        print("Number of Angles: %-5d" % self.nAngles)

    def genDihedrals(self):
        self.dihedrals = self.getNext(self.angles)
        self.dihedrals.sort()
        self.nDihedrals = len(self.dihedrals)

        print("Number of Dihedrals: %-5d" % self.nDihedrals)

    def genImDihedrals(self):
        if self.coordObj:
            self.imDihedrals = ImproperDihedralGenerator(self.coordObj).gen()
            self.nImDihedrals = len(self.imDihedrals)

            print("Number of Improper dihedrals: %-5d" % self.nImDihedrals)
        else:
            self.imDihedrals = []  # can not generate improper from bond list
            self.nImDihedrals = len(self.imDihedrals)

            print("Number of Improper dihedrals: %-5d" % self.nImDihedrals)

    def genOneFours(self, onlyCyclic=False):
        ringUtil = RingUtil(self.bonds)

        if not onlyCyclic:
            dihedrals = self.dihedrals
        else:

            dihedrals = []
            for dihedral in self.dihedrals:
                _, mid1, mid2, _ = dihedral

                if (
                    ringUtil.isRingMember(mid1)
                    and ringUtil.isRingMember(mid2)
                    and ringUtil.isFormRing(mid1, mid2)
                ):
                    dihedrals.append(dihedral)

        oneFours = [(i, j) for (i, _, _, j) in dihedrals]
        excludes = self.bonds + [(i, j) for (i, _, j) in self.angles]

        self.oneFours = []
        for oneFour in oneFours:
            if oneFour not in excludes:
                self.oneFours.append(oneFour)

        self.oneFours = list(set(self.oneFours))
        self.oneFours.sort()
        self.nOneFours = len(self.oneFours)

        print("Number of 1-4s: %-5d" % self.nOneFours)

    def getNext(self, currentList):
        """
        It generates next internal coordinate,
        for example, if angles is provided as currentList, it will
        return dihedral.
        """

        nextList = []

        n = len(currentList[0]) - 1

        for clist in currentList:
            for bond in self.bonds:

                if (bond[0] not in clist) and (bond[1] not in clist):
                    continue
                elif (bond[0] in clist) and (bond[1] in clist):
                    continue

                tempList = None
                if bond[0] == clist[n]:
                    tempList = clist + (bond[1],)

                elif bond[1] == clist[n]:
                    tempList = clist + (bond[0],)

                elif bond[0] == clist[0]:
                    tempList = (bond[1],) + clist

                elif bond[1] == clist[0]:
                    tempList = (bond[0],) + clist

                if tempList:
                    if tempList[0] > tempList[-1]:
                        tempList = tempList[::-1]
                    nextList.append(tempList)

        return list(set(nextList))

    def write(self, file_name):

        FH = open(file_name, "w")

        FH.write("#nBonds: %d\n" % self.nBonds)
        for (i, j) in self.bonds:
            FH.write("%6d  %6d\n" % (i, j))

        FH.write("\n#nAngles: %d\n" % self.nAngles)
        for (i, j, k) in self.angles:
            FH.write("%6d  %6d  %6d\n" % (i, j, k))

        FH.write("\n#nDihedrals: %d\n" % self.nDihedrals)
        for (i, j, k, l) in self.dihedrals:
            FH.write("%6d  %6d  %6d  %6d\n" % (i, j, k, l))

        if self.coordObj:
            FH.write("\n#nImDihedrals: %d\n" % self.nImDihedrals)
            for (i, j, k, l) in self.imDihedrals:
                FH.write("%6d  %6d  %6d  %6d\n" % (i, j, k, l))

        FH.write("#n14s: %d\n" % self.nOneFours)
        for (i, j) in self.oneFours:
            FH.write("%6d  %6d\n" % (i, j))

        FH.close()
