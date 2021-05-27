from GenTopo.PeriodicTable import elements as elements
from GenTopo.Warning import non_orthogonal_box, connectivity_missing


class PDBobj:
    """
    This a molecule class, which is initialized with a pdb file.
    This object contains coordinates and connectivity.
    """

    def __init__(self, coordFile, box=None, lpbc=(True, True, True)):
        self.coordFile = coordFile
        self.box = box
        self.lpbc = lpbc
        self.ffPresent = True

        if self.box and len(self.box) != 3:
            raise RuntimeError(non_orthogonal_box)

        self.process()

    def process(self):
        self._symbols = []
        self._x = []
        self._y = []
        self._z = []
        self._resNames = []
        self._resIDs = []
        self.bonds = []
        self.nBonds = 0
        self.nAtoms = 0

        self.atomTypes = []
        self.atomQQs = []

        coordFH = open(self.coordFile, "r")
        self.lines = coordFH.readlines()
        coordFH.close()

        self.readCoords()
        self.readBonds()
        self.readFF()

    def readBonds(self):

        foundCONECT = False
        for line in self.lines:
            if line.startswith("CONECT"):
                foundCONECT = True
                keys = line.split()
                iatom = int(keys[1])
                for i in range(2, len(keys)):
                    jatom = int(keys[i])
                    if iatom < jatom:
                        self.bonds.append((iatom, jatom))
                    else:
                        self.bonds.append((jatom, iatom))

        self.bonds = list(set(self.bonds))
        self.bonds.sort()
        self.nBonds = len(self.bonds)

        if not foundCONECT:
            print(connectivity_missing)
            self.genBonds()

    def genBonds(self):

        self.radii = []
        self.bonds = []

        # assign vdw radius
        for symbol in self._symbols:
            radius = elements[symbol[0]]["vdw_radius"]
            self.radii.append(radius)

        for iatom in range(self.nAtoms - 1):
            iradius = self.radii[iatom]
            for jatom in range(iatom + 1, self.nAtoms):
                jradius = self.radii[jatom]
                rcut = 0.6 * (iradius + jradius)
                dx = self._x[iatom] - self._x[jatom]
                dy = self._y[iatom] - self._y[jatom]
                dz = self._z[iatom] - self._z[jatom]
                r = dx * dx + dy * dy + dz * dz
                if r < rcut ** 2:
                    if iatom < jatom:
                        bond = (iatom + 1, jatom + 1)
                    else:
                        bond = (jatom + 1, iatom + 1)

                    if bond not in self.bonds:
                        self.bonds.append(bond)
        self.nBonds = len(self.bonds)

    def readCoords(self):
        for line in self.lines:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                self._symbols.append(line[12:16].strip())
                self._resNames.append(line[17:20])
                self._resIDs.append(int(line[22:26]))
                self._x.append(float(line[30:38]))
                self._y.append(float(line[38:46]))
                self._z.append(float(line[46:54]))

                self.nAtoms += 1

    def readFF(self):
        # check if FF information is here

        for line in self.lines:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                try:
                    atomType, qq = self.splitLine(line)
                except:
                    self.ffPresent = False
                    return

        for line in self.lines:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                atomType, qq = self.splitLine(line)
                self.atomTypes.append(atomType)
                self.atomQQs.append(qq)

    @staticmethod
    def splitLine(line):
        line = line[80:]
        keys = line.split()

        atomType = keys[0]
        if len(keys) > 1:
            qq = float(keys[1])
        else:
            qq = 0.000

        return atomType, qq

    def applyPBC(self, dx, dy, dz):

        if self.lpbc[0]:
            dx -= self.box[0] * round(dx / self.box[0])

        if self.lpbc[1]:
            dy -= self.box[1] * round(dy / self.box[1])

        if self.lpbc[2]:
            dz -= self.box[2] * round(dy / self.box[2])

        return dx, dy, dz

    # Property decorators are used to avoid accidental overwrite of Coord Object
    @property
    def symbols(self):
        return self._symbols

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def z(self):
        return self._z

    @property
    def resIDs(self):
        return self._resIDs

    @property
    def resNames(self):
        return self._resNames
