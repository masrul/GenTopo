import copy


class Topo:
    def __init__(self, mol, molGraph, itpFile="topol.top"):
        self.itpFile = itpFile
        self.molGraph = molGraph
        self.mol = mol
        self.atomTypes = mol.atomTypes
        self.atomQQs = mol.atomQQs
        self.assignTypes()
        self.setFuncID()

    def assignTypes(self):

        # bond types
        bondTypes = []
        for bond in self.molGraph.bonds:
            iatom = bond[0] - 1
            jatom = bond[1] - 1

            iType = self.atomTypes[iatom]
            jType = self.atomTypes[jatom]

            bondType = (iType, jType)

            if bondType not in bondTypes and bondType[::-1] not in bondTypes:
                bondTypes.append(bondType)
        self.bondTypes = bondTypes

        # angle types
        angleTypes = []
        for angle in self.molGraph.angles:
            iatom = angle[0] - 1
            jatom = angle[1] - 1
            katom = angle[2] - 1

            iType = self.atomTypes[iatom]
            jType = self.atomTypes[jatom]
            kType = self.atomTypes[katom]

            angleType = (iType, jType, kType)

            if angleType not in angleTypes and angleType[::-1] not in angleTypes:
                angleTypes.append(angleType)
        self.angleTypes = angleTypes

        # dihedral types
        dihedralTypes = []
        for dihedral in self.molGraph.dihedrals:
            iatom = dihedral[0] - 1
            jatom = dihedral[1] - 1
            katom = dihedral[2] - 1
            latom = dihedral[3] - 1

            iType = self.atomTypes[iatom]
            jType = self.atomTypes[jatom]
            kType = self.atomTypes[katom]
            lType = self.atomTypes[latom]

            dihedralType = (iType, jType, kType, lType)

            if (
                dihedralType not in dihedralTypes
                and dihedralType[::-1] not in dihedralTypes
            ):
                dihedralTypes.append(dihedralType)
        self.dihedralTypes = dihedralTypes

    def setFuncID(self):
        self.setDefaults()
        self.setBondFuncID()
        self.setAngleFuncID()
        self.setDihedralFuncID()
        self.setImDihedralFuncID()
        self.setOneFourFuncID()

    def setBondFuncID(self, bondFuncID=None):
        self.bondFuncID = bondFuncID

    def setAngleFuncID(self, angleFuncID=None):
        self.angleFuncID = angleFuncID

    def setDihedralFuncID(self, dihedralFuncID=None):
        self.dihedralFuncID = dihedralFuncID

    def setImDihedralFuncID(self, imDihedralFuncID=None):
        self.imDihedralFuncID = imDihedralFuncID

    def setOneFourFuncID(self, oneFourFunID=None):
        self.oneFourFunID = oneFourFunID

    def setDefaults(self, NBFunc=1, CombRule=2, GenPairs=True, FudgeFactors=(0.5, 0.5)):

        self.nbFunc = NBFunc
        self.combRule = CombRule

        if GenPairs:
            self.genPairs = "yes"
        else:
            self.genPairs = "no"
        self.fudgeFactors = FudgeFactors

    def write(self):
        self.topFH = open(self.itpFile, "w")
        self.writeDefaults()
        self.writeAtomTypes()
        self.writeBondTypes()
        self.writeAngleTypes()
        self.writeDihedralTypes()
        self.writeHeader()
        self.writeAtoms()
        self.writeBonds()
        self.writeAngles()
        self.writeDihedrals()
        self.writePairs()

        self.topFH.close()

    def writeDefaults(self):

        self.topFH.write("[ defaults ]\n")
        self.topFH.write(
            ";%9s  %10s  %10s  %10s  %10s\n"
            % ("nbfunc", "comb-rule", "gen-pairs", "fudgeLJ", "fudgeQQ")
        )
        self.topFH.write(
            "%10d  %10d  %10s  %10.4f  %10.4f\n"
            % (
                self.nbFunc,
                self.combRule,
                self.genPairs,
                self.fudgeFactors[0],
                self.fudgeFactors[1],
            )
        )

    def writeHeader(self):
        molName = "MOL"
        self.topFH.write("\n[ moleculetype ]\n")
        self.topFH.write(";name    nrexcl\n")
        self.topFH.write("%-s       3  ; Note: Adjust nrexcl\n\n" % (molName))

    def writeAtomTypes(self):

        _atomTypes = copy.deepcopy(self.atomTypes)
        _atomTypes = list(set(_atomTypes))
        _atomTypes = sorted(_atomTypes)

        self.topFH.write("\n[ atomtypes ]  ;nAtomTypes:%3d\n" % len(_atomTypes))

        self.topFH.write(
            "; name  at.num      mass     charge   ptype     sigma     epsilon\n"
        )
        for atype in _atomTypes:
            self.topFH.write(
                "%6s       0   0.00000    0.00000       A   0.00000     0.00000\n"
                % atype
            )

    def writeBondTypes(self):
        self.topFH.write("\n")
        self.topFH.write("[ bondtypes]   ; nBondTypes: %d\n" % len(self.bondTypes))

        if self.bondFuncID:
            self.topFH.write(";%5s  %6s  %6s\n" % ("atom1", "atom2", "func"))
            for (i, j) in self.bondTypes:
                self.topFH.write("%6s  %6s  %6d\n" % (i, j, self.bondFuncID))
        else:
            self.topFH.write(";%5s  %6s\n" % ("atom1", "atom2"))
            for (i, j) in self.bondTypes:
                self.topFH.write("%6s  %6s\n" % (i, j))

    def writeAngleTypes(self):
        self.topFH.write("\n")
        self.topFH.write("[ angletypes ]   ; nAngleTypes: %d\n" % len(self.angleTypes))

        if self.angleFuncID:
            self.topFH.write(
                ";%5s  %6s  %6s  %6s\n" % ("atom1", "atom2", "atom3", "func")
            )
            for (i, j, k) in self.angleTypes:
                self.topFH.write("%6s  %6s  %6s  %6d\n" % (i, j, k, self.angleFuncID))
        else:
            self.topFH.write(";%5s  %6s  %6s\n" % ("atom1", "atom2", "atom3"))
            for (i, j, k) in self.angleTypes:
                self.topFH.write("%6s  %6s  %6s\n" % (i, j, k))

    def writeDihedralTypes(self):
        self.topFH.write("\n")
        self.topFH.write(
            "[ dihedraltypes ]   ; nDihedrals: %d\n" % len(self.dihedralTypes)
        )

        if self.dihedralFuncID:
            self.topFH.write(
                ";%5s  %6s  %6s  %6s  %6s\n"
                % ("atom1", "atom2", "atom3", "atom4", "func")
            )
            for (i, j, k, l) in self.dihedralTypes:
                self.topFH.write(
                    "%6s  %6s  %6s  %6s  %6d\n" % (i, j, k, l, self.dihedralFuncID)
                )
        else:
            self.topFH.write(
                ";%5s  %6s  %6s  %6s\n" % ("atom1", "atom2", "atom3", "atom4")
            )
            for (i, j, k, l) in self.dihedralTypes:
                self.topFH.write("%6s  %6s  %6s  %6s\n" % (i, j, k, l))

    def writeAtoms(self):
        self.topFH.write("[ atoms ]   ; nAtoms: %d\n" % self.molGraph.nAtoms)
        self.topFH.write("; nr  type  resnr residue atom cgnr charge\n")

        for i in range(self.molGraph.nAtoms):
            id = i + 1
            self.topFH.write(
                "%6d %10s %3d %8s %8s %6d %14.8f\n"
                % (
                    id,
                    self.atomTypes[i],
                    self.mol.resIDs[i],
                    self.mol.resNames[i],
                    self.mol.symbols[i],
                    id,
                    self.atomQQs[i],
                )
            )

    def writeBonds(self):
        self.topFH.write("\n")
        self.topFH.write("[ bonds ]   ; nBonds: %d\n" % self.molGraph.nBonds)

        if self.bondFuncID:
            self.topFH.write(";%5s  %6s  %6s\n" % ("atom1", "atom2", "func"))
            for (i, j) in self.molGraph.bonds:
                self.topFH.write("%6d  %6d  %6d\n" % (i, j, self.bondFuncID))
        else:
            self.topFH.write(";%5s  %6s\n" % ("atom1", "atom2"))
            for (i, j) in self.molGraph.bonds:
                self.topFH.write("%6d  %6d\n" % (i, j))

    def writeAngles(self):
        self.topFH.write("\n")
        self.topFH.write("[ angles ]   ; nAngles: %d\n" % self.molGraph.nAngles)
        if self.angleFuncID:
            self.topFH.write(
                ";%5s  %6s  %6s  %6s\n" % ("atom1", "atom2", "atom3", "func")
            )
            for (i, j, k) in self.molGraph.angles:
                self.topFH.write("%6d  %6d  %6d  %6d\n" % (i, j, k, self.angleFuncID))
        else:
            self.topFH.write(";%5s  %6s  %6s\n" % ("atom1", "atom2", "atom3"))
            for (i, j, k) in self.molGraph.angles:
                self.topFH.write("%6d  %6d  %6d\n" % (i, j, k))

    def writeDihedrals(self):
        self.topFH.write("\n")
        self.topFH.write(
            "[ dihedrals ]   ; nDihedrals: %d\n" % self.molGraph.nDihedrals
        )
        if self.dihedralFuncID:
            self.topFH.write(
                ";%5s  %6s  %6s  %6s  %6s\n"
                % ("atom1", "atom2", "atom3", "atom4", "func")
            )
            for (i, j, k, l) in self.molGraph.dihedrals:
                self.topFH.write(
                    "%6d  %6d  %6d  %6d  %6d\n" % (i, j, k, l, self.dihedralFuncID)
                )
        else:
            self.topFH.write(
                ";%5s  %6s  %6s  %6s\n" % ("atom1", "atom2", "atom3", "atom4")
            )
            for (i, j, k, l) in self.molGraph.dihedrals:
                self.topFH.write("%6d  %6d  %6d  %6d\n" % (i, j, k, l))

        if self.molGraph.nImDihedrals == 0:
            return

        self.topFH.write("\n;Followings are imporper dihedrals\n")
        self.topFH.write(
            "[ dihedrals ]   ; nDihedrals: %d\n" % self.molGraph.nImDihedrals
        )
        if self.imDihedralFuncID:
            self.topFH.write(
                ";%5s  %6s  %6s  %6s  %6s\n"
                % ("atom1", "atom2", "atom3", "atom4", "func")
            )
            for (i, j, k, l) in self.molGraph.imDihedrals:
                self.topFH.write(
                    "%6d  %6d  %6d  %6d  %6d\n" % (i, j, k, l, self.imDihedralFuncID)
                )
        else:
            self.topFH.write(
                ";%5s  %6s  %6s  %6s\n" % ("atom1", "atom2", "atom3", "atom4")
            )
            for (i, j, k, l) in self.molGraph.imDihedrals:
                self.topFH.write("%6d  %6d  %6d  %6d\n" % (i, j, k, l))

    def writePairs(self):
        self.topFH.write("\n")
        self.topFH.write("[ paris ]   ; nPairs: %d\n" % self.molGraph.nOneFours)

        if self.oneFourFunID:
            self.topFH.write(";%5s  %6s  %6s\n" % ("atom1", "atom2", "func"))
            for (i, j) in self.molGraph.oneFours:
                self.topFH.write("%6d  %6d  %6d\n" % (i, j, self.oneFourFunID))
        else:
            self.topFH.write(";%5s  %6s\n" % ("atom1", "atom2"))
            for (i, j) in self.molGraph.oneFours:
                self.topFH.write("%6d  %6d\n" % (i, j))
