###############################################################################
# @todo add Pilot2-splash-app disclaimer
# @authors  Helgi I. Ingolfsson <ingolfsson1@llnl.gov>
#           Changhoan Kim <kimchang@us.ibm.com>
###############################################################################

"""Contains definitions of all lipids both Martini CG and CHARMM AA."""

import logging
import MDAnalysis
import numpy
import os

LOGGER = logging.getLogger(__name__)

# Define all lipids
lipidTypeList = []
# Structure of lipid type definition list
#  0 Martini name (3-4 char)
#  1 Flippes -- 0 if none of these lipid flip-flop / 1 lipid flip-flop
#    (Note then we need to update the count each frame and can not use for
#     leaflet finding)
#  2 Lipid count upper/outer leaflet (in current frame) --  -1 is not set
#  3 Lipid count lower/inner leaflet (in current frame) --  -1 is not set
#  4 Lipid charge
#  5 Headgroup group ID (see lipidGroup list below)
#  6 Headgroup bond count (this is how many of the bond list should excluded
#    from tail order and to find headgroup vs tail bead names
#  7 Lipid tails (on letter code)
#  8 Bond list (used for order param calculations and to find
#  9 Representative headgroup bead - use for lipid tilt and to define leaflets

# WARNING lipid definitions have changed over time so a few lipids have different versions - that can be selected here
version = 0  # This is default current version included normally
# version = 1 # This is the PM initial version e.g. used in the normal (not new) complex PM simulations

# GroupID         0     1     2     3     4        5     6     7       8      9      10     11
lipidGroup = ["CHOL", "PC", "PE", "SM", "PS", "Glyco", "PI", "PA", "PIPs", "CER", "Lyso", "DAG"]

if version == 1:
    LIPID_DEFS = \
        os.path.join(Naming.dir_res('martini'), "lipiddefs/lipiddefs_v1.0.2.csv")
else:
    LIPID_DEFS = \
        os.path.join(Naming.dir_res('martini'), "lipiddefs/lipiddefs_v0.0.2.csv")

class LipidType():
    """Lipid type object - parses the lipid type string and knows what the lipid type is."""

    def __init__(self, defList):
        """
        Initialize the LipidType.

        :param defList: lipid type definition list see below,
        """
        self.ldef = defList

    def setCountUpper(self, count):
        """Set number of lipid in upper leaflet."""
        self.ldef[2] = count

    def setCountLower(self, count):
        """Set number of lipid in lower leaflet."""
        self.ldef[3] = count

    def __str__(self):
        """Print lipid info."""
        return self.ldef[0] + " - upper count " + str(self.ldef[2]) + " - lower count " + str(self.ldef[3])

    def getLeafletSelection(self, syst):
        """Get representative headgroup bead for lipid types that define the bilayer (i.e. none flip-flopping lipids)."""
        if self.ldef[1] == 0:
            # return syst.select_atoms("resname %s and name %s" % (self.ldef[0], self.getHeadBeads()))
            return syst.select_atoms("resname %s and name %s" % (self.ldef[0], self.ldef[9]))
        else:
            return syst.atoms[[]]

    def getNoneLeafletSelection(self, syst):
        """Get representative headgroup bead for lipid types don't define the bilayer (i.e. flip-flopping lipids)."""
        if self.ldef[1] == 1:
            return syst.select_atoms("resname %s and name %s" % (self.ldef[0], self.ldef[9]))
        else:
            return syst.atoms[[]]

    def getTailBeads(self):
        """Get all lipid tail bead names."""
        # ldef[8] is the bonded list
        # 1) Get all bonds for tails (Linker should be included so exclude first ldef[6] - 1
        cList = (self.ldef[8].split())[self.ldef[6] - 1:]
        # 2) Split the bonds into beads and flatten into 1D array
        # cList = numpy.array(map(lambda x:x.split("-"), cList)).flatten()
        cList = [item for sublist in [x.split("-") for x in cList] for item in sublist]
        # 3) make unique and return as a space separated string
        return " ".join(numpy.unique(cList))

    def getHeadBeads(self):
        """Get all lipid head bead names."""
        # ldef[8] is the bonded list
        # 1) Get all bonds for head (Linker should be included so include ldef[6]
        cList = (self.ldef[8].split())[:self.ldef[6]]
        # 2) Split the bonds into beads and flatten into 1D array
        # cList = numpy.array(map(lambda x:x.split("-"), cList)).flatten()
        cList = [
            item for sublist in [x.split("-") for x in cList]
            for item in sublist
        ]
        # 3) make unique and return as a space separated string
        return " ".join(numpy.unique(cList))


# NOTE: This really should be passed as a command line parameter. We're hard
# coding a lot of things in this code.
with open(LIPID_DEFS, "r") as lipidfile:
    lipiddefs = lipidfile.readlines()
    for lipiddef in lipiddefs:
        #  Only parse none comment lines (lines starting with #) and none empty lines
        if lipiddef[0:1] != "#" and len(lipiddef.strip()) > 0:
            lipid = lipiddef.rstrip().split(",")
            # Convert some of the definition entries to ints
            # Strip spaces off the ends
            lipid[0] = lipid[0].strip().strip("\"").strip()
            lipid[1] = int(lipid[1])
            lipid[2] = int(lipid[2])
            lipid[3] = int(lipid[3])
            lipid[4] = int(lipid[4])
            lipid[5] = int(lipid[5])
            lipid[6] = int(lipid[6])
            lipid[7] = lipid[7].strip().strip("\"").strip()
            lipid[8] = lipid[8].strip().strip("\"").strip()
            lipid[9] = lipid[9].strip().strip("\"").strip()
            lipidTypeList.append(lipid)
