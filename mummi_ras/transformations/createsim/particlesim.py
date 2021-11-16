###############################################################################
# @todo add Pilot2-splash-app disclaimer
# @authors  Helgi I. Ingolfsson <ingolfsson1@llnl.gov>
#           Changhoan Kim <kimchang@us.ibm.com>
###############################################################################

"""A particle simulation instance - keeps track of what is in the simulation and all the sub-simulations."""

import logging
import os
import re
import datetime
import pickle as pickle
import numpy as np
import glob

import mummi_core
import mummi_ras
from mummi_core.utils import Naming
import mummi_core.utils.utilities as utils

# Logger has to be initialized the first thing in the script
LOGGER = logging.getLogger(__name__)

# Innitilize MuMMI if it has not been done before
MUMMI_ROOT = mummi_core.init(True)

# ------------------------------------------------------------------------------

pSimFileName = "psim.npz"
pSimSubPatchFileUpper = "psimSubPatch_upper.pat"
pSimSubPatchFileLower = "psimSubPatch_lower.pat"
lGridName = "lipid_counts.npz"

dirddcMD = Naming.dir_res('ddcMD')  # dirResources+"ddcMD/"
dirMartini = Naming.dir_res('martini')  # dirResources+"Martini/"
dirCHARMM = Naming.dir_res('charmm')  # dirResources+"CHARMM/"
dirKRASStructures = Naming.dir_res('structures')  # dirResources+"KRAS-structures/"


class ParticleSim:
    """Defines a specific particle simulation instance - It knows what went into the simulations and their status."""
    class_version = '3.0'  # Update if you add new variables


    def __init__(self, lipidTypes, lipidUpperRatioArray, lipidLowerRatioArray, lipidSuPatches=[1, 1],
                 proteinUnits=[], asym=0, simName="", lipidFullAsignment=False):
        """
        Initialize the ParticleSimulation - only modify the default values where nessisary.

        :param lipidTypes: An array of lipid resnames, needs to be "supported" and listed in lipiddef.lipidTypeList,
                length n, string
        :param lipidUpperRatioArray: upper leaflet subpatch consentrations for all subpatches, dim (m,n), int
        :param lipidLowerRatioArray: lower leaflet subpatch consentrations for all subpatches, dim (m,n), int
               The sum of each subpatch should be modular to the number of lipids in that subpatch
               e.g. if 64 lipids in each subpatch the sum mod 64 = 0, else we can have rounding errors.
        :param lipidSuPatches: a list of x,y number of subpatches (e.g. [5,5] for x25 subpatches)
        :param proteinUnits: List of protein units - each one has the structures .gro and posision of that protein units
        :param asym: If bilayer as diffrent number of lipids in out/inner leaflet - passed to insane.
               Int number of lipids to be randomly removed in the upper (positive number) or lower (negative numbers) leaflets
        """
        self.version = self.__class__.class_version   # Opject know their version and can be updated if of older version

        # Variables for lipid assingmnet
        self.lipidFullAsignment = lipidFullAsignment # True if using full assignment of each lipid or False for probabilistic assignment
        self.lipidTypes = np.asarray(lipidTypes)
        self.lipidUpperRatioArray = np.asarray(lipidUpperRatioArray)
        self.lipidLowerRatioArray = np.asarray(lipidLowerRatioArray)
        self.lipidSuPatches = np.asarray(lipidSuPatches)

        self.rasList = [] # not used anymore, np.asarray(rasList)
        self.rasPos = [] # not used anymore, np.asarray(rasPos)
        self.proteinUnits = proteinUnits
        self.asym = asym

        if len(lipidTypes) != self.lipidUpperRatioArray.shape[1] or len(lipidTypes) != self.lipidLowerRatioArray.shape[1]:
            error = "lipidTypes and lipidRatios (upper and lower) need to be of same length"
            LOGGER.error(error)
            raise AttributeError(error)

        if self.lipidUpperRatioArray.shape[0] != self.lipidLowerRatioArray.shape[0] or \
           lipidSuPatches[0] * lipidSuPatches[1] != self.lipidLowerRatioArray.shape[0]:
            error = "Wrong number of subpatches"
            LOGGER.error(error)
            raise AttributeError(error)

        # Set default params
        # Sim size - for production use  31x31 nm or 3200 lipids (1600,1600 - 40x40 per leaflet
        #  - this will support x25 8x8 lipid subpatches)
        self.sizeX = 31  # For small 450 lipid sims use 12 (some old test sims used 15.5)
        self.sizeY = 31  # For small 450 lipid sims use 12 (some old test sims used 15.5)
        self.sizeZ = 16  # use 13 for Martini 3 (less gap between protein and bilayer) used 12 for bilayer only sims
        self.ioncons = 0.15  # 150 mM NaCl
        self.water = "W"  # Use PW for polirizable and WN for Martini3 normal water
        self.defTopFile = "system-cg-default-M2.top"  # Alt system-cg-default-M2PW.top or system-cg-default-M3.top

        # Default names for RAS and RAS-RAF Martini parameters - 4B and 4A
        self.defrasItp = "{}/RAS_scfix-CYFpos".format(dirMartini)  # New scfix RAS for C3
        self.defrasrafItp = "{}/RAS_RAF_TERNARY_scfix-CYFpos".format(dirMartini)  # New scfix RAS-RAF for C3
        self.defras4aItp = "{}/RAS4A_scfix-CYFpos".format(dirMartini)
        self.defras4arafItp = "{}/RAS4A_RAF_scfix-CYFpos".format(dirMartini)

        # @TODO add RAS4A to feedback
        self.rasItp = self.defrasItp
        self.rasrafItp = self.defrasrafItp
        self.ras4aItp = self.defras4aItp
        self.ras4arafItp = self.defras4arafItp

        # Default is empty list - that is no feedback has been done and use standard referance structure
        # This is a list as there can be 1 or 2 structure per protein if changes was made in HVR or CRD or both
        self.prefRefPdbs = []
        self.ddcMDParams = dirddcMD # The default - no feedback is form the mummi_resources
        self.ddcMDParams_time = "-1" # Keeps the time of the current ddcMD folder (so we can select the latest)
        dirFeedback = Naming.dir_root('feedback-aa')

        # Check if AA to CG feedback for all protein types - use the newest of each
        # For ddcMD params use the newest of all
        # RAS-ONLY .itp (4B)
        self.rasItp = self.getProteinFeedback(dirFeedback, "RAS-ONLY", self.rasItp, "RAS_scfix-CYFpos")
        # RAS-RAF .itp (4B)
        self.rasrafItp = self.getProteinFeedback(dirFeedback, "RAS-RAF", self.rasrafItp, "RAS_RAF_TERNARY_scfix-CYFpos")
        # RAS4A-ONLY .itp (4A)
        self.ras4aItp = self.getProteinFeedback(dirFeedback, "RAS4A-ONLY", self.ras4aItp, "RAS4A_scfix-CYFpos")
        # RAS4A-RAF .itp (4A)
        self.ras4arafItp = self.getProteinFeedback(dirFeedback, "RAS4A-RAF", self.ras4arafItp, "RAS4A_RAF_scfix-CYFpos")

        LOGGER.info("Final itp's after feedback check: RAS4B {}, RAS4B-RAF {}, RAS4A {}, RAS4A-RAF {}, ddcMD folder {}, ref pdbs {}".format(self.rasItp, self.rasrafItp, self.ras4aItp, self.ras4arafItp, self.ddcMDParams, str(self.prefRefPdbs)))

        self.cgnt = 0  # If > 0 this sets a limit to the domain max size (e.g. will use gmx mdrun -nt "nt")
        self.descriptiveName = "NA"  # Only used for testing simulations
        # Generate name
        if len(simName) > 0:
            self.name = simName
        else:
            self.name = "sim_%s_KRAS_%i_h%s" % (datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S"), len(rasList), hash(self))


    def getProteinFeedback(self, dirFeedback, proteinName, proteinOldItp, feedbackName):
        """Get the newest feedback .itp and reference .pdb files if available - saving relative paths only - use getItp """
        proteinNewItp = ""
        #feedbackItpListFull = glob.glob("{}/{}*[0-9].itp".format(dirFeedback, feedbackName))  # Force a number [0-9] after . so not the include PULL .itp's
        #feedbackItpList = [os.path.basename(x) for x in feedbackItpListFull] # get file name only - not full path
        # @WARNING this is still full path - need to restructure to make relative paths - i.e. add all latest files to one folder so default vs feedback is not in diffrent location
        feedbackItpList =glob.glob("{}/{}*[0-9].itp".format(dirFeedback, feedbackName))  # Force a number [0-9] after . so not the include PULL .itp's
        if feedbackItpList == []:
            LOGGER.info("Feedback AA to CG - no new .itp's found in {} use inital/default {} .itp - {}".format(dirFeedback, proteinName, proteinOldItp))
            return proteinOldItp
        else:
            sortList = sorted(feedbackItpList, key = lambda x: x[-19:-4])
            proteinNewItp = sortList[-1][0:-4] # last is newest file - and store file without the .itp ending (this is added in .top file)
            LOGGER.info("Feedback AA to CG - found these new .itp's {} use newest as {} .itp - {}".format(str(feedbackItpList), proteinName, proteinNewItp))

            # check on PULL file
            fullFile = "{}-PULL.itp".format(proteinNewItp)
            if not os.path.isfile(fullFile):
                error = "Feedback AA to CG - in folder {} could not find {} -PULL.itp assosiated with {} .itp".format(dirFeedback, proteinName, proteinNewItp)
                LOGGER.error(error)
                raise NameError(error)

            # Find correct ref .pdb files and add to prefRefPdbs list.
            pdb_ref_files_full = glob.glob("{}/*{}.pdb".format(dirFeedback, proteinNewItp[-15:]))
            pdb_ref_files = [os.path.basename(x) for x in pdb_ref_files_full] # get file name only - not full path
            if len(pdb_ref_files) == 0:
                error = "Feedback AA to CG - in folder {} could not find referance .pdb(s) for new .itp {}".format(dirFeedback, proteinNewItp)
                LOGGER.error(error)
                raise NameError(error)
            else:
                self.prefRefPdbs = self.prefRefPdbs + pdb_ref_files
                LOGGER.info("Feedback AA to CG - new ref protein AA structure(s) for updated .itp {} are {}".format(proteinNewItp, str(pdb_ref_files)))

            # get dir with correct ddcMD parameter files, we need to point to the newest feedback for the most up to date ddcMD params
            # @TODO do this simpler
            latest_feedback_time = sorted([self.ddcMDParams_time, proteinNewItp[-15:]])
            if latest_feedback_time[-1] == proteinNewItp[-15:]:
                # @WARNING this is full path not relaive as it's use right away in creatsims and default is also full path
                self.ddcMDParams = "{}/ddcMDParams_{}".format(dirFeedback, proteinNewItp[-15:])
                self.ddcMDParams_time = proteinNewItp[-15:]
                LOGGER.info("Feedback AA to CG - newer ddcMD parameters found for updated .itp {} currently use {}".format(proteinNewItp, str(self.ddcMDParams)))

            return proteinNewItp


    def create(self, simpath):
        """Creates Particle simulation directory and as well as inner/outer sub-grid concentration files."""
        # @Note Workflow should have created this dir, but if calling in local mode create
        if not os.path.isdir(simpath):
            utils.sys_call('mkdir -p ' + simpath, test=False)
        os.chdir(simpath)
        # Make subpatch files
        with open(pSimSubPatchFileUpper, 'w') as OUT:
            OUT.write("%i %i\n" % (self.lipidSuPatches[0], self.lipidSuPatches[1]))
            for i in self.lipidUpperRatioArray:
                current = ""
                for j in range(len(i)):
                    if i[j] > 0:
                        current += str(self.lipidTypes[j]) + ":" + str(i[j]) + " "
                OUT.write("%s\n" % current.strip())

        with open(pSimSubPatchFileLower, 'w') as OUT:
            OUT.write("%i %i\n" % (self.lipidSuPatches[0], self.lipidSuPatches[1]))
            for i in self.lipidLowerRatioArray:
                current = ""
                for j in range(len(i)):
                    if i[j] > 0:
                        current += str(self.lipidTypes[j]) + ":" + str(i[j]) + " "
                OUT.write("%s\n" % current.strip())
        # Done
        LOGGER.info("-> Particle simulation %s has been created", self.name)


    def save(self):
        """Save Particle simulation to disk - NOTE has to be called to make any changes presistant."""
        pickle.dump(self, open(pSimFileName, 'wb'), pickle.HIGHEST_PROTOCOL)


    def upgrade(self):
        """Upgrade object to be compatible with nevest vesrion, adding default values to all new params."""
        version = float(self.version)
        if version < 1.:
            self.version = '1.0'
        if version < 2.:
            self.version = '2.0'
            self.water = "?"  # For old sims not defined check system.top file
            self.defTopFile = "?"  # For old sims not defined check system.top file
            self.rasItp = "?"  # For old sims not defined check system.top file
            self.cgnt = 0
            self.descriptiveName = "NA"
        if version < 3.:
            self.rasrafItp = "?"
            self.proteinUnits = [] # used self.rasList and self.rasPos - if really important can reconstruct here?


def load_sim(simName):
    """Load simulation with simName.

    :returns: a ParticleSim object
    """

    LOGGER.info("sname= {}, pSim= {}, cDir= {} ".format(simName, pSimFileName, os.getcwd()))
    pSim = pickle.load(open(simName+'/'+pSimFileName, 'rb'), encoding='bytes')
    if not hasattr(pSim, 'version'):
        pSim.version = '0.0'
    if pSim.version != pSim.class_version:
        pSim.upgrade()
    return pSim


def load_sims():
    """
    Load all previously run simulations - in current directory.

    :returns: A list of ParticleSim in the current directory
    """
    simDirs = os.listdir("./")
    simList = []
    for i in simDirs:
        if i[0:3] == "sim" or i[0:7] == "pfpatch":  # Loead only sim* or pfpatch directories
            simList.append(load_sim(i))
    return simList


''' Old code not used after refactor
def get_sim_dir():
    """Get current dir."""
    return Naming.dir_sim('cg', simname=self.name)


def get_splash_root_sim_dir():
    """Get current SPLASH ROOT dir."""
    # @Helgi: Harsh changed this on Jul 24, 2020
    return Naming.dir_root('cg')
    #return Naming.dir_alt_org(MUMMI_ROOT, 'sims')


def get_alt_sim_dir(root):
    """Get current dir."""
    # @Helgi: Harsh changed this on Jul 24, 2020
    return Naming.dir_sim('cg', root)
    #return Naming.dir_alt_org(root, 'sims')


def set_sim_dir():
    """Set current dir to simulation directory."""
    set_spec_sim_dir(get_sim_dir())


def set_spec_sim_dir(simdir):
    """Set current dir to provided simulation directory."""
    if not os.path.isdir(simdir):
        error = "-> Can't find simulation dir {}".format(simdir)
        LOGGER.error(error)
        raise ValueError(error)
    os.chdir(simdir)


def set_sub_dir(simdir):
    """For a simulaiton create and enter sub dir - works on relative paths"""
    if os.path.isdir(simdir):
        error = "---> Sim already has subdir {}".format(simdir)
        LOGGER.error(error)
        raise ValueError(error)
    utils.sys_call('mkdir -p ' + simdir, test=False)
    os.chdir(simdir)


def get_sim_sub_dir(pSim, subDirRef):
    """For a simulaiton get full sub dir name"""
    # @Helgi: this is where your code will break unfortunately.
    # you're passing subDirRef as "cg-eq" or "cg-prod" but neither exists now
    # Instead, we need to pass "cg" (please look at createsim.py and change as needed)
    return Naming.dir_sim(subDirRef, pSim.name)


def get_sim_eq_pdb_name(pSimName):
    """For a simulaiton get full sub dir name"""
    return Naming.eq_pdb(pSimName, "cg000")


def add_simpath(filename):
    """Add CURRENT path to filename"""
    return os.path.join(os.getcwd(), filename)
'''

# All the testing and test function that used to be below here have been removed but are in backup/particlesim.py in C1 git rebo
