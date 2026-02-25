"""Python module for cgSelections class"""

import mummi_ras
import mummi_ras.online.cg.cnlibs as cnlibs

from logging import getLogger
LOGGER = getLogger(__name__)

# ----------------------------------------------------------------------------
class CgSelections(object):
    """This class encapuslates a set of MDAnalysis selections we can use for processing CG simulaitons."""

    def __init__(self, cg):

        # @TODO this should probably also be in a yaml file?
        # Constants for analysis
        # cutoff for protein-protein contacts (in Angstroms)
        self.contactcutoff = 6.0
        # cutoff for protein-lipid contacts (in Angstroms)
        self.contactcutoffLong = 12.0
        self.numResPerProtein = 184
        self.numResPerProteinGTPMg = 186
        # self.numAtomPerProtein = 418
        # self.numAtomPerProteinGTPMg = 428
        self.numResGdomain = 165
        self.proteinIdentResidueString = "THR"
        self.numProteinIdentResiduePerProtein = 26

        # Can fail on first call
        try:
            self.pdbRef = \
                mummi_ras.online.cg.get_tiltrot_z_state_multiple \
                .kras_ref_universe
            # This is a MDAnalysis universe with a reference KRAS only loaded
            # x1 time
        except ImportError as e:
            print(e)
            print(
                "> ERROR: Splash installation not found. Please handle "
                "pdbRef explicitly!")

        # Define the number of proteins in the system, has been taken out of the getKRASstates function so only has to be done once
        self.kras_indices = mummi_ras.online.cg.get_tiltrot_z_state_multiple.get_protein_info(cg.syst,'resname ACE1 and name BB')

        # find the first protein residue and number of proteins
        self.firstProteinStartResidue, self.numProteins = \
            cnlibs.findFirstProteinResidue(
                cg.syst, self.proteinIdentResidueString,
                self.numProteinIdentResiduePerProtein)

        # groupContactProt, groupPlContactLipid = \
        #   cnlibs.getSel_findProtLipidContacts(
        #       cg.syst, numProteins, firstProteinStartResidue,
        #       numResPerProteinGTPMg)
        # groupContactProt, groupContactProtBB, groupPlContactLipid = \
        #   cnlibs.getSel_findProtLipidContacts(
        #        cg.syst, numProteins, firstProteinStartResidue,
        #        numResPerProteinGTPMg)
        self.groupContactProt, self.groupContactG, self.groupContactHVR, \
            self.groupContactProtBB, self.groupPlContactLipid = \
            cnlibs.getSel_findProtLipidContacts(
                cg.syst, self.numProteins, self.firstProteinStartResidue,
                self.numResPerProteinGTPMg, self.numResPerProtein,
                self.numResGdomain)
        # could have getSel_findProtLipidContacts() return
        # groupContactAllProt but Chris' other branch is ahead of develop so
        # will leave this as a separate routine for now
        # groupContactAllProt = cnlibs.getSel_allProteinFromAllRas(
        #    pdbPatch, numProteins, firstProteinStartResidue,
        #    numResPerProteinGTPMg)

        self.groupLipidHeadSelection, self.groupLipidTailASelection, \
            self.groupLipidTailBSelection = \
            cnlibs.getSel_findLipidLeaflet(cg.syst)

        self.groupGdomainBB = cnlibs.getSel_findGdomainBB(
            cg.syst, self.numProteins, self.firstProteinStartResidue,
            self.numResPerProteinGTPMg, self.numResGdomain)
        self.groupRefGdomainBB = cnlibs.getSel_findGdomainBB(
            self.pdbRef, 1, 1, self.numResPerProteinGTPMg, self.numResGdomain)

        # options are 'serial' and 'OpenMP'
        self.distanceBackendType = 'Serial'

        return
