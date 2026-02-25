import logging
import numpy as np
import MDAnalysis as mda
from collections import defaultdict

LOGGER = logging.getLogger(__name__)

# hardcoded anchor points 
RESIDXS_GDOMAIN = np.arange(0, 165, dtype=int)

# for the 4B 
RESIDXS_FARNESYL_4B = np.array([183])
RESIDXS_CRD_4B = np.array([276,277,278,279, 289,290,291,292])

RESIDXS_FARNESYL_4A = np.array([183]) + 1 
RESIDXS_CRD_4A = np.array([276,277,278,279, 289,290,291,292]) + 1

RAS_BRAF_GDOMAIN = (0, 380)
RAS_BRAF_CRD_4B = (649, 660, 685, 694)
RAS_BRAF_FARYNESYL = (430, 431)


alignment_anchors = {"RAS4B" : {"origin" : RESIDXS_FARNESYL_4B,
                                "xaxis" : RESIDXS_GDOMAIN,
                                "xyaxis" : np.array([0])},

                     "RAS4B_RAF" : {"origin" : RESIDXS_GDOMAIN,
                                    "xaxis" : RESIDXS_CRD_4B,
                                    "xyaxis" : np.array([0])},

                     "RAS4A" : {"origin" : RESIDXS_FARNESYL_4A,
                                "xaxis" : RESIDXS_GDOMAIN,
                                "xyaxis" : np.array([0])},

                     "RAS4A_RAF" : {"origin" : RESIDXS_GDOMAIN,
                                    "xaxis" : RESIDXS_CRD_4A,
                                    "xyaxis" : np.array([0])},

                     "RAS-BRAF" : {"origin" : f"index {RAS_BRAF_GDOMAIN[0]} to {RAS_BRAF_GDOMAIN[1]}",
                                   "xaxis" : f"index {RAS_BRAF_CRD_4B[0]} to {RAS_BRAF_CRD_4B[1]} or index {RAS_BRAF_CRD_4B[2]} to {RAS_BRAF_CRD_4B[3]}",
                                   "xyaxis" : None}
                    }

# ------------------------------------------------------------------------------
# Class to realign a mdanalysis trajectory
# ------------------------------------------------------------------------------
class CGFrameAligner():

    """
    There are two types of orientations we can ensure.

    'membrane' orientation:
        here, membrane orientation stays consistent (unrotated in z)
        i.e., we leave the protein tilted wrt the membrane

        we use two anchors and make two transformations
        1. translate to move anchor0 to the center of the box
        2. rotate within xy plane to ensure anchor1 is (+x,0,z) wrt anchor0

    'protein' orientation:
        here, protein orientation stays consistent
        i.e., we "lay down" the protein along +x axis

        we use two/three anchors to make these transformations
        1. translate to move anchor0 to (0,0,0)
        2. rotate the anchor1 to make it (+x,0,0) wrt anchor0
        3. if a third anchor is used, we move anchor3 to make sure it is
           in xy plane, i.e., (+x,+y,0)
    """

    def __init__(self, orient_type='protein', origin_at_anchor0=True):

        assert orient_type in ['protein', 'membrane']
        if orient_type == 'protein':
            assert origin_at_anchor0

        self.orient_type = orient_type
        self.origin_at_anchor0 = origin_at_anchor0

        self.mdu = None
        self.protein = None

        # these are the MDAnalysis atom selection strings
        self.anchor_residxs_origin = None
        self.anchor_residxs_xaxis = None
        self.anchor_redixs_xyaxis = None

        # these the MDAnalysi atomGroups for the above selections
        self.anchor_origin = None
        self.anchor_xaxis = None
        self.anchor_xyaxis = None

        self.angles = []
        self.xydisplacements = []


        # for membrane orientation, need another translation at the end
        # to make this work!
        assert self.origin_at_anchor0, 'Currently forcing origin at anchor0'

    # --------------------------------------------------------------------------
    def set_mdu(self, mdu):
        assert isinstance(mdu, mda.Universe)
        self.mdu = mdu

    def set_anchors(self, residxs_origin, residxs_xaxis, residxs_xyaxis=None):

        if self.orient_type == 'protein':
            if residxs_xyaxis is not None:
                assert self.orient_type == 'protein', \
                    "Third anchor is needed only for orient_type = 'protein'"

        LOGGER.debug(f'Setting anchors: ({residxs_origin}), '
                     f'({residxs_xaxis}), ({residxs_xyaxis})')
        self.anchor_residxs_origin = residxs_origin
        self.anchor_residxs_xaxis = residxs_xaxis
        self.anchor_residxs_xyaxis = residxs_xyaxis


    def get_angles(self):
        return np.array(self.angles)

    def get_xy_displacements(self):
        return np.array(self.xydisplacements)

    # --------------------------------------------------------------------------
    def init(self):
        #self.protein = self.mdu.select_atoms('protein')
        self.protein = self.mdu.select_atoms('name BB')
            
        def get_residues(residxs):
            resids = ' '.join([str(self.protein.resids[0]+i) for i in residxs])
            return self.protein.select_atoms(f'resid {resids}')

        # convert from a selection string to an atomGroup
        def get_selection(selection):
            return self.protein.select_atoms(selection)

        #self.anchor_origin = get_residues(self.anchor_residxs_origin)
        #self.anchor_xaxis = get_residues(self.anchor_residxs_xaxis)
        self.anchor_origin = get_selection(self.anchor_residxs_origin)
        self.anchor_xaxis = get_selection(self.anchor_residxs_xaxis)

        LOGGER.info(f'origin: {self.anchor_origin.resnames.shape}: {list(self.anchor_origin.resnames)}')
        LOGGER.info(f'xaxis:  {self.anchor_xaxis.resnames.shape}: {list(self.anchor_xaxis.resnames)}')

        if self.anchor_residxs_xyaxis is not None:
            #self.anchor_xyaxis = get_residues(self.anchor_residxs_xyaxis)
            self.anchor_xyaxis = get_selection(self.anchor_residxs_xyaxis)
            LOGGER.info(f'xyaxis: {self.anchor_xyaxis.resnames.shape}: {list(self.anchor_xyaxis.resnames)}')

    # --------------------------------------------------------------------------
    # anchors for this alignment
    def get_anchor_residxs(self):
        return self.anchor_residxs_origin, \
               self.anchor_residxs_xaxis, \
               self.anchor_residxs_xyaxis

    def get_anchors(self):
        xy = None
        if self.anchor_xyaxis is not None:
            xy = self.anchor_xyaxis.center_of_mass(unwrap=True)

        return self.anchor_origin.center_of_mass(unwrap=True), \
               self.anchor_xaxis.center_of_mass(unwrap=True), \
               xy

    # --------------------------------------------------------------------------
    # fix the periodic boundary
    # --------------------------------------------------------------------------
    def fix_periodic_bfs(self, node_pos,  bbox_min, bbox_max):

        print (f'   bbox:  min={bbox_min}; max={bbox_max}')

        bbox_x = bbox_max[0]-bbox_min[0]
        bbox_y = bbox_max[1]-bbox_min[1]
        bbox_z = bbox_max[2]-bbox_min[2]

        seed_idxs = np.array([0]) #np.array([183]) # RESIDXS_FARNESYL self.anchor_residxs_origin

        # use all pairs instead edges to capture the not "connected" beads
        protein_ids = np.arange(0, node_pos.shape[0], dtype=int)

        all_pairs = []
        for e in protein_ids:
            for f in protein_ids:
                if e < f:
                    all_pairs.append([e,f])

        # create a set of outgoing edges from each node
        outgoing = defaultdict(list)
        for (a,b) in all_pairs:
            outgoing[a].append(b)
            outgoing[b].append(a)

        print (f'Fixing periodic-xyz artifacts using BFS (seed = {seed_idxs}, bbox_x = {bbox_x:.3f}, bbox_y = {bbox_y:.3f}, bbox_z = {bbox_z:.3f})')

        # perform a bfs on the graph!
        visited = [s for s in seed_idxs]
        queue = [s for s in seed_idxs]

        while queue:
            s = queue.pop(0)
            for nbr in outgoing[s]:

                # check and fix periodic boundary artifact
                d = node_pos[nbr] - node_pos[s]
                # in x
                if d[0] > 0.5*bbox_x:
                    node_pos[nbr, 0] -= bbox_x
                elif d[0] < -0.5*bbox_x:
                    node_pos[nbr, 0] += bbox_x

                # in y
                if d[1] > 0.5*bbox_y:
                    node_pos[nbr, 1] -= bbox_y
                elif d[1] < -0.5*bbox_y:
                    node_pos[nbr, 1] += bbox_y

                # in z
                if d[2] > 0.5*bbox_z:
                    node_pos[nbr, 2] -= bbox_z
                elif d[2] < -0.5*bbox_z:
                    node_pos[nbr, 2] += bbox_z

                if nbr not in visited:
                    visited.append(nbr)
                    queue.append(nbr)
        return node_pos
    # --------------------------------------------------------------------------
    def align(self):

        xaxis = np.array([1.,0.,0.])
        zaxis = np.array([0.,0.,1.])

        # ----------------------------------------------------------------------
        # normalize a vector
        do_norm = lambda x: x / np.linalg.norm(x)

        # get periodic displacement vector
        def get_displacement_2d(a,b):
            d = b-a
            for i in range(2):
               if d[i] > 0.5*self.mdu.dimensions[i]:
                  d[i] -= self.mdu.dimensions[i]
               if d[i] < -0.5*self.mdu.dimensions[i]:
                  d[i] += self.mdu.dimensions[i]
            return d


        def minmax_dist():
            import scipy
            x = self.protein.positions
            pd = scipy.spatial.distance.pdist(x)
            return pd.min(), pd.max()

        # start with an unwrapped protein
        # self.protein.unwrap()
        # ----------------------------------------------------------------------
        # protein alignment
        # ----------------------------------------------------------------------
        if self.orient_type == 'protein':

            # ------------------------------------------------------------------
            # step 1. bring anchor 1 to the center
            x0, x1, x2 = self.get_anchors()
            dx = -x0

            self.mdu.atoms.translate(dx)

            # ------------------------------------------------------------------
            # step 2. bring anchor 2 to +x axis
            x0, x1, x2 = self.get_anchors()

            # find the angle between displacement vector (current) and x axis (target)
            # and rotate along the normal vector
            dx = do_norm(x1-x0)
            norm = do_norm(np.cross(dx, xaxis))
            angle = np.rad2deg(np.arccos(np.dot(dx, xaxis)))

            self.mdu.atoms.rotateby(angle, norm, point=x0)

            # ------------------------------------------------------------------
            # step 3. bring anchor3 to +y axis and in yz plane
            # i.e., z-component should be zero (wrt anchor2)
            if self.anchor_xyaxis is not None:

                x0, x1, x2 = self.get_anchors()

                # operate in yz plane
                dx = x2-x0

                # angle in (y,z) plane
                angle = np.rad2deg(np.arctan2(dx[2], dx[1]))

                self.mdu.atoms.rotateby(-angle, xaxis, point=x0)

        # ----------------------------------------------------------------------
        # membrane alignment
        # ----------------------------------------------------------------------
        else:
            # ------------------------------------------------------------------
            # step 1. bring anchor 1 to the center
            x0, x1, _ = self.get_anchors()
            ref_crd = self.mdu.dimensions[0:3] / 2
            ref_crd[2] += 22.5
            dx = ref_crd - x0
            #dx[2] = 0    # we dont want to move z
            self.mdu.atoms.translate(dx)
            self.mdu.atoms.wrap()

            # add addition z translation to correct for PDC
            #x0, x1, _ = self.get_anchors()
            #dz = -x0
            #dz[0] = 0
            #dz[1] = 0
            #self.mdu.atoms.translate(dz)
            #self.mdu.atoms.wrap()

            # ------------------------------------------------------------------
            # step 2. bring anchor 2 to +x axis
            x0, x1, _ = self.get_anchors()
            dvec = get_displacement_2d(x0, x1)
            angle = np.rad2deg(np.arctan2(dvec[1], dvec[0]))
            self.angles.append(angle)

            self.mdu.atoms.rotateby(-angle, zaxis, point=x0)
            #self.mdu.atoms.wrap()

            # now move 
            x0, x1, _ = self.get_anchors()
            dx = -x0
            self.mdu.atoms.translate(dx)

            # step 3. bring  m2 Gdomain to +x axis
            x0, x1, _ = self.get_anchors()
            a = x0
            b = self.mdu.atoms.select_atoms("index 2897 to 3277").center_of_geometry() # G-domain m2
 
            # get displacement vector
            d = b-a
            # d = get_displacement_2d(a, b)
            zaxis = np.array([0.,0.,1.])
            # Calculate the angle to rotate d to the positive x-axis
            angle_to_x = np.rad2deg(np.arctan2(d[1], d[0]))  # Current angle of d
            # Calculate the angle needed to rotate d to align with the positive x-axis
            rotation_angle = -angle_to_x  # We want to rotate in the opposite direction
            self.mdu.atoms.rotateby(rotation_angle, zaxis, point=a)

            dx[2] = -self.mdu.atoms.positions[430, 2] # translate z to zero farnesyl

        # ----------------------------------------------------------------------
        # now validate!
        # ----------------------------------------------------------------------
        atol = 1e-04
        x0, x1, x2 = self.get_anchors()

        # x1.x > x0.x and x1.y == x0.y
        # assert x1[0] > x0[0] and np.abs(x1[1]-x0[1]) < atol, \
        #        f'Cannot force x-axis: x0={x0}; x1={x1}'
                # :: {x1[0]>x0[0]} {np.isclose(x1[1], x0[1], atol)} '

        # x0 == 0
        if self.origin_at_anchor0:
            assert np.abs(x0[0]) < atol and np.abs(x0[1]) < atol, \
                   f'Cannot force origin: {x0}'

        if self.orient_type == 'protein':
            # x2.x > 0, x2.y > 0, x2.z = 0
            if self.anchor_xyaxis is not None:
                assert np.abs(x2[2]) < atol, f'Cannot force anchor2: {x2}'

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
