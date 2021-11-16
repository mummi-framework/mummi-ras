"""Patch Creator module."""

import os
import numpy as np
from logging import getLogger

from mummi_core.utils import Naming
from mummi_core.utils.timer import Timer
from mummi_ras.datastructures.patch import PatchConfig, Patch
from mummi_ras.datastructures.patch import Patch as pfPatch
from mummi_ras.datastructures.names_of_lipids_and_proteins import MACRO_RAS_ALL_NAMES
import mummi_ras
from mummi_ras.transformations.patch_creator.macroParser import MacroParser, show_macro_data

LOGGER = getLogger(__name__)


DO_DISABLE_MU18 = False
KNOWN_MUS = [-24,0,12,18]


# ------------------------------------------------------------------------------
# This module creates patches
# ------------------------------------------------------------------------------
class MacroPatchCreator (object):
    """Macro Patch Creator."""

    @staticmethod
    def get_ras_idxs(protein_states):
        return np.where(np.isin(protein_states, MACRO_RAS_ALL_NAMES))[0]

    # --------------------------------------------------------------------------
    # create patches from a macro snapshot *at native resolution*
    # --------------------------------------------------------------------------
    @classmethod
    def _initialize(cls, macro_gsz, macro_psz, patch_psz):
        ''' Intitialize the class to extract patches.
                macro_gsz:  size of macro domain in grid points
                macro_psz:  size of macro domain in physical units
                patch_psz:  size of patches in (same) physical units
        '''
        assert isinstance(macro_gsz, int)
        assert isinstance(macro_psz, float)
        assert isinstance(patch_psz, float)

        cls.macro_gsz = np.array([macro_gsz, macro_gsz])
        cls.macro_psz = np.array([macro_psz, macro_psz])
        cls.patch_psz = np.array([patch_psz, patch_psz])

        cls.dx = cls.macro_psz[0] / cls.macro_gsz[0]
        cls.patch_gsz = (cls.patch_psz / cls.dx).astype(int) + 1

        assert cls.patch_gsz[0] == cls.patch_gsz[1]
        assert cls.patch_gsz[0] % 2 == 1

        cls.pixel_area = cls.dx * cls.dx

        npixels = cls.patch_gsz[0] - 1
        total_area = cls.pixel_area * npixels * npixels
        patch_area = cls.patch_psz[0] * cls.patch_psz[0]
        assert np.isclose(total_area, patch_area)

        # half radius of macro and patch
        cls.half_macro_psz = 0.5 * cls.macro_psz[0]
        cls.half_patch_psz = 0.5 * cls.patch_psz[0]
        cls.half_patch_gsz = int(cls.patch_gsz[0] / 2)

        s = 'Initialized PatchCreator: macro = p:{} g:{}, patch = p:{} g:{}'
        LOGGER.info(s.format(cls.macro_psz, cls.macro_gsz,
                             cls.patch_psz, cls.patch_gsz))

    # --------------------------------------------------------------------------
    @classmethod
    def _extract_patches(cls, data, all_positions, patch_positions, npixels=None):

        LOGGER.debug('Extracting patches: data = {}, patch_positions = {}, all_positions = {}'.format(data.shape, patch_positions.shape, all_positions.shape))

        assert isinstance(data, np.ndarray)
        assert isinstance(patch_positions, np.ndarray)
        assert isinstance(all_positions, np.ndarray)

        # already assumed that this macro psz is a square
        assert patch_positions.min() >= 0.0 and patch_positions.max() <= cls.macro_psz[0]
        assert all_positions.min() >= 0.0 and all_positions.max() <= cls.macro_psz[0]

        if 1 == len(patch_positions.shape):
            patch_positions = patch_positions[np.newaxis, :]

        if 1 == len(all_positions.shape):
            all_positions = all_positions[np.newaxis, :]

        assert 3 == len(data.shape)
        assert 2 == patch_positions.shape[1]
        assert 2 == all_positions.shape[1]

        pshape = (cls.patch_gsz[0], cls.patch_gsz[1], data.shape[-1])
        if npixels is None:
            npixels = cls.patch_gsz[0]

        # output data: for each patch
        pdata = []          # list of data points
        pneighbors = []     # list of neighbor ids
        pcenters = []       # list of the physical coords of central pixel
        ppositions = []     # list of the patch_positions wrt the central pixel

        #show_macro_data(data, all_positions[0:1], 1000, False)

        # ----------------------------------------------------------------------
        i = 0
        for p in patch_positions:

            #LOGGER.debug('doing {}: {}'.format(i, p))
            i += 1

            # get the grid location of this particle
            gpos = np.array(p / cls.dx)
            gpix = np.round(gpos).astype(int)
            #print('\np =', p, 'gpos =', gpos, 'gpix =', gpix)

            # the bounds of the patch (both inclusive)
            gp0 = gpix - cls.half_patch_gsz
            gp1 = gpix + cls.half_patch_gsz
            #print('  gbounds =', gp0, gp1)

            pp0 = np.array(gp0 * cls.dx)
            pp1 = np.array(gp1 * cls.dx)
            #print('  pbounds =', pp0, pp1)

            # make sure we got the correct sizes
            gdiff = gp1-gp0
            pdiff = pp1-pp0
            #print(' gdiff, pdff =', gdiff, pdiff)
            #print(' patch_sz =', cls.patch_gsz, cls.patch_psz)

            assert np.allclose(gdiff, cls.patch_gsz - [1,1])
            assert np.allclose(pdiff, cls.patch_psz)

            # make sure the center RAS is within half pixel from the center
            pcent = pp0+cls.half_patch_psz
            p_wrt_cent = np.abs(p-pcent)
            #print(' pcent =', pcent)
            #print(' p_wrt_cent =',p_wrt_cent)

            # by checking for 0.51*dx instead of 0.5*dx accounts
            # for precision errors
            assert p_wrt_cent[0] <= 0.51*cls.dx and p_wrt_cent[1] <= 0.51*cls.dx

            # ------------------------------------------------------------------
            # extract the ids to be copied
            ix = np.arange(gp0[0], gp1[0]+1)
            iy = np.arange(gp0[1], gp1[1]+1)

            assert ix.shape[0] == cls.patch_gsz[0]
            assert iy.shape[0] == cls.patch_gsz[1]

            # account for periodic boundary
            ix[ix < 0] += cls.macro_gsz[0]
            iy[iy < 0] += cls.macro_gsz[0]
            ix[ix >= cls.macro_gsz[1]] -= cls.macro_gsz[1]
            iy[iy >= cls.macro_gsz[1]] -= cls.macro_gsz[1]

            # ------------------------------------------------------------------
            # now extract the concentrations
            patch = np.zeros(pshape, dtype=data.dtype)
            for y in range(pshape[1]):
                for x in range(pshape[0]):
                    patch[y,x,:] = data[iy[y], ix[x], :]

            # subsample!
            if cls.patch_gsz[0] > npixels:
                patch = pfPatch._subsample_mean(patch, cls.patch_psz[0], npixels)

            # ------------------------------------------------------------------
            # find the other positions near this patch's center
            pos_wrt_cent = all_positions - pcent

            # handle periodic domain
            pos_wrt_cent[pos_wrt_cent[:,0] < 0, 0] += cls.macro_psz[0]
            pos_wrt_cent[pos_wrt_cent[:,1] < 0, 1]  += cls.macro_psz[1]
            pos_wrt_cent[pos_wrt_cent[:,0] > cls.half_macro_psz, 0] -= cls.macro_psz[0]
            pos_wrt_cent[pos_wrt_cent[:,1] > cls.half_macro_psz, 1] -= cls.macro_psz[1]

            # get the neighbors
            pnbrs = np.abs(pos_wrt_cent) < cls.half_patch_psz
            pnbrs = np.where(2 == pnbrs.sum(axis=1))[0]
            nnbrs = pnbrs.shape[0]

            if nnbrs > 0:
                # these are the patch_positions of neighbors
                pos_wrt_cent = pos_wrt_cent[pnbrs]

                # rearrange the neighbors with respect to dist from center_pixel
                sidx = np.argsort(np.linalg.norm(pos_wrt_cent, axis=1))
                pnbrs = pnbrs[sidx]
                pos_wrt_cent = pos_wrt_cent[sidx]

                # make sure that all relative neighbor patch_positions
                # are within the half patch radius
                assert np.abs(pos_wrt_cent).max() <= cls.half_patch_psz

                # repeat the test that the first RAS is indeed near the center pixel
                assert pos_wrt_cent[0,0] <= 0.5*cls.dx and pos_wrt_cent[0,1] <= 0.5*cls.dx

            else:
                pos_wrt_cent = []

            # ------------------------------------------------------------------
            # add to the list
            pdata.append(patch)
            pcenters.append(pcent)
            pneighbors.append(pnbrs)
            ppositions.append(pos_wrt_cent)

        # ----------------------------------------------------------------------
        LOGGER.debug('Extracted patches: {}'.format(len(pdata)))

        # ----------------------------------------------------------------------
        return pdata, pcenters, pneighbors, ppositions

    # --------------------------------------------------------------------------
    @classmethod
    def _get_random_position_far_from(cls, avoid_positions, min_distance):

        # do 100 trials to find a position
        # if not, simply return the farthest
        ntrials = 100

        farthest_pos, farthest_dist = None, 0.
        for i in range(ntrials):

            pos = np.random.random((1,2)) * cls.macro_psz
            diff = np.abs(avoid_positions - pos)    # absolute distance to avoid_positions

            fx = diff[:,0] >= 0.5*cls.macro_psz[0]
            fy = diff[:,1] >= 0.5*cls.macro_psz[1]

            diff[fx,0] = cls.macro_psz[0] - diff[fx,0]
            diff[fy,1] = cls.macro_psz[1] - diff[fy,1]

            diff = np.linalg.norm(diff, axis=1)     # compute the distance
            dist = diff.min()

            if dist > min_distance:
                return pos

            if farthest_pos is None:
                farthest_pos, farthest_dist = pos, dist

            elif dist > farthest_dist:
                farthest_pos, farthest_dist = pos, dist

        s = f'Could not find a position at least {min_distance} away ' \
            f'from {avoid_positions.shape[0]} points. ' \
            f'Returning one that is {farthest_dist} away!'

        LOGGER.warning(s)
        return farthest_pos

    @classmethod
    def _create_random_positions(cls, n, avoid_positions, min_distance):

        assert isinstance(avoid_positions, np.ndarray)
        LOGGER.debug('creating {} random positions, avoiding {} at min_dist {}'.format(n, avoid_positions.shape, min_distance))

        # create random positions, sort them based on nearest ras distance,
        # and then pick the farthest

        avoid_positions2 = np.copy(avoid_positions)
        positions = np.zeros((n,2), dtype=avoid_positions.dtype)

        for i in range(n):
            new_pos = cls._get_random_position_far_from(avoid_positions2, min_distance)

            positions[i,:] = new_pos
            avoid_positions2 = np.concatenate((avoid_positions2, new_pos))

        return np.array(positions)

    # --------------------------------------------------------------------------
    @staticmethod
    def _create_patches(patch_config, lipid_concentrations,
                        protein_ids, protein_states, protein_positions,
                        complex_ids, pwidth, idcounter, npatches_random=0):

        assert isinstance(patch_config, PatchConfig)
        assert isinstance(lipid_concentrations, np.ndarray)
        assert isinstance(protein_ids, np.ndarray)
        assert isinstance(protein_states, np.ndarray)
        assert isinstance(protein_positions, np.ndarray)
        assert isinstance(complex_ids, np.ndarray)
        assert isinstance(pwidth, float)
        assert isinstance(npatches_random, int)
        assert npatches_random >= 0
        assert len(protein_positions.shape) == 2
        assert protein_positions.shape[1] == 2

        # ----------------------------------------------------------------------
        patches = []

        # ----------------------------------------------------------------------
        # first, create ras-cented patches
        # ----------------------------------------------------------------------
        if True:
            # the indices of the ras particles "in the list"
            ras_idxs = MacroPatchCreator.get_ras_idxs(protein_states)
            nras = ras_idxs.shape[0]

            # the indices of the ras particles "given by macro data"
            ras_ids = protein_ids[ras_idxs]
            assert 1 == ras_ids.min() and nras == ras_ids.max()

            ras_positions = protein_positions[ras_idxs]
            npatches = ras_positions.shape[0]

            pconcs, pcenters, pneighbor_idxs, pneighbor_positions = \
                MacroPatchCreator._extract_patches(lipid_concentrations,
                                                   protein_positions,
                                                   ras_positions)

            # now create the patch objects
            for pidx in range(npatches):

                pnbrs = pneighbor_idxs[pidx]            # indices in the list
                p = Patch(Naming.pfpatch(idcounter+pidx),
                          patch_config, pconcs[pidx],
                          pcenters[pidx],
                          protein_ids[pnbrs],
                          protein_states[pnbrs],
                          pneighbor_positions[pidx])

                patches.append(p)
            LOGGER.warning('Need to add back code to compute complexes')
            idcounter += npatches

        # ----------------------------------------------------------------------
        # now, create lipid only patches
        # ----------------------------------------------------------------------
        if npatches_random > 0:
            # a patch is 30 nm wide.
            # the farthest point from the patch center = 15 * sqrt(2) = 21.21
            # we enforce (the center of) a lipid only patch to be at least 50 nm
            # away from any protein
            random_positions = \
                MacroPatchCreator._create_random_positions(npatches_random,
                                                           protein_positions[:,:2],
                                                           50.)

            pconcs, pcenters, pneighbor_idxs, pneighbor_positions = \
                MacroPatchCreator._extract_patches(lipid_concentrations,
                                                   protein_positions,
                                                   random_positions)

            # now create the patch objects
            for pidx in range(npatches_random):

                pnbrs = pneighbor_idxs[pidx]
                if len(pnbrs) > 0:
                    LOGGER.warning('Discard a lipid-only patch due to presence of a protein')
                    continue

                p = Patch(Naming.pfpatch(idcounter+pidx),
                          patch_config, pconcs[pidx], pcenters[pidx])
                patches.append(p)

            idcounter += npatches

        # ----------------------------------------------------------------------
        return patches

    # --------------------------------------------------------------------------
    @staticmethod
    def create_gc_patches(inpath, sim_name, frame_ids, iointerface, fail_fast=False):

        assert isinstance(inpath, str)
        assert isinstance(sim_name, str)
        assert isinstance(frame_ids, list)

        inpath = os.path.join(inpath, sim_name)
        LOGGER.info(f'Extracting {len(frame_ids)} patches from ({inpath})')
        patches = {}

        if not iointerface.namespace_exists(inpath):
            raise AttributeError(f'Namespace ({inpath}) does not exist')
            return {}

        for fidx in frame_ids:

            snapshot_name, filenames = MacroParser.get_filenames(fidx)
            tag = f'({inpath})/({snapshot_name})'

            try:
                st = Timer()
                data = iointerface.load_files(inpath, filenames)
                if data is None:
                    LOGGER.debug(f'Did not find {tag}')
                    if fail_fast:
                        break
                    else:
                        continue

                assert isinstance(data, list)
                assert len(data) == len(filenames)
                LOGGER.debug(f'Read {tag} {st}')

                # parse these bytestreams into the actual data
                st = Timer()
                simtime, tunit, lunit, concs, prot_ids, prot_states, prot_pos, complexes = \
                    MacroParser.parse_snapshot_from_streams(data[0], data[1], data[2], sim_name=sim_name)

                LOGGER.debug(f'Parsed {tag} {st}')

                patch_id = Naming.pfpatch(fidx, sim_name)
                config = PatchConfig(sim_name, simtime, tunit, lunit,
                                 {'snapshot': snapshot_name})

                patch_cent = None
                p = Patch(patch_id, config, concs,
                          patch_cent, prot_ids, prot_states, prot_pos,
                          complexes)
                LOGGER.debug(f'Created ({patch_id})')

                # center around first RAS
                if prot_states is not None:
                    assert prot_states.shape[0] > 0, f'Found no proteins'
                    ras_idxs = MacroPatchCreator.get_ras_idxs(prot_states)
                    assert ras_idxs.shape[0] > 0, f'Found no RAS (states = {prot_states})'
                    first_ras_pos = np.copy(prot_pos[ras_idxs[0]][:2])

                    #show_patch(p.concentrations, p.protein_positions)
                    p.wrap_to_center(first_ras_pos)
                    #show_patch(p.concentrations, p.protein_positions)

                patches[patch_id] = p

            except Exception as e:
                LOGGER.warning(f'Failed to read and parse {tag}. moving on! error=({e})')
                break

        return patches

    # --------------------------------------------------------------------------
    # now, enable using this class as an object with state
    # --------------------------------------------------------------------------
    def __init__(self, config, in_interface, out_interface):

        self.config = config
        self.in_interface = in_interface
        self.out_interface = out_interface

        self.inpath = config['params'].get('inpath', Naming.dir_root('macro'))
        self.outpath = config['params'].get('outpath', Naming.dir_root('patches'))
        self.chkpt = os.path.join(config['params'].get('chkpath', self.outpath), 'patch_creator.chk')

        self.max_frames_per_read = config['params']['max_frames_per_read']
        self.sim2frame = config['params']['sim2frame']

        self.mgsize = config['params']['macro_gsz']
        self.mpsize = config['params']['macro_psz']
        self.ppsize = config['params']['patch_psz']
        self.is_gc = bool(config['params']['is_gc'])

        self.npatches = 0
        self.nlipid_only = config['params'].get('nlipid_only', 0)

        # macro simulations
        if not self.is_gc:
            assert len(self.sim2frame) == 1
            MacroPatchCreator._initialize(self.mgsize, self.mpsize, self.ppsize)

    # --------------------------------------------------------------------------
    # run the next iteration of patch selector
    def run_for_macro(self, max_frames_per_read=1000):

        max_frames_per_read = min(max_frames_per_read, self.max_frames_per_read)
        patches = []

        for sim_name, frame_idx in self.sim2frame.items():
            for fidx in range(frame_idx, frame_idx+max_frames_per_read):

                snapshot_name, filenames = MacroParser.get_filenames_c3(fidx)
                tag = f'{sim_name}/{snapshot_name}'

                LOGGER.info(f'Looking for ({tag}) in ({self.inpath})')

                st = Timer()
                data = self.in_interface.load_files(self.inpath, filenames)
                if data is None:
                    break

                assert isinstance(data, list)
                assert len(data) == len(filenames)
                LOGGER.info(f'Read ({tag}) {st}')

                # parse these bytestreams into the actual data
                st = Timer()
                simtime, tunit, lunit, concs, prot_ids, prot_states, prot_pos, complexes = \
                    MacroParser.parse_snapshot_from_streams(data[0], data[1], data[2])

                LOGGER.info(f'Parsed ({snapshot_name}) {st}')

                assert len(concs.shape) == 3
                assert concs.shape[0] == MacroPatchCreator.macro_gsz[0]
                assert concs.shape[1] == MacroPatchCreator.macro_gsz[0]

                config = PatchConfig(sim_name, simtime, tunit, lunit,
                                     {'snapshot': snapshot_name})

                p = MacroPatchCreator._create_patches(config, concs,
                                                      prot_ids, prot_states,
                                                      prot_pos[:,:2], complexes,
                                                      self.ppsize,
                                                      self.npatches,
                                                      self.nlipid_only)

                LOGGER.info(f'Created {len(p)} patches ({tag}) {st}')
                patches.extend(p)
                #p[0].render()

                self.sim2frame[sim_name] = fidx+1
                self.npatches += len(p)

        return patches

    # --------------------------------------------------------------------------
    def run_for_gc(self, max_frames_per_read=1000):

        max_frames_per_read = min(max_frames_per_read, self.max_frames_per_read)
        patches = []

        # go over each simulation, and read a few frames
        for sim_name, frame_idx in self.sim2frame.items():

            if DO_DISABLE_MU18:
                mu = Naming.get_mu_from_simname(sim_name)
                assert mu in KNOWN_MUS, f'Found invalid mu = {mu}'
                if mu == 18:
                    continue

            frame_ids = [frame_idx+i for i in range(max_frames_per_read)]
            sim_patches = \
                MacroPatchCreator.create_gc_patches(self.inpath, sim_name,
                                                    frame_ids, self.in_interface,
                                                    fail_fast=True)

            patches.extend(sim_patches.values())
            n = len(sim_patches)
            self.sim2frame[sim_name] += n
            self.npatches += n

        LOGGER.info(f'Created {len(patches)} patches')
        return patches

    # --------------------------------------------------------------------------
    def save_patches(self, _):
        self.out_interface.save_patches(self.outpath, _)

    # --------------------------------------------------------------------------
    def checkpoint(self):

        chkpt = {
            'npatches': self.npatches,
            'is_gc': self.is_gc,
            'inpath': self.inpath,
            'outpath': self.outpath,
            'mgsize': self.mgsize,
            'mpsize': self.mpsize,
            'ppsize': self.ppsize,
            'nlipid_only': self.nlipid_only,
        }
        if self.is_gc:
            chkpt['sim2frame'] = self.sim2frame
        else:
            chkpt['frame'] = self.frame

        self.out_interface.save_checkpoint(self.chkpt, chkpt)

        LOGGER.info( f"Checkpoint'd MacroPatchCreator:  npatches = {self.npatches}")

    # --------------------------------------------------------------------------
    def restore(self):
        chkpt = self.in_interface.load_checkpoint(self.chkpt)
        LOGGER.debug(f'read checkpoint {chkpt}')
        if chkpt:
            assert chkpt['is_gc'] == self.is_gc, "RESTORE PATCH: is_gc doesn't match"
            assert chkpt['inpath'] == self.inpath, "RESTORE PATCH: inpath doesn't match"
            assert chkpt['outpath'] == self.outpath, "RESTORE PATCH: outhpath doesn't match"
            assert chkpt['mgsize'] == self.mgsize, "RESTORE PATCH: mgsize doesn't match"
            assert chkpt['mpsize'] == self.mpsize, "RESTORE PATCH: mpsize doesn't match"
            assert chkpt['ppsize'] == self.ppsize, "RESTORE PATCH: ppsize doesn't match"
            assert chkpt['nlipid_only'] == self.nlipid_only, "RESTORE PATCH: nlipid_only doesn't match"

            self.npatches = chkpt['npatches']
            if self.is_gc:
                sim2frame = chkpt['sim2frame']
                k1 = np.array(list(sim2frame.keys()))
                k2 = np.array(list(self.sim2frame.keys()))
                assert np.array_equal(np.sort(k1), np.sort(k2)), 'Inconsistent simulation list'

                # need to take the one from checkpoint (it has the correct values)
                self.sim2frame = sim2frame.copy()

            else:
                self.frame = chkpt['frame']

            LOGGER.info(f"Restored MacroPatchCreator:  npatches = {self.npatches}: {self.sim2frame.keys()}")
        else:
            LOGGER.info(f"No {self.chkpt} file found")

    # --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# create an abstraction to grab the patches out of gc simulations
# --------------------------------------------------------------------------
def fetch_gc_patches(inpath, patch_ids):

    iointerface = mummi_ras.get_io('mummi')

    # --------------------------------------------------------------------------
    # read a list of patches and return a dictionary
    if isinstance(patch_ids, list):

        # split the requested ids based on sim names
        sim2frames = {}
        for p in patch_ids:
            sim, frame = Naming.pfpatch_parse(p)
            if sim not in sim2frames:
                sim2frames[sim] = [frame]
            else:
                sim2frames[sim].append(frame)

        # now, go over these frames
        patches = {}
        for sim, frames in sim2frames.items():

            sim_patches = MacroPatchCreator.create_gc_patches(inpath, sim, frames, iointerface, fail_fast=False)
            for patch_id, patch in sim_patches.items():
                patches[patch_id] = patch

        LOGGER.info(f'Read {len(patches)} patches (out of {len(patch_ids)} requested)')
        return patches

    # --------------------------------------------------------------------------
    # read a single patch
    else:
        sim, frame = Naming.pfpatch_parse(patch_ids)
        sim_patches = MacroPatchCreator.create_gc_patches(inpath, sim, [frame], iointerface, fail_fast=False)
        return None if len(sim_patches) == 0 else sim_patches[patch_ids]

    # --------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if __name__ == '__main__':

    mummi_ras.init_logger()

    # --------------------------------------------------------------------------
    inpath = './data'
    inpath = '/p/gpfs1/splash/campaign4/bootstrap-macro/good-sign'
    sim_name = 'mu-10ras0raf0'
    sim_name = 'mu-24-0ras0raf-instance1'
    patch_ids = [Naming.pfpatch(i, sim_name) for i in range(1,11)]
    # patches = fetch_gc_patches(inpath, patch_ids)
    # print (patches)

    patch = fetch_gc_patches(inpath, patch_ids)

    exit()
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    if False:   # test reading a c3 patch
        inpath = '/Users/bhatia4/work/data/pilot2/3/mummi_c3_20201223/patches/tomaso'
        p = Patch.load_npz(os.path.join(inpath, 'pfpatch_000000000027.npz'))
        print (p)
        exit()

    # --------------------------------------------------------------------------
    if False:   # test reading c4 patches for both old and new data
        inpath = '/Users/bhatia4/work/data/pilot2/campaign4/macro_from_c3/patches'
        p = Patch.load_npz(os.path.join(inpath, 'pfpatch_000000000330.npz'))
        print(p)

        inpath = '/Users/bhatia4/work/data/pilot2/campaign4/macro-20210413/patches'
        p = Patch.load_npz(os.path.join(inpath, 'gc-sim2_00000001.npz'))
        print(p)
        exit()

    # --------------------------------------------------------------------------
    if 0:   # test fetching a gc patch
        inpath = '/Users/bhatia4/work/data/pilot2/campaign4/macro-20210413'
        sim_name = 'gc-sim'
        patch_ids = [Naming.pfpatch(i, sim_name) for i in range(2)]
        #patches = fetch_gc_patches(inpath, patch_ids)
        #print (patches)

        patch = fetch_gc_patches(inpath, patch_ids[1])
        patch.render(False)

        #patch.concentrations = patch.subsample_intg()
        #patch.render(False)

        import matplotlib.pyplot as plt
        plt.show()
        exit()

    # --------------------------------------------------------------------------
    # test creating and saving/rendering c4 patches
    # --------------------------------------------------------------------------
    DO_GC = True

    config = {}
    if DO_GC:
        in_int = mummi_ras.get_io('mummi')
        out_int = mummi_ras.get_io('simple')
        config['inpath'] = '/Users/bhatia4/work/data/pilot2/campaign4/macro-20210413'

        config['sim2frame'] = {'gc-sim': 1} #, 'gc-sim2': 1}

        config['macro_gsz'] = 72
        config['macro_psz'] = 30.
        config['nlipid_only'] = 0

    else:
        in_int = mummi_ras.get_io('simple')
        out_int = mummi_ras.get_io('simple')
        config['inpath'] = '/Users/bhatia4/work/data/pilot2/campaign4/macro_from_c3'

        config['sim2frame'] = {"macro_c3": 1000}

        config['macro_gsz'] = 2400
        config['macro_psz'] = 1000.
        config['nlipid_only'] = 33

    config['patch_psz'] = 30.
    config['max_frames_per_read'] = 2000
    config['outpath'] = os.path.join(config['inpath'], 'patches')
    config['chkpath'] = config['inpath']

    pc = MacroPatchCreator(config, in_int, out_int)
    if DO_GC:
        patches = pc.run_for_gc()
    else:
        patches = pc.run_for_macro()

    if True:
        print (f'---- rendering {len(patches)} patches')
        import matplotlib.pyplot as plt
        for p in patches:

            filename = os.path.join(config['outpath'], f'{p.id}.png')
            p.render(False)
            plt.savefig(filename, bbox_inches='tight')
            plt.close()

    # --------------------------------------------------------------------------
