# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
import os
import time
import shutil
import numpy as np
import logging

from mummi_core.utils import Naming
from mummi_ras.datastructures.patch import pfPatchType, gcPatchType
from mummi_ras.datastructures.patch import Patch as pfPatch
from mummi_ras.interfaces import get_io

if __package__ or "." in __name__:
    from . import dynim
else:
    import dynim

from .model_pretrained import PretrainedModel as PatchEncoder

LOGGER = logging.getLogger(__name__)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
DO_JOINT_QUEUES = False
DO_PROTEIN_SCALING = False

DO_DISABLE_MU18 = False
KNOWN_MUS = [-24,0,12,18]

marker = '------------------------------'


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
class PatchSelector(object):

    # --------------------------------------------------------------------------
    def __init__(self, stype: str, workspace: str, config: dict):

        assert stype in ['random', 'importance']
        assert isinstance(workspace, str)
        assert isinstance(config, dict)

        buffer_size = int(config['buffer_sz'])
        min_cands_b4_sel = int(config['min_cands_b4_sel'])
        min_rand_b4_importance = int(config['min_rand_b4_importance'])
        assert min_cands_b4_sel >= 0
        assert min_rand_b4_importance >= 0
        assert buffer_size >= 0

        LOGGER.info(f'{marker} Initializing Patch Selector')
        self.type = stype
        self.samplers = []
        self.recent_selections = []
        self.workspace = workspace
        os.makedirs(self.workspace, exist_ok=True)

        # ----------------------------------------------------------------------
        # random sampler
        if self.type == 'random':
            self.Nsamplers = 1
            self.sampler_names = ['macro-sampler']
            self.samplers = [dynim.SamplerRandom('macro-sampler', workspace,
                                                 buffer_size=buffer_size,
                                                 min_cands_b4_sel=min_cands_b4_sel)]
            return

        # ----------------------------------------------------------------------
        # ml-based importance sampler for C4!
        mpath = os.path.join(config['path'], config['encoder'])
        lspace_fname = os.path.join(mpath, 'lspace.idx')

        self.patchencoder = PatchEncoder(mpath)
        self.Pshape = self.patchencoder.input_shp

        # create these types of samplers!
        ''' * Q1 — Lipid only
            * Q2 — 4x monomers (RASa, RASb, RASa-RAFa, RASb-RAFb)
            * Q3 — 4+3+2+1 == dimers
            * May be: 10%, 50%, 40% (estimated: 20-30 K)
        '''
        self.Nsamplers = 4
        self.sampler_names = ['lipid-only', 'monomers', 'dimers', 'invalid']
        self.sampler_ratios = np.array([0.1, 0.5, 0.4, 0.0])

        # mapping between patch type and queue type (notes from ML Meeting, Jun 03)
        self.patch2sampler = {
            'ras0raf0': 0,
            'ras1raf0': 1,
            'ras1raf1': 1,
            'ras2raf0': 2,
            'ras2raf1': 2,
            'ras2raf2': 2,
            'invalid': 3
        }

        self.sampler_names = np.array(self.sampler_names)
        assert self.sampler_names.shape == self.sampler_ratios.shape
        assert np.isclose(sum(self.sampler_ratios), 1.0)

        # latent space (load the default one from resources)
        if DO_JOINT_QUEUES:
            self.lspace = dynim.HDSpace()
            self.lspace.restore(lspace_fname)

        for i, sname in enumerate(self.sampler_names):
            if sname == 'dimers':
                s = dynim.SamplerBinned(sname, workspace,
                                        novelty_factor = 0.2,
                                        binning_type='1d-nonperiodic',
                                        min_cands_b4_sel=min_cands_b4_sel)
            else:
                s = dynim.SamplerImportance(self.sampler_names[i], workspace,
                                            buffer_size=buffer_size,
                                            min_cands_b4_sel = min_cands_b4_sel,
                                            min_rand_b4_importance = min_rand_b4_importance)
            if DO_JOINT_QUEUES:
                s.set_hdspace(self.lspace)
            else:
                s.set_hdspace(lspace_fname)

            self.samplers.append(s)

        # ----------------------------------------------------------------------
        # additional info to update ranks in "other" queue based on protein counts
        if DO_PROTEIN_SCALING:
            self.PROT_CNT_QUEUE = 3
            self.samplers[self.PROT_CNT_QUEUE].set_rank_scaler(self.rank_scalar_protein_count)

            # a dictionary of the scaling of ranks, wrt protein counts
            #self.pcnts_scales = {0: 1., 1: 1., 2: 1., 3: 0.98, 4: 0.97, 5: 0.96, 6: 0.95}
            self.pcnts_scales = {0: 1., 1: 1., 2: 1., 3: 0.95, 4: 0.9, 5: 0.85, 6: 0.8}

            # a dictionary of all patches in the "other" queue
            # for each patch id, stores a tuple (nras, nraf) [[ncomplexes = nraf]]
            self.patch2pcnts = {}

        # ----------------------------------------------------------------------
        if DO_DISABLE_MU18:
            for i in range(len(self.samplers)):
                self.samplers[i].set_rank_scaler(self.rank_scalar_mu)

        # ----------------------------------------------------------------------
        LOGGER.info(f'{marker} Initialized Patch Selector')
        for s in self.samplers:
            LOGGER.debug(f'created {s}')

    # --------------------------------------------------------------------------
    def rank_scalar_protein_count(self, patch_id, patch_rank):
        assert DO_PROTEIN_SCALING
        assert 0
        #LOGGER.debug(f'Scaling rank ({patch_rank}) for id ({patch_id})')

        pcounts = self.patch2pcnts.get(patch_id, None)
        if pcounts is None:
             LOGGER.error(f'Failed to find the protein counts for ({patch_id})')
             return patch_rank

        nraf = pcounts[1]
        s = self.pcnts_scales.get(nraf, None)
        if s is None:
             LOGGER.error(f'Failed to find the scaling factor for nraf={nraf}')
             return patch_rank

        return np.float32(patch_rank * s)

    # --------------------------------------------------------------------------
    def rank_scalar_mu(self, patch_id, patch_rank):
        assert DO_DISABLE_MU18
        mu = Naming.get_mu_from_simname(patch_id)
        assert mu in KNOWN_MUS, f'Found invalid mu = {mu}'
        return patch_rank if mu != 18 else 0.

    # --------------------------------------------------------------------------
    def add_candidates(self, patches):

        assert isinstance(patches, list)
        assert all([isinstance(p, pfPatch) for p in patches])

        if len(patches) == 0:
            return

        LOGGER.profile(f'Adding {len(patches)} candidate patches')

        # ----------------------------------------------------------------------
        # random sampler
        if self.type == 'random':
            HDsamples = [dynim.HDPoint(p.id, np.empty(0)) for p in patches]
            self.samplers[0].add_candidates(HDsamples)
            self.checkpoint()
            return

        # ----------------------------------------------------------------------
        # ml-based importance sampler

        # for each patch, figure out which sampler they will go to
        sampler_idxs = [self.add_criteria(_) for _ in patches]
        assert min(sampler_idxs) >= 0 and max(sampler_idxs) <= self.Nsamplers

        # make sure all the shapes are correct
        pconcs = [_.concentrations for _ in patches]
        pconcs = np.array(pconcs)

        #if pconcs.shape[1:] != self.Pshape:
        #    LOGGER.error(f'Forcing C4 patches {pconcs.shape[1:]} to C3 shape {self.Pshape} for ML')
        #    pconcs = pconcs[:,:self.Pshape[0],:self.Pshape[1],:self.Pshape[2]]
        assert pconcs.shape[1:] == self.Pshape

        lcoords = self.patchencoder.encode(pconcs)

        # for each patch, get the latent space coordinates
        # and figure out which sampler they will go to
        encoded_by_ids = [[] for _ in range(self.Nsamplers)]
        for p,s,l in zip(patches, sampler_idxs, lcoords):
            if isinstance(self.samplers[s], dynim.SamplerBinned):
                d = p.distance_between_dimers()
                d = np.float32(d)   # all latent coords need to be 32-bit

                LOGGER.debug(f'Attaching dimer distance (={d}) for patch {p.id}')
                l = np.concatenate((l, [d]))

            encoded_by_ids[s].append(dynim.HDPoint(p.id, l))

        # ----------------------------------------------------------------------
        # protein-count based scaling requires storing the protein counts
        if DO_PROTEIN_SCALING:
            for p,s in zip(patches, sampler_idxs):
                if s == self.PROT_CNT_QUEUE:        # but, only for others queue
                    self.patch2pcnts[p.id] = p.protein_counts()

        # ----------------------------------------------------------------------
        # now, add these patches to the respective samplers
        for sidx in range(self.Nsamplers):
            encoded = encoded_by_ids[sidx]
            n = len(encoded)
            if n == 0:
                continue

            # LOGGER.debug('Adding {} patches to sampler {}'.format(n, sidx))
            self.samplers[sidx].add_candidates(encoded)

        LOGGER.profile(f'Added {len(patches)} candidate patches')
        self.checkpoint()

    # --------------------------------------------------------------------------
    def select(self, k, do_confirm=True):

        assert isinstance(k, int)
        assert k >= 0

        LOGGER.info(f'Selecting {k} patches (confirm = {do_confirm})')

        # ----------------------------------------------------------------------
        # random sampler
        if self.type == 'random':
            selections = self.samplers[0].select(k, do_confirm)
            self.recent_selections = [(s.id, 0) for s in selections]

        # ----------------------------------------------------------------------
        # ml-based importance sampler
        else:
            nsamplers = len(self.samplers)
            num_cands = np.array([s.num_candidates() for s in self.samplers]).astype(int)
            num_selected = np.array([s.num_selected() for s in self.samplers]).astype(int)

            sampler_cnts = self.select_criteria(k, self.sampler_ratios, num_cands, num_selected)
            assert len(sampler_cnts) == nsamplers

            self.recent_selections = []
            for sidx in range(nsamplers):

                cnt = sampler_cnts[sidx]
                selection = self.samplers[sidx].select(int(float(cnt)), do_confirm)
                selection = [(s.id, sidx) for s in selection]

                self.recent_selections.extend(selection)

        # ----------------------------------------------------------------------
        selections = [s[0] for s in self.recent_selections]

        # no need to maintain this if i am confirming
        if do_confirm:
            self.recent_selections = []
            self.checkpoint()

        LOGGER.info(f'Selected {len(selections)} patches: {selections}')
        return selections

    # --------------------------------------------------------------------------
    def confirm_selection(self, patch_ids):

        assert isinstance(patch_ids, list)
        assert all([isinstance(p, str) for p in patch_ids])

        n = len(patch_ids)
        if n == 0:
            return

        LOGGER.info(f'Confirming selection of {n} patches')

        assert n == len(self.recent_selections)
        assert all([patch_ids[i] == self.recent_selections[i][0] for i in range(n)])

        # ----------------------------------------------------------------------
        # random sampler
        if self.type == 'random':
            self.samplers[0].confirm_selections(patch_ids)
            self.checkpoint()
            return

        # ----------------------------------------------------------------------
        # ml-based importance sampler

        # split these patches based on the selectors they were used in
        selections = [[] for i in range(self.Nsamplers)]
        for i in range(n):
            selections[self.recent_selections[i][1]].append(patch_ids[i])

        # now, confirm the selections in dynim
        for sidx in range(self.Nsamplers):
            if len(selections[sidx]) == 0:
                continue

            self.samplers[sidx].confirm_selections(selections[sidx])

        self.recent_selections = []
        self.checkpoint()

    # --------------------------------------------------------------------------
    def get_weights(self, patches = None):

        if patches is None:
            LOGGER.info('getting weights for all patches')
        else:
            LOGGER.info(f'getting weights for {patches.shape}, patches: {patches}')

        do_norm = True
        do_history = True
        try:
            # compute weights for all patches
            sweights = [s.get_weights(do_norm) for s in self.samplers]
            LOGGER.debug(f'sweights = {sweights}')

            sim_ids = np.concatenate([_[0] for _ in sweights])
            sim_wgts = np.concatenate([_[1] for _ in sweights])
            LOGGER.debug(f'sim_ids = {sim_ids}')
            LOGGER.debug(f'sim_wgts = {sim_wgts}')

            if sim_ids.shape[0] == 0:
                raise Exception('Failed to find any selected patches!')

            # filter on requested patches only
            if patches is not None:
                LOGGER.debug(f'filtering on requested patches: {patches}')
                patches = np.array(patches)

                # sim_ids is likely a large list
                # remove all the ones that are not requested
                indices = np.where(np.in1d(sim_ids, patches))[0]
                LOGGER.debug(f'filter for = {indices}')

                if indices.shape[0] == 0:
                    raise Exception('Could not find any requested patches in sims from samplers')

                sim_ids = sim_ids[indices]
                sim_wgts = sim_wgts[indices]

                # now, sim_ids contains only the patches that are required
                # we want to reorder them wrt the order of patches
                indices = np.array([np.where(sim_ids==_)[0][0] for _ in patches])
                LOGGER.debug(f'reorder using {indices}')

                sim_ids = sim_ids[indices]
                sim_wgts = sim_wgts[indices]

                # the ids should be exactly like input patches
                if not np.array_equal(sim_ids, patches):
                    raise Exception(f'sim_ids != patches. {sim_ids} != {patches}')

        except Exception as e:
            LOGGER.error(f'Weight computation failed due to error: {e}')
            return patches, np.ones(patches.shape[0], dtype=np.float32)

        LOGGER.info(f'Successfully calculated weights: {sim_ids} {sim_wgts}')
        if do_history:
            tfile = os.path.join(self.workspace, "weights_history.tar")
            wfile = 'weights_{}.npz'.format(time.strftime("%Y%m%d-%H%M%S"))
            LOGGER.debug(f'Writing weights to ({tfile}) / ({wfile})')
            get_io('taridx').save_npz(tfile, wfile,
                                      {'sims': sim_ids, 'weights': sim_wgts})
        return sim_ids, sim_wgts

    # --------------------------------------------------------------------------
    def test(self):
        for sampler in self.samplers:
            sampler.test()

    # --------------------------------------------------------------------------
    def checkpoint(self, do_lspace = True):

        # ----------------------------------------------------------------------
        # checkpoint the patch counts
        if DO_PROTEIN_SCALING:
            n = len(self.patch2pcnts)
            fname = os.path.join(self.workspace, 'PatchSelector-pcounts.npy')
            LOGGER.debug(f'Checkpointing ({fname}) with protein counts for {n} patches')
            if os.path.isfile(fname):
                shutil.move(fname, fname[:-3]+'bak.npy')
            np.save(fname, self.patch2pcnts)

        # ----------------------------------------------------------------------
        for sampler in self.samplers:
            sampler.checkpoint()

    def restore(self):

        LOGGER.info(f'{marker} Restoring Patch Selector')
        # ----------------------------------------------------------------------
        # restore the patch counts
        self.patch2pcnts = {}
        if DO_PROTEIN_SCALING:
           fname = os.path.join(self.workspace, 'PatchSelector-pcounts.npy')
           LOGGER.debug(f'Restoring {fname}')
           self.patch2pcnts = np.load(fname, allow_pickle=True).item()
           assert isinstance(self.patch2pcnts, dict)
           n = len(self.patch2pcnts)
           LOGGER.debug(f'Found protein counts for {n} patches')

        # ----------------------------------------------------------------------
        # now, restore the rest of the dynim files
        restore_flags = [s.restore() for s in self.samplers]
        nrestored = sum(restore_flags)
        if nrestored == 0:
            LOGGER.error(f'Failed to Restore Patch Selector')
            return False

        if nrestored == len(self.samplers):
            LOGGER.info(f'{marker} Successfully Restored Patch Selector')
            for s in self.samplers:
                LOGGER.debug(s)
            return True

        LOGGER.error(f'Restore flags: {restore_flags}')
        assert 0, 'Restored some samplers but not all'

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    def add_criteria(self, patch):
        assert isinstance(patch, pfPatch)
        return self.patch2sampler[patch.type()]

    @staticmethod
    def add_criteria_c3(patch):

        assert isinstance(patch, pfPatch)

        # assign queues based on the patch type
        ptype = patch.type()
        if ptype == pfPatchType.Lipid_Only:
            return 0
        if ptype == pfPatchType.RAS_Only_Single:
            return 1
        if ptype == pfPatchType.RAS_RAF_Single:
            return 2
        if ptype == pfPatchType.Invalid_Complex:
            return 4
        return 3

    @staticmethod
    def select_criteria(k, ratios, num_cands, num_selected):

        # helper function to adjust the ratios based on counts
        def _adjust_ratios(orig_ratios, counts):
            ratios = np.copy(orig_ratios)
            ratios[np.where(counts <= 0)[0]] = 0
            return np.divide(ratios, ratios.sum())

        assert isinstance(num_cands, np.ndarray)
        assert isinstance(num_selected, np.ndarray)
        assert isinstance(ratios, np.ndarray)
        assert num_cands.shape == num_selected.shape
        assert num_cands.shape == ratios.shape
        assert isinstance(k, int)
        assert k >= 0
        assert ratios.min() >= 0. and ratios.max() <= 1.

        nsamplers = num_cands.shape[0]
        tot_cands = num_cands.sum()
        tot_selected = num_selected.sum()
        selected = np.zeros(nsamplers, dtype=int)

        LOGGER.info(f'Computing select_criteria for k={k} with candidates={num_cands}, selected={num_selected}')

        if k == 0 or tot_cands == 0:
            return selected

        # can only select these many at the most!
        krequested = min(k, tot_cands)
        k = krequested

        # initial target based on num_selected (but using adjusted ratios)
        ratios = _adjust_ratios(ratios, num_cands)

        target = (tot_selected + k) * ratios
        target = np.maximum(0., target - num_selected)
        target *= (k/target.sum())

        # should not require more than n attempts!
        for i in range(nsamplers):

            LOGGER.debug(f'attempt {i}: k = {k}')
            LOGGER.debug(f'num_cands = {num_cands}, selected = {selected}')
            LOGGER.debug(f'ratios = {ratios}, target = {target}')

            # use ceil(target) to always select a bit more that needed
            # then use min(target, num_cands) to ensure that target can be selected
            s = np.minimum(np.ceil(target), num_cands).astype(int)
            LOGGER.debug(f'selected in this turn: {s}')

            # if we overselected (due to taking ceil)?
            overselection = s.sum() - k
            if overselection > 0:
                LOGGER.debug('fixing overselection of {} elements'.format(overselection))
                assert overselection <= nsamplers
                # reduce from the ones with largest error
                idxs_of_max_rounding_error = np.argsort(s-target)[::-1]
                s[idxs_of_max_rounding_error[:overselection]] -= 1
                LOGGER.info(f'updated selection: {s}')

            # now, update the overall counts
            selected += s
            num_cands -= s
            k -= s.sum()

            # no more to select
            if k == 0:
                break

            # create new target for the next attempt
            ratios = _adjust_ratios(ratios, num_cands)
            target = k * ratios

        assert krequested == selected.sum()
        assert selected.min() >= 0
        LOGGER.info(f'select_criteria = {selected}')
        return selected

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
