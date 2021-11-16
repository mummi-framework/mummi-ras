# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
"""Patch module."""

import io
import enum
import numpy as np
import logging
LOGGER = logging.getLogger(__name__)


from .names_of_lipids_and_proteins import *

# ------------------------------------------------------------------------------
# Patch Configuration
# ------------------------------------------------------------------------------
class PatchConfig:
    """A Generic Configuration for a Patch."""

    def __init__(self, simname, simtime, tunit, lunit, filenames):
        assert isinstance(simname, str)
        assert isinstance(simtime, (int,float))
        assert isinstance(tunit, str)
        assert isinstance(lunit, str)
        assert isinstance(filenames, dict)

        self.filenames = filenames  # filenames that create the patch
        self.simname = simname      # name of the simulation
        self.simtime = simtime      # simulation time
        self.tunit = tunit          # time unit
        self.lunit = lunit          # length unit


# ------------------------------------------------------------------------------
@enum.unique
class pfPatchType(enum.IntEnum):
    Other = 0
    Lipid_Only = 1
    RAS_Only_Single = 2
    RAS_RAF_Single = 3
    Invalid_Complex = 4


@enum.unique
class gcPatchType(enum.IntEnum):
    ras0raf0 = 0,
    ras1raf0 = 1,
    ras1raf1 = 2,
    ras2raf0 = 3,
    ras2raf1 = 4,
    ras2raf2 = 5,
    invalid = 6

    @staticmethod
    def get_name(_):
        return str(_).split('.')[1]


# ------------------------------------------------------------------------------
# Patch Base Class
# ------------------------------------------------------------------------------
class Patch:
    """Base class for a patch."""

    def __init__(self, patch_id, patch_config, lipid_concentrations,
                 patch_center=None,
                 protein_ids=None, protein_states=None, protein_positions=None,
                 complex_ids=None):

        assert isinstance(patch_id, str)
        assert isinstance(patch_config, PatchConfig)
        assert isinstance(lipid_concentrations, np.ndarray)

        if patch_center is not None:
            assert isinstance(patch_center, np.ndarray)
            assert patch_center.shape == (2,)

        if protein_ids is not None:
            assert (protein_states is not None) and (protein_positions is not None)
            assert isinstance(protein_ids, np.ndarray)
            assert isinstance(protein_states, np.ndarray)
            assert isinstance(protein_positions, np.ndarray)

            assert protein_ids.shape[0] == protein_states.shape[0]
            assert protein_ids.shape[0] == protein_positions.shape[0]

            if complex_ids is not None:
                assert isinstance(complex_ids, np.ndarray)
                assert len(complex_ids.shape) == 2
                assert complex_ids.shape[1] == 2

        else:
            assert (protein_states is None) and (protein_positions is None)
            assert complex_ids is None

        self.id = patch_id
        self.center = patch_center
        self.config = patch_config

        self.concentrations = lipid_concentrations
        self.protein_ids = protein_ids
        self.protein_states = protein_states
        self.protein_positions = protein_positions
        self.complex_ids = complex_ids

    # --------------------------------------------------------------------------
    def __str__(self):

        s = f'{self.id}: '
        s += f'time = ({self.config.simtime} {self.config.tunit}), '
        s += f'concs = {self.concentrations.shape}'
        if self.protein_ids is not None:
            s += f', ids = {list(self.protein_ids)}'
            s += f', states = {list(self.protein_states)}'
        if self.complex_ids is not None:
            s += f', complexes = {[tuple(_) for _ in self.complex_ids]}'
        return s

    def __repr__(self):
        return self.__str__()

    # --------------------------------------------------------------------------
    @staticmethod
    def protein_counts_c3(pstates):
        nras = sum([int('RAS' in _) for _ in pstates])
        nraf = sum([int('RAF' in _) for _ in pstates])
        return nras, nraf

    def type_c3(self):
        nproteins = self.protein_ids.shape[0]
        ncomplexes = self.complex_ids.shape[0]

        # no proteins present!
        if nproteins == 0:
            assert ncomplexes == 0
            return pfPatchType.Lipid_Only

        # center has to be RAS or zRAF (zRAF and zRAS overlap)
        assert ('RAS' in self.protein_states[0]) or ('zRAFa' == self.protein_states[0])

        # check if all complexes are fully represented
        is_found = np.in1d(self.complex_ids.reshape(-1), self.protein_ids)
        if not np.all(is_found):
            LOGGER.debug(f'Found an invalid type-1 patch: '
                         f'patch = {self.id}, ids = {self.protein_ids}, '
                         f'complexes = {self.complex_ids}')
            return pfPatchType.Invalid_Complex

        # decided to ignore the patches with high aggregation
        nras = sum([int('RAS' in _) for _ in self.protein_states])
        if nras > 6:
            LOGGER.debug(f'Found an invalid type-2 patch: '
                         f'patch = {self.id}, nras = {nras}')
            return pfPatchType.Invalid_Complex

        # 2021.02.19. additional condition to ignore patches
        nraf = sum([int('RAF' in _) for _ in self.protein_states])
        if nraf + nras > 4:
            LOGGER.info(f'Found an invalid type-3 patch: '
                        f'patch = {self.id}, nraf = {nraf}, nras = {nras}')
            return pfPatchType.Invalid_Complex

        # Only one bead present - has to be only x1 RAS-only
        if nproteins == 1:
            assert ncomplexes == 0
            return pfPatchType.RAS_Only_Single

        # one complex present!
        if ncomplexes == 1 and nproteins == 2:
            assert ('RAF' in self.protein_states[1]) or (self.protein_states[1] == 'zRASa')
            return pfPatchType.RAS_RAF_Single

        # otherwise, all is good.
        return pfPatchType.Other

    # --------------------------------------------------------------------------
    def protein_counts(self):

        if self.protein_states is None:
            return 0, 0, 0
        if self.protein_states.shape[0] == 0:
            return 0, 0, 0

        nras = np.sum(np.in1d(self.protein_states, MACRO_RAS_NAMES))
        nrasf = np.sum(np.in1d(self.protein_states, MACRO_RASF_NAMES))
        nraf = np.sum(np.in1d(self.protein_states, MACRO_RAF_NAMES))

        assert self.protein_states.shape[0] == nras + nrasf + nraf
        assert nrasf == nraf
        return nras, nrasf, nraf

    def type(self):

        nras, nrasf, nraf = self.protein_counts()

        # no RAF
        if nrasf == 0 and nraf == 0:
            if nras == 0:
                return gcPatchType.get_name(gcPatchType.ras0raf0)
            if nras == 1:
                return gcPatchType.get_name(gcPatchType.ras1raf0)
            if nras == 2:
                return gcPatchType.get_name(gcPatchType.ras2raf0)

        # 1 RAS-RAF
        elif nrasf == 1 and nraf == 1:
            if nras == 0:
                return gcPatchType.get_name(gcPatchType.ras1raf1)
            if nras == 1:
                return gcPatchType.get_name(gcPatchType.ras2raf1)

        # 2 RAS-RAF
        elif nrasf == 2 and nraf == 2:
            if nras == 0:
                return gcPatchType.get_name(gcPatchType.ras2raf2)

        # else, did not find
        return gcPatchType.get_name(gcPatchType.invalid)

    # --------------------------------------------------------------------------
    def complex_idxs(self):

        if self.complex_ids is None:
            return None

        _, pid2pidx = np.unique(self.protein_ids, return_inverse=True)
        assert self.protein_ids.shape == _.shape
        assert 1 == _.min() and self.protein_ids.shape[0] == _.max()
        return pid2pidx[self.complex_ids-1]

    def complex_midpoints(self, force_center_to_first_ras=False):

        # no protein
        if self.protein_positions is None:
            assert self.protein_states is None
            assert self.complex_ids is None
            return None, None

        nprots = self.protein_positions.shape[0]
        assert nprots == self.protein_states.shape[0]

        # no protein
        if 0 == nprots:
            assert self.complex_ids is None
            return None, None

        # first protein is already at the center
        dx = 30./self.concentrations.shape[0]
        d = np.abs(self.protein_positions[0, :2] - np.array([15., 15.]))
        assert d[0] < 0.5*dx and d[1] < 0.5*dx, 'Should be centered around first protein'

        # ----------------------------------------------------------------------
        def mid_point(idxa, idxb):
            if force_center_to_first_ras and idxa == 0:
                return self.protein_positions[idxa]
            elif force_center_to_first_ras and idxb == 0:
                return self.protein_positions[idxb]
            else:
                return 0.5 * (self.protein_positions[idxa] + self.protein_positions[idxb])

        # ----------------------------------------------------------------------
        # 1 protein
        if 1 == nprots:
            assert self.complex_ids is None
            return self.protein_positions[0], None

        # 2 proteins that are both RAS
        if 2 == nprots and self.complex_ids is None:
            return self.protein_positions[0], self.protein_positions[1]

        # 2 proteins that are a complex
        if 2 == nprots and self.complex_ids is not None:
            assert (1,2) == self.complex_ids.shape
            return mid_point(0, 1), None

        # 3 proteins (1 RAS and 1 RAS-RAF complex)
        if 3 == nprots:
            assert (1,2) == self.complex_ids.shape

            cidxs = self.complex_idxs()
            pcid = cidxs[0]
            ncid = np.setdiff1d(np.array([0, 1, 2]), pcid)
            assert ncid.shape[0] == 1

            posnc = self.protein_positions[ncid[0]]
            posc = mid_point(pcid[0], pcid[1])
            return (posnc, posc) if ncid[0] == 0 else (posc, posnc)

        # 4 proteins (2 RAS-RAF complex)
        if 4 == nprots:
            assert (2,2) == self.complex_ids.shape

            cidxs = self.complex_idxs()
            pa,pb = cidxs[0], cidxs[1]
            posa = mid_point(pa[0], pa[1])
            posb = mid_point(pb[0], pb[1])
            return (posa, posb) if 0 in pa else (posb, posa)

        assert False, f'Unknown configuration: {self.protein_states} {self.complex_ids}'
        return None, None

    # --------------------------------------------------------------------------
    # distance between two dimers
    def distance_between_dimers(self):
        a,b = self.complex_midpoints()
        return np.linalg.norm(b-a) if (a is not None and b is not None) else np.nan

    # --------------------------------------------------------------------------
    def wrap_to_center(self, from_p, psize=30.):

        assert isinstance(from_p, np.ndarray)
        assert from_p.shape == (2,)

        LOGGER.debug(f'Wrapping {from_p} to center')
        npix = self.concentrations.shape[0]

        dx = psize / npix
        hnpix, hpsize = npix // 2, 0.5 * psize

        to_p = np.array([hpsize, hpsize])
        to_g = np.array([hnpix, hnpix])

        # ----------------------------------------------------------------------
        from_g = np.round(from_p / dx).astype(int)

        # compute the displacement in grid indices
        for d in range(2):
            if from_g[d] == npix:
                from_g[d] = 0

        disg = to_g - from_g
        for d in range(2):
            if disg[d] < 0:
                disg[d] += npix

        # compute the displacement in physical coordinates
        disp = disg * dx

        # ----------------------------------------------------------------------
        # perform the displacements
        self.protein_positions[:,:2] += disp
        for d in range(2):
            self.protein_positions[self.protein_positions[:,d] < 0, d] += psize
            self.protein_positions[self.protein_positions[:,d] > psize, d] -= psize

        self.concentrations = np.roll(self.concentrations, shift=(disg[1], disg[0]), axis=(0, 1))

        # ----------------------------------------------------------------------
        # make sure that the protein closest to the center is no more
        # than half a pixel away!
        pnorm = np.abs(self.protein_positions[:,:2]-to_p).max(axis=1)
        pnorm = pnorm.min()
        assert pnorm <= 0.5*dx,\
               f'Protein closest to center is more than half a pixel away ({pnorm} > {0.5*dx})'

        LOGGER.debug(f'Centered ({self.id}) at {from_p}; disg = {disg}, disp = {disp}')

    # --------------------------------------------------------------------------
    def subsample_mean(self, psize=30., npixels=5):
        return Patch._subsample_mean(self.concentrations, psize, npixels)

    def subsample_intg(self, psize=30., npixels=5):
        return Patch._subsample_intg(self.concentrations, psize, npixels)

    # --------------------------------------------------------------------------
    def reorient(self):
        """This function reorients a patch with respect to protein positions.
        -  0 proteins: do nothing
        - >0 proteins: center around the first RAS
        - >1 proteins: center around the first RAS and
                       put RAS-RAS vector on [0,45] degrees
        """
        # ----------------------------------------------------------------------
        # no proteins. nothing to do
        if self.protein_positions is None:
            return self.concentrations, self.protein_positions

        if len(self.protein_positions.shape) == 0:
            return self.concentrations, self.protein_positions

        if self.protein_positions.shape[0] == 0:
            return self.concentrations, self.protein_positions

        # ----------------------------------------------------------------------
        # start working
        psize, gsize = 30., self.concentrations.shape[0]
        is_centered_at_0 = False
        assert not is_centered_at_0  # need to debug is_centered_at_0

        if is_centered_at_0:
            pcent = np.array([0., 0., 0.])
            pext = np.array([-0.5 * psize, 0.5 * psize, 0.])
        else:
            pcent = np.array([0.5 * psize, 0.5 * psize, 0.])
            pext = np.array([0., psize, 0.])

        p2g = psize / gsize
        g2p = gsize / psize

        # ----------------------------------------------------------------------
        def _center_at(_img, _pos, _p):
            disg = np.round(g2p * (pcent - _p)).astype(np.int)
            if disg[0] == 0 and disg[1] == 0:
                return _img, _pos
            _pos += p2g * disg
            for i in range(2):
                _pos[_pos[:, i] < pext[0], i] += psize
                _pos[_pos[:, i] > pext[1], i] -= psize
            _img = np.roll(_img, shift=(disg[1], disg[0]), axis=(0, 1))
            return _img, _pos

        def _flip(_img, _pos, _dir):
            if _dir == 'x':
                _pos[:, 0] = psize - _pos[:, 0]
                _img = np.flip(_img, axis=1)
            elif _dir == 'y':
                _pos[:, 1] = psize - _pos[:, 1]
                _img = np.flip(_img, axis=0)
            else:
                assert 0
            return _img, _pos

        def _rot(_img, _pos, _dir):
            _pos[:, [1, 0]] = _pos[:, [0, 1]] - pcent[:2]
            if _dir == 'ccw':
                _pos[:, 0] *= -1
                _img = np.rot90(_img, axes=(1, 0))
            elif _dir == 'cw':
                _pos[:, 1] *= -1
                _img = np.rot90(_img, axes=(0, 1))
            else:
                assert 0
            _pos += pcent
            return _img, _pos

        # ----------------------------------------------------------------------
        # create a copy
        img = np.copy(self.concentrations)
        ppositions = np.copy(self.protein_positions)

        pa, pb = self.complex_midpoints(force_center_to_first_ras=True)

        # ----------------------------------------------------------------------
        # 1 protein (center around first RAS)
        if pb is None:
            img, ppositions = _center_at(img, ppositions, pa)
            return img, ppositions

        # ----------------------------------------------------------------------
        # 2 proteins

        # center the first
        img, ppositions = _center_at(img, ppositions, pa)

        # compute where the vector falls
        disp = pb - pa
        ang = np.rad2deg(np.arctan2(disp[1], disp[0]))
        if ang < 0:
            ang += 360
        q = int(np.floor(ang / 45))
        assert 0 <= q < 8

        # based on the half-quadrant, different transformations are needed
        if q == 1:
            img, ppositions = _flip(img, ppositions, 'y')
            img, ppositions = _rot(img, ppositions, 'ccw')

        elif q == 2:
            img, ppositions = _rot(img, ppositions, 'cw')

        elif q == 3:
            img, ppositions = _flip(img, ppositions, 'x')

        elif q == 4:
            img, ppositions = _flip(img, ppositions, 'x')
            img, ppositions = _flip(img, ppositions, 'y')

        elif q == 5:
            img, ppositions = _flip(img, ppositions, 'y')
            img, ppositions = _rot(img, ppositions, 'cw')

        elif q == 6:
            img, ppositions = _rot(img, ppositions, 'ccw')

        elif q == 7:
            img, ppositions = _flip(img, ppositions, 'y')

        # ----------------------------------------------------------------------
        return img, ppositions

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    @staticmethod
    def _subsample_mean(data, psize, npixels):
        '''
        Subsample a patch using mean values of grid points
        * patch is a (n, n, c) numpy array (e.g., 37x37x14)
        * psize is the physical span of the patch (i.e., 30 nm)
        * npixels is the size of the requested subsampling (e.g., 5)
        '''
        assert len(data.shape) == 3
        assert data.shape[0] == data.shape[1]

        ngridpts = data.shape[0]
        nlipids = data.shape[-1]

        if ngridpts == npixels:
            return data

        # LOGGER.debug('Subsampling using mean {} to {}'.format(patch.shape, npixels))
        # ----------------------------------------------------------------------
        # decide how to divide the patch
        if npixels == 5 and ngridpts == 37:
            w = np.array([7, 7, 9, 7, 7], dtype=np.uint32)

        elif npixels == 37 and ngridpts == 73:
            w = 2 * np.ones(npixels, dtype=np.uint32)
            w[18] = 1

        else:
            raise Exception('{} to {} is not implemented'.format(ngridpts, npixels))

        assert w.sum() == ngridpts and w.shape[0] == npixels
        subdata = np.zeros((npixels, npixels, nlipids), dtype=data.dtype)

        w = np.cumsum(np.concatenate((np.array([0]), w)))
        for y in range(npixels):
            for x in range(npixels):
                subgrid = data[w[x]:w[x + 1], w[y]:w[y + 1], :]
                subdata[x, y] = subgrid.mean(axis=(0, 1))

        return subdata

    # --------------------------------------------------------------------------
    @staticmethod
    def _subsample_intg(data, psize, npixels):
        '''
        Subsample a patch using trapezoidal integration.
        * patch is a (n, n, c) numpy array (e.g., 37x37x14)
        * psize is the physical span of the patch (i.e., 30 nm)
        * npixels is the size of the requested subsampling (e.g., 5)
        '''
        # tested with these values, but not really need to assume!
        assert (npixels == 5 and data.shape[0] == 37) or \
               (npixels == 5 and data.shape[0] == 73) or \
               (npixels == 37 and data.shape[0] == 73) or \
               (npixels == 5 and data.shape[0] == 72)

        # other assertions
        assert len(data.shape) == 3
        assert data.shape[0] == data.shape[1]

        from scipy.interpolate import interp2d

        LOGGER.debug('Subsampling using interpolation {} to {}'.format(data.shape, npixels))

        # -----------------------------------------------------
        # in this function, I am calling
        # "grid" or "grid point" to refer to the points and spacing
        # of the macro grid (also, of the input patch)
        # and "pixel" as one of the pixels in the subsampled patch

        ngridpts = data.shape[0]
        nlipids = data.shape[-1]

        ngrid_per_pixel = (ngridpts - 1) / npixels

        grid_dx = psize / (ngridpts - 1)
        pixel_dx = psize / npixels
        pixel_area = pixel_dx * pixel_dx
        # print (grid_dx, pixel_dx, pixel_area)

        # the "end points" of the pixels of small grid
        pixel_endpoints = ngrid_per_pixel * np.arange(npixels + 1)

        # -----------------------------------------------------
        # the grid points for each pixel (in a single dimension)
        gridpoints = []
        physpoints = []
        npts = int(np.ceil(ngrid_per_pixel)) + 1
        for i in range(npixels):
            # bounds of this pixel
            l, r = pixel_endpoints[i], pixel_endpoints[i + 1]

            # create a small grid for this pixel
            p = np.arange(np.ceil(l), np.ceil(r))
            p = np.unique(np.append(p, [l, r]))

            assert np.abs(p.shape[0] - npts) <= 1

            gridpoints.append(p)
            physpoints.append(grid_dx * p)

        # print ('gridpoints', gridpoints)

        # -----------------------------------------------------
        # now, actually subsample
        subdata = np.zeros((npixels, npixels, nlipids), dtype=data.dtype)

        # remove negative values
        data = np.maximum(data, 0.0)

        # create a grid axis for the input patch
        grid_axs = np.arange(0, ngridpts)

        for c in range(nlipids):

            # create an interpoltor for this lipid concentration
            f = interp2d(grid_axs, grid_axs, data[:, :, c],
                         kind='linear', copy=False,
                         bounds_error=True, fill_value=np.nan)

            # for each of the pixels
            for iy in range(npixels):
                for ix in range(npixels):
                    # extract the subgrid for this pixel
                    # not optimal: i should interpolate only when needed
                    subgrid = f(gridpoints[ix], gridpoints[iy])

                    # now, we do the 2d trapezoidal integration
                    subgrid = np.trapz(subgrid, x=physpoints[iy], axis=0)
                    subgrid = np.trapz(subgrid, x=physpoints[ix])

                    # need to normalize the area under the curve
                    # by the pixel area
                    subdata[iy, ix, c] = subgrid / pixel_area

        # -----------------------------------------------------
        return subdata

    # --------------------------------------------------------------------------
    @staticmethod
    def save_npz(fptr, p):

        # write the patch to either a filename or an in-memory byte stream
        assert isinstance(fptr, (str, io.BytesIO))
        assert isinstance(p, Patch)
        np.savez_compressed(fptr,
                            id = p.id, config = p.config,
                            concentrations = p.concentrations,
                            center = p.center,
                            ids = p.protein_ids,
                            states = p.protein_states,
                            positions = p.protein_positions,
                            complex_ids = p.complex_ids)
        return fptr

    # --------------------------------------------------------------------------
    @staticmethod
    def load_npz(fptr):
        assert isinstance(fptr, (str, io.BytesIO))
        data = np.load(fptr, allow_pickle=True)

        # ----------------------------------------------------------------------
        # Campaign 1 patches
        if 'rasIds' in data.files:

            config = PatchConfig(str(data['simname']),
                                 dict(data['filenames'][()]),
                                 float(data['simtime']),
                                 str(data['tunit']), str(data['lunit']))

            p = Patch(str(data['id']), config, data['concentrations'],
                      data['center'], data['rasIds'], data['rasStates'],
                      data['rasPositions'])

        # ----------------------------------------------------------------------
        # Campaign 3 patches
        elif 'simname' in data.files:

            config = PatchConfig(str(data['simname']),
                                 float(data['simtime']),
                                 str(data['tunit']), str(data['lunit']),
                                 data['filenames'][()])

            p = Patch(str(data['id']), config, data['concentrations'],
                      data['center'], data['ids'], data['states'],
                      data['positions'], data['complex_ids'][()])

        # ----------------------------------------------------------------------
        # Campaign 4 patches
        else:
            p = Patch(str(data['id']), data['config'][()],
                      data['concentrations'], data['center'][()],
                      data['ids'][()], data['states'][()],
                      data['positions'][()], data['complex_ids'][()])

        # ----------------------------------------------------------------------

        data.close()
        return p

    # --------------------------------------------------------------------------
    @staticmethod
    def write_combined(filename, patches, simname, simpath):
        np.savez_compressed(filename,
                            simname = simname, simpath = simpath,
                            ids = [p.id for p in patches],
                            types = [p.type() for p in patches],
                            protein_ids=[p.protein_ids for p in patches],
                            protein_states=[p.protein_states for p in patches],
                            protein_positions = [p.protein_positions for p in patches],
                            complex_ids = [p.complex_ids for p in patches],
                            concentrations = [p.concentrations for p in patches])

    @staticmethod
    def write_combined_c3(filename, patches, simname, time, fgrid, fhycop):
        np.savez_compressed(filename,
                            fgrid = fgrid, fhycop = fhycop,
                            simname = simname, time = time,
                            types = [int(p.type()) for p in patches],
                            ids = [p.ids for p in patches],
                            states = [p.states for p in patches],
                            positions = [p.positions for p in patches],
                            centers = [p.center for p in patches],
                            complex_ids = [p.complex_ids for p in patches],
                            complex_dists = [p.complex_dists for p in patches],
                            concentrations = [p.concentrations for p in patches])

    # --------------------------------------------------------------------------
    def render(self, show=True):

        # in C3 patches, positions were stored with respect to the center
        is_c3_patch = 'snapshot' not in self.config.filenames
        if self.protein_positions is not None:
            positions = np.copy(self.protein_positions)
            if is_c3_patch:
                positions += np.array([15., 15.])
        else:
            positions = None

        # ----------------------------------------------------------------------
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        psz, ssz = 30., 15.
        fig, axs = plt.subplots(nrows=2, ncols=8, figsize=(16, 4))
        for i in range(14):

            r, c = i // 8, i % 8
            if i == 13:
                c = 7
            ax = axs[r][c]

            l = self.concentrations[:, :, i]
            im = ax.imshow(l, origin='bottom', extent=[0, psz, 0, psz])

            if positions is not None:
                ax.scatter(positions[:, 0], positions[:, 1], color='r', s=ssz)

            ax.set_xlim([0, psz])
            ax.set_ylim([0, psz])
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(MACRO_FIELDS_DENSITIES[i], fontsize=9)

            if True:
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                fig.colorbar(im, cax=cax, ticks=[l.min(), l.mean(), l.max()])

        for c in [5, 6]:
            axs[1][c].axis('off')

        fig.suptitle(self.id)
        plt.tight_layout()
        plt.subplots_adjust(top=0.85)

        if show:
            plt.show()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
