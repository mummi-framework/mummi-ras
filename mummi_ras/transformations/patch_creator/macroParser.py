"""Macro Parser module."""

import os
import io
import re
import numpy as np

import logging
LOGGER = logging.getLogger(__name__)

from mummi_ras.datastructures.names_of_lipids_and_proteins import *


# ------------------------------------------------------------------------------
# This module parses the macro data
#   densities, atoms, and bonds files
# ------------------------------------------------------------------------------
def show_macro_data(lipids, positions, psz, title='', show=True):

    assert isinstance(lipids, np.ndarray)
    assert isinstance(positions, np.ndarray)
    assert isinstance(psz, (float,int))
    assert isinstance(show, bool)

    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    is_patch = np.isclose(psz, 30.)
    ssz = 15 if is_patch else 5

    fig, axs = plt.subplots(nrows=2, ncols=8, figsize=(16,4))
    for i in range(14):

        r,c = i//8, i%8
        if i == 13:
            c = 7
        ax = axs[r][c]

        l = lipids[:,:,i]
        im = ax.imshow(l, origin='bottom',extent=[0,psz,0,psz])
        ax.scatter(positions[:,0], positions[:,1], color='r', s=ssz)

        if 1: #is_patch:
           ax.set_xlim([0,psz])
           ax.set_ylim([0,psz])
           ax.set_xticks([])
           ax.set_yticks([])
        else:
            p = positions[0]
            ax.set_xlim([p[0] - 15, p[0] + 15])
            ax.set_ylim([p[1] - 15, p[1] + 15])

        ax.set_title(MACRO_FIELDS_DENSITIES[i], fontsize=9)

        if True:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(im, cax=cax, ticks=[l.min(), l.mean(), l.max()])

    for c in [5,6]:
        axs[1][c].axis('off')

    fig.suptitle(title)
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)

    if show:
        plt.show()


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
class MacroParser:

    # --------------------------------------------------------------------------
    @staticmethod
    def get_filenames(snapshot_idx, path=''):

        snapshot_name = f'snapshot.{snapshot_idx:08d}'
        filenames = [os.path.join(snapshot_name, 'densities#000000'),
                     os.path.join(snapshot_name, 'atoms#000000'),
                     os.path.join(snapshot_name, 'bonds#000000')]
        return snapshot_name, [os.path.join(path, f) for f in filenames]

    @staticmethod
    def get_filenames_c3(snapshot_idx, path=''):

        snapshot_name = f'snapshot.{snapshot_idx:08d}'
        filenames = [f'cfield.iter{snapshot_idx:d}#000000',
                     os.path.join(snapshot_name,'atoms#000000'),
                     os.path.join(snapshot_name,'bonds#000000')]
        return snapshot_name, [os.path.join(path, f) for f in filenames]

    # --------------------------------------------------------------------------
    @staticmethod
    def _process_header_atoms(header):
        #for k, v in header.items():
        #    print(f'[{k}] ==> {v}')

        for k in ['time']:
            header[k][0] = float(header[k][0])
        for k in ['loop','nrecord','nfields','nfiles']:
            header[k] = int(header[k])

        assert header['time'][-1] == 'fs'

        h = header['h']
        assert h[-1] == 'Ang'
        for i in range(9):
            h[i] = float(h[i])
        assert np.isclose(h[0], h[4])
        header['h'] = h

        for i, f in enumerate(MACRO_FIELDS_ATOMS):
            if f in ['rx', 'ry', 'rz']:
                assert header['field_units'][MACRO_FIDXS_ATOMS[i]] == 'Ang'

        species = header['species']
        assert np.all(np.in1d(species, MACRO_PROTEIN_NAMES)), f'Found unknown proteins {species}'

        if 'complexes' not in header:
            LOGGER.warning(f'Did not find possbile complex list in file')

        complexes = header.get('complexes', [])
        complexes = [_ for _ in complexes if _ != '<->' and len(_.split('_')) == 2]
        assert np.all(np.in1d(complexes, MACRO_COMPLEX_NAMES)), f'Found unknown complexes {complexes}'

        return header

    @staticmethod
    def _process_header_bonds(header):
        # for k, v in header.items():
        #    print(f'[{k}] ==> {v}')
        for k in ['time']:
            header[k][0] = float(header[k][0])
        for k in ['loop','nrecord','nfields','nfiles']:
            header[k] = int(header[k])

        return header

    @staticmethod
    def _process_header_densities(header):
        #for k, v in header.items():
        #    print(f'[{k}] ==> {v}')

        for k in ['Lx','Ly','delta','time']:
            header[k][0] = float(header[k][0])
        for k in ['nx','ny','nrecord','nLipids','nfields','nfiles','nl']:
            header[k] = int(header[k])

        assert header['nx']*header['ny'] == header['nrecord']
        assert header['Lx'][-1] == 'nm'
        assert header['Ly'][-1] == 'nm'
        assert header['delta'][-1] == 'nm'
        assert header['time'][-1] == 'fs'
        return header

    @staticmethod
    def _process_data_atoms(data, header):

        ids = data[:, 0].astype(np.uint16)
        states = data[:, 1]
        assert np.all(np.in1d(states, MACRO_PROTEIN_NAMES))

        assert np.isclose(header['h'][0], header['h'][4])
        PSIZE = header['h'][0]

        positions = data[:, 2:].astype(np.float)
        if positions[:,:2].min() < 0 or positions[:,:2].max() > PSIZE:
            positions[:,:2] -= (PSIZE * np.floor(positions[:,:2] / PSIZE))

        assert positions[:,:2].min() >= 0 and positions[:,:2].max() <= PSIZE, f'Atom positions out of bounds: {positions}'
        return ids, states, positions

    @staticmethod
    def _process_data_bonds(data, header):
        # for k, v in header.items():
        #    print(f'[{k}] ==> {v}')
        ids = data[:,:2].astype(np.uint16)
        dists = data[:, 2].astype(np.float32)
        return ids, dists

    @staticmethod
    def _process_data_densities(data, header):
        return data

    # --------------------------------------------------------------------------
    @staticmethod
    def _parse_stream(stream, read_fields, read_field_idxs,
                      process_header, process_data, name=''):

        assert isinstance(stream, bytes)
        assert isinstance(read_fields, np.ndarray)
        assert isinstance(read_field_idxs, np.ndarray)
        assert len(read_fields) == len(read_field_idxs)
        assert callable(process_header) and callable(process_data)

        # ----------------------------------------------------------------------
        # extract header
        header = {}

        # go through the lines to extract header
        lines = stream.split(b'\n')
        for lidx, line in enumerate(lines):

            line = line.decode('ascii')
            if line[-1] == '{':
                continue
            if line == '}':
                break

            toks = [_.strip() for _ in re.split('[=;]', line)]
            assert len(toks) == 3 and toks[2] == '', f'unfamiliar header line in stream {name}: {lidx}: ({type(line)}, {len(line)})'

            k, v = toks[0], [_ for _ in toks[1].split(' ') if len(_) > 0]
            header[k] = v[0] if len(v) == 1 else v

        # ----------------------------------------------------------------------
        # validate and process the header
        #LOGGER.debug(header['field_names'])
        fields = np.array(header['field_names'])[read_field_idxs]
        assert np.array_equal(read_fields, fields)

        header['nrecord'] = int(header['nrecord'])
        header = process_header(header)

        # ----------------------------------------------------------------------
        # extract data
        ndata = header['nrecord']
        if ndata == 0:
            return header, None

        # for easy processing, re-join the lines into a single stream
        data = lines[len(header) + 3:]
        data = b'\n'.join(data)

        if header['datatype'] == 'FIXRECORDBINARY':
            data = np.frombuffer(data, dtype=np.float32).reshape(ndata,-1)
            gidx = data[:, 1] * header['nx'] + data[:, 0]
            aidx = np.argsort(gidx)
            data = data[aidx,:]
            data = data[:,read_field_idxs]

        elif header['datatype'] == 'VARRECORDASCII':
            data = np.genfromtxt(io.StringIO(data.decode('ascii')),
                                 usecols=read_field_idxs, dtype=str)

        else:
            assert 0

        # ----------------------------------------------------------------------
        # validate and process the data
        if 1 == len(data.shape):
            data = data[np.newaxis, :]
        assert data.shape == (ndata, len(read_field_idxs))
        data = process_data(data, header)

        # ----------------------------------------------------------------------
        return header, data

    # --------------------------------------------------------------------------
    @staticmethod
    def parse_snapshot_from_directory(snapshot_idx, dirpath=None):
        assert isinstance(snapshot_idx, int)
        assert isinstance(dirpath, str)

        dirpath, filepaths = MacroParser.get_filenames(snapshot_idx, dirpath)
        for filepath in filepaths:
            if not os.path.isfile(filepath):
                LOGGER.debug(f'File ({filepath}) does not exist!')
                return None

        with open(filepaths[0], "rb") as fp:
            densities = fp.read()

        with open(filepaths[1], "rb") as fp:
            atoms = fp.read()

        with open(filepaths[2], "rb") as fp:
            bonds = fp.read()

        return MacroParser.parse_snapshot_from_streams(densities, atoms, bonds)

    # --------------------------------------------------------------------------
    @staticmethod
    def parse_snapshot_from_streams(densities, atoms, bonds, sim_name='simname'):

        assert isinstance(densities, bytes)
        assert isinstance(atoms, bytes)
        assert isinstance(bonds, bytes)

        # ----------------------------------------------------------------------
        #LOGGER.debug('parsing atoms')
        ah,ad = MacroParser._parse_stream(atoms,
                                          MACRO_FIELDS_ATOMS, MACRO_FIDXS_ATOMS,
                                          MacroParser._process_header_atoms,
                                          MacroParser._process_data_atoms, f'{sim_name}-atoms')
        #LOGGER.debug('parsing bonds')
        bh,bd = MacroParser._parse_stream(bonds,
                                          MACRO_FIELDS_BONDS, MACRO_FIDXS_BONDS,
                                          MacroParser._process_header_bonds,
                                          MacroParser._process_data_bonds, f'{sim_name}-bonds')
        #LOGGER.debug('parsing densities')
        dh,dd = MacroParser._parse_stream(densities,
                                          MACRO_FIELDS_DENSITIES, MACRO_FIDXS_DENSITIES,
                                          MacroParser._process_header_densities,
                                          MacroParser._process_data_densities, f'{sim_name}-densities')

        # ----------------------------------------------------------------------
        # extract the data
        atom_hdr, bond_hdr, density_hdr = ah,bh,dh
        if ad is not None:
            prot_ids, prot_states, prot_pos = ad
        else:
            prot_ids, prot_states, prot_pos = None, None, None
        if bd is not None:
            bond_ids, bond_dists = bd
        else:
            bond_ids, bond_dists = None, None
        if dd is not None:
            density_data = dd
        else:
            density_data = None

        # now combine the data
        t,tu = atom_hdr['time'][0], atom_hdr['time'][1]
        lu = density_hdr['Lx'][1]

        assert np.isclose(t, density_hdr['time'][0]) and np.isclose(t, bond_hdr['time'][0])
        assert tu == density_hdr['time'][1] and tu == bond_hdr['time'][1]

        x, y, c = density_hdr['nx'], density_hdr['ny'], len(MACRO_FIELDS_DENSITIES)
        density_data = density_data.reshape((y,x,-1), order='C')
        assert atom_hdr['h'][-1] == 'Ang'

        if prot_pos is not None:
            prot_pos *= 0.1  # Angs to nm

        #show_macro_data(density_data, pplot, density_hdr['Lx'][0], False)

        # print('\tatoms:', prot_ids, prot_states, prot_pos)
        # print('\tbonds:', bond_ids)
        # print('\tdensities:', density_data.shape)

        # ----------------------------------------------------------------------
        return t,tu,lu, density_data, prot_ids, prot_states, prot_pos, bond_ids

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
