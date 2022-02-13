# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
import os
import io
import yaml
import logging
import traceback
import time
import random
import numpy as np

from multiprocessing import Pool
from logging import getLogger
from pathlib import Path
LOGGER = getLogger(__name__)

import mummi_core
from mummi_core.interfaces import get_io
from mummi_core.interfaces.base import IO_Base, check_extn, check_filename
from mummi_core.interfaces.default_functions import write_npz, read_npz
from mummi_core.utils.utilities import sig_ign_and_rename_proc
from mummi_ras.datastructures.patch import Patch
from mummi_ras.datastructures.rdf import RDF

LOGGER = logging.getLogger(__name__)


# To add a new function to IO_Mummi:
# 1. Define the function in IO_Simple, IO_DBR, or IO_Tar
# 2. Function must have argument 'interfaces' with default IO
# 3. Call appropriate abstraction with the function name and arguments
class IO_Mummi (IO_Base):

    io = {
        'simple': mummi_core.get_io('simple'),
        'taridx': mummi_core.get_io('taridx'),
        'redis': mummi_core.get_io('redis'),
    }

    MACRO_FB_MAX_KEYS_TO_READ = 100000
    MACRO_FB_PARALLEL_CHUNK_SIZE = 1000
    MACRO_FB_PARALLEL_NPROCS = 100

    @classmethod
    def get_type(cls):
        return 'mummi'

    @classmethod
    def check_environment(cls):
        return True

    @classmethod
    def list_keys(cls, namespace, keypattern, interfaces=['simple']):
        assert(len(interfaces) == 1), 'Only one interface allowed'
        return cls.io[interfaces[0]].list_keys(namespace, keypattern)

    @classmethod
    def file_exists(cls, namespace, key, interfaces=['taridx']):
        assert(len(interfaces) == 1), 'Only one interface allowed'
        return cls.io[interfaces[0]].file_exists(namespace, key)

    @classmethod
    def namespace_exists(cls, namespace, interfaces=['taridx']):
        assert(len(interfaces) == 1), 'Only one interface allowed'
        return cls.io[interfaces[0]].namespace_exists(namespace)

    @classmethod
    def move_key(cls, namespace, key, prefix='done', suffix='.npz', interfaces=['simple']):
        for interface in interfaces:
            cls.io[interface].move_key(namespace, key, prefix, suffix)

    @classmethod
    def load_files(cls, namespace, filenames, interfaces=['taridx']):
        assert(len(interfaces) == 1), 'Only one interface allowed'
        return cls.io[interfaces[0]].load_files(namespace, filenames)

    @classmethod
    def save_files(cls, namespace, keys, data, interfaces=['taridx']):
        for interface in interfaces:
            cls.io[interface].save_files(namespace, keys, data)

    @classmethod
    def remove_files(cls, namespace, keys, interfaces=['simple']):
        for interface in interfaces:
            cls.io[interface].remove_files(namespace, keys)

    @classmethod
    def load_npz(cls, namespace, keys, reader_func=read_npz, interfaces=['taridx']):
        assert(len(interfaces) == 1), 'Only one interface allowed'
        return cls.io[interfaces[0]].load_npz(namespace, keys, reader_func)

    @classmethod
    def save_npz(cls, namespace, keys, data, writer_func=write_npz, interfaces=['taridx']):
        for interface in interfaces:
            cls.io[interface].save_npz(namespace, keys, data, writer_func)

    # --------------------------------------------------------------------------
    # App specific layer
    # --------------------------------------------------------------------------

    @classmethod
    def load_patches(cls, inpath, patch_ids, filename='patches', interfaces=['taridx']):

        assert(len(interfaces) == 1), 'Only one interface allowed'
        intfc = cls.io[interfaces[0]]

        if intfc.get_type() == 'taridx':
            inpath = check_filename(inpath, filename)

        LOGGER.debug(f'Loading {len(patch_ids)} patches from ({inpath})')

        patch_ids = [check_extn(p, '.npz') for p in patch_ids]
        patches = intfc.load_npz(inpath, patch_ids, Patch.load_npz)

        LOGGER.info('Loaded {} patches from ({})'.format(len(patches), inpath))
        return patches

    @classmethod
    def save_patches(cls, outpath, patches, filename='patches', interfaces=['taridx']):

        for interface in interfaces:
            intfc = cls.io[interface]

            if intfc.get_type() == 'taridx':
                outpath = check_filename(outpath, filename)

            LOGGER.debug(f'Saving {len(patches)} patches to ({outpath})')

            patch_ids = [check_extn(p.id, '.npz') for p in patches]
            intfc.save_npz(outpath, patch_ids, patches, Patch.save_npz)

            LOGGER.info(f'Saved {len(patches)} patches to ({outpath})')

    @classmethod
    def load_rdfs(cls, inpath, rdfkey, interfaces=['simple']):

        assert isinstance(inpath, str) and isinstance(rdfkey, str)

        assert(len(interfaces) == 1), 'Only one interface allowed'
        intfc = cls.io[interfaces[0]]

        LOGGER.debug(f'Loading RDF ({rdfkey}) from ({inpath})')
        rdf = intfc.load_npz(inpath, rdfkey, RDF.load_npz)
        LOGGER.debug(f'Loading RDF ({rdfkey}) from ({inpath})')
        return rdf

    @classmethod
    def save_rdfs(cls, outpath, rdfkey, rdf, interfaces=['simple']):

        assert isinstance(outpath, str) and isinstance(rdfkey, str) and isinstance(rdf, RDF)

        for interface in interfaces:
            intfc = cls.io[interface]
            if intfc.save_npz(outpath, rdfkey, rdf, RDF.save_npz):
                LOGGER.info(f'Saved RDF ({rdfkey}) to ({outpath})')
            else:
                LOGGER.info(f'Failed to save RDF ({rdfkey}) to ({outpath})')

    # try to match pattern of load_patches
    # move this to mummi-ras
    @classmethod
    def load_snapshots(cls, inpath, frames, extensions, filename='positions', interfaces=['taridx']):

        # 07/01/21 waiting for Xiouhua's function
        raise Exception('Not implemented yet')

        # from mummi_ras.parsers.cg import CGParser

        # fname_gro = 'lipids-water-eq4.gro'
        # fname_tpr = 'topol.tpr'
        # data = np.empty( (len(extensions), len(frames)), dtype=object)
        # cg = CGParser(os.path.join(inpath,fname_gro), os.path.join(inpath,fname_tpr), empty=False)

        # for cFrameIndex in range(len(frames)):
        #     cFrame = frames[cFrameIndex]
        #     for cKeyIndex in range(len(extensions)):
        #         cKey = extensions[cKeyIndex]
        #         fname = "pos.{:012d}".format(cFrame)

        #         try:
        #             buff = cls.load_files(inpath, fname)
        #             with open(os.path.join(inpath,'pos#000000'),'wb') as fp:
        #                 fp.write(buff)
        #             cg.update_frame(os.path.join(inpath,'pos#000000'))
        #             x = cg.syst.select_atoms(cKey)
        #             data[cKeyIndex, cFrameIndex] = x.positions.copy()
        #         except KeyError:
        #             LOGGER.warning("Load analysis - subfile {} not found".format(fname))
        #             data[cKeyIndex, cFrameIndex] = 'empty'

        # LOGGER.info('Load snapshots data from {} for frames {} and extensions {}'.format(path, frames, extensions))
        # return data

    @classmethod
    def save_snapshots(cls, keys, inpath, outpath, filename='positions', interfaces=['taridx']):

        # for API purposes, we want to faciliate simple as well
        # but the function will need more work and we will probably never use it
        assert len(interfaces) == 1 and interfaces[0] == 'taridx'
        intfc = cls.io[interfaces[0]]

        outpath = check_filename(outpath, filename)

        if len(keys) == 0:
            LOGGER.warning('Nothing to do!')
            return

        # read in from simple, and write to tar
        LOGGER.info(f'Saving snapshots {keys} '
                    f'from ({inpath})[simple] to ({outpath})[{interfaces[0]}]')

        fnames = [f'snapshot.{key:012d}/subset#000000' for key in keys]
        outnames = [f'pos.{key:012d}' for key in keys]

        data = cls.io['simple'].load_files(inpath, fnames)

        intfc.save_files(outpath, outnames, data)
        LOGGER.info(f'Saved all snapshots for frames {keys}')

    @classmethod
    def list_snapshots(cls, outpath, filename='positions', interfaces=['taridx']):

        assert len(interfaces) == 1 and interfaces[0] == 'taridx'
        intfc = cls.io[interfaces[0]]

        outpath = check_filename(outpath, filename)

        data = intfc.load_index(outpath, idx_start=0, idx_end=-1)[:, 0]
        return [p for p in data if Path(p).match('pos.????????????')]

    @classmethod
    def save_analysis_data(cls, cOutputList, outpath, filename='analysis', interfaces=['taridx']):

        # for API purposes, we want to faciliate simple as well
        # but the function will need more work and we will probably never use it
        assert len(interfaces) == 1 and interfaces[0] == 'taridx'
        intfc = cls.io[interfaces[0]]

        outpath = check_filename(outpath, filename)

        # @TODO make this handle GPFS errors
        for cDic in cOutputList:
            cFrame = cDic['frameNum']

            # For each output set write all keys except frameNum
            for cKey, cVal in cDic.items():
                if cKey != 'frameNum':
                    cFileName = f'analysis.{cFrame:012d}_{cKey}.npz'
                    intfc.save_npz(outpath, cFileName, {cKey: cVal})

        clist = list([x['frameNum'] for x in cOutputList[:]])
        LOGGER.info(f'Saved all analysis data for frames {clist}')

    @classmethod
    def load_analysis_data(cls, inpath, frames, extensions, filename='analysis', interfaces=['taridx']):
        assert len(interfaces) == 1 and interfaces[0] == 'taridx'

        inpath = check_filename(inpath, filename)
        
        cFileNames = [None] * len(frames) * len(extensions)
        for cFrameIndex, cFrame in enumerate(frames):
            for cKeyIndex, cKey in enumerate(extensions):
                idx = cFrameIndex * len(extensions) + cKeyIndex
                cFileNames[idx] = "analysis.{:012d}_{}.npz".format(int(cFrame), cKey)
        
        if None in cFileNames:
            LOGGER.debug('Failed to load all analysis data!')

        data = np.empty( (len(extensions), len(frames)), dtype=object)
        data_array = cls.io[interfaces[0]].load_files(inpath, cFileNames)
        for cFrameIndex, cFrame in enumerate(frames):
            for cKeyIndex, cKey in enumerate(extensions):
                idx = cFrameIndex * len(extensions) + cKeyIndex
                data[cKeyIndex, cFrameIndex] = data_array[idx]
        
        LOGGER.info(f'Loaded analysis data from {inpath} for {len(frames)} frames and {len(extensions)} extensions {extensions}')
        return data

    # 8/3/21: archiving old version of this function
    # @classmethod
    # def load_analysis_data(cls, inpath, frames, extensions, filename='analysis', interfaces=['taridx']):
        
    #     assert len(interfaces) == 1 and interfaces[0] == 'taridx'
    #     intfc = cls.io[interfaces[0]]

    #     inpath = check_filename(inpath, filename)

    #     data = np.empty( (len(extensions), len(frames)), dtype=object)
    #     for cFrameIndex in range(len(frames)):
    #         cFrame = frames[cFrameIndex]
    #         for cKeyIndex in range(len(extensions)):
    #             cKey = extensions[cKeyIndex]
    #             cFileName = "analysis.{:012d}_{}.npz".format(cFrame, cKey)
    #             try:
    #                 data[cKeyIndex, cFrameIndex] = intfc.load_npz(inpath, cFileName)
    #             except KeyError:
    #                 LOGGER.warning("Load analysis - subfile {} not found".format(cFileName))
    #                 data[cKeyIndex, cFrameIndex] = 'empty'

    #     LOGGER.info('Load analysis data from {} for frames {} and extensions {}'.format(path, frames, extensions))
    #     return data

    @classmethod
    def save_pdb_data_aa(cls, frame_idx, selection, path, interfaces=['taridx']):

        import MDAnalysis as mda

        # for API purposes, we want to faciliate simple as well
        # but the function will need more work and we will probably never use it
        assert len(interfaces) == 1 and interfaces[0] == 'taridx'
        intfc = cls.io[interfaces[0]]

        namespace = os.path.join(path, "pdb.tar")
        LOGGER.info(f'saving pdb data for frame {frame_idx} to ({namespace})')

        # TODO: CS, requires import MDAnalysis; may want to pass this in from AA
        memfile = io.StringIO()
        f = mda.lib.util.NamedStream(memfile, "output.pdb")  # pdb name is placeholder to get extension
        #memfile.close()
        selection.write(f)
        #memfile.close()  # HB moved this close after the write
        b = ''.join(f.readlines())
        memfile.close()  # CS moved this close to after the readlines

        LOGGER.debug(f'created string stream of len ({len(b)})')
        intfc.save_npz(namespace, f"protein.{frame_idx:012d}.pdb", b.encode('utf-8'))

    @classmethod
    def save_rmsfs_all(cls, inpath, keys, rmsfs, coords, interfaces=['taridx']):

        assert isinstance(inpath, str)
        assert isinstance(keys, list) and isinstance(rmsfs, list) and isinstance(coords, list)
        assert len(keys) == len(rmsfs) and len(keys) == len(coords)

        for interface in interfaces:
            intfc = cls.io[interface]
            data = [{'rmsf': rmsfs[i], 'coord': coords[i]} for i in range(len(keys))]
            intfc.save_npz(inpath, keys, data)

    @classmethod
    def load_cg(cls, patch_counter, max_patches, patch_frame, path, interfaces=['taridx']):
        
        assert(len(interfaces) == 1), 'Only one interface allowed'
        intfc = cls.io[interfaces[0]]

        cg_files = {}
        files = intfc.list_keys(path, '*.npz')
        for file in files:
            try:
                idx = re.findall(r'\d{12}', filename)[0]
                idx = int(idx)
                assert (idx not in cg_files.keys()), "CG data has a duplicate patch index!"
                cg_files[idx] = os.path.basename(file)
            except IndexError:
                continue

        last_idx = 0
        idxs = list(cg_files.keys())
        idxs.sort()
        cg_data = {}
        for idx in idxs:
            if len(cg_data) >= max_patches:
                break
            if idx >= patch_counter:
                npz_file = intfc.load_npz(path, cg_files[idx])
                cg = npz_file['ras_alldist']
                if patch_frame >= cg.shape[0]:
                    continue
                else:
                    cg = cg[[patch_frame],:]
                pid = "{:06d}_{:012d}".format(patch_frame, idx) # ensure pid for each frame in each patch is unique
                cg_data[pid] = cg
                last_idx = idx

        LOGGER.info('Loaded {} CG patches from {}, idxs {} to {}'.format(len(cg_data), path, patch_counter, last_idx))
        return cg_data, last_idx+1

    # --------------------------------------------------------------------------
    # functions for macro feedback
    # --------------------------------------------------------------------------
    @classmethod
    def get_keys_for_macro_feedback(cls, namespace, keypattern, interfaces=['redis']):
        assert(len(interfaces) == 1), 'Only one interface allowed'
        intf = interfaces[0]

        t0 = time.time()

        if intf == 'redis':
            servers_to_keys = cls.io[intf].list_servers_to_keys(namespace, keypattern)

            # inverse dict
            key_to_server = {}
            for server in servers_to_keys:
                for key in servers_to_keys[server]:
                    key_to_server[key] = server
            keys = list(key_to_server.keys())

            # cap with max_keys
            # shuffle so we don't exclude any servers
            if len(keys) > cls.MACRO_FB_MAX_KEYS_TO_READ:
                random.shuffle(keys)
                keys = keys[:cls.MACRO_FB_MAX_KEYS_TO_READ]

            # re-inverse dict
            servers_to_keys = {}
            for key in keys:
                server = key_to_server[key]
                servers_to_keys.setdefault(server, []).append(key)

        else:
            keys = cls.io[intf].list_keys(namespace, keypattern)
            keys = keys[:cls.MACRO_FB_MAX_KEYS_TO_READ]
            servers_to_keys = {'local': keys}

        LOGGER.info(f'Found {len(keys)} keys using ({intf}) interface. '
                    f'Took {time.time() - t0:.5f} sec')
        return servers_to_keys

    @classmethod
    def load_rdfs_for_macro_feedback(cls, namespace, servers_to_keys, interfaces=['redis']):
        assert(len(interfaces) == 1), 'Only one interface allowed'
        intf = interfaces[0]

        t0 = time.time()

        keys = list(set([k for d in servers_to_keys.values() for k in d]))
        nkeys = len(keys)
        rdfs = [None] * nkeys

        if intf == 'redis':
            # create pool arguments for loading
            pool_args = [[namespace, server_keys[k:k + cls.MACRO_FB_PARALLEL_CHUNK_SIZE], server, RDF.load_npz]
                         for server, server_keys in servers_to_keys.items()
                         for k in range(0, len(keys), cls.MACRO_FB_PARALLEL_CHUNK_SIZE)
                         ]
            pool = Pool(processes=cls.MACRO_FB_PARALLEL_NPROCS, initializer=sig_ign_and_rename_proc, initargs=('pool_load_rdfs_at_server',))
            key_to_rdfs_list = pool.starmap(cls.io[intf].load_npz_at_server, pool_args)
            pool.close()
            pool.join()

            # iterate keys because we need to preserve order in rdfs list
            for i, key in enumerate(keys):
                for key_to_rdfs in key_to_rdfs_list:
                    if key in key_to_rdfs:
                        rdfs[i] = key_to_rdfs[key].data
                        break

        else:
            # create pool arguments for loading
            pool_args = [[namespace, keys[k:k + cls.MACRO_FB_PARALLEL_CHUNK_SIZE], RDF.load_npz]
                         for k in range(0, len(keys), cls.MACRO_FB_PARALLEL_CHUNK_SIZE)
                         ]
            pool = Pool(processes=cls.MACRO_FB_PARALLEL_NPROCS, initializer=sig_ign_and_rename_proc, initargs=('pool_load_rdfs_macro',))
            rdfs_list = pool.starmap(cls.io[intf].load_npz, pool_args)
            rdfs_list = [_ for chunk in rdfs_list for _ in chunk]
            pool.close()
            pool.join()
            for i, rdf in enumerate(rdfs_list):
                if rdf is None:
                    LOGGER.warning("Found null RDF -- likely file is '%s'", keys[i])
                    rdfs[i] = None
                else:
                    rdfs[i] = rdf.data

            # 2021.09.30: HB removed this again. instead, push a None for the corrupt file!
            #rdfs = [rdf.data for rdf in rdfs_list if rdf is not None]

        LOGGER.info(f'Fetched {len(rdfs)} rdfs using ({intf}) interface. '
                    f'Took {time.time() - t0:.5f} sec')
        return rdfs

    @classmethod
    def remove_keys_for_macro_feedback(cls, namespace, servers_to_keys, interfaces=['redis']):
        assert(len(interfaces) == 1), 'Only one interface allowed'
        intf = interfaces[0]

        keys = list(set([k for d in servers_to_keys.values() for k in d]))
        nkeys = len(keys)

        archive_namespace = os.path.join(namespace, 'processed_rdfs')

        if intf == 'redis':
            t0 = time.time()
            pool_args = []
            for server in servers_to_keys:
                _keys = servers_to_keys[server]
                for chunk_start in range(0, len(_keys), cls.MACRO_FB_PARALLEL_CHUNK_SIZE):
                    keys_chunk = _keys[chunk_start:chunk_start + cls.MACRO_FB_PARALLEL_CHUNK_SIZE]
                    keys_chunk = [key for key in keys_chunk if key in keys]
                    pool_args.append([namespace, archive_namespace, keys_chunk, server])
            pool = Pool(processes=cls.MACRO_FB_PARALLEL_NPROCS, initializer=sig_ign_and_rename_proc, initargs=('pool_remove_keys_macro_redis',))
            pool.starmap(cls.io[intf].rename_files_at_server, pool_args)
            pool.close()
            pool.join()
            LOGGER.info(f'Archived {nkeys} keys using {intf} interface. '
                        f'Took {time.time() - t0} sec.')

        else:
            ts = time.strftime("%Y%m%d-%H%M%S")
            archive_namespace = f'{archive_namespace}.{ts}'

            # Load files from interface
            t0 = time.time()
            pool_args_load = [[namespace, keys[k:k + cls.MACRO_FB_PARALLEL_CHUNK_SIZE]]
                              for k in range(0, len(keys), cls.MACRO_FB_PARALLEL_CHUNK_SIZE)
                              ]
            pool_load = Pool(processes=cls.MACRO_FB_PARALLEL_NPROCS, initializer=sig_ign_and_rename_proc, initargs=('pool_remove_keys_macro',))
            rdfs = pool_load.starmap(cls.io[intf].load_files, pool_args_load)
            rdfs = [_ for chunk in rdfs for _ in chunk]
            pool_load.close()
            pool_load.join()

            LOGGER.info(f'Loaded {len(rdfs)} rdfs using {intf} interface. '
                        f'Took  {time.time() - t0} sec.')


            # Save files in taridx archive. 
            t0 = time.time()
            if True:
                for k in range(0, len(keys), cls.MACRO_FB_PARALLEL_CHUNK_SIZE):
                    cls.io['taridx'].save_files(archive_namespace, keys[k:k + cls.MACRO_FB_PARALLEL_CHUNK_SIZE], rdfs[k:k + cls.MACRO_FB_PARALLEL_CHUNK_SIZE])
            else:
                # JM 7/16/21: Don't use parallel pool because pytaridx will crash.
                # mummi_core.interfaces.tar:_save_files:130 - ERROR - Failed to save files: Master block of file /g/g20/moon15/pilot2/mummi_root/feedback-cg2macro/processed_rdfs.20210716-103124.tar.pytree_ invalid!
                # mummi_core.interfaces.tar:_save_files:130 - ERROR - Failed to save files: [Errno 2] No such file or directory: '/g/g20/moon15/pilot2/mummi_root/feedback-cg2macro/processed_rdfs.20210716-103124.tar.pylst_' -> '/g/g20/moon15/pilot2/mummi_root/feedback-cg2macro/processed_rdfs.20210716-103124.tar.pylst'
                pool_args_save = [[archive_namespace, keys[k:k + cls.MACRO_FB_PARALLEL_CHUNK_SIZE], rdfs[k:k + cls.MACRO_FB_PARALLEL_CHUNK_SIZE]]
                                  for k in range(0, len(keys), cls.MACRO_FB_PARALLEL_CHUNK_SIZE)
                                  ]
                pool_save = Pool(processes=cls.MACRO_FB_PARALLEL_NPROCS, initializer=sig_ign_and_rename_proc, initargs=('pool_remove_keys_pytaridx',))
                pool_save.starmap(cls.io['taridx'].save_files, pool_args_save)
                pool_save.close()
                pool_save.join()

            LOGGER.info(f'Archived {len(rdfs)} rdfs using taridx interface. '
                        f'Took  {time.time() - t0} sec.')

            # Delete files
            t0 = time.time()
            pool_args_remove = [[namespace, keys[k:k + cls.MACRO_FB_PARALLEL_CHUNK_SIZE]]
                                for k in range(0, len(keys), cls.MACRO_FB_PARALLEL_CHUNK_SIZE)
                                ]
            pool_remove = Pool(processes=cls.MACRO_FB_PARALLEL_NPROCS, initializer=sig_ign_and_rename_proc, initargs=('pool_remove_keys',))
            pool_remove.starmap(cls.io[intf].remove_files, pool_args_remove)
            pool_remove.close()
            pool_remove.join()
            LOGGER.info(f'Deleted {len(rdfs)} rdfs using {intf} interface. '
                        f'Took  {time.time() - t0} sec.')

    # --------------------------------------------------------------------------
    # Unused functions
    # --------------------------------------------------------------------------
    # @classmethod
    # def load_itp(cls, inpath, itpkey):
    #     raise Exception('Empty function!')

    #     assert isinstance(inpath, str) and isinstance(itpkey, str)
    #     filename = os.path.join(inpath, itpkey)
    #     LOGGER.debug('Loaded ITP from ({})'.format(filename))
    #     return filename

    # @classmethod
    # def save_itp(cls, outpath, itpkey, itpid):
    #     raise Exception('Empty function!')

    #     assert isinstance(outpath, str) and isinstance(itpkey, str)
    #     itpkey = itpkey[:-4] + '_' + str(itpid) + itpkey[-4:]
    #     filename = os.path.join(outpath, itpkey)
    #     LOGGER.info('Saved ITP ({}) to ({})'.format(itpkey, outpath))
    #     return filename

    # @classmethod
    # def load_pdb(cls, inpath, pdbkey):
    #     raise Exception('Empty function!')

    #     assert isinstance(inpath, str) and isinstance(pdbkey, str)
    #     filename = os.path.join(inpath, pdbkey)
    #     LOGGER.debug('Loading PDB from ({})'.format(filename))
    #     return filename

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
