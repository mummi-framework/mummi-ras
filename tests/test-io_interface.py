#!/usr/bin/env python3

# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------

import os
import numpy as np
import shutil, logging, sys, time, pickle

import mummi_core, mummi_ras
from mummi_ras.datastructures.patch import PatchConfig, Patch
from mummi_ras.datastructures.rdf import RDF
from mummi_core.utils import timeout, Naming

LOGGER = logging.getLogger(__name__)


# ------------------------------------------------------------------------------
def print_separator():
    print('\n-----------------\n')


# ------------------------------------------------------------------------------
def test_keys(iointerface):
    LOGGER.debug('TEST IO: keys')

    iointerface.save_files('_test_io/dir', 'testkey', 'testdata')
    iointerface.save_files('_test_io/dir', ['testkey2', 'testkey3'], ['testdata2', 'testdata3'])
    LOGGER.debug(iointerface.load_files('_test_io/dir', 'testkey'))
    LOGGER.debug(iointerface.load_files('_test_io/dir', ['testkey2', 'testkey3']))
    LOGGER.debug(f'Good namespace exists: {iointerface.namespace_exists("_test_io/dir")}')
    LOGGER.debug(f'Good key exists: {iointerface.file_exists("_test_io/dir", "testkey")}')
    LOGGER.debug(f'Bad namespace exists: {iointerface.namespace_exists("_test_io/baddir")}')
    LOGGER.debug(f'Bad key exists: {iointerface.file_exists("_test_io/dir", "badkey")}')
    LOGGER.debug(f'Key in bad namespace exists: {iointerface.file_exists("_test_io/baddir", "testkey")}')
    print_separator()


def test_npz(iointerface):
    LOGGER.debug('TEST IO: npz')
    arrays = {'a':np.random.rand(4, 6), 'b':np.random.rand(3, 5)}
    iointerface.save_npz('_test_io/file.tar', 'key', arrays)
    loaded = iointerface.load_npz('_test_io/file.tar', 'key')
    maxVal = 0.0
    for key in arrays:
        difference = np.subtract(arrays[key], loaded[key])
        maxVal = max(maxVal, np.max(difference))
    LOGGER.debug("Maximum difference read: {}".format(maxVal))
    print_separator()


def test_patches(iointerface):
    LOGGER.debug('TEST IO: patches')
    patches = []
    patch_ids = []

    for i in range(3):
        config = PatchConfig('dummy_simname', 1.0, 'dummy_unit', 'dummy_unit', {'grid':'dummy', 'hycop':'dummy'})
        pid = 'dummy ' + str(i)
        p = Patch(
            patch_id=pid,
            patch_config=config, 
            lipid_concentrations=np.ndarray((2,)), 
            patch_center=np.ndarray((2,)),
            protein_ids=np.array([['RAS', 'RAS'], []], dtype=object), 
            protein_states=np.ndarray((2,)), 
            protein_positions=np.ndarray((2,)),
        )
        patches.append(p)
        patch_ids.append(pid)

    iointerface.save_patches('_test_io', patches)
    loaded = iointerface.load_patches('_test_io', patch_ids)
    LOGGER.debug(loaded)
    print_separator()


def test_rdfs(iointerface):
    LOGGER.debug('TEST IO: rdfs')
    rdf = RDF()
    rdf.set_zero()
    rdfkey = 'rdfkey'

    iointerface.save_rdfs('_test_io/test_path', rdfkey, rdf)
    loaded = iointerface.load_rdfs('_test_io/test_path', rdfkey)
    LOGGER.debug(loaded)
    print_separator()


def test_checkpoint(iointerface):
    if iointerface.get_type() != 'simple':
        return
    LOGGER.debug('TEST IO: checkpoints')
    iointerface.save_checkpoint('_test_io/test_checkpoint', {'a':1, 'b':2})
    data = iointerface.load_checkpoint('_test_io/test_checkpoint')
    LOGGER.debug(data)
    iointerface.send_signal('_test_io/', 'test_signal.txt')
    LOGGER.debug(iointerface.test_signal('_test_io/', 'test_signal.txt'))
    iointerface.take_backup('_test_io/test_signal.txt')
    print_separator()


def test_saveload(iointerface):
    LOGGER.debug('TEST IO: files')
    iointerface.save_files('_test_io/test_namespace', 'test_key', 'blahblahblah')
    LOGGER.debug(iointerface.load_files('_test_io/test_namespace', 'test_key'))
    print_separator()


def test_performance(iointerface):
    LOGGER.debug('TEST IO: performance')

    N = 500
    keys = [f'k_{i}' for i in range(N)]
    data = [f'd_{i}' for i in range(N)]

    t0 = time.time()
    iointerface.save_files('_test_io/dir', keys, data)
    iointerface.load_files('_test_io/dir', keys)
    print(f'Time: {time.time() - t0} sec')
    print_separator()


def test_heterogenous(iointerface):
    LOGGER.debug('TEST IO: heterogenous')

    iointerface.save_files('_test_io/dir', ['testkey2', 'testkey3', 'testkey4'], 
        [pickle.dumps(5), 'testdata3', pickle.dumps({'a':5})])
    print(iointerface.load_files('_test_io/dir', ['testkey2', 'testkey3', 'testkey4']))
    print_separator()


def test_macro_keys(iointerface):
    wspace_fb = Naming.dir_root('feedback-cg')
    print(iointerface.get_keys_for_macro_feedback(wspace_fb, Naming.fb_macro_key('worker')))
    print_separator()


def test_snapshots(iointerface, path):
    LOGGER.debug('TEST IO: snapshots')
    NDUMMIES = 3
    frames = list(range(NDUMMIES))
    rawpath = '_test_io/raw_pfpatch'
    outpath = '_test_io/saved_pfpatch'
    dummy_txt = \
"""particle FILEHEADER {type=MULTILINE; datatype=FIXRECORDASCII; checksum=CRC32; create_time=Sun May 23 06:31:19 2021; run_id=0x60a82ea3;
code_version=TARGET rGITVERSION (Compiled on May 20 2021 17:24:28); srcpath=PATH;
loop=82500000; time=1650000000.000000 fs;
nfiles=1; nrecord=136286; lrec=216; nfields=11; endian_key=875770417;
field_names=checksum id class type group rx ry rz vx vy vz;
field_types=u u s s s f f f f f f;
field_units=1 1 1 1 1 Ang Ang Ang Ang/fs Ang/fs Ang/fs;
field_format=%08x  %16.16lx  %s %s %s %21.13e %21.13e %21.13e %21.13e %21.13e %21.13e;
reducedcorner=    -0.50000000000000     -0.50000000000000     -0.50000000000000;
h=   293.55874029120287      0.00000000000000      0.00000000000000
       0.00000000000000    293.55874029120287      0.00000000000000
       0.00000000000000      0.00000000000000    180.09629035409318 Ang;
random = lcg64 ;
randomFieldSize = 27 ;
groups = group free ;
}"""
    io_simple = mummi_core.get_io('simple')
    for frame in frames:
        namespace = os.path.join(rawpath, f"snapshot.{frame:012d}")
        io_simple.save_files(namespace, "subset#000000", dummy_txt)

    iointerface.save_snapshots(frames, rawpath, outpath)
    print(iointerface.list_snapshots(outpath))

    shutil.copy(os.path.join(path, 'lipids-water-eq4.gro'), outpath)
    shutil.copy(os.path.join(path, 'topol.tpr'), outpath)
    extensions = ['name SC1'] * NDUMMIES
    print_separator()


def test_analysis(iointerface, path):
    LOGGER.debug('TEST IO: analysis')
    io_tar = mummi_core.get_io('taridx')
    keys = io_tar.list_keys(path, '*')
    frames = []
    extensions = []
    for key in keys:
        chunks = key.split('.')
        frame, ext = chunks[1].split('_', 1)
        if not frame in frames:
            frames.append(frame)
        if not ext in extensions:
            extensions.append(ext)
    data = iointerface.load_analysis_data(path, frames, extensions)
    print_separator()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if __name__ == '__main__':
    mummi_core.init_logger(level=1)

    Naming.init()
    iointerface = mummi_ras.get_io('mummi')

    test_keys(iointerface)
    test_npz(iointerface)
    test_patches(iointerface)
    test_rdfs(iointerface)
    test_checkpoint(iointerface)
    test_saveload(iointerface)
    test_performance(iointerface)
    test_heterogenous(iointerface)
    # test_snapshots(iointerface, '/p/gpfs1/splash/moon15/')
    # test_analysis(iointerface, '/g/g20/moon15/pilot2/analysis_test/analysis.tar')

    shutil.rmtree('_test_io', ignore_errors=True)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
