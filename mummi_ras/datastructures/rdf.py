# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
"""RDF module."""

import io
import os
import numpy as np


from logging import getLogger
LOGGER = getLogger(__name__)

# ------------------------------------------------------------------------------
# RDF settings

range_max = 120.0    # in Angstroms, now was 70 in campain 1
#range_max = 70.0
range_min = 0.1
bin_size = 0.2
nbins = int((range_max-range_min)/bin_size)

count, edges = np.histogram([-1], bins=nbins, range=(range_min, range_max))


LIPID_TYPES = ['POPC', 'PAPC', 'POPE', 'DIPE', 'DPSM', 'PAPS', 'PAP6', 'CHOL']
#STATE_NAMES = ['RASa', 'RASb', 'RASc', 'zRASa', 'mRASa', 'mRASb', 'zRAFa', 'mRAFa', 'mRAFb']
STATE_NAMES = ['RASa', 'RASb', 'RASc', 'zRASa', 'mRASa', 'mRASb', 'zRAFa', 'mRAFa', 'mRAFb',
               'RAS4Aa', 'RAS4Ab', 'RAS4Ac', 'zRAS4Aa', 'mRAS4Aa', 'mRAS4Ab', 'zRAF4Aa', 'mRAF4Aa', 'mRAF4Ab']

nstates = len(STATE_NAMES)
nlipids = len(LIPID_TYPES)


# ------------------------------------------------------------------------------
# RDF module
# ------------------------------------------------------------------------------
class RDF:
    """Utility class to store an RDF"""

    DTYPE = np.int64
    SHAPE = (nstates, nlipids, nbins)

    RNG_MAX = range_max
    RNG_MIN = range_min
    SZ_BIN = bin_size
    NBINS = nbins
    BINS = 0.5 * (edges[:-1] + edges[1:])

    def __init__(self):
        self.naggregate = 0
        self.data = None

    def sum(self):
        return self.data.sum()

    def set_zero(self):
        self.naggregate = 0
        self.data = np.zeros(RDF.SHAPE, dtype=RDF.DTYPE)

    def set(self, data, naggregate=0):
        if data.shape != RDF.SHAPE:
            s = 'Trying to create invalid RDF. expect {} array!'
            raise ValueError(s.format(RDF.SHAPE))

        self.data = data.astype(RDF.DTYPE)
        self.naggregate = naggregate

    def aggregate(self, rdf):
        self.naggregate = self.naggregate+1

        if isinstance(rdf, RDF):
            np.add(self.data, rdf.data, self.data)
        elif isinstance(rdf, np.ndarray):
            np.add(self.data, rdf, self.data)
        else:
            assert False

    def scale(self, weight):
        self.data = (weight * self.data).astype(RDF.DTYPE)

    def scale_and_aggregate_many(self, rdfs, weights=None):
        assert isinstance(rdfs, np.ndarray)
        assert len(rdfs.shape) == 4
        assert rdfs.shape[1:] == RDF.SHAPE
        if weights is not None:
            assert isinstance(weights, np.ndarray)
            assert rdfs.shape[0] == weights.shape[0]
            rdfs = np.multiply(rdfs.astype(np.float64),
                               weights[:, None, None, None].astype(np.float64))
            rdfs = rdfs.astype(RDF.DTYPE)

        np.add(self.data, rdfs.sum(axis=0), self.data)
        self.naggregate = self.naggregate+rdfs.shape[0]

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # HB copied this from mummi_ras/online/cg/counts2rdf.py
    @staticmethod
    def counts2rdf(numpyRDF):

        assert isinstance(numpyRDF, np.ndarray)
        assert numpyRDF.shape == (8, RDF.NBINS)

        ncols, nbins = numpyRDF.shape
        actualbin = (RDF.RNG_MAX -  RDF.RNG_MIN) / RDF.NBINS

        density = np.zeros(ncols)
        for i in range(ncols):
            summing_bins = 10
            density[i] = (np.sum(numpyRDF[i, -(summing_bins+1):-1])) /\
                (((actualbin * RDF.NBINS)**2 - (actualbin * (RDF.NBINS-summing_bins))**2)*np.pi)

        RDFcurves = np.zeros(numpyRDF.shape)
        for i in range(ncols):
            for j in range(nbins):
                RDFcurves[i][j] = numpyRDF[i][j] /\
                    (density[i]*(((actualbin * (j+2))**2 - (actualbin * (j+1))**2)*np.pi))

        return RDFcurves

    # --------------------------------------------------------------------------
    @staticmethod
    def wnpz(fname, dh):
        np.savez_compressed(fname, data=dh[0], header=dh[1])

    @staticmethod
    def write_rdfs_individually(rdf, outpath, iointerface):

        assert isinstance(rdf, np.ndarray)
        assert rdf.shape == RDF.SHAPE
        assert iointerface.get_type() == 'simple'

        LOGGER.debug(f'Writing RDFs individually to ({outpath})')
        zeros = np.zeros(rdf.shape[1], dtype=rdf.dtype)
        bins = np.arange(0.0, 0.1*range_max, 0.1*bin_size)
        bins = np.expand_dims(bins, 1)

        header = ['bin'] + LIPID_TYPES
        fmt = ['%.2f'] + ['%12d' for l in LIPID_TYPES]

        hdr = ' '.join(header)
        for i in range(nstates):
            outfname = os.path.join(outpath, f'rdf_{STATE_NAMES[i]}')
            d = np.hstack((bins, np.vstack((zeros, rdf[i].T))))

            np.savetxt(outfname+'.txt', d, header=hdr, fmt=fmt)
            iointerface.save_npz(outpath, outfname+'.npz', {'data': d, 'header': hdr})

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # save and load rdfs as numpy archives (npz)
    @staticmethod
    def save_npz(fptr, rdf):
        assert isinstance(fptr, (str, io.BytesIO))
        assert isinstance(rdf, RDF)
        np.savez_compressed(fptr, naggregates = rdf.naggregate, data = rdf.data)
        return fptr

    @staticmethod
    def load_npz(fptr):


        assert isinstance(fptr, (str, io.BytesIO))
        try:
            data = np.load(fptr, allow_pickle=True)
        except Exception as _excp:
            LOGGER.error('Failed to load RDF from ({}): {}'.format(fptr, _excp))
            return

        rdf = RDF()
        rdf.naggregate = data['naggregates']
        rdf.data = data['data']

        data.close()
        return rdf

    # --------------------------------------------------------------------------
    # RW of text files for RDF data
    @staticmethod
    def parse(lines):
        """Parse the lines read from a *single* file into a numpy array."""

        # i may get an empty line at the end!
        if len(lines) == RDF.SHAPE[2]+1 and len(lines[-1]) == 0:
            lines = lines[:-1]

        if len(lines) != RDF.SHAPE[2]:
            raise ValueError('RDF.parse expects {} lines, but got {}'
                             .format(RDF.SHAPE[2], len(lines)))

        rdf = np.genfromtxt(lines, comments='%', delimiter='')

        # match the bins!
        if not np.allclose(rdf[:, 0], RDF.BINS, atol=1e-03):
            raise ValueError('RDF.parse found mismatching bin')

        rdf = rdf[:, 1:]
        rdf = np.transpose(rdf.astype(np.float).astype(RDF.DTYPE))
        if rdf.shape != (RDF.SHAPE[1], RDF.SHAPE[2]):
            s = "RDF.parse expected to create {} array from {} lines, "
            s += "but created {} array"
            raise ValueError(s.format(RDF.SHAPE[1:], len(lines), rdf.shape))
        return rdf

    @staticmethod
    def format(rdf):
        """Format the rdf into a set of lines into a *single* string."""

        if rdf.shape != (RDF.SHAPE[1], RDF.SHAPE[2]):
            s = "RDF.format expected a {} array, but got {}"
            raise ValueError(s.format(RDF.SHAPE[1:], rdf.shape))

        line_width = 8*RDF.SHAPE[1] + (RDF.SHAPE[1]-1) + 2

        # create as many lines!
        text = ''
        for j in range(rdf.shape[1]):
            s = np.array2string(rdf[:, j],
                                max_line_width=line_width,
                                formatter={'int_kind': lambda x: "%d" % x})
            text = text + ('%.3f' % RDF.BINS[j]) + ' ' + s[1:-1] + '\n'
        return text

    @staticmethod
    def load(filenames):
        """Read and parse several files to create rdf."""

        if len(filenames) != RDF.SHAPE[0]:
            s = 'RDF.load expects {} rdf files, one each for a state'
            raise ValueError(s.format(RDF.SHAPE[0]))

        rdfs = []
        for f in filenames:
            lines = RDF.read(f)
            rdfs.append(RDF.parse(lines))

        rdfs = np.asarray(rdfs)
        if rdfs.shape != RDF.SHAPE:
            s = 'RDF.load loaded incorrect shape {}'
            raise ValueError(s.format(rdfs.shape))

        return rdfs

    @staticmethod
    def save(rdfs, filenames):
        # format and save rdfs into files

        if len(filenames) != RDF.SHAPE[0] or len(filenames) != rdfs.shape[0]:
            s = 'RDF.save expects {} rdf files, one each for a state, '
            s += 'but found {} rdfs and {} files'
            raise ValueError(s.format(RDF.SHAPE[0], rdfs.shape[0],
                                      len(filenames)))

        for i in range(len(filenames)):
            txt = RDF.format(rdfs[i])
            RDF.write(filenames[i], txt)

    @staticmethod
    def read(filename):
        text_file = open(filename, 'r')
        lines = text_file.read()
        text_file.close()
        return lines

    @staticmethod
    def write(filename, txt):
        text_file = open(filename, 'w')
        text_file.write(txt)
        text_file.close()

# ------------------------------------------------------------------------------
