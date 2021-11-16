
import numpy as np
from MDAnalysis.lib.util import blocks_of
from MDAnalysis.lib import distances

from mummi_ras.online.cg.baseFastSingleFrame import AnalysisBase

"""
WARNING this function is modified from MDAnalysis

1) during campaign 1, change to x,y only distance + maybe somethign more
2) during campaign 3, now g1 is not all atoms but center of mass of those atoms.

"""

class InterRDF(AnalysisBase):
    """Intermolecular pair distribution function

    InterRDF(g1, g2, nbins=75, range=(0.0, 15.0))

    Arguments
    ---------
    g1 : AtomGroup
      First AtomGroup
    g2 : AtomGroup
      Second AtomGroup
    nbins : int (optional)
          Number of bins in the histogram [75]
    range : tuple or list (optional)
          The size of the RDF [0.0, 15.0]
    exclusion_block : tuple (optional)
          A tuple representing the tile to exclude from the distance
          array. [None]
    start : int (optional)
          The frame to start at (default is first)
    stop : int (optional)
          The frame to end at (default is last)
    step : int (optional)
          The step size through the trajectory in frames (default is
          every frame)

    Example
    -------
    First create the :class:`InterRDF` object, by supplying two
    AtomGroups then use the :meth:`run` method ::

      rdf = InterRDF(ag1, ag2)
      rdf.run()

    Results are available through the :attr:`bins` and :attr:`rdf`
    attributes::

      plt.plot(rdf.bins, rdf.rdf)

    The `exclusion_block` keyword allows the masking of pairs from
    within the same molecule.  For example, if there are 7 of each
    atom in each molecule, the exclusion mask `(7, 7)` can be used.

    .. versionadded:: 0.13.0

    """
    def __init__(self, g1, g2,
                 nbins=75, range=(0.0, 15.0), exclusion_block=None,
                 **kwargs):
        super(InterRDF, self).__init__(g1.universe.trajectory, **kwargs)
        self.g1 = g1
        self.g2 = g2
        self.u = g1.universe

        self.rdf_settings = {'bins': nbins,
                             'range': range}
        self._exclusion_block = exclusion_block

    def _prepare(self):
        # Empty histogram to store the RDF
        count, edges = np.histogram([-1], **self.rdf_settings)
        count = count.astype(np.float64)
        count *= 0.0
        self.count = count
        self.edges = edges
        self.bins = 0.5 * (edges[:-1] + edges[1:])

        # Need to know average volume
        self.volume = 0.0

        # Allocate a results array which we will reuse
        self._result = np.zeros((1, len(self.g2)), dtype=np.float64)
        #self._result = np.zeros((len(self.g1), len(self.g2)), dtype=np.float64)
        # If provided exclusions, create a mask of _result which
        # lets us take these out
        if self._exclusion_block is not None:
            self._exclusion_mask = blocks_of(self._result,
                                             *self._exclusion_block)
            self._maxrange = self.rdf_settings['range'][1] + 1.0
        else:
            self._exclusion_mask = None

    def _single_frame(self):
        distances.distance_array(np.float32(self.g1.center_of_mass()*[1, 1, 0]), np.float32(self.g2.positions*[1, 1, 0]),
                                 box=self.u.dimensions, result=self._result)
        #distances.distance_array(np.float32(self.g1.positions*[1, 1, 0]), np.float32(self.g2.positions*[1, 1, 0]),
        #                         box=self.u.dimensions, result=self._result)
        # Maybe exclude same molecule distances
        if self._exclusion_mask is not None:
            self._exclusion_mask[:] = self._maxrange

        count = np.histogram(self._result, **self.rdf_settings)[0]
        self.count += count

        self.volume += (self.u.dimensions[0]*self.u.dimensions[1])

    def _conclude(self):
        # Number of each selection
        nA = 1
        #nA = len(self.g1)
        nB = len(self.g2)
        N = nA * nB

        # If we had exclusions, take these into account
        if self._exclusion_block:
            xA, xB = self._exclusion_block
            nblocks = nA / xA
            N -= xA * xB * nblocks

        # Volume in each radial shell
        vol = np.power(self.edges[1:], 2) - np.power(self.edges[:-1], 2)
        vol *= np.pi

        # Average number density
        box_vol = self.volume / self.n_frames
        density = N / box_vol

        rdf = self.count / (density * vol * self.n_frames)

        self.rdf = rdf
