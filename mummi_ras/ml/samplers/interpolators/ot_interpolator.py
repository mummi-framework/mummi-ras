"""
Optimal transport based interpolator.


NOTES:

to find the LS points to initialize the 
"""

import ot
import itertools
import numpy as np
from typing import List
from mummi_ras.ml.ls_point import LSPoint
from logging import getLogger

LOGGER = getLogger(__name__)

default = {
                "k" : 10,            # number of neighbors (used to calculate the k nearest neighbors)
               "ep" : 0.0001,        # perturbation value
               "lb" : 0.0,           # lower bound on lambda 
               "ub" : 1.0,           # upper bound on lambda 
      "int_samples" : None,           # # of internal samples
 "int_samples_mult" : 10.0,           # # of internal samples
      "edge_weight" : 'edge_size',      # selection frequency of the edge
           "method" : 'kdist_min_rank'
            }

class OTInterpolator:
    """
    Interpolator using an optimal transport solution.
    """

    def __init__(self,
                 states: List[np.ndarray],
                 kneigh=10,
                 lowerbound=0.4,
                 upperbound=0.6,
                 int_samples=None,
                 edge_weight='edge_size',
                 method='kdist_min_rank',
                 num_iter_max=1e6,
                 safe_mode=True):
        """Constructor for this interpolator"""

        # concatenate arrays and make list of index mappings
        assert len(states) >= 2, "Number of states must be equal or larger to 2."

        # make one large state 
        self.states = states  # list of LSPoint lists
        _all_states = []
        LOGGER.info(f"creating all states")

        for st in self.states:
            _all_states += st
        self.all_points = self._get_lsp_coordinate_array(_all_states)  # list of lsp coordinate ndarray

        self.index_offsets = []
        _offset = 0
        for si in self.states:
            # add offset to list
            self.index_offsets.append(_offset)
            # move counter
            _offset += len(si)

        # That loop can be long so we do not check if we are sure
        if safe_mode:
            LOGGER.info(f"Doing assert")
            for _pairs in itertools.product(states, states):
                assert len(_pairs[0]) == len(_pairs[1]), "All sets of LS states need to be the same shape."

        self.kneigh = kneigh
        self.lowerbound = lowerbound
        self.upperbound = upperbound
        self.int_samples = int_samples
        self.edge_weight = edge_weight
        self.method = method
        self.num_iter_max = num_iter_max
        self.edges = {}  # dictionary holding the OT edges 

        self.ot_schedule = self.compute_optimal_transport()

        return

    def compute_optimal_transport(self):
        """Solve optimal transport schedule for arbitrary list of source and targets"""

        # reset edges list
        self.edges = []

        for pairs in itertools.product(range(len(self.states)), range(len(self.states))):
            i,j = pairs

            # avoid self and double counting
            if i >= j:
                continue

            i_offset = self.index_offsets[i]
            j_offset = self.index_offsets[j]

            # compute pairwise schedule, will return a dictionary of i to j interal index mappings
            _edges = self.compute_optimal_transport_pairwise(self._get_lsp_coordinate_array(self.states[i]),
                                                             self._get_lsp_coordinate_array(self.states[j]))

            # add new edges in, offset indexes by the number of indexes already in. 
            _edge_index_global = [ (k[0] + i_offset, k[1] + j_offset) for k in _edges ]

            self.edges.extend(_edge_index_global)

        return

    def _get_lsp_coordinate_array(self, points):
        """Return a single ndarray from a list of lspoints"""

        _array =np.array([lsp.get_coordinates() for lsp in points])

        return _array

    def compute_optimal_transport_pairwise(self, source, target):
        """
        Solve optimal transport solution using POT package.

        expects a nparray of source and target points.

        this function solves the OT problem using the POT package
        between two sets A and B and stores the solution in
        the object dto.

        """
        LOGGER.info(f"compute_optimal_transport_pairwise with source={source.shape} / target={target.shape}")

        n    = source.shape[0]
        a, b = np.ones((n,)) / n, np.ones((n,)) / n  # uniform distribution on samples
        M    = ot.dist(source, target)
        LOGGER.info(f"ot.emd with a={a.shape} / b={b.shape} and M={M.shape}")
        G0 = ot.emd(a, b, M, numItermax = self.num_iter_max)

        #-- calculates the solution (i.e. the edge with the highest weight)
        edges   = [ (i, np.argmax(G0[i])) for i in range(G0.shape[0]) ]

        return edges

    def draw_k_points(self, npoints):

        # arcs is index in edges dictionary
        arcs = np.random.randint(len(self.edges), size=npoints)
        lamb = (self.upperbound - self.lowerbound) * np.random.random(size=npoints) + self.lowerbound

        _edges = [ self.edges[arc] for arc in arcs ]

        _npts = [self.all_points[ed[0]] + l * (self.all_points[ed[1]] - self.all_points[ed[0]])  for ed,l in zip(_edges, lamb)]

        return _npts

    def size(self) -> int:
        """ Return the size of the interpolator"""
        return len(self.all_points)