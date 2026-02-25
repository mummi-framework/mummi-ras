"""
Optimal transport based interpolator with a mechanism for dynamic feedback
"""

import ot
import itertools
import numpy as np
import random
from typing import List
from mummi_ras.ml.ls_point import LSPoint
from logging import getLogger

LOGGER = getLogger(__name__)

default = {"k": 10,            # number of neighbors (used to calculate the k nearest neighbors)
           "ep": 0.0001,        # perturbation value
           "lb": 0.0,           # lower bound on lambda
           "ub": 1.0,           # upper bound on lambda
           "int_samples": None,           # # of internal samples
           "int_samples_mult": 10.0,           # # of internal samples
           "edge_weight": 'edge_size',      # selection frequency of the edge
           "method": 'kdist_min_rank'}

class FeedbackOTInterpolator:
    """
    Interpolator using an optimal transport solution with dynamic feedback from
    sampled trajectories.
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
                 feedback_edge_probability=0.10,
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


        # That loop can be long so we do not check if we are sure
        if safe_mode:
            LOGGER.debug(f"Doing assert")
            for _pairs in itertools.product(states, states):
                assert len(_pairs[0]) == len(_pairs[1]), "All sets of LS states need to be the same shape."

        self.kneigh = kneigh
        self.lowerbound = lowerbound
        self.upperbound = upperbound
        self.int_samples = int_samples
        self.edge_weight = edge_weight
        self.feedback_edge_probability = feedback_edge_probability
        self.method = method
        self.num_iter_max = num_iter_max
        self.ot_edges = {}  # dictionary holding the OT edges 
        self.feedback_edges = {}  # dictionary holding the OT edges 

        self.compute_optimal_transport()

        return

    def compute_optimal_transport(self):
        """Solve optimal transport schedule for arbitrary list of source and targets"""

        # reset edges list
        self.reset_ot_edges()

        for pairs in itertools.product(range(len(self.states)), range(len(self.states))):
            i,j = pairs

            # avoid self and double counting
            if i >= j:
                continue

            # compute pairwise schedule, will return a dictionary of i to j interal index mappings
            _edges = self.compute_optimal_transport_pairwise(self._get_lsp_coordinate_array(self.states[i]),
                                                             self._get_lsp_coordinate_array(self.states[j]))

            # add new edges in, storing them as an lsp_coordiante array 
            for k in _edges:
                ename = f"ot_edge_{len(self.ot_edges)+1}"
                lse = np.array([k[0],
                                k[1]])

                self.add_new_ot_edge(ename, lse)

            # _edge_index_global = [ (k[0] + i_offset, k[1] + j_offset) for k in _edges ]

            # self.edges.extend(_edge_index_global)

        return

    def add_new_ot_edge(self, ename, ecrd):
        """add new entry to the list of sampled trajectories"""

        assert isinstance(ename, str), "trajectory label must be a string."

        if ename in self.ot_edges:
            return
        self.ot_edges[ename] = ecrd

        return

    def add_new_fb_edge(self, ename, ecrd):
        """add new entry to the list of sampled trajectories"""

        assert isinstance(ename, str), "trajectory label must be a string."

        if ename in self.feedback_edges:
            return
        self.feedback_edges[ename] = ecrd

        return

    def _get_lsp_coordinate_array(self, points):
        """Return a single ndarray from a list of lspoints"""

        _array = np.array([lsp.get_coordinates() for lsp in points])

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
        edges   = [ (source[i], target[np.argmax(G0[i])]) for i in range(G0.shape[0]) ]

        return edges


    def draw_one_point(self, edge, lmb):
        """
        Return one LS point.
        """

        pwdist = self._compute_piecewise_distances(edge)

        normalized_pwd = np.cumsum(pwdist) / pwdist.sum()

        pw_lambda = np.concatenate((np.zeros(1), normalized_pwd))

        rint = random.random()

        # index to get edge from
        segid = (normalized_pwd < rint).sum()

        _gamma = (rint - pw_lambda[0:segid].sum()) / (pw_lambda[segid+1] - pw_lambda[segid])

        _lsp = edge[segid] + _gamma * (edge[segid+1] - edge[segid])

        return _lsp

    def _compute_piecewise_distances(self, edge):
        """
        Return one LS point.
        """

        pwld = [np.linalg.norm(edge[i+1] - edge[i]) for i in range(0, edge.shape[0]-1)]

        return np.array(pwld)


    def draw_k_points(self, npoints):

        # flip weighted coin whether to draw from OT or from feedback.
        # feedback_coin array is boolean where true means darw from feedback_edges
        # and false means draw from ot_edges.
        if len(self.feedback_edges) > 0:
            feedback_coin = np.random.random(size=npoints) < self.feedback_edge_probability
        #if we don't have any feedback edges, draw with prob 1 from OT
        else:
            feedback_coin = np.zeros(npoints, dtype=bool)

        # arcs is index in edges dictionary
        num_ot_arcs = npoints-feedback_coin.sum()
        num_feedback_arcs = feedback_coin.sum()

        # get OT draws
        arcs = np.random.randint(len(self.ot_edges), size=num_ot_arcs)
        lamb = (self.upperbound - self.lowerbound) * np.random.random(size=num_ot_arcs) + self.lowerbound
        _all_lamb = list(lamb)

        _ot_keys = list(self.ot_edges.keys())
        _ot_edges = [ self.ot_edges[_ot_keys[arc]] for arc in arcs ]

        _npts = [ self.draw_one_point(ed,l) for ed,l in zip(_ot_edges, lamb) ]

        # get feedback draws
        if num_feedback_arcs > 0:
            arcs = np.random.randint(len(self.feedback_edges), size=num_feedback_arcs)
            lamb = (self.upperbound - self.lowerbound) * np.random.random(size=num_feedback_arcs) + self.lowerbound
            _fb_keys = list(self.feedback_edges.keys())
            _fb_edges = [ self.feedback_edges[_fb_keys[arc]] for arc in arcs ]

            _npts.extend([ self.draw_one_point(ed,l) for ed,l in zip(_fb_edges, lamb) ])
            _all_lamb.extend(lamb)

        return _npts, _all_lamb

    def size(self) -> int:
        """ Return the size of the interpolator"""
        return len(self.all_points)

    def remove_ot_edge(self, edge_key):
        """Removes an edge from the LSSampler"""

        del self.ot_edges[edge_key]

        return

    def remove_fb_edge(self, edge_key):
        """Removes an edge from the LSSampler"""

        del self.feedback_edges[edge_key]

        return

    def reset_ot_edges(self):
        """remove all ot_edges from the edge dictionary"""

        self.ot_edges = {}

        return

    def reset_feedback_edges(self):
        """remove all ot_edges from the edge dictionary"""

        self.feedback_edges = {}

        return

    def reset_all_edges(self):
        """remove all entries from LSSampler, OT and feedback database entries."""

        self.reset_feedback_edges()
        self.reset_ot_edges()

        return
