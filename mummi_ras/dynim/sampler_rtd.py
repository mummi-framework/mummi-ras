import numpy as np
import scipy

import logging
LOGGER = logging.getLogger(__name__)

from .hdpoint import HDPoint
from .hdspace import HDSpace
from .sampler import Sampler
from .utils import take_backup

# ------------------------------------------------------------------------------
# An importance sampler of hd points
# ------------------------------------------------------------------------------

class SamplerRTD (Sampler):

    # currently developed for 3d binning
    ldim = 3
    def __init__(self, name: str, workspace: str,
                       BufferSize: int, nMinSamples: int):

        self.type = 'SamplerRTD'
        super(SamplerRTD, self).__init__(name, workspace, BufferSize)

        self.nMinSamples = nMinSamples
        self.rdmax = 180.
        self.ranges       = np.array([[-180.0,180.0], [0.,180.0], [0.,10.]])
        self.dist_factors = np.array([1./(r[1]-r[0]) for r in self.ranges])

    # --------------------------------------------------------------------------
    def _add_selections(self, points):
        pass

    # --------------------------------------------------------------------------
    def _rtd_distance(self, u, v):
        d = u-v
        if d[0] < -self.rdmax:  d[0] += self.rdmax
        if d[0] >  self.rdmax:  d[0] -= self.rdmax
        d *= self.dist_factors
        return (d**2).sum()

    def _update_ranks(self):

        LOGGER.info('Updating ranks ({})'.format(self.__str__()))
        ncandidates = len(self.candidates)
        nselected = len(self.selected)

        # if no selection has been made, compare against other candidates
        if nselected == 0:
            compare_coords = np.array([s.coords for s in self.candidates])
            knn = min(10, ncandidates)

        # otherwise, compare against previous selections
        else:
            compare_coords = np.array([s.coords for s in self.selected])
            knn = min(10, nselected)

        # now, assign ranks to all candidates
        cand_coords = np.array([s.coords for s in self.candidates])
        dists = scipy.spatial.distance.cdist(cand_coords, compare_coords,
                                                metric=self._rtd_distance)

        for cid in range(ncandidates):
            d = dists[cid,:]
            self.candidates[cid].rank = np.sort(dists[cid,:])[:knn].mean()

        '''
        for cand in self.candidates:
            cand_coords = np.expand_dims(cand.coords, axis=0)
            d = scipy.spatial.distance.cdist(cand_coords, compare_coords,
                                                  metric=self._rtd_distance)

            # rank = average of k nearest neighbors
            cand.rank = np.sort(d[0])[:knn].mean()
        '''
        LOGGER.info('Updated ranks ({})'.format(self.__str__()))

    # --------------------------------------------------------------------------
    def checkpoint(self):
        super(SamplerRTD, self).checkpoint()

    def restore(self):
        super(SamplerRTD, self).restore()

    def test(self):
        pass
# ------------------------------------------------------------------------------
