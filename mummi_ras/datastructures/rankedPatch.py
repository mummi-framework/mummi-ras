# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
"""Patch module."""

import numpy as np


# ------------------------------------------------------------------------------
# Ranked patch is used by patch selector
# ------------------------------------------------------------------------------
class RankedPatch:
    """A Patch with its lspace coordinates and rank."""

    def __str__(self):
        return '<' + str(self.id) + ':' + str(self.rank) + '>'

    def __repr__(self):
        return self.__str__()

    def __init__(self, id, lcoords, rank = float('nan')):
        """Initialize a RankedPatch.

        Args:
            id (str):           Patch ID
            lcoords (ndarray):  encoded coordinates of this patch
            rank (float):       rank of this patch

        """
        assert isinstance(id, (int, str))
        assert isinstance(lcoords, np.ndarray)
        assert isinstance(rank, float)

        self.id = id                # id
        self.lcoords = lcoords      # latent space coordinates
        self.rank = rank            # ranking

'''
# ------------------------------------------------------------------------------
# A coarse grained patch
# ------------------------------------------------------------------------------

class cgPatch(Patch):
    """A Phase Field Patch."""

    def __str__(self):

        s = '({}) at time ({} {}) was created from ({})'
        return s.format(self.id, self.config.simtime, self.config.tunit,
                        self.pfpatch)

    def __repr__(self):
        return self.__str__()

    def __init__(self, config, pfpatchid, id, ngrid, nspecies):
        """Initialize a cgPatch."""

        self.config = config
        self.pfpatch = pfpatchid    # id of the pfpatch it was created from
        Patch.__init__(self, id, (ngrid, ngrid, nspecies))
'''
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
