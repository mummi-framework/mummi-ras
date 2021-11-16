# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
import io
import numpy as np

from logging import getLogger
LOGGER = getLogger(__name__)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
class CGSnapshot (object):

    @staticmethod
    def is_valid(tilt, rot, depth):
        return (0. <= tilt <= np.pi) and \
               (-np.pi <= rot <= np.pi) and not np.isnan(rot) and \
               (depth >= 0.)

    def __init__(self, simname, frame_id, tilt, rot, depth):
        self.simname = simname
        self.frame_id = frame_id

        if np.isnan(depth) or depth < 0.:
            LOGGER.warning(f'Clamping depth={depth} to 0')
            depth = 0.

        if CGSnapshot.is_valid(tilt, rot, depth):
            self.data = np.array([tilt, rot, depth])
        else:
            self.data = None

    def __str__(self):
        return '<{}:{}  {}>'.format(self.simname, self.frame_id, self.data)

    def __repr__(self):
        return self.__str__()

    # --------------------------------------------------------------------------
    # save and load patches as numpy archives (npz)
    @staticmethod
    def save_npz(fptr, p):
        assert isinstance(fptr, (str, io.BytesIO))
        assert isinstance(p, CGSnapshot)
        np.savez_compressed(fptr, type='cgSnapshot', data = p.data,
                            simname = p.simname, frame_id = p.frame_id)
        return fptr

    @staticmethod
    def load_npz(fptr):
        assert isinstance(fptr, (str, io.BytesIO))
        try:
            data = np.load(fptr, allow_pickle=True)
            assert data['type'] == 'cgSnapshot'
            d = data['data']
            cg = CGSnapshot(data['simname'], data['frame_id'], d[0], d[1], d[2])
            data.close()

        except Exception as e:
            LOGGER.error(f'Failed to load CGSnapshot from ({fptr}): {e}')
            cg = None

        return cg

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
