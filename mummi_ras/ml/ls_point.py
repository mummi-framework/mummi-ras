#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2023, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------

import h5py

class LSPoint:
    """
    Class for representing a point in latent space.
    """
    def __init__(self, pid, parent=None, ls_coord=None, valid=False, initial_lamda=None):
        """Construct a LS point"""

        assert type(pid) is int, "LS point ID must be an integer." 

        self._pid = pid # a unique ID, we enforce to be immutable  

        # data structures describing LS point 
        self.parent = parent # which set thsi point belongs to (we assume can only be "A", "B" or None)
        self.validated = valid # whether this point has been validated 
        self.ls_coord = ls_coord  # LS coordinate
        self.initial_lamda = initial_lamda
        return

    """
    # override setattr and delattr in order to protect the LS point ID from being modified after initialization.  
    def __setattr__(self, name, value):
        if name == "_id":
            raise AttributeError("Cannot modify LSPoint ID.")
        else:
            return type.__setattr__(self, name, value)

    def __delattr__(self, name):
        if name == "_id":
            raise AttributeError("Cannot delete LSPoint ID.")
        else:
            return type.__delattr__(self, name)
    """

    def __repr__(self):
        """String representation of LSPoint."""

        return f"lspoint_{self._pid:012}" 

    def __eq__(self, other):
        """Equivalence testing for LS points. Only checks PID, not LS coordinates."""

        return self._pid == other._pid

    def get_validation_status(self):
        """Get whether this point has been validated."""

        return self.validated

    def set_validation_status(self, validation_status):
        """Set the validation status."""

        self.validated = bool(validation_status)

        return 

    def set_coordinates(self, new_coord):
        """Set the LS coordinates."""

        #assert new_coord.shape == self.ls_coord.shape, "New coordinate is not in the expected shape." 

        # assign numpy array. is this the best way to do this? we want to avoid shallow copies 
        self.ls_coord = new_coord.copy()

        return 

    def get_coordinates(self):
        """Return the ls_coordinates"""
        return self.ls_coord

    def get_lamda(self):
        """Return the initial_lamda"""
        return self.initial_lamda

    def get_parent(self):
        """Return a string representing the parent set for this point."""

        return self.parent

    def set_parent(self, new_parent):
        """Set the parent."""

        self.parent = new_parent

        return 

    def get_id(self):
        """Return the pid"""

        return self._pid
 
def write_ls_points(filename, points):
    """Write an interable of LSPoints to disk in .npz format.
    
    For now we are storing this as an HDF5 file.

    The format in the hdf5 file is to store the coordinates of the point as a dataset of name "lspoint_**" with the pid label in the name.

    The PID, parent and validation status fields are stored as meta-data via the attributes feature of the dataset. 
    """

    with h5py.File(filename, "w") as fh:
        for p in points:
            ds = fh.create_dataset(f"lspoint_{p._pid:012}", data=p.get_coordinates())

            # store these criterion as attributes to the coordinate dataset
            ds.attrs.create("pid", p._pid)
            ds.attrs.create("validation_status", p.get_validation_status())
            ds.attrs.create("parent", p.get_parent())

    return


def read_ls_points(filename):
    """Read an HDF5 file of LSPoints."""

    points = []

    with h5py.File(filename, "r") as fh:

        # read in all the points from disk
        for pname in fh.keys():

            p_ds = fh[pname]
            pnt = LSPoint(int(p_ds.attrs.get("pid")))

            pnt.set_validation_status(p_ds.attrs.get("validation_status"))
            pnt.set_parent(p_ds.attrs.get("parent"))

            pnt.set_coordinates(p_ds[:])

            points.append(pnt)

    return points

