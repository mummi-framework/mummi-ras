# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
from abc import ABC, abstractmethod, abstractproperty
import os
from shutil import copy2, which
from subprocess import PIPE, Popen
import signal
from time import sleep

from logging import getLogger
LOGGER = getLogger(__name__)


class SuccessCriteria(ABC):

    @abstractmethod
    def __init__(self, sim_dir, **kwargs):
        pass

    @abstractmethod
    def __bool__(self):
        pass


class ddcMDLengthCriteria(SuccessCriteria):
    def __init__(self, sim_dir, length=1.0e9, snapshot_freq=25000, dt=20):
        self.length = length
        self.snap_freq = snapshot_freq
        self.dt = dt
        self.sim_dir = sim_dir

    def get_snapshot(self):
        return f"snapshot.{int(self.length / self.dt):012d}"

    def __bool__(self):
        snapshot_dir = \
            os.path.join(self.sim_dir, self.get_snapshot(), "subset#000000")
        return os.path.exists(snapshot_dir)


class Simulation(ABC):

    def __init__(self, sim_name, input_dir, out_dir, bin_path, r_limit,
                 **kwargs):
        # Back end process for simulation management
        self.sim_name = sim_name
        self.restart_limit = r_limit  # Limit to number of restarts
        self._setup = False           # Is the simulation setup?
        self._process = None          # Initialize the process to None
        self._success = {}            # Label success criteria

        # Pathing for inputs and outputs
        self.sim_dir = os.path.abspath(out_dir)     # Path for simulation run
        self.input_dir = os.path.abspath(input_dir) # Path containing inputs
        self.bin_path = which(bin_path)   # Path to simulation code
        LOGGER.info("Was ddcMD in PATH? [bin_path=%s]: %s", bin_path, self.bin_path)
        if not self.bin_path:
            self.bin_path = os.path.abspath(bin_path)
            LOGGER.info("ddcMD was not in PATH: abspath=%s", self.bin_path)

    def add_success_criteria(self, criteria):
        if not isinstance(SuccessCriteria, criteria):
            raise TypeError(
                "A simulation success criteria must be of type "
                "'SuccessCriteria'. Received an object of type "
                f"'{type(criteria)}'.")

        criteria_name = criteria.__class__.__name__
        if criteria_name in self._success:
            raise ValueError(
                f"A SuccessCriteria of type '{criteria_name}' has already "
                "been added."
            )

        self._success[criteria_name] = criteria

    @abstractmethod
    def initialize(self):
        pass

    @abstractmethod
    def setup_restart(self):
        pass

    @abstractmethod
    def equilibrate(self):
        pass

    @abstractproperty
    def command(self):
        pass

    def run(self):
        if not self._setup:
            raise Exception("Simulation not setup prior to running. Aborting.")

        #stdio_log = os.path.join(self.sim_dir, "ddcmd.log")
        #stdio_log = open(stdio_log, "w")
        LOGGER.info("Starting ddcMD....")

        kwargs = {
            "shell":                True,
            "universal_newlines":   True,
            "preexec_fn":           os.setpgrp,
            "cwd":                  self.sim_dir,
            "env":                  os.environ,
            "stdout":               PIPE,
            "stderr":               PIPE,
            "close_fds":            False,
        }

        cmd = self.command
        LOGGER.info("cmd =       %s", cmd)
        LOGGER.info("os.pwd =    %s", os.getcwd())
        LOGGER.info("simdir =    %s", self.sim_dir)
        self._process = Popen(cmd, **kwargs)

        sleep(5)
        LOGGER.info("Process Running? %s", self._process.poll())
        cuda_visible_devices = os.environ.get("CUDA_VISIBLE_DEVICES", "UNSET")
        LOGGER.info("CUDA_VISIBLE_DEVICES=%s", cuda_visible_devices)

        if self._process.poll():
            out, err = self._process.communicate()
            LOGGER.info(
                "---------------- %s stdout --------------\n%s",
                self.sim_name, out)
            LOGGER.error(
                "---------------- %s stderr --------------\n%s",
                self.sim_name, err)

    def __bool__(self):
        return bool(self._process is not None)

    def stop(self):
        if self._process and self._process.poll() is None:
            # self._process.terminate()
            os.killpg(os.getpgid(self._process.pid), signal.SIGTERM)

    @property
    def running(self):
        return bool(
            self._process is not None and
            self._process.poll() is None
        )

    @property
    def successful(self):
        # TODO: Come back to this and implement a more flexible success
        # check.
        pass


class ddcMD(Simulation):

    def initialize(self):
        # If we've already run, we shouldn't setup again.
        if self._process is not None:
            return

        # Create the simulation workspace.
        if not os.path.exists(self.sim_dir):
            os.mkdir(self.sim_dir)

        # Copy input files to the directory.
        # This path will usually be the createsim directory as it is
        # supposed to
        # inputs = [
        #     'ConsAtom.data', 'martini.data', 'molecule.data', 'object.data',
        #     'restraint.data', 'topol.tpr', 'lipids-water-eq4.gro',
        #     'lipids-water-eq4hvr.gro', 'POPX_Martini_v2.0_lipid.itp',
        #     'resItpList'
        # ]
        inputs = []

        try:
            for in_file in inputs:
                src = os.path.join(self.input_dir)
                copy2(src, self.sim_dir)
        except Exception as err:
            LOGGER.error(err.message)
            raise err

        # TODO: Commented for now, but will be added in the future.
        # Set up restart file (copy over and symlink)
        # restart_path = os.path.join(self.input_dir, self.restart_snap)
        # copytree(restart_path, self.sim_dir)
        # src = os.path.join(self.sim_dir, self.restart_snap, "restart")
        # dst = os.path.join(self.sim_dir, "restart")
        # os.symlink(src, dst)

        self._setup = True

    @property
    def command(self):
        return f"{self.bin_path} -o object.data molecule.data >> ddcmd.out"

    def equilibrate(self):
        pass

    def setup_restart(self):
        pass

    @property
    def successful(self):
        success = ddcMDLengthCriteria(self.sim_dir)
        return bool(success)
