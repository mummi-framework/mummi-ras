#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2025, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------
# @authors
#           Loic Pottier <pottier1@llnl.gov>
# ------------------------------------------------------------------------------

import os
import signal
import traceback
import numpy as np
from typing import List, Union, Any
import time
import datetime
import tempfile
import sys
from filelock import FileLock
from timeit import default_timer as timer

from logging import getLogger, Logger
LOGGER = getLogger(__name__)

class LSSamplingDB_C3:
    def __init__(self, database_path: str, model_name: List[str], ls_encoder = None, lockfile: str = None, logger: Logger = LOGGER):
        self.database_path = database_path
        self.model_name = model_name
        self.lockfile = lockfile if lockfile else self.create_tempfile(prefix="mummi-sampling-", suffix=".lock")
        self.logger = logger
        self.ls_encoder = ls_encoder

        self.logger.info(f"Initialization of {self.__class__.__name__} with:")
        self.logger.info(f" > database_path = {database_path}")
        self.logger.info(f" > model_name = {model_name}")
        if self.ls_encoder:
            self.logger.info(f" > LS Encoder = {self.ls_encoder.get_model_name()}")

        with FileLock(self.lockfile):
            if self.exist_db():
                self.logger.info(f"The DB {self.database_path} exists")
                if not self.checking_db():
                    self.logger.info(f"The DB {self.database_path} is corrupted")
                    # Rename corrupted DB
                    self.rename_file(self.database_path, f"{self.database_path}.orig")
                    # Recreate new DB
                    self.create_db()
                else:
                    self.logger.info(f"The DB {self.database_path} is valid")
            else:
                self.create_db()

    def create_timestamp(self):
        """
        Return a timestamp with a predefined format
        """
        return datetime.datetime.now().strftime("%d.%m.%Y-%H:%M:%S")

    def create_tempfile(self, prefix: str = None, suffix: str = None, dir: str=None):
        """
        Return the path for a temporary file
        """
        return tempfile.mkstemp(prefix=prefix, suffix=suffix, dir=dir)[1]

    def print_db(self) -> str:
        """
        Return the path for a temporary file
        """
        try:
            with np.load(self.database_path, allow_pickle = True) as sampling_db:
                created_at = sampling_db["created_at"]
                updated_at = sampling_db["updated_at"]
                structs = sampling_db["structure_names"]
                validation = sampling_db["validation_status"]
                createsims = sampling_db["createsims_status"]

                return f"Created at {created_at} / updated at {updated_at} {len(structs) / len(validation)} structures with createsims = {createsims}"
        except Exception as e:
            return ""

    def rename_file(self, src: str, dst: str):
        """
        Rename a file in an atomic way (safe)
        """
        self.logger.debug(f"Renaming {src} -> {dst} ")
        try:
            os.replace(src, dst)
        except OSError as e:
            self.logger.error(f"{e}")

    def create_db(self):
        """
        Create an empty sampling feeback DB.
        """
        dirname = os.path.dirname(self.database_path)
        ts = self.create_timestamp()
        if dirname != '':
            os.makedirs(dirname, exist_ok = True)
        data = {
            "model_name": self.model_name,
            "model_encoder": self.ls_encoder.get_model_name(),
            "created_at": ts,
            "updated_at": ts,
            "structure_names": np.empty(0, dtype = str),
            "v_coords": [],
            "v2ucg_coords": [],
            "ucg2cg_coords": [],
            "ls_coords": [],
            "validation_status": np.empty(0, dtype = bool),
            "createsims_status": np.empty(0, dtype = int)
        }
        np.savez_compressed(self.database_path, **data)
        self.logger.info(f"Created fresh DB for C3 at {self.database_path}")

    def exist_db(self) -> bool:
        """
        Check if the DB exists
        """
        return os.path.isfile(self.database_path)

    def checking_db(self) -> bool:
        """
        Check if the sampling database is valid. A sampling database
        is only valid for a given latent space.
        """
        if not self.exist_db():
            return None
        try:
            sampling_db = np.load(self.database_path, allow_pickle=True)
        except Exception as e:
            self.logger.error(f"{self.database_path} {e}")
            return False

        existing_model_name = sampling_db['model_name'] # Here we have two ML models
        if self.ls_encoder is not None and "model_encoder" in sampling_db.files:
            existing_encoding_name = sampling_db['model_encoder'] # UCG encoder name
            if existing_encoding_name != self.ls_encoder.get_model_name():
                self.logger.warning(f"Sampling DB cannot be used with different UCG Encoder models. {existing_encoding_name} != {self.ls_encoder.get_model_name()}")
                return False

        if self.ls_encoder is None and "model_encoder" in sampling_db.files:
            model_name_encoder = sampling_db['model_encoder']
            self.logger.warning(f"Sampling DB requires a valid UCG Encoder models as the original DB has an encoder ({model_name_encoder}).")
            return False

        if self.ls_encoder is not None and "model_encoder" not in sampling_db.files:
            self.logger.warning(f"Sampling DB cannot be used with an UCG Encoder models when the original DB has no encoder.")
            return False

        if self.model_name[0] != existing_model_name[0] or self.model_name[1] != existing_model_name[1]:
            self.logger.warning(f"Sampling DB cannot be used with different ML model.")
            self.logger.warning(f"This DB {self.database_path} has been created with models:")
            self.logger.warning(f" - {existing_model_name[0]}")
            self.logger.warning(f" - {existing_model_name[1]}")
            self.logger.warning(f"Your are using:")
            self.logger.warning(f" - {self.model_name[0]}")
            self.logger.warning(f" - {self.model_name[1]}")        
            self.logger.warning(f"Feedback is deactivated.")
            return False
        return True

    def update_field(self, data: dict, key: str, existing_data, new_data):
        if new_data != []:
            data[key] = np.append(existing_data, new_data)
        else:
            data[key] = existing_data
        return data

    def read_data(self):
        """
        Return the underlying data from DB
        """
        data = {}
        with FileLock(self.lockfile):
            try:
                with np.load(self.database_path, allow_pickle = True) as np_data:
                    # We cannot return directly np_data as it will be closed after the return
                    for k,v in np_data.items():
                        data[k] = v
            except Exception as e:
                self.logger.error(f"Error: could not read {self.database_path}: {e}")
                return None
            return data

    def is_db_corrupted(self, database: str) -> bool:
        """
        Test if a database is corrupted or not
        """
        try:
            with np.load(database, allow_pickle = True) as test:
                self.logger.debug(f"{database} is valid {test.files}")
                return False
        except Exception as e:
            self.logger.error(f"{database} seems to be corrupted. We keep the old {database} intact")
            return True

    def _update_array(self, data, key: str, prev_values: List[float], new_values: List[List[float]]):
        """
        Update an array that originated from a NPZ file.

        new_values are batched and on the form (10, X, 3) where X varies and 10 is batch number
        """
        if len(new_values) == 0:
            data[key] = prev_values
            return data
        if len(prev_values) != 0:
            data[key] = np.append(prev_values, new_values, axis=0)
        else:
            data[key] = new_values
        return data

    def encode(self, cg_positions):
        """
        Return UCG encodings given one CG position of shape (4716, 3)

        The UCG encodings are of shape (20,)
        """
        if len(cg_positions) == 0 or self.ls_encoder is None:
            return np.array([])
        return self.ls_encoder.encode(cg_positions)

    def update_db(self,
        structure_names: List[str],
        num_new_sample: int = 0,
        v_coords: List[List[float]] = [],
        v2ucg_coords: List[List[float]] = [],
        ucg2cg_coords: List[List[float]] = [],
        validation_status: List[bool] = [],
        createsims_status: List[bool] = []
        ) -> int:
        """
        Update sampling feeback DB with new samples.
        """
        self.logger.debug(f"Starting updating sampling DB")

        with FileLock(self.lockfile):
            # To avoid corrupting the file in case of interruption
            database_tmp = self.create_tempfile(prefix=f"{os.path.splitext(self.database_path)[0]}-", suffix=".npz")

            if self.exist_db():
                if not self.checking_db():
                    return None
            else:
                self.create_db()

            try:
                with np.load(self.database_path, allow_pickle=True) as sampling_db:
                    created_at = sampling_db["created_at"]
                    prev_struct = sampling_db["structure_names"]
                    prev_v_coords = sampling_db["v_coords"]
                    prev_v2ucg_coords = sampling_db["v2ucg_coords"]
                    prev_ucg2cg_coords = sampling_db["ucg2cg_coords"]
                    prev_ls_coords = sampling_db["ls_coords"]
                    prev_validation = sampling_db["validation_status"]
                    prev_createsims = sampling_db["createsims_status"]
            except Exception as e:
                self.logger.error(f"{self.database_path} {e}")
                return None

            new_records = len(structure_names)

            nb_updated_structures = 0
            for struct, status in zip(structure_names, prev_createsims):
                if struct not in createsims_status:
                    nb_updated_structures += 1
                if struct in createsims_status and status != createsims_status[struct]:
                    nb_updated_structures += 1

            data = {
                "model_name": self.model_name,
                "model_encoder": self.ls_encoder.get_model_name(),
                "created_at": created_at,
                "updated_at": self.create_timestamp()
            }

            if self.ls_encoder is not None:
                ls_coords = []
                for cg_position in ucg2cg_coords:
                    ls_coords.append(self.encode(cg_position))
                ls_coords = np.array(ls_coords)

            # Update sampling DB for feedback
            data = self.update_field(data, "structure_names", prev_struct, structure_names)
            data = self.update_field(data, "validation_status", prev_validation, validation_status)

            data = self._update_array(data, "v_coords", prev_v_coords, v_coords)
            data = self._update_array(data, "v2ucg_coords", prev_v2ucg_coords, v2ucg_coords)
            data = self._update_array(data, "ucg2cg_coords", prev_ucg2cg_coords, ucg2cg_coords)
            
            if self.ls_encoder is not None:
                # self.logger.info(f"Writting ls_coords = {ls_coords}")
                # self.logger.info(f"Writting prev_ls_coords = {prev_ls_coords}")
                data = self._update_array(data, "ls_coords", prev_ls_coords, ls_coords)
                # self.logger.info(f"Writting after data = {data}")

            # 0 means False, 1 True and 2 is unknown state for createsims
            # As no createsims for these structures are running, they are all unknowns
            if createsims_status != []:
                data["createsims_status"] = createsims_status
            else:
                new_stat = np.full(num_new_sample, 2, dtype=int)
                data["createsims_status"] = np.append(prev_createsims, new_stat)

            try:
                np.savez_compressed(database_tmp, **data)
            except Exception as e:
                self.logger.error(f"{e}: data = {data}")
                raise e

            # We make sure the new file is not corrupted somehow
            # and if it is corrupted we do not rename it to keep the old 
            # file intact
            if not self.is_db_corrupted(database_tmp):
                self.rename_file(database_tmp, self.database_path)

            self.logger.debug(f"Added new structures: {structure_names}")
            self.logger.info(f"Updated DB at {self.database_path} with {new_records} new records (updated {nb_updated_structures} createsims status)")

        return nb_updated_structures
