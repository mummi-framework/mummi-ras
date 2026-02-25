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

class LSSamplingDB_C2:
    def __init__(self, database_path: str, model_name: str, lockfile: str = None, logger: Logger = LOGGER):
        self.database_path = database_path
        self.model_name = model_name
        self.lockfile = lockfile if lockfile else self.create_tempfile(prefix="mummi-sampling-", suffix=".lock")
        self.logger = logger

        if self.exist_db():
            if not self.checking_db():
                # Rename corrupted DB
                self.rename_file(self.database_path, f"{self.database_path}.orig")
                # Recreate new DB
                self.create_db()
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
                createsims = sampling_db["createsims_status"]

                return f"Created at {created_at} / updated at {updated_at} {len(structs)} structures with createsims = {createsims}"
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
            "created_at": ts,
            "updated_at": ts,
            "structure_names": np.empty(0, dtype = str),
            "ls_coords": np.empty(0, dtype = float),
            "lambda_values": np.empty(0, dtype = bool),
            "validation_status": np.empty(0, dtype = bool),
            "createsims_status": np.empty(0, dtype = int)
        }
        np.savez_compressed(self.database_path, **data)
        self.logger.info(f"Created fresh DB at {self.database_path}")

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

        existing_model_name = sampling_db['model_name']

        if self.model_name != existing_model_name:
            self.logger.error(f"Old sampling DB cannot be used with different ML model.")
            self.logger.error(f"This DB {self.database_path} has been created with ML model {existing_model_name}")
            self.logger.error(f"You are currently used Latent Space {self.model_name}")
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

    def update_db(self,
            structure_names: List[str],
            ls_coords: List[List[float]] = [],
            lambda_values: List[float] = [],
            validation_status: List[bool] = [],
            createsims_status: List[bool] = [],
            num_new_sample: int = 0
        ) -> int:
        """
        Update sampling feeback DB with new samples.
        """
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
                    prev_coords = sampling_db["ls_coords"]
                    prev_lambda = sampling_db["lambda_values"]
                    prev_validation = sampling_db["validation_status"]
                    prev_createsims = sampling_db["createsims_status"]
            except Exception as e:
                self.logger.error(f"{self.database_path} {e}")
                return None
            
            nb_updated_structures = 0
            for struct, status in zip(structure_names, prev_createsims):
                if struct not in createsims_status:
                    nb_updated_structures += 1
                if struct in createsims_status and status != createsims_status[struct]:
                    nb_updated_structures += 1

            data = {
                "model_name": self.model_name,
                "created_at": created_at,
                "updated_at": self.create_timestamp()
            }

            # Update sampling DB for feedback
            data = self.update_field(data, "structure_names", prev_struct, structure_names)

            if ls_coords != []:
                if len(prev_coords) > 0:
                    data["ls_coords"] = np.append(prev_coords, ls_coords, axis=0)
                else:
                    data["ls_coords"] = ls_coords
            else:
                data["ls_coords"] = prev_coords

            data = self.update_field(data, "lambda_values", prev_lambda, lambda_values)
            data = self.update_field(data, "validation_status", prev_validation, validation_status)

            # 0 means False, 1 True and 2 is unknown state for createsims
            # As no createsims for these structures are running, they are all unknowns
            if createsims_status != []:
                data["createsims_status"] = createsims_status
            else:
                data["createsims_status"] = np.append(prev_createsims, np.full(num_new_sample, 2, dtype=int))

            np.savez_compressed(database_tmp, **data)

            # We make sure the new file is not corrupted somehow
            # and if it is corrupted we do not rename it to keep the old 
            # file intact
            if not self.is_db_corrupted(database_tmp):
                self.rename_file(database_tmp, self.database_path)

            self.logger.debug(f"Added new structures: {structure_names}")
            self.logger.info(f"Upgraded DB at {self.database_path} with {nb_updated_structures} new records")

        return nb_updated_structures
