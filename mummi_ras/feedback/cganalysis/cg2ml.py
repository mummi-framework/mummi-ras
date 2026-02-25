#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2023, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------

from datetime import datetime
import glob
from logging import getLogger, Logger
import multiprocessing
from multiprocessing.pool import ThreadPool as Pool
import os
import tarfile
import tempfile
import threading
import time
from timeit import default_timer as timer
import traceback
from typing import List, Tuple
import signal
import sys

from filelock import FileLock

import numpy.typing as npt
import numpy as np

from mummi_ras.ml.models.autoencoders import CGAutoencoder

LOGGER = getLogger(__name__)

class FeedbackFrames(object):
    """
    Class that watches the directory where cganalysis writes feedback frames.
    Tar up the frames and encode the frames into LS coordinates for the sampler.
    It goes CG -> ML
    """
    def __init__(self,
            database: str,
            path: str,
            model_name: str = None,
            model: CGAutoencoder = None,
            refresh_time: int = 60,
            max_size_mb: int = 50,
            nthreads: int = 1,
            lockfile: str = None,
            logger: Logger = getLogger(__name__)
        ):
        self.path = path
        os.makedirs(self.path, exist_ok=True)
        self.logger = logger 
        self.model = model
        self.model_name = model_name
        self.database = database
        self.max_size = max_size_mb # Size in MB
        self.lastdb_id = 0
        self.refresh_time = refresh_time # in seconds
        self.nthreads = nthreads

        if self.model_name is None:
            self.model = None
        if self.model is None:
            self.model_name = None

        self.lockfile = lockfile or tempfile.mkstemp(prefix="mummi-sampling-", suffix=".lock", dir=None)[1]

        p = multiprocessing.current_process()
        signal.signal(signal.SIGTERM, self.signal_wrapper("Feedback_cg2ml", p.pid))
        signal.signal(signal.SIGINT, self.signal_wrapper("Feedback_cg2ml", p.pid))

        self.logger.info(f"Starting to observe feedback files")
        self.logger.info(f" > Directory               {self.path}")
        if self.model_name is None:
            self.logger.info(f" > No ML model provided")
        else:
            self.logger.info(f" > ML model                {self.model_name}")
        self.logger.info(f" > Feeback Frame DB        {self.database}")
        self.logger.info(f" > Every (seconds)         {self.refresh_time}")
        self.logger.info(f" > Max size (MB)           {self.max_size}")
        self.logger.info(f" > Threads                 {self.nthreads}")

        if self.model is not None and not os.path.exists(database):
            self.init_db(self.database, self.model_name)
        if self.model is not None and not self.check_db(self.database, self.model_name):
            raise ValueError("Wrong frame feedback DB")

    def signal_wrapper(self, name, pid):
        def handler(signum, frame):
            print("")
            self.logger.warning(f"Received (SIGNUM={signum}) for {name}[pid={pid}], stopping the process...")
            try:
                self.logger.debug(f'Removing {self.lockfile}')
                os.remove(self.lockfile)
            except FileNotFoundError as _:
                self.logger.warning(f'Lock file not found {self.lockfile}')
                pass
            sys.exit(1)
        return handler


    def start(self, clean_up: bool = True):
        self.running = True
        self.logger.info(f"Starting feedback frame process")

        def print_failed(caught):
            traceback.print_exc(file = sys.stdout)

        while self.running:
            futures = {}
            frame_processed = {}
            with Pool(self.nthreads) as pool:
                list_folders = list(os.scandir(self.path))
                for folder in list_folders:
                    struct_name = os.path.basename(folder.path)
                    futures[struct_name] = pool.apply_async(self.work, (folder, clean_up), error_callback = print_failed)

                for k,v in futures.items():
                    frame_processed[k] = v.get()

                time.sleep(self.refresh_time)

    def stop(self):
        self.running = False

    def check_db(self, database: str, model_name: str) -> bool:
        """
        Check if DB exists and is valid.
        """
        if not os.path.isfile(database):
            return True

        with np.load(database, allow_pickle=True) as db:
            prev_model_name = db['model_name']
            if model_name != prev_model_name:
                self.logger.error(f"Old feedback DB cannot be used with different ML model.")
                self.logger.error(f"This frame DB {database} has been created with ML model {prev_model_name}")
                self.logger.error(f"You are currently used Latent Space {model_name}")
                return False
        return True

    def init_db(self, database: str, model_name: str):
        """
        Create an empty frame feeback DB.
        """
        with FileLock(self.lockfile):
            dirname = os.path.dirname(database)
            ts = datetime.now().strftime('%d.%m.%Y-%H:%M:%S')
            data = {"model_name" : model_name, "created_at": ts, "updated_at": ts, "structures" : {}}
            if dirname != '':
                os.makedirs(dirname, exist_ok=True)
            np.savez_compressed(database, **data)
            self.logger.info(f"Created new frame DB at {database}")

    def update_db(self, batch_data: dict):
        """
        Update feedback frame feeback DB with new frames.

        batch_data: dict of dict containing LS positions.
            => {
                'structure_iter00_000000000017': [frame_1, frame_2, ..., frame_k],
                'structure_iter00_000000000032': [frame_1, frame_2, ..., frame_k]
            }

            where each frame_k is a numpy array of shape (1, 32).
        """
        start = timer()
        with FileLock(self.lockfile):
            # If the DB is larger than max_size, we move the current database to 
            # self.database_k where k is an index that we increase
            # The new DB self.database is now empty and ready to get new frames
            if os.path.getsize(self.database) >= self.max_size * 1024 * 1024:
                db_name, ext = os.path.splitext(self.database)
                old_db = f"{db_name}_{self.lastdb_id:012d}{ext}"
                os.rename(self.database, old_db)
                self.logger.info(f"Moving current DB ({self.database}) to {old_db}")
                self.lastdb_id += 1
                # We init a new DB with the old name
                self.init_db(self.database, self.model_name)

            # We write into a temporary file to avoid possible file corruption
            timestr = time.strftime("%Y%m%d-%H%M%S")
            db_path = os.path.splitext(self.database)
            database_tmp = f"{db_path[0]}-{timestr}{db_path[1]}"

            new_data = {}

            with np.load(self.database, allow_pickle=True) as db:
                new_data["model_name"] = db["model_name"]
                new_data["created_at"] = db["created_at"]
                new_data["updated_at"] = datetime.now().strftime('%d.%m.%Y-%H:%M:%S')
                new_data["structures"] = db["structures"].tolist() # .tolist() is needed it turns a 0-d array into dict

            self.logger.info(f"Updating {len(batch_data)} structures in {self.database}")
            # Copying new frames
            for struct, all_frames in batch_data.items():
                if not struct in new_data["structures"]:
                    new_data["structures"][struct] = []
                for frame in sorted(all_frames):
                    new_data["structures"][struct].append(all_frames[frame])

            new_elements = sum([len(x) for x in batch_data.values()])

            np.savez_compressed(database_tmp, **new_data)
            # We make sure the new file is not corrupted somehow
            try:
                test = np.load(database_tmp, allow_pickle=True)
                self.logger.debug(f"{database_tmp} is valid {test.files}")
            except Exception as e:
                self.logger.warning(f"{database_tmp} seems to be corrupted. We keep the old {self.database} intact")
                return

            # Now we can replace the temporary file with the real DB. On Unix systems this operation is atomic
            os.replace(database_tmp, self.database)
            end = timer() - start
            self.logger.info(f"Updated {self.database} database with {new_elements} frames in {end:.3f} seconds ({new_elements/end} frame/sec)")

    def parse_fname(self, fname: str) -> Tuple[str, str]:
        """
        Parse the name of feedback frame:
            example: structure_iter00_000000000803_f000000625000.npz
        return [structure_iter00_000000000803, 000000625000]
        """
        fname = os.path.basename(fname)
        fname = os.path.splitext(fname)[0]
        try:
            structure_name, fname_id = fname.split("f")
        except ValueError as _:
            return ["", ""]
        structure_name = structure_name[:-1] # Remove final '_'
        return [structure_name, fname_id]

    def read_frame(self, frame_path: str) -> npt.NDArray:
        """
        Return the content a of frame as a npz file buffer.
        """
        if not os.path.exists(frame_path):
            self.logger.error(f"{frame_path} does not exist exists, aborting.")
            return
        try:
            with np.load(frame_path, allow_pickle = True) as data:
                if "positions" in data.files:
                    return data["positions"]
                elif "data" in data.files:
                    return data["data"]
                else:
                    return np.array([])
        except EOFError as e:
            return np.array([])

    def convert_frame(self, positions: npt.NDArray) -> npt.NDArray:
        """
        Convert physical positions into Latent Space coordinates
        """
        if positions.size == 0:
            return np.array([])
        return self.model.encode(positions)

    def process_frame(self, frame_path: str) -> npt.NDArray:
        """
        Convert a feedback frame into Latent Space coordinates
        """
        return self.convert_frame(self.read_frame(frame_path))

    def updates_all_tar(self, root_path: str, frames_to_tar: dict, clean_up: bool):
        """
        Update a collection of tar archives with fresh frames
        """
        start = timer()
        for struct in frames_to_tar:
            list_files = [frames_to_tar[struct][f] for f in frames_to_tar[struct]]
            name_file = os.path.join(root_path, struct)+".tar"
            self.update_tar(name_file, list_files, clean_up)
            self.logger.info(f"Updated {name_file} with {len(list_files)} frames")

        end = timer() - start
        self.logger.info(f"Tar-ed up {len(frames_to_tar)} structures frames in {end:.3f} seconds ({len(frames_to_tar)/end} struct/sec)")

    def update_tar(self, name: str, list_files: List[str], clean_up: bool):
        """
        Update an existing archive (uncompressed) with a list of files.
        """
        with FileLock(self.lockfile):
            if list_files == []:
                self.logger.error(f"List of files = {list_files}, aborting.")
                return

            opening_mode = 'a' if os.path.exists(name) else 'w'

            self.logger.debug(f"Updating {name} with {len(list_files)} files.")
            try:
                with tarfile.open(name, opening_mode) as archive:
                    for file in list_files:
                        if os.path.isfile(file):
                            archive.add(file, arcname = os.path.basename(file))
                        else:
                            self.logger.warning(f"Skipping {file}, it's not file.")
                        if clean_up:
                            try:
                                os.remove(file)
                            except OSError as e:
                                self.logger.error(f"{e.filename}: {e.strerror}")
            except Exception as e:
                self.logger.error(f"Error while opening {name}: {e}")

    def work(self, folder: os.DirEntry, clean_up: bool = True) -> int:
        # List of top-level directories
        start = timer()
        buffer = glob.glob(folder.path+'/*.npz')
        # True if a model is provided
        encode_frames = self.model is not None
        # structure that contains structure_name => ls_coords
        batch_data = {}
        # Keep track of every structures/frames processed until they are tared up
        frames_to_tar = {}
        tid = threading.current_thread().name

        if buffer == []:
            self.logger.info(f"[{tid}] No new feedbacks files detected in {folder.path}")
            return 0
        self.logger.debug(f"[{tid}] Processing {len(buffer)} new feedbacks files in {folder.path}")

        frame_processed = 0
        struct_name = ""
        for frame in buffer:
            struct_name, fname_id_long = self.parse_fname(frame)
            if struct_name == "" or fname_id_long == "":
                self.logger.debug(f"[{tid}] Skipped ill-formed frame {frame}")
                continue
            fname_id = int(fname_id_long)
            # We encode the frame in LS
            if encode_frames:
                ls_coords = self.process_frame(frame)
                if ls_coords.size == 0:
                    self.logger.warning(f"[{tid}] empty coordinates: {folder.path}/{struct_name}_f{fname_id_long}.npz")
                    continue
                if not struct_name in batch_data:
                    batch_data[struct_name] = {}
                batch_data[struct_name][fname_id] = ls_coords

            if not struct_name in frames_to_tar:
                frames_to_tar[struct_name] = {}
            frames_to_tar[struct_name][fname_id] = frame

            frame_processed += 1
        if struct_name != "" and frame_processed > 0:
            self.logger.debug(f"[{tid}] Found {frame_processed} valid frames for {struct_name}")

        if encode_frames and len(batch_data) > 0:
            # We upgrade the DB
            self.update_db(batch_data)
        if len(frames_to_tar) > 0:
            # We tar everything and deleted the files after
            self.updates_all_tar(folder, frames_to_tar, clean_up = clean_up)

        end = timer() - start
        if frame_processed > 0:
            self.logger.info(f"[{tid}] Added {frame_processed} frames to {self.database} in {end:.3f} seconds ({frame_processed/end} frame/sec)")
        return frame_processed