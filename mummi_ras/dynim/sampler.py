# Copyright 2020 Lawrence Livermore National Security, LLC and other
# DynIm Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: MIT

################################################################################

import os
from typing import List, Union
import logging

import numpy as np
from .hdpoint import HDPoint
from .utils import format_time, take_backup, filter_list, write_history

LOGGER = logging.getLogger(__name__)


# ------------------------------------------------------------------------------
# A base class for sampling functionality
# ------------------------------------------------------------------------------
class Sampler(object):

    # --------------------------------------------------------------------------
    def __init__(self, name: str, workspace: str,
                 buffer_size: int = 0,
                 min_cands_b4_sel: int = 0):

        """ Initialize a Sampler.
        Args:
            workspace (str):    a workspace directory where
                                history and checkpoints will be written
            buffer_size (int):  size of the candidate buffer_size
            min_cands_b4_sel (int):  min number of samples before any selection
        """
        assert isinstance(name, str)
        assert isinstance(workspace, str)
        assert isinstance(buffer_size, int) and int(buffer_size) >= 0
        assert isinstance(min_cands_b4_sel, int) and int(min_cands_b4_sel) >= 0

        self.name = name
        self.workspace = workspace
        self.buffer_size = buffer_size
        self.min_cands_b4_sel = min_cands_b4_sel
        os.makedirs(self.workspace, exist_ok=True)

        # candidates for the next set of selection (sorted list of samples)
        self.candidates = np.array([])

        # samples received since last selection (list of samples)
        self.cached = []

        # samples that have been selected (in the order of selection)
        self.selected = []

        # samples that have been dropped (maintain for recalculation of weights)
        self.discarded = []

        # function to explicitly scale the rank
        self.rank_scaler = None

        # prepare to write the history and checkpoints
        self.schkpt = os.path.join(self.workspace, f'{self.__tag__()}.data.npz')
        self.hists = os.path.join(self.workspace, f'{self.__tag__()}.history.csv')

    # --------------------------------------------------------------------------
    def __type__(self):
        return type(self).__name__

    def __tag__(self):
        return f'{self.__type__()}-{self.name}'

    def __str__(self):
        return f'[{self.__type__()}<{self.name}>:' \
               f' {len(self.cached)} cached,' \
               f' {len(self.candidates)} candidates,' \
               f' {len(self.selected)} selected,' \
               f' {len(self.discarded)} discarded]'

    def __repr__(self):
        return self.__str__()

    def _state_dict(self):
        return {'ncached': len(self.cached),
                'ncandidates': len(self.candidates),
                'nselected': len(self.selected),
                'ndiscarded': len(self.discarded)}

    # --------------------------------------------------------------------------
    def num_candidates(self):
        return len(self.cached) + len(self.candidates)

    def num_selected(self):
        return len(self.selected)

    def set_rank_scaler(self, rank_scaler: callable) -> None:
        assert callable(rank_scaler)
        LOGGER.info(f'Setting rank scalar ({rank_scaler}) for ({self.__tag__()})')
        self.rank_scaler = rank_scaler

    # --------------------------------------------------------------------------
    def _apply_ranks_to_candidates(self, ranks):
        assert isinstance(ranks, np.ndarray)
        assert ranks.dtype == np.float32
        assert ranks.shape[0] == len(self.candidates)

        if self.rank_scaler is None:
            LOGGER.info(f'Applying {len(self.candidates)} ranks to candidates')
            for i, sample in enumerate(self.candidates):
                sample.rank = ranks[i]

        else:
            LOGGER.info(f'Scaling {len(self.candidates)} ranks using ({self.rank_scaler})')
            for i, sample in enumerate(self.candidates):
                sample.rank = np.float32(self.rank_scaler(sample.id, ranks[i]))

    # --------------------------------------------------------------------------
    # add new candidates to the sampler
    # --------------------------------------------------------------------------
    def add_candidates(self, points: List[HDPoint]) -> None:

        assert isinstance(points, list)
        assert all([isinstance(p, HDPoint) for p in points])

        n = len(points)
        if n == 0:
            return

        # ----------------------------------------------------------------------
        def _filter(_points, _check_against, _tag):
            _points, _discarded = filter_list(_points, _check_against)
            _n = len(_discarded)
            if _n == 0:
                return _points
            LOGGER.warning(f'Rejecting {_n} already selected points: {_discarded}')
            write_history(self.hists, f'rejected:add_cands:already_{_tag}',
                          _discarded, self._state_dict())
            return _points

        # ----------------------------------------------------------------------
        LOGGER.debug(f'Adding {n} candidates to ({self})')

        # remove the ones we have already seen
        points = _filter(points, self.selected, 'selected')
        points = _filter(points, self.discarded, 'discarded')
        points = _filter(points, self.cached, 'cached')
        points = _filter(points, self.candidates, 'candidates')

        # finally, add the ones we have not seen before
        self.cached.extend(points)
        LOGGER.info(f'Added {n} candidates to ({self})')
        write_history(self.hists, 'added', points, self._state_dict())

    # --------------------------------------------------------------------------
    # update the state of the sampler
    # --------------------------------------------------------------------------
    def update(self):

        LOGGER.debug(f'Updating sampler ({self})')

        # move the cached candidates to actual candidate list
        self.candidates = np.append(self.candidates, self.cached)
        self.cached = []

        n = len(self.candidates)
        if n > 0:
            # update the rankings (to be implemented by the actual selector!)
            self._update_ranks()

            # sort the candidates based on ranking!
            ranks = np.array([p.rank for p in self.candidates])
            order = np.argsort(ranks)[::-1]
            self.candidates = self.candidates[order]

            # discard the tail if needed!
            if 0 < self.buffer_size < n:
                discard = list(self.candidates[self.buffer_size:])
                self.candidates = self.candidates[:self.buffer_size]
                self.discarded.extend(discard)

                # write discards
                write_history(self.hists, 'discarded', discard, self._state_dict())

        LOGGER.debug(f'Updated sampler ({self})')

    # --------------------------------------------------------------------------
    # select k new candidates
    # --------------------------------------------------------------------------
    def select(self, k: int, confirm_selection: bool = True) -> List[HDPoint]:

        assert isinstance(k, int)
        assert k >= 0
        if k == 0:
            return []

        if (self.min_cands_b4_sel > 0) and \
           (self.num_selected() == 0) and \
           (self.num_candidates() < self.min_cands_b4_sel):
            LOGGER.debug(f'Not selecting due to too few candidates from {self}')
            return []

        self.update()
        LOGGER.debug(f'Selecting {k} samples from {self}')

        # pick the top k candidates
        k = min(k, len(self.candidates))
        selection = list(self.candidates[:k])

        LOGGER.info(f'Selected {k} samples from {self}')

        if confirm_selection:
            self._confirm_selections(selection)
            LOGGER.info(f'Confirmed selection of {k} samples from {self}')

        return selection

    # --------------------------------------------------------------------------
    # confirm the selection of these points
    # --------------------------------------------------------------------------
    def confirm_selections(self, selection: Union[List[HDPoint], List[str]]) -> bool:

        assert isinstance(selection, list)

        is_ids = False
        is_samples = all([isinstance(p, HDPoint) for p in selection])
        if not is_samples:
            is_ids = all([isinstance(p, str) for p in selection])

        if not is_ids and not is_samples:
            raise AttributeError('confirm_selections expects a list of '
                                 'ids or list of samples')

        k = len(selection)
        if k == 0:
            return False

        LOGGER.debug(f'Confirming selection of {k} samples from {self}')

        # ----------------------------------------------------------------------
        # make sure that the given selection are the top candidates
        if is_samples:
            for i in range(k):
                if self.candidates[i] != selection[i]:
                    raise AttributeError(f'Cannot confirm selection: '
                                         f'[{i} = {selection[i]}] not a candidate')
        else:
            for i in range(k):
                if self.candidates[i].id != selection[i]:
                    raise AttributeError(f'Cannot confirm selection: '
                                         f'[{i} = {selection[i]}] not a candidate')

        # ----------------------------------------------------------------------
        selection = list(self.candidates[:k])

        # now, check that these are not previously selected
        for i in range(k):
            if selection[i] in self.selected:
                raise AttributeError(f'Cannot confirm selection: '
                                     f'[{i} = {selection[i]}] already selected')

            if selection[i] in self.discarded:
                raise AttributeError(f'Cannot confirm selection: '
                                     f'[{i} = {selection[i]}] already discarded')

        # ----------------------------------------------------------------------
        self._confirm_selections(selection)
        LOGGER.info(f'Confirmed selection of {k} samples from {self}')

    # --------------------------------------------------------------------------
    def _confirm_selections(self, selection: List[HDPoint]) -> bool:

        assert isinstance(selection, list)
        assert all([isinstance(p, HDPoint) for p in selection])

        k = len(selection)
        if k == 0:
            return False

        LOGGER.debug(f'Confirming selection of {k} samples')

        # mark the selections as sampled (to be implemented by the child class!)
        self._add_selections(selection)

        # remove selected from the candidate list
        self.candidates = self.candidates[k:]

        # and add to the selected list
        self.selected.extend(selection)

        # write selections
        write_history(self.hists, 'selected', selection, self._state_dict())

    # --------------------------------------------------------------------------
    def invalidate_patches(self, validator):
        assert callable(validator)

        def _filter_valid(_data, tag):
            _ids = [_.id for _ in _data]
            _flags = [_ for _ in validator(_ids)]
            _valids = [_data[i] for i, _ in enumerate(_flags) if _]
            _invalids = [_data[i] for i, _ in enumerate(_flags) if not _]

            LOGGER.info(f'{tag} = {len(_data)}: '
                        f'valid = {len(_valids)}, invalid = {len(_invalids)}')
            return _valids, _invalids

        LOGGER.info(f'testing for invalid patches for {self}')

        self.cached, _invalid_cached = _filter_valid(self.cached, 'cached')
        self.candidates, _invalid_cands = _filter_valid(self.candidates, 'candidates')
        LOGGER.info(f'after invaliding patches: {self}')

        if len(_invalid_cands) > 0:
            _hfile = self.hists
            write_history(_hfile, 'invalidated:cached', _invalid_cached, self._state_dict())
            write_history(_hfile, 'invalidated:candidates', _invalid_cands, self._state_dict())

        return len(_invalid_cands)

    # --------------------------------------------------------------------------
    # checkpoint and restore
    # --------------------------------------------------------------------------
    def _checkpoint(self):

        st = format_time()
        LOGGER.info(f'Checkpointing Sampler data {self} at {st}')

        # take backup of previous checkpoint
        take_backup(self.schkpt)
        np.savez_compressed(self.schkpt,
                            t=st, type=self.__type__(), name=self.name,
                            cached=self.cached, candidates=self.candidates,
                            selected=self.selected, discarded=self.discarded)

    # --------------------------------------------------------------------------
    def _restore(self, filename=None):

        if filename is None:
            filename = self.schkpt

        try:
            LOGGER.debug(f'Restoring from ({filename})')
            data = np.load(filename, allow_pickle=True)

            assert self.name == data['name']
            assert self.__type__() == data['type']

            self.cached = list(data['cached'])
            self.selected = list(data['selected'])
            self.candidates = np.array(data['candidates'])
            self.discarded = list(data['discarded'])
            data.close()

        except Exception as e:
            raise Exception(f'{type(e).__name__}: {e} ({filename})')

        LOGGER.info(f'Successfully restored {self}')
        return True

    # --------------------------------------------------------------------------
    '''
    def _del(self):
        def restore_v1(self):
    
            # ----------------------------------------------------------------------
            def _read_chkpts(fstate, fdata):
                try:
                    with open(fstate, 'r') as infile:
                        state = yaml.load(infile, Loader=yaml.Loader)
    
                    # restore data!
                    LOGGER.debug(f"Restoring data from {state['t']}")
                    sys.stdout.flush()
    
                    data = np.load(fdata, allow_pickle=True)
                    assert self.name == data['name']
                    assert self.__type__() == data['type']
                    cached = list(data['cached'])
                    selected = list(data['selected'])
                    candidates = np.array(data['candidates'])
                    discarded = list(data['discarded'])
    
                except Exception as e:
                    LOGGER.error(f'Failed to restore {fdata}. Error = {e}')
                    raise e
                return state, data
    
            # ----------------------------------------------------------------------
            _s,_d = self.chkpts['state'], self.chkpts['data']
    
            # if the primary checkpoint file does not exist, just go back!
            for f in [_s,_d]:
                if not os.path.isfile(f):
                    LOGGER.debug(f'Checkpoint file ({f}) does not exist!')
                    return False
    
            # prepare to read primary and backup files (if needed)
            _files = [(_s,_d), (backup_fname(_s), backup_fname(_d))]
            for i in range(2):
                try:
                    state, data = _read_chkpts(_files[i][0], _files[i][1])
                    break
    
                except Exception as e:
                    if i == 1:      # can't try again!
                        return False
    
            # otherwise, copy over
            self.cached = list(data['cached'])
            self.selected = list(data['selected'])
            self.candidates = np.array(data['candidates'])
            self.discarded = list(data['discarded'])
            data.close()
            LOGGER.info(f'Restored {self}')
            return True
    
        # --------------------------------------------------------------------------
        def restore_v0(self):
    
            for k in self.chkpts.keys():
                if not os.path.isfile(self.chkpts[k]):
                    LOGGER.debug(f'Checkpoint file ({self.chkpts[k]}) does not exist!')
                    return False
    
            # load yaml
            with open(self.chkpts['state'], 'r') as infile:
                state = yaml.load(infile, Loader=yaml.Loader)
    
            # restore data!
            LOGGER.debug(f"Restoring data from {state['t']}")
            sys.stdout.flush()
    
            try:
                data = np.load(self.chkpts['data'], allow_pickle=True)
                assert self.name == data['name']
                assert self.__type__() == data['type']
    
            except Exception as e:
                LOGGER.error(f'Failed to restore {self.chkpts["data"]}. Error = {e}')
                raise e
    
            self.cached = list(data['cached'])
            self.selected = list(data['selected'])
            self.candidates = np.array(data['candidates'])
            self.discarded = list(data['discarded'])
    
            assert self.test()
            data.close()
            LOGGER.info(f'Restored {self}')
            return True
    '''
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
