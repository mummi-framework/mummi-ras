#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2023, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------

import os
import yaml
import numpy as np
import logging
from collections import defaultdict

import torch
from torch.utils.data import DataLoader, TensorDataset
# from torchsummary import summary

USE_HOROVOD = bool(os.getenv('USE_HOROVOD', False))
if USE_HOROVOD:
    import horovod.torch as hvd

USE_DDP = bool(os.getenv('USE_DDP', False))
if USE_DDP:
    import torch.distributed as dist
LOGGER = logging.getLogger(__name__)

# ------------------------------------------------------------------------------
# A base class for a pytorch model
# allows managing some simple state variables so we can quickly change
# the model type and architecture
# ------------------------------------------------------------------------------
class BaseModel (torch.nn.Module):
    """An Abstract Base Model Class."""

    # --------------------------------------------------------------------------
    def __init__(self, config, tstamp, create_subdirs, **kwargs):
        super().__init__(**kwargs)

        self.config = config
        self.type = self.config['model']['type']
        self.arch_name = self.config['model']['arch']
        self.input_shp = self.config['model']['inshape']
        self.output_shp_ae = self.config['model']['outshape']
        self.z_dim = self.config['model']['z_dim']

        if isinstance(self.z_dim, list):
            sign = "\u00D7"
            self.model_name = f'model-{self.type}-{self.arch_name}-' \
                              f'c{self.input_shp[-1]}_z{(sign.join(str(x) for x in self.z_dim))}-{tstamp}'
        else:
            self.model_name = f'model-{self.type}-{self.arch_name}-' \
                              f'c{self.input_shp[-1]}_z{self.z_dim}-{tstamp}'

        ## self.model_name = 'test_model'
        self.outpath = os.path.join(self.config['outpath'], self.model_name)
        LOGGER.info(f'Initializing ({self.model_name})')

        self.encoder = None
        self.decoder = None         # could be None

        # ----------------------------------------------------------------------
        # training settings
        self.nvalidation_vis = config['data']['nvalidation_vis']
        self.batch_size = config['training']['batch_size']
        self.num_epochs = config['training']['nepochs']
        self.learning_rate_fn = config['training']['learning_rate'] #0.001

        # both are lists: use constant decay!
        lrate_vals = self.config['training'].get('learning_rate_vals', None)
        lrate_bdrs = self.config['training'].get('learning_rate_bdrs', None)

        if isinstance(lrate_vals, list) and isinstance(lrate_bdrs, list):
            assert len(lrate_vals) == 1 + len(lrate_bdrs)
            LOGGER.info(f'Using Learning rate: {lrate_vals}')
            assert False, 'how to use this for torch?'
            #self.learning_rate_fn = tf.keras.optimizers.schedules.PiecewiseConstantDecay(lrate_bdrs, lrate_vals)

        elif isinstance(lrate_vals, (int, float)):
            self.learning_rate_fn = float(lrate_vals)

        # need to define optimizer *after* the model
        self.optimizer = None

        # should define the loss function here!
        self.loss_fn = None

        # default device is "cpu"
        self.device = torch.device("cpu")

        # ----------------------------------------------------------------------
        # the model will also save intermediate data
        self.tpath = os.path.join(self.outpath, 'training')
        self.cpath = os.path.join(self.outpath, 'checkpoints')

        if create_subdirs:
            os.makedirs(self.outpath, exist_ok=True)
            os.makedirs(self.tpath, exist_ok=True)
            os.makedirs(self.cpath, exist_ok=True)

        # ----------------------------------------------------------------------

    # --------------------------------------------------------------------------
    def display(self):
        return
        print (self.encoder, self.input_shp)
        # summary(self.encoder, self.input_shp)
        if self.decoder is not None:
            summary(self.decoder, self.z_dim)

    # --------------------------------------------------------------------------
    def save_model(self, save_configs=True, save_weights=True):
        LOGGER.info(f'saving the model to ({self.outpath})')
        # with open(os.path.join(self.outpath, 'model.pth'), 'w') as fh:
        fh = os.path.join(self.outpath, 'model.pth')
        torch.save(self.state_dict(), fh)

        # ----------------------------------------------------------------------
        if save_configs:
            # write the model config (created by me)
            if self.config is not None:
                with open(os.path.join(self.outpath, 'config.yaml'), 'w') as fh:
                    yaml.dump(self.config, fh)

            # write the model summary as a string
            with open(os.path.join(self.outpath, 'model_summary.txt'), 'w') as fh:
                print(self, file=fh)
                fh.write('\n Total params:')
                print(sum([param.nelement() for param in self.parameters()]), file=fh)

    def load_model(self, model_path, outpath, std_path = "CG_pos_data_summary_pos_dis_C2.npz", map_location = torch.device('cpu'), compile = False):
        LOGGER.info(f'loading model ({model_path})')
        fname = os.path.join(outpath, std_path)
        checkpoint = torch.load(model_path, map_location = map_location)
        self.load_state_dict(checkpoint['model_state_dict'])

        LOGGER.info(f'loading model standardization file ({fname})')
        with np.load(fname) as npz:
            mean_mt = npz['mean_pos']
            std_mt = npz['std_pos']
            self.datamean = mean_mt
            self.datastd = std_mt

    def standardize_zero_std_2d(self, data, label, ind, dmean=None, dstd=None):
        if dmean is None or dstd is None:
            dmean, dstd = self.compute_standardization(data, label)
        # LOGGER.info(f'Standardizing {data.shape} {label}')
        data -= dmean
        for j in range(data.shape[0]):
            if j == ind:
                data[j,0:2] = data[j,0:2] / dstd[j,0:2]
            else:
                data[j,:] = data[j,:] / dstd[j,:]
        return data, dmean, dstd

    def standardize_zero_std_3d(self, data, label, ind, dmean=None, dstd=None):
        if dmean is None or dstd is None:
            dmean, dstd = self.compute_standardization(data, label)
        # LOGGER.info(f'Standardizing {data.shape} {label}')
        data -= dmean
        for j in range(data.shape[1]):
            if j == ind:
                data[:,j,0:2] = data[:,j,0:2] / dstd[j,0:2]
            else:
                data[:,j,:] = data[:,j,:] / dstd[j,:]
        return data, dmean, dstd

    def standardize(self, positions) -> np.ndarray:
        # index 430 is farnesyl (CYF BB) which is at z=0
        if positions.ndim == 2:
            positions, dmean, dstd = self.standardize_zero_std_2d(positions, 'positions', 430, self.datamean, self.datastd)
        elif positions.ndim == 3:
            positions, dmean, dstd = self.standardize_zero_std_3d(positions, 'positions', 430, self.datamean, self.datastd)
        return positions

    # --------------------------------------------------------------------------
    def predict(self, x, max_chunk=2500):

        self.eval()
        
        n, d = x.shape[0], self.z_dim
        if isinstance(d, list):
            d = tuple(d)
            codes = np.zeros((n,) + d, dtype=np.float32)
        else:
            codes = np.zeros((n, d), dtype=np.float32)
        idxs = list(range(0, n, max_chunk)) + [n]
        for i in range(len(idxs) - 1):
            i0, i1 = idxs[i], idxs[i + 1]
            codes[i0:i1] = self._encode(x[i0:i1]).cpu().numpy()
        return codes

    def reconstruct_ae(self, x, max_chunk=2500):

        self.eval()

        n, d = x.shape[0], self.output_shp_ae[0]
        if isinstance(d, list):
            d = tuple(d)
            codes = np.zeros((n,) + d, dtype=np.float32)
        else:
            codes = np.zeros((n, d), dtype=np.float32)
        idxs = list(range(0, n, max_chunk)) + [n]
        for i in range(len(idxs) - 1):
            i0, i1 = idxs[i], idxs[i + 1]
            # codes[i0:i1] = self._encode(x[i0:i1]).cpu().numpy()
            # pos, dis, stats = self._reconstruct(x[i0:i1])
            pos, dis = self._reconstruct(x[i0:i1])
            codes[i0:i1] = pos.cpu().numpy()
        return codes

    def reconstruct_ae_joint(self, x, max_chunk=2500):

        self.eval()

        n, d = x.shape[0], self.output_shp_ae[0]
        if isinstance(d, list):
            d = tuple(d)
            codes = np.zeros((n,) + d, dtype=np.float32)
        else:
            codes = np.zeros((n, d), dtype=np.float32)
        idxs = list(range(0, n, max_chunk)) + [n]
        for i in range(len(idxs) - 1):
            i0, i1 = idxs[i], idxs[i + 1]
            codes[i0:i1] = self._reconstruct(x[i0:i1]).cpu().numpy()
        return codes

    def encode_ae(self, x, max_chunk=2500):

        self.eval()

        n, d = x.shape[0], self.z_dim
        if isinstance(d, list):
            d = tuple(d)
            codes = np.zeros((n,) + d, dtype=np.float32)
        else:
            codes = np.zeros((n, d), dtype=np.float32)
        idxs = list(range(0, n, max_chunk)) + [n]
        for i in range(len(idxs) - 1):
            i0, i1 = idxs[i], idxs[i + 1]
            xt = torch.Tensor(x[i0:i1])
            codes[i0:i1] = self._encode(xt.to(self.device)).detach().cpu().numpy()
        return codes

    def reconstruct_from_ls_ae(self, x, max_chunk=2500):

        self.eval()

        n, d = x.shape[0], self.input_shp[0]
        if isinstance(d, list):
            d = tuple(d)
            codes = np.zeros((n,) + d, dtype=np.float32)
        else:
            codes = np.zeros((n, d), dtype=np.float32)
        idxs = list(range(0, n, max_chunk)) + [n]
        for i in range(len(idxs) - 1):
            i0, i1 = idxs[i], idxs[i + 1]
            xt = torch.Tensor(x[i0:i1])
            codes[i0:i1] = self.decoder(xt.to(self.device))[0].detach().cpu().numpy()
        return codes

    def decode(self, z):
        return self._decode(z).cpu().numpy()

    def reconstruct(self, x):
        return self._decode(self._encode(x)).cpu().numpy()

    def predict_labels(self, x, max_chunk=2500):
        n = x.shape[0]

        labels = np.zeros(n, dtype=np.uint8)
        idxs = list(range(0, n, max_chunk)) + [n]
        for i in range(len(idxs) - 1):
            i0, i1 = idxs[i], idxs[i + 1]
            labels[i0:i1] = self._predict_labels(x[i0:i1]).cpu().numpy()
        return labels

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    def train_step(self, data):
        self.train() # prep model for training

        # zero the parameter gradients
        self.optimizer.zero_grad()

        # forward + backward + optimize
        total_loss, loss_metrics = self._compute_losses(data, training=True)

        total_loss.backward()  # compute the gradient
        self.optimizer.step()  # update the parameters
        for k,v in loss_metrics.items():
            loss_metrics[k] = v.item()
        return loss_metrics

    def test_step(self, data):

        self.eval()    # prep model for evaluation

        # since we're not training, we don't need to calculate the gradients for our outputs
        with torch.no_grad():
            total_loss, loss_metrics = self._compute_losses(data, training=False)

        val_loss_metrics = {}
        for k,v in loss_metrics.items():
            val_loss_metrics['val_'+k] = v.item()
        return val_loss_metrics

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    def train_model(self, train_data, validation_data, **kwargs):

        callback = kwargs.get('callback', None)
        self.model_ddp = kwargs.get('model_ddp', None)
        epochs = kwargs.get('epochs', self.num_epochs)
        batch_size = kwargs.get('batch_size', self.batch_size)
        learning_rate = kwargs.get('learning_rate', self.learning_rate_fn)

        metrics = {}

        # ----------------------------------------------------------------------
        # input assertions
        assert isinstance(train_data, tuple)
        assert isinstance(validation_data, tuple)
        assert len(train_data) == len(validation_data)

        ncomps = len(train_data)
        ntrain = train_data[0].shape[0]
        nvalid = validation_data[0].shape[0]
        for i in range(1, ncomps):
            assert ntrain == train_data[i].shape[0]
            assert nvalid == validation_data[i].shape[0]

        assert 2 == ncomps, 'Forcing only two components for now!'

        # ----------------------------------------------------------------------
        self.optimizer = torch.optim.Adam(self.parameters(), lr=learning_rate)

        if USE_HOROVOD:
            self.optimizer = hvd.DistributedOptimizer(self.optimizer,
                                                      named_parameters=self.named_parameters())
            # Broadcast parameters from rank 0 to all other processes.
            hvd.broadcast_parameters(self.state_dict(), root_rank=0)

        # if USE_DDP:
        #     # Broadcast parameters from rank 0 to all other processes.
        #     dist.broadcast_object_list(self.state_dict(), 0)


        # ----------------------------------------------------------------------
        data_kwargs = {'batch_size': batch_size}
        if 'cuda' == self.device.type:
            cuda_kwargs = {'num_workers': 1,
                           'pin_memory': True,
                           'shuffle': True}
            data_kwargs.update(cuda_kwargs)

        dataset_train = TensorDataset(torch.Tensor(train_data[0]), torch.Tensor(train_data[1]))
        dataset_val = TensorDataset(torch.Tensor(validation_data[0]), torch.Tensor(validation_data[1]))

        train_loader = torch.utils.data.DataLoader(dataset_train, **data_kwargs)
        test_loader = torch.utils.data.DataLoader(dataset_val, **data_kwargs)

        # ----------------------------------------------------------------------
        if callback is not None:
            callback.on_train_begin({})

        for epoch in range(1, int(epochs) + 1):

            if callback is not None:
                callback.on_epoch_begin(epoch)
                callback.on_epoch_train_begin(epoch)

            m1_dict = defaultdict(list)
            for batch_idx, (data, target) in enumerate(train_loader):

                # get the inputs
                data, target = data.to(self.device), target.to(self.device)

                m1 = self.train_step((data, target))
                for key, val in m1.items():
                    m1_dict[key].append(np.float32(val))

            if callback is not None:
                callback.on_epoch_train_end(epoch)
                callback.on_epoch_validate_begin(epoch)

            m2_dict = defaultdict(list)
            for batch_idx, (data, target) in enumerate(test_loader):
                # get the inputs
                data, target = data.to(self.device), target.to(self.device)
                m2 = self.test_step((data, target))
                for key, val in m2.items():
                    m2_dict[key].append(np.float32(val))
            if callback is not None:
                callback.on_epoch_validate_end(epoch)

            # ------------------------------------------------------------------
            # collect the metrics from training and validation
            for key, val in m1_dict.items():
                metrics[key] = np.asarray(val)

            for key, val in m2_dict.items():
                metrics[key] = np.asarray(val)

            for key, val in metrics.items():
                 if "confusion" not in key:
                     metrics[key] = np.nanmean(val, axis=0)
                 else:
                     metrics[key] = np.sum(val, axis=0)

            # ------------------------------------------------------------------
            if callback is not None:
                callback.on_epoch_end(epoch, metrics)

            # # barrier at each epoch!
            if USE_DDP:
                dist.barrier()

        if callback is not None:
            callback.on_train_end({})

        # return all the metrics captured during training
        return metrics



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# --------------------------------------------------------------------------
    def train_model_ae(self, train_data, **kwargs):

        callback = kwargs.get('callback', None)
        epochs = kwargs.get('epochs', self.num_epochs)
        batch_size = kwargs.get('batch_size', self.batch_size)
        learning_rate = kwargs.get('learning_rate', self.learning_rate_fn)

        metrics = {}

        # ----------------------------------------------------------------------
        # input assertions
        assert isinstance(train_data, tuple)

        ncomps = len(train_data)
        ntrain = train_data[0].shape[0]
        for i in range(1, ncomps):
            assert ntrain == train_data[i].shape[0]

        assert 2 == ncomps, 'Forcing only two components for now!'

        # ----------------------------------------------------------------------
        self.optimizer = torch.optim.Adam(self.parameters(), lr=learning_rate)

        if USE_HOROVOD:
            self.optimizer = hvd.DistributedOptimizer(self.optimizer,
                                                      named_parameters=self.named_parameters())

            # Broadcast parameters from rank 0 to all other processes.
            hvd.broadcast_parameters(self.state_dict(), root_rank=0)

        # ----------------------------------------------------------------------
        data_kwargs = {'batch_size': batch_size}
        if 'cuda' == self.device.type:
            cuda_kwargs = {'num_workers': 1,
                           'pin_memory': True,
                           'shuffle': True}
            data_kwargs.update(cuda_kwargs)

        dataset_train = TensorDataset(torch.Tensor(train_data[0]), torch.Tensor(train_data[1]))

        train_loader = torch.utils.data.DataLoader(dataset_train, **data_kwargs)

        # ----------------------------------------------------------------------
        if callback is not None:
            callback.on_train_begin({})

        for epoch in range(1, int(epochs) + 1):

            if callback is not None:
                callback.on_epoch_begin(epoch)
                callback.on_epoch_train_begin(epoch)

            m1_dict = defaultdict(list)
            for batch_idx, (data, target) in enumerate(train_loader):

                # get the inputs
                data = data.to(self.device)

                m1 = self.train_step((data))
                for key, val in m1.items():
                    m1_dict[key].append(np.float32(val))

            if callback is not None:
                callback.on_epoch_train_end(epoch)

            # ------------------------------------------------------------------
            # collect the metrics from training and validation
            for key, val in m1_dict.items():
                metrics[key] = np.asarray(val)

            for key, val in metrics.items():
                 if "confusion" not in key:
                     metrics[key] = np.nanmean(val, axis=0)
                 else:
                     metrics[key] = np.sum(val, axis=0)

            # ------------------------------------------------------------------
            if callback is not None:
                callback.on_epoch_end(epoch, metrics)

        if callback is not None:
            callback.on_train_end({})

        # return all the metrics captured during training
        return metrics

# ------------------------------------------------------------------------------
