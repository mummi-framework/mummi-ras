# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
import os
import yaml
import numpy as np
import warnings

os.environ['KERAS_BACKEND'] = 'theano'
import keras

from logging import getLogger
LOGGER = getLogger(__name__)

# ------------------------------------------------------------------------------
# Harsh copied this architecture manually from the training repo
#from architectures import model6_metric
def model6_metric(encoder_in, encoder_out, decoder_in):

    encoder = keras.Sequential([
        keras.layers.InputLayer(input_shape=encoder_in),

        # across channels
        keras.layers.SeparableConv2D(filters=6, depth_multiplier=6, kernel_size=1, strides=1, activation='relu'),
        keras.layers.BatchNormalization(),
        keras.layers.Conv2D(filters=16, kernel_size=3, strides=(2, 2), activation='relu'),
        keras.layers.BatchNormalization(),
        keras.layers.Conv2D(filters=16, kernel_size=3, strides=(2, 2), activation='relu'),
        keras.layers.BatchNormalization(),
        keras.layers.Flatten(),

        # No activation
        keras.layers.Dense(encoder_out),
        keras.layers.BatchNormalization(),
    ])
    return encoder, None
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
class PretrainedModel:
    """Load a saved model for inference."""

    def __init__(self, mpath):

        LOGGER.info('> Loading pretrained model ({})'.format(mpath))
        self.path = mpath

        # read the config file
        config_fname = os.path.join(self.path, 'config.yaml')
        with open(config_fname, 'r') as data:
            config = yaml.load(data, Loader=yaml.Loader)

        # C3 had model_name key
        self.name = config.get('model_name', None)
        if self.name is None:
            # C4 changed the key :-/
            self.name = config['model'].get('name', None)
        assert self.name is not None
        self.input_shp = config['model']['inshape']
        self.z_dim = config['model']['z_dim']

        # create the keras architectire
        assert config['model']['arch'] == 'model6_metric'
        self.encoder, _ = model6_metric(self.input_shp, self.z_dim, self.z_dim)
        assert self.input_shp == self.encoder.layers[0].input_shape[1:]
        assert self.z_dim == self.encoder.layers[-1].output_shape[1]

        # load the weights
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning)
            self.encoder.load_weights(os.path.join(self.path, 'encoder_weights.h5'))

        LOGGER.info(f'Successfully loaded pretrained model ({mpath})')

        # load the data standardizer
        norm_type = config['data'].get('norm_type', None)
        assert 'standardize' == norm_type
        self.do_std = True
        if self.do_std:
            std_file = os.path.join(self.path, 'data_summary_standardized.npz')
            with np.load(std_file) as std_data:
                self.data_mean = std_data['mean']
                self.data_std = std_data['std']
            assert self.data_mean.shape == self.data_std.shape
            assert self.data_mean.shape[0] == self.input_shp[-1]
            LOGGER.info(f'Successfully loaded data standardizer: {list(self.data_mean)} {list(self.data_std)}')

        self.channel_weights = config['data'].get('channel_weights', None)
        if self.channel_weights is not None:
            self.channel_weights = np.array(self.channel_weights)
            assert self.channel_weights.shape[0] == 14
            LOGGER.info(f'Successfully loaded channel_weights: {list(self.channel_weights)}')

    def display(self):
        self.encoder.summary()

    def encode(self, x):
        if self.do_std:
            for i in range(14):
                x[:,:,:,i] -= self.data_mean[i]
                x[:,:,:,i] /= self.data_std[i]

        if self.channel_weights is not None:
            for i in range(14):
                x[:,:,:,i] *= self.channel_weights[i]

        return self.encoder.predict(x)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
