#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------

import logging

import torch
from .model_torch import BaseModel

LOGGER = logging.getLogger(__name__)

# ------------------------------------------------------------------------------
# A regression model
# ------------------------------------------------------------------------------
class RegressionModel(BaseModel):
    """An Abstract Model."""

    # --------------------------------------------------------------------------
    def __init__(self, inshape, outshape, config, tstamp, create_subdirs, **kwargs):
        super().__init__(config, tstamp, create_subdirs, **kwargs)

        assert self.type == 'xyz_reconstruction', f'Invalid type = {self.type}'

        self.input_shp = inshape
        self.output_shp = outshape
        self.show_plots = False

        # architectures can be defined in architectures.py
        # and chosen using config file
        #arch = getattr(architectures, self.arch_name)
        #self.encoder = arch(self.input_shp, self.z_dim)
        self.encoder = None
        self.decoder = None

        self.loss_fn = torch.nn.MSELoss()

    def _compute_losses(self, data, training=True):

        x, y = data[0], data[1]
        y_pred = self._encode(x, training=training)
        y_pred = y_pred.reshape(y_pred.shape[0], y_pred.shape[1]*y_pred.shape[2])
        # loss_value = torch.sqrt(self.loss_fn(y, y_pred))
        loss_value = self.loss_fn(y, y_pred)
        loss_report = {'loss': loss_value}
        return loss_value, loss_report

    def _encode(self, data, training=False):
        _ = self.forward(data)
        return _ if training else _.detach()

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# a custom regression model
# ------------------------------------------------------------------------------
class CustomRegressionModel(RegressionModel):
    # def __init__(self, inshape, outshape, config, tstamp, create_subdirs):
    #     super(CustomRegressionModel,self).__init__(inshape, outshape, config, tstamp, create_subdirs)

    #     self.linear1 = torch.nn.Linear(inshape[0], 5612)
    #     self.bn1 = torch.nn.BatchNorm1d(5612)
    #     self.linear2 = torch.nn.Linear(5612, 1870)
    #     self.bn2 = torch.nn.BatchNorm1d(1870)
    #     self.linear3 = torch.nn.Linear(1870, 552)
    #     self.bn3 = torch.nn.BatchNorm1d(552)
    #     self.linear4 = torch.nn.Linear(552, 552)
    #     self.drop = torch.nn.Dropout(0.2)
    #     self.act = torch.nn.ReLU()
    #     self.reshape_dim = outshape

    # def forward(self, x):
    #     if not isinstance(x, torch.Tensor):
    #         x = torch.Tensor(x)

    #     x = x.to(self.device)
    #     x = self.drop(self.bn1(self.act(self.linear1(x))))
    #     x = self.drop(self.bn2(self.act(self.linear2(x))))
    #     x = self.drop(self.bn3(self.act(self.linear3(x))))
    #     x = self.linear4(x)
    #     x = torch.reshape(x, (-1, self.reshape_dim[0], self.reshape_dim[1]))
    #     return x


    def __init__(self, inshape, outshape, config, tstamp, create_subdirs):
        super(CustomRegressionModel,self).__init__(inshape, outshape, config, tstamp, create_subdirs)

        # hierarchical design

        self.linear1 = torch.nn.Linear(inshape[0], 7724)
        self.bn1 = torch.nn.BatchNorm1d(7724)
        self.linear2 = torch.nn.Linear(7724, 3862)
        self.bn2 = torch.nn.BatchNorm1d(3862)
        self.linear3 = torch.nn.Linear(3862, 2253)
        self.bn3 = torch.nn.BatchNorm1d(2253)
        self.linear4 = torch.nn.Linear(2253, 2253)

        self.drop1 = torch.nn.Dropout(0.2)
        self.drop2 = torch.nn.Dropout(0.2)
        self.drop3 = torch.nn.Dropout(0.2)

        self.act = torch.nn.ReLU()
        self.reshape_dim = outshape

    def forward(self, x):
        if not isinstance(x, torch.Tensor):
            x = torch.Tensor(x)

        x = x.to(self.device)

        # x = self.drop1(self.bn1(self.act(self.linear1(x))))
        # x = self.drop2(self.bn2(self.act(self.linear2(x))))
        # x = self.drop3(self.bn3(self.act(self.linear3(x))))
        x = (self.bn1(self.act(self.linear1(x))))
        x = (self.bn2(self.act(self.linear2(x))))
        x = (self.bn3(self.act(self.linear3(x))))

        x = self.linear4(x)

        x = torch.reshape(x, (-1, self.reshape_dim[0], self.reshape_dim[1]))
        return x
    # def __init__(self, inshape, outshape, config, tstamp, create_subdirs):
    #     super(CustomRegressionModel,self).__init__(inshape, outshape, config, tstamp, create_subdirs)

    #     # representative CG clusters

    #     self.linear1 = torch.nn.Linear(inshape[0], 2380)
    #     self.bn1 = torch.nn.BatchNorm1d(2380)
    #     self.linear2 = torch.nn.Linear(2380, 793)
    #     self.bn2 = torch.nn.BatchNorm1d(793)
    #     self.linear3 = torch.nn.Linear(793, 360)
    #     self.bn3 = torch.nn.BatchNorm1d(360)
    #     self.linear4 = torch.nn.Linear(360, 360)
    #     self.drop1 = torch.nn.Dropout(0.2)
    #     self.drop2 = torch.nn.Dropout(0.2)
    #     self.drop3 = torch.nn.Dropout(0.2)

    #     self.act = torch.nn.ReLU()
    #     self.reshape_dim = outshape

    # def forward(self, x):
    #     if not isinstance(x, torch.Tensor):
    #         x = torch.Tensor(x)

    #     x = x.to(self.device)

    #     # x = self.drop1(self.bn1(self.act(self.linear1(x))))
    #     # x = self.drop2(self.bn2(self.act(self.linear2(x))))
    #     # x = self.drop3(self.bn3(self.act(self.linear3(x))))
    #     x = (self.bn1(self.act(self.linear1(x))))
    #     x = (self.bn2(self.act(self.linear2(x))))
    #     x = (self.bn3(self.act(self.linear3(x))))

    #     x = self.linear4(x)

    #     x = torch.reshape(x, (-1, self.reshape_dim[0], self.reshape_dim[1]))
    #     return x
    # def __init__(self, inshape, outshape, config, tstamp, create_subdirs):
    #     super(CustomRegressionModel,self).__init__(inshape, outshape, config, tstamp, create_subdirs)

    #     # good for all distances
    #     self.linear1 = torch.nn.Linear(56325, 6000)
    #     self.bn1 = torch.nn.BatchNorm1d(6000)
    #     self.linear2 = torch.nn.Linear(6000, 3000)
    #     self.bn2 = torch.nn.BatchNorm1d(3000)
    #     self.linear3 = torch.nn.Linear(15000, 5000)
    #     self.bn3 = torch.nn.BatchNorm1d(5000)
    #     self.linear4 = torch.nn.Linear(5000, 2253)
    #     self.bn4 = torch.nn.BatchNorm1d(2253)
    #     self.linear5 = torch.nn.Linear(2253, 2253)
    #     self.drop1 = torch.nn.Dropout(0.2)
    #     self.drop2 = torch.nn.Dropout(0.2)
    #     self.drop3 = torch.nn.Dropout(0.2)
    #     self.drop4 = torch.nn.Dropout(0.2)
    #     self.act = torch.nn.ReLU()
    #     self.reshape_dim = outshape

    # def forward(self, x):
    #     if not isinstance(x, torch.Tensor):
    #         x = torch.Tensor(x)

    #     x = x.to(self.device)
    #     x1, x2, x3, x4, x5= torch.chunk(x, 5, 1)
    #     # print(x1.shape)
    #     # print(x2.shape)
    #     # print(x3.shape)
    #     x1 = self.drop1(self.bn1(self.act(self.linear1(x1))))
    #     x1 = self.drop2(self.bn2(self.act(self.linear2(x1))))

    #     x2 = self.drop1(self.bn1(self.act(self.linear1(x2))))
    #     x2 = self.drop2(self.bn2(self.act(self.linear2(x2))))

    #     x3 = self.drop1(self.bn1(self.act(self.linear1(x3))))
    #     x3 = self.drop2(self.bn2(self.act(self.linear2(x3))))

    #     x4 = self.drop1(self.bn1(self.act(self.linear1(x4))))
    #     x4 = self.drop2(self.bn2(self.act(self.linear2(x4))))

    #     x5 = self.drop1(self.bn1(self.act(self.linear1(x5))))
    #     x5 = self.drop2(self.bn2(self.act(self.linear2(x5))))

    #     x = torch.cat((x1, x2, x3, x4, x5),1)
    #     x = self.drop3(self.bn3(self.act(self.linear3(x))))
    #     x = self.drop4(self.bn4(self.act(self.linear4(x))))
    #     x = self.linear5(x)
    #     x = torch.reshape(x, (-1, self.reshape_dim[0], self.reshape_dim[1]))
    #     return x
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# a custom regression model
# ------------------------------------------------------------------------------
class CustomRegressionModelBackmapping(RegressionModel):
    def __init__(self, inshape, outshape, config, tstamp, create_subdirs):
        super(CustomRegressionModelBackmapping,self).__init__(inshape, outshape, config, tstamp, create_subdirs)

        # self.linear1 = torch.nn.Linear(inshape[0], 5000)
        # self.bn1 = torch.nn.BatchNorm1d(5000)
        # self.linear2 = torch.nn.Linear(5000, 7500)
        # self.bn2 = torch.nn.BatchNorm1d(7500)
        # self.linear3 = torch.nn.Linear(7500, 2253)
        # self.bn3 = torch.nn.BatchNorm1d(2253)
        # self.linear4 = torch.nn.Linear(2253, 2253)

        # self.linear1 = torch.nn.Linear(inshape[0], 2700)
        # self.bn1 = torch.nn.BatchNorm1d(2700)
        # self.linear2 = torch.nn.Linear(2700, 8100)
        # self.bn2 = torch.nn.BatchNorm1d(8100)
        # self.linear3 = torch.nn.Linear(8100, 2253)
        # self.bn3 = torch.nn.BatchNorm1d(2253)
        # self.linear4 = torch.nn.Linear(2253, 2253)

        # self.linear1 = torch.nn.Linear(inshape[0], 2700)
        # self.bn1 = torch.nn.BatchNorm1d(2700)
        # self.linear2 = torch.nn.Linear(2700, 2253)
        # self.bn2 = torch.nn.BatchNorm1d(2253)
        # self.linear3 = torch.nn.Linear(2253, 2253)


        self.linear1 = torch.nn.Linear(inshape[0], 10000)
        self.bn1 = torch.nn.BatchNorm1d(10000)
        self.linear2 = torch.nn.Linear(10000, 2253)
        self.bn2 = torch.nn.BatchNorm1d(2253)
        self.linear3 = torch.nn.Linear(2253, 2253)

        # self.drop1 = torch.nn.Dropout(0.2)
        # self.drop2 = torch.nn.Dropout(0.2)
        # self.drop3 = torch.nn.Dropout(0.2)

        self.act = torch.nn.ReLU()
        self.reshape_dim = outshape

    def forward(self, x):
        if not isinstance(x, torch.Tensor):
            x = torch.Tensor(x)

        x = x.to(self.device)

        # x = self.drop1(self.bn1(self.act(self.linear1(x))))
        # x = self.drop2(self.bn2(self.act(self.linear2(x))))
        # x = self.drop3(self.bn3(self.act(self.linear3(x))))

        # x = (self.bn1(self.act(self.linear1(x))))
        # x = (self.bn2(self.act(self.linear2(x))))
        # x = (self.bn3(self.act(self.linear3(x))))

        # x = self.linear4(x)

        x = (self.bn1(self.act(self.linear1(x))))
        x = (self.bn2(self.act(self.linear2(x))))

        x = self.linear3(x)


        x = torch.reshape(x, (-1, self.reshape_dim[0], self.reshape_dim[1]))
        return x