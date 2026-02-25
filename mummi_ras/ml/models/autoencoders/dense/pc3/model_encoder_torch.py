#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------

import logging

import torch
import torch.nn as nn

LOGGER = logging.getLogger(__name__)


# ------------------------------------------------------------------------------
# An encoder class
# ------------------------------------------------------------------------------
class Encoder(nn.Module):

    # --------------------------------------------------------------------------
    def __init__(self, inshape, lantent_dims):
        super(Encoder, self).__init__()

        self.input_shp = inshape
        self.output_shp = lantent_dims
        self.show_plots = False
        # 17581
        self.linear1 = torch.nn.Linear(inshape[0], 8790)
        self.bn1 = torch.nn.BatchNorm1d(8790, momentum=0.01)
        self.linear2 = torch.nn.Linear(8790, 4395)
        self.bn2 = torch.nn.BatchNorm1d(4395, momentum=0.01)
        self.linear3 = torch.nn.Linear(4395, lantent_dims[0])
        self.bn3 = torch.nn.BatchNorm1d(lantent_dims[0], momentum=0.01)
        
        self.linear4 = torch.nn.Linear(lantent_dims[0], lantent_dims[0])

        # self.linear3 = torch.nn.Linear(4395, 2197)
        # self.bn3 = torch.nn.BatchNorm1d(2197, momentum=0.01)
        # self.linear4 = torch.nn.Linear(2197, lantent_dims[0])
        # self.bn4 = torch.nn.BatchNorm1d(lantent_dims[0], momentum=0.01)
        # self.linear4 = torch.nn.Linear(lantent_dims[0], lantent_dims[0])

        # self.drop1 = torch.nn.Dropout(0.2)
        # self.drop2 = torch.nn.Dropout(0.2)
        # self.drop3 = torch.nn.Dropout(0.2)

        self.act = torch.nn.ReLU()
        # self.reshape_dim = lantent_dims

    def forward(self, x):

        # x = self.drop1(self.bn1(self.act(self.linear1(x))))
        # x = self.drop2(self.bn2(self.act(self.linear2(x))))
        # x = self.drop3(self.bn3(self.act(self.linear3(x))))
        x = (self.bn1(self.act(self.linear1(x))))
        x = (self.bn2(self.act(self.linear2(x))))
        x = (self.bn3(self.act(self.linear3(x))))

        x = self.linear4(x)

        # x = torch.reshape(x, (-1, self.reshape_dim[0], self.reshape_dim[1]))
        return x
    