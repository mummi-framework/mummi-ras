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
import torch.nn as nn

LOGGER = logging.getLogger(__name__)


# ------------------------------------------------------------------------------
# A decoder class
# ------------------------------------------------------------------------------
class Decoder(nn.Module):

    # --------------------------------------------------------------------------
    def __init__(self, inshape, lantent_dims):
        super(Decoder, self).__init__()

        self.input_shp = lantent_dims
        self.output_shp = inshape
        self.show_plots = False

        self.linear1 = torch.nn.Linear(lantent_dims[0], 4395)
        self.bn1 = torch.nn.BatchNorm1d(4395)
        self.linear2 = torch.nn.Linear(4395, 8790)
        self.bn2 = torch.nn.BatchNorm1d(8790)
        self.linear3 = torch.nn.Linear(8790, inshape[0])
        self.bn3 = torch.nn.BatchNorm1d(inshape[0])
        # self.linear4 = torch.nn.Linear(inshape[0], inshape[0])

        # self.drop1 = torch.nn.Dropout(0.2)
        # self.drop2 = torch.nn.Dropout(0.2)
        # self.drop3 = torch.nn.Dropout(0.2)

        self.act = torch.nn.ReLU()
        # self.reshape_dim = outshape

    def forward(self, x):

        # x = self.drop1(self.bn1(self.act(self.linear1(x))))
        # x = self.drop2(self.bn2(self.act(self.linear2(x))))
        # x = self.drop3(self.bn3(self.act(self.linear3(x))))
        x = (self.bn1(self.act(self.linear1(x))))
        x = (self.bn2(self.act(self.linear2(x))))
        x = (self.bn3(self.act(self.linear3(x))))

        # x = self.linear4(x)

        # x = torch.reshape(x, (-1, self.reshape_dim[0], self.reshape_dim[1]))
        return x
    