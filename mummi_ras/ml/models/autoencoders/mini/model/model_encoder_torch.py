#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2023, Lawrence Livermore National Security, LLC. All rights
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
        # 140513
        self.linear1 = torch.nn.Linear(inshape[0], 1024) #512 is also good
        self.bn1 = torch.nn.BatchNorm1d(1024)
        self.linear2 = torch.nn.Linear(1024, lantent_dims[0])
        self.bn2 = torch.nn.BatchNorm1d(lantent_dims[0])

        self.linear3 = torch.nn.Linear(lantent_dims[0], lantent_dims[0])

        # self.drop1 = torch.nn.Dropout(0.2)
        # self.drop2 = torch.nn.Dropout(0.2)
        # self.drop3 = torch.nn.Dropout(0.2)

        self.act = torch.nn.ReLU()

    def forward(self, x):

        x = (self.act(self.bn1(self.linear1(x))))
        x = (self.act(self.bn2(self.linear2(x))))

        x = self.linear3(x)

        return x
