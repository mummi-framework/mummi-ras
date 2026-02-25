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
# A decoder class
# ------------------------------------------------------------------------------
class Decoder_pos(nn.Module):

    # --------------------------------------------------------------------------
    def __init__(self, inshape, outshape):
        super(Decoder_pos, self).__init__()

        self.input_shp = inshape
        self.output_shp = outshape
        self.show_plots = False


        self.linear1 = torch.nn.Linear(inshape[0], 1024) #1250
        self.bn1 = torch.nn.BatchNorm1d(1024) # momentum=0.01
        self.linear2 = torch.nn.Linear(1024,outshape)
        self.bn2 = torch.nn.BatchNorm1d(outshape)


        self.linear3 = torch.nn.Linear(outshape, outshape)

        # self.drop1 = torch.nn.Dropout(0.2)
        # self.drop2 = torch.nn.Dropout(0.2)

        self.act = torch.nn.ReLU()

    def _encode(self, data, training=False):
        pos, dis = self.forward(data)
        # return _ if training else _.detach()
        if training:
            return pos, dis
        else:
            pos = pos.detach()
            dis = dis.detach()
            return pos, dis

    def forward(self, x):

        dist = (self.act(self.bn1(self.linear1(x))))
        dist = (self.act(self.bn2(self.linear2(dist))))

        pos = self.linear3(dist)

        return pos, dist
