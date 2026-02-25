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
    def __init__(self, inshape, outshape, encoder_size = 1024):
        super(Decoder_pos, self).__init__()

        self.input_shp = inshape
        self.output_shp = outshape
        self.show_plots = False
        self.encoder_size = encoder_size

        self.linear1 = torch.nn.Linear(inshape[0], self.encoder_size)
        self.bn1 = torch.nn.BatchNorm1d(self.encoder_size)
        self.linear2 = torch.nn.Linear(self.encoder_size, outshape[0])
        self.bn2 = torch.nn.BatchNorm1d(outshape[0])
        self.linear3 = torch.nn.Linear(outshape[0], outshape[0])

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
