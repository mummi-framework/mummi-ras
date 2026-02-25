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
# A decoder class
# ------------------------------------------------------------------------------
class Decoder_pos(nn.Module):

    # --------------------------------------------------------------------------
    def __init__(self, inshape, outshape):
        super(Decoder_pos, self).__init__()

        self.input_shp = inshape
        self.output_shp = outshape
        self.show_plots = False


        self.linear1 = torch.nn.Linear(inshape[0], 4395)
        self.bn1 = torch.nn.BatchNorm1d(4395, momentum=0.01)
        self.linear2 = torch.nn.Linear(4395, 8790)
        self.bn2 = torch.nn.BatchNorm1d(8790, momentum=0.01)
        self.linear3 = torch.nn.Linear(8790, 15328)
        self.bn3 = torch.nn.BatchNorm1d(15328, momentum=0.01)

        self.linear4 = torch.nn.Linear(15328+inshape[0], 8790)
        self.bn4 = torch.nn.BatchNorm1d(8790, momentum=0.01)
        self.linear5 = torch.nn.Linear(8790, 4395)
        self.bn5 = torch.nn.BatchNorm1d(4395, momentum=0.01)
        self.linear6 = torch.nn.Linear(4395, outshape)
        self.bn6 = torch.nn.BatchNorm1d(outshape, momentum=0.01)
        # self.linear4 = torch.nn.Linear(inshape[0], inshape[0])

        self.linear7 = torch.nn.Linear(outshape, outshape)

        # self.drop1 = torch.nn.Dropout(0.2)
        # self.drop2 = torch.nn.Dropout(0.2)
        # self.drop3 = torch.nn.Dropout(0.2)

        self.act = torch.nn.ReLU()
        # self.reshape_dim = outshape

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

        # x = self.drop1(self.bn1(self.act(self.linear1(x))))
        # x = self.drop2(self.bn2(self.act(self.linear2(x))))
        # x = self.drop3(self.bn3(self.act(self.linear3(x))))
        dist = (self.bn1(self.act(self.linear1(x))))
        dist = (self.bn2(self.act(self.linear2(dist))))
        dist = (self.bn3(self.act(self.linear3(dist))))

        new_x = torch.cat((x, dist),1)

        pos = (self.bn4(self.act(self.linear4(new_x))))
        pos = (self.bn5(self.act(self.linear5(pos))))
        pos = (self.bn6(self.act(self.linear6(pos))))

        # x = self.linear4(x)
        pos = self.linear7(pos)

        # x = torch.reshape(x, (-1, self.reshape_dim[0], self.reshape_dim[1]))
        return pos, dist
    