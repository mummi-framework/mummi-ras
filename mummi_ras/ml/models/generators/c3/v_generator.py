#!/usr/bin/env python3
#
# Copyright (c) 2025, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# @authors
#           Tim Hsu <hsu16@llnl.gov>
# ------------------------------------------------------------------------------

import torch
import numpy as np


class V_Generator:
    def __init__(self, seed: int | None = None):
        if seed is not None:
            self.rng_gen = torch.Generator()
            self.rng_gen.manual_seed(seed)
        else:
            self.rng_gen = None

    def generate(self, num_interp_points: int, num_samples_per_interp_point: int) -> torch.Tensor:
        stddev = torch.tensor([4, 4, 4, 10]).view(1, 1, 4)
        gaussian_noises = torch.randn(1, num_samples_per_interp_point, 4, generator=self.rng_gen) * stddev

        v_start = torch.tensor([9.0749, 8.0228, 14.9984, 14.0154])
        v_interp = 10 ** torch.linspace(0, np.log10(80), num_interp_points).view(-1, 1, 1)

        v = v_start.view(1, 1, 4) + v_interp + gaussian_noises
        v = v.view(-1, 4)
        return v[(v > v_start).all(dim=1)]
