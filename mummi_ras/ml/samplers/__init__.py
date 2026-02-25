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
from .cg_sampler import CGSampler, plot_pca, plot_tsne, destandardize

# For C3
from .ls_sampler_c3 import LSSampler_C3 as LSSampler

KNOWN_INTERPOLATOR_INTERFACES = ["naive", "ot", "ot_feedback"]

def get_interpolator_interfaces():
    return KNOWN_INTERPOLATOR_INTERFACES

def get_interpolator(_):
    if _ == "naive":
        from .interpolators.naive_interpolator import NaiveInterpolator
        interface = NaiveInterpolator
    elif _ == "ot":
        from .interpolators.ot_interpolator import OTInterpolator
        interface = OTInterpolator
    elif _ == "ot_feedback":
        from .interpolators.feedback_interpolator import FeedbackOTInterpolator
        interface = FeedbackOTInterpolator
    else:
        raise ValueError(f"Invalid interpolator requested ({_})")

    return interface

