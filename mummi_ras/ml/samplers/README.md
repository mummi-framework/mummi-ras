# ML sampler module

This directory contains code for sampling the latent space.

## Overview

The sampler module provides two primary classes: 

* LSSampler
* LSInterpolator

The LSSampler implements the sampling algorithm by first drawing points uniformly from an interpolation schedule connecting sets A and B and then accepting or rejecting that new point based on a dynamically changing acceptance criterion.

The LSInterpolator module implements the interpolation scheme between sets A and B. 

Feedback to the importance sampling routine is implemented in the accpet_point function of LSSampler, which accepts or rejects a proposed point (draw uniformly over the interpolator set) based on a criterion which can be updated on the fly.  

If a user wants to alter the core algorithm used in of the default importance sampler, override the accept_point routine via a subclass.

## Usage Example

The following pieces of code show how the sampler can be used. In the code below, a path to a directory containing latent space training data is required. 


```python

import numpy as np
import pathlib
from mummi_ras.ml import ls_point as lsp
from mummi_ras.ml.samplers import ls_sampler
from mummi_ras.ml.samplers.interpolators import naive_interpolator as ni 

# first we need to read in the validation points from disk
training_data_files = pathlib.Path("/p/gpfs1/tempkin1/ls_sampling/training_data").glob("*")

setA = []
setB = []

for fn in training_data_files:
    data = np.load(fn)
    _a = data['labels'].astype(bool)
    _b = np.logical_not(_a)
    if len(setA) == 0:
        setA = data['x_test_encoded'][_a] 
    else: 
        setA = np.concatenate((setA, data['x_test_encoded'][_a]))

    if len(setB) == 0:
        setB = data['x_test_encoded'][_b] 
    else: 
        setB = np.concatenate((setB, data['x_test_encoded'][_b]))
    

source = [lsp.LSPoint(i, 
                      parent="A",
                      ls_coord=setA[i], 
                      valid=True) for i in range(setA.shape[0])]

target = [lsp.LSPoint(i + setA.shape[0], 
                      parent="B",
                      ls_coord=setB[i], 
                      valid=True) for i in range(setB.shape[0])]

# instantiate a naive interpolator
interpolator = ni.NaiveInterpolator(source, target)

# construct an LSSampelr by specifying the source and target sets and an interpolator
lss = ls_sampler.LSSampler(source, target, interpolator)

# draw 10 new points from the latent space
new_ls_points = lss.get_new_ls_points(10)

print("Done!")
```
