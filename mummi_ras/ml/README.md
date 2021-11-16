## ML-driven Sampling Module

To facilitate multiscale simulations, MuMMI uses a ML-driven sampling process 
to analyze given configurations and select the ones that are deemed of 
"most interest".

As the macro simulation proceeds, the resulting patches are provided for 
analysis to the `PatchSelector`, which projects these patches into a lower-dimensional 
space using a deep neural network and records the results. Once MuMMI workflow 
decides a new CG simulation is to be spawned, it makes a request to the 
`PatchSelector`, which responds with the patches of highest interest score.

A similar process is used for analyzing frames from CG trajectories to select 
the "most interesting" ones to spawn corresponding AA simulations. One key 
difference is that the `CGSelector` uses reaction coordinates for sampling, 
instead of a learned metric/latent space using deep learning.

The two types of selectors, therefore, provide two examples of sampling 
functionality in MuMMI --- deep learning as well as traditional statistics.

#### What is "most important"?

Toward the goal of sampling "most important" configurations, we use the notion 
of diversity sampling, which focuses on selecting the configurations that are 
most dissimilar from the ones selected previously. 

Here, we make use of the dynamic-importance sampling framework (DynIm)
([github.com/LLNL/dynim](https://github.com/LLNL/dynim)), which implements 
diversity sampling. Since initial release, we have also added hybrid sampling 
to provide additional features to users' toolbox.

Currently, the `PatchSelector` uses a hybrid sampling approach through 
multiple queues (with respect to different protein constellations), where 
diversity is maintained individually within each queue. The `CGSelector` offers 
a different form of hybrid sampling by allowing to add a "random" commponent 
to the sampling to also maintain the effect of the data distribution.

#### Publications

1. Bhatia et al. <b>Machine Learning Based Dynamic-Importance Sampling for Adaptive Multiscale Simulations</b>.
   <i>Nature Machine Intelligence</i>, vol. 3, pp. 401–409, May 2021. 
   [doi:10.1038/s42256-021-00327-w](https://doi.org/10.1038/s42256-021-00327-w).