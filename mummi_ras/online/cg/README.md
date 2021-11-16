## Coarse-grained (CG) Module

The coarse-grained (CG) online module is responsible for the management of both
simulating and analyzing lipid-only, RAS, and RAS+RAF simulations run in
Martini using the ddcMD molecular dynamics code. The online analysis module
(`cganalysis`) is responsible for the following tasks:

1. Copy over all relevant simulation inputs to the local filesystem
2. Initialize a ddcMD simulation and monitor time-series output
3. In-situ process each timestepped frame of the running simulation collecting
raw data, calculating radial distribution functions (RDFs) for feedback, and
reporting success/failure of the simulation.

### Execution Requirements

The execution of the `cganalysis` requires at least 2 CPU-cores and a single
GPU in order to execute. The ddcMD code executes a lipid-protein system using
at least a single CPU-core and a single GPU, with the remaining processes
running on the remaining cores.