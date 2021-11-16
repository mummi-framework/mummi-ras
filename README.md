## MuMMI RAS v0.1
#### Released: Nov 16, 2021

<b>MuMMI RAS</b> is the application component of the <b>MuMMI framework</b> 
developed to create large-scale ML-driven multiscale simulation ensembles to 
study the interactions of RAS proteins and RAS-RAF protein complexes with lipid 
plasma membranes.

MuMMI framework was developed as part of the <b><i>Pilot2</b></i> project of the 
[Joint Design of Advanced Computing Solutions for Cancer](https://cbiit.cancer.gov/ncip/hpc/jdacs4c)
funded jointly by the [Department of Energy](http://www.doe.gov) (DOE) and the 
[National Cancer Institute](http://www.cancer.gov) (NCI).
 
The Pilot 2 project focuses on developing multiscale simulation models for 
understanding the interactions of the lipid plasma membrane with the 
[RAS and RAF](https://www.cancer.gov/research/key-initiatives/ras) proteins. The broad computational tool development 
aims of this pilot are: 
* Developing scalable multi-scale molecular dynamics code that will automatically switch between phase field, coarse-grained and all-atom simulations. 
* Developing scalable machine learning and predictive models of molecular simulations to: 
    * identify and quantify states from simulations 
    * identify events from simulations that can automatically signal change of resolution between phase field, coarse-grained and all-atom simulations 
    * aggregate information from the multi-resolution simulations to efficiently feedback to/from machine learning tools
* Integrate sparse information from experiments with simulation data

MuMMI RAS defines the specific functionalities needed for the various components 
and scales of a target multiscale simulation. The application components need to 
define the scales, how to read the corresponding data, how to perform ML-based 
selection, how to run the simulations, how to perform analysis, and how to perform 
feedback. This code uses several utilities made available through "MuMMI Core".

#### Publications

MuMMI framework is described in the following publications.

1. Bhatia et al. <b>Generalizable Coordination of Large Multiscale Ensembles: Challenges and Learnings at Scale</b>.
   In Proceedings of the ACM/IEEE International Conference for High Performance Computing, Networking, Storage and Analysis, SC '21, 
   Article No. 10, November 2021.
   [doi:10.1145/3458817.3476210](https://doi.org/10.1145/3458817.3476210).

2. Di Natale et al. <b>A Massively Parallel Infrastructure for Adaptive Multiscale Simulations: Modeling RAS Initiation Pathway for Cancer</b>.
   In Proceedings of the ACM/IEEE International Conference for High Performance Computing, Networking, Storage and Analysis, SC '19, Article No. 57, November 2019. 
   [doi:10.1145/3295500.3356197](https://doi.org/10.1145/3295500.3356197).
   <br/><b><i>Best Paper at SC 2019</i></b>.

3. Ingólfsson et al. <b>Machine Learning-driven Multiscale Modeling Reveals Lipid-Dependent Dynamics of RAS Signaling Protein</b>.
   Proceedings of the National Academy of Sciences (PNAS), accepted, 2021. [preprint](https://www.researchsquare.com/article/rs-50842/v1).

4. <b>Reciprocal Coupling of Coarse-Grained and All-Atom scales</b>. In preparation.


#### Installation
```
git clone https://github.com/mummi-framework/mummi-ras
cd mummi-ras
pip3 install .

export MUMMI_ROOT=/path/to/outputs
export MUMMI_CORE=/path/to/core/repo
export MUMMI_APP=/path/to/app/repo
export MUMMI_RESOURCES=/path/to/resources
```

##### The installaton process as described above installs the MuMMI framework. The simulation codes (`gridsim2d`, `ddcMD`, `AMBER`, `GROMACS`) are not included and are to be installed separately.

##### Spack installation. We are also working towards releasing the option of installing MuMMI and its dependencies through [spack](https://github.com/spack/spack).


#### Authors and Acknowledgements
MuMMI was developed at Lawrence Livermore National Laboratory, in collaboration 
with Los Alamos National Laboratory, Oak Ridge National Laboratory, and International Business Machines. A 
list of main contributors is given below.

* <b>LLNL</b>:
Harsh Bhatia, Francesco Di Natale, Helgi I Ingólfsson, Joseph Y Moon, Xiaohua Zhang, 
Joseph R Chavez, Fikret Aydin, Tomas Oppelstrup, Timothy S Carpenter, 
Shiv Sundaram (previously LLNL), Gautham Dharuman (previously LLNL), 
Dong H Ahn, Stephen Herbein, Tom Scogland, Peer-Timo Bremer, and James N Glosli.   

* <b>LANL</b>:
Chris Neale and Cesar Lopez

* <b>ORNL</b>:
Chris Stanley

* <b>IBM</b>:
Sara K Schumacher


MuMMI was funded by the Pilot2 project led by Dr. Fred Streitz (DOE) and 
Dr. Dwight Nissley (NIH). We acknowledge contributions from the entire
Pilot 2 team.

This work was performed under the auspices of the U.S. Department of Energy (DOE) by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344, Los Alamos National Laboratory (LANL) under Contract DE-AC5206NA25396, and Oak Ridge National Laboratory under Contract DE-AC05-00OR22725.

Contact: Lawrence Livermore National Laboratory, 7000 East Avenue, Livermore, CA 94550.

#### Contributing

Contributions may be made through pull requests and/or issues on github.

### License

MuMMI RAS is distributed under the terms of the MIT License.

Livermore Release Number: LLNL-CODE-827655
