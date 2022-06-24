# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------

import os
import sys
from setuptools import setup, find_packages

# ------------------------------------------------------------------------------
# relevant paths
src_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
src_path = os.path.join(src_path, 'mummi_ras')

# ------------------------------------------------------------------------------
setup(
    name='mummi_ras',
    version='1.0.0',
    description='MuMMI for Pilot2.',
    url='https://github.com/mummi-framework/mummi-ras',

    author='Harsh Bhatia, Helgi I Ingólfsson, Francesco Di Natale, Joseph Moon, Joseph R Chavez',
    author_email='hbhatia@llnl.gov, ingolfsson1@llnl.gov, dinatale3@llnl.gov, moon15@llnl.gov, chavez35@llnl.gov',

    packages=find_packages(),
    entry_points={
      'console_scripts': [
          'mummi_createsim = mummi_ras.transformations.createsim.createsims:main',
          'mummi_backmapping = mummi_ras.transformations.backmapping.backmapping:main',
          'mummi_cganalysis = mummi_ras.online.cg.cganalysis:main',
          'mummi_aaanalysis = mummi_ras.online.aa.aa_analysis:main',
          'mummi_workflow = mummi_ras.scripts.run_wfmanager:main',
      ]
    },

    install_requires=['pip>=21.2.4',
                      'pytest>=6.2.4',
                      'numpy>=1.20.2',
                      'scipy>=1.6.3',
                      'parmed>=3.2.0',
                      'h5py>=2.8.0',
                      'keras==2.2.4',
                      'theano==1.0.4',
                      'matplotlib>=3.4.1',
                      'dynim @ git+https://github.com/LLNL/dynim@main',
                      'mummi-core @ git+https://github.com/mummi-framework/mummi-core.git@main'
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: System :: Distributed Computing",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.7",
    ],
)
