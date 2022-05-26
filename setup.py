# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
import sys
from os.path import abspath, split, join
from setuptools import setup, find_packages

'''
import platform
dir = '{}-{}-{}.{}'.format(platform.system().lower(), platform.machine(), platform.python_version_tuple()[0], platform.python_version_tuple()[1])
print dir
'''

# ------------------------------------------------------------------------------
# relevant paths
src_path = split(abspath(sys.argv[0]))[0]
src_path = join(src_path, 'mummi_ras')

# ------------------------------------------------------------------------------
setup(
    name='mummi_ras',
    description='MuMMI for Pilot2.',
    version='0.0.1',
    author='Harsh Bhatia, Helgi I Ingólfsson, Francesco Di Natale, Joseph Moon, Joseph R Chavez',
    author_email='hbhatia@llnl.gov, ingolfsson1@llnl.gov, dinatale3@llnl.gov, moon15@llnl.gov, chavez35@llnl.gov',
    url='https://github.com/mummi-framework/mummi-ras',
    entry_points={
      'console_scripts': [
          'mummi_createsim = mummi_ras.transformations.createsim.createsims:main',
          'mummi_cganalysis = mummi_ras.online.cg.cganalysis:main',
          'mummi_backmapping = mummi_ras.transformations.backmapping.backmapping:main',
          'mummi_aaanalysis = mummi_ras.online.aa.aa_analysis:main',
          'mummi_workflow = mummi_ras.scripts.run_wfmanager:main',
      ]
    },
    packages=find_packages(),
    install_requires=['pip>=21.2.4',
                      'pytest>=6.2.4',
                      'numpy>=1.20.2',
                      'scipy>=1.6.3',
                      #'scikit-learn>=1.0.1',
                      'parmed>=3.2.0',
                      'keras==2.2.4',
                      #'theano>=1.0.4',
                      #'h5py>=3.5.0',
                      'matplotlib>=3.4.3',
                      'dynim @ git+ssh://git@github.com/LLNL/dynim@main',
                      'mummi-core @ git+ssh://git@github.com/mummi-framework/mummi-core@main',
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
