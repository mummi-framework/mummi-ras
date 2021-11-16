# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# This module simply lists the names of lipid and protein states
# ------------------------------------------------------------------------------
import numpy as np

# ------------------------------------------------------------------------------
# LIPID NAMES
# ------------------------------------------------------------------------------
MACRO_FIELDS_DENSITIES = np.array(['nInner_POPC', 'nInner_PAPC', 'nInner_POPE',
                                   'nInner_DIPE', 'nInner_DPSM', 'nInner_PAPS',
                                   'nInner_PAP6', 'nInner_CHOL', 'nOuter_POPC',
                                   'nOuter_PAPC', 'nOuter_POPE', 'nOuter_DIPE',
                                   'nOuter_DPSM', 'nOuter_CHOL'])
MACRO_FIELDS_ATOMS = np.array(['id', 'type', 'rx', 'ry', 'rz'])
MACRO_FIELDS_BONDS = np.array(['a', 'b', 'length'])

# for performance, assume (and check) the order of these fields in the file
MACRO_FIDXS_DENSITIES = np.array([3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
MACRO_FIDXS_ATOMS = np.array([0, 2, 4, 5, 6])
MACRO_FIDXS_BONDS = np.array([0, 1, 2])

# ------------------------------------------------------------------------------
# PROTEIN NAMES
# ------------------------------------------------------------------------------
MACRO_RAS4b_NAMES =  np.array(['RASa',   'RASb',   'RASc'])
MACRO_RAS4bF_NAMES = np.array(['mRASa',  'mRASb',  'zRASa'])
MACRO_RAF4b_NAMES =  np.array(['mRAFa',  'mRAFb',  'zRAFa'])

MACRO_RAS4a_NAMES =  np.array(['RAS4Aa', 'RAS4Ab', 'RAS4Ac'])
MACRO_RAF4a_NAMES =  np.array(['mRAF4Aa','mRAF4Ab','zRAF4Aa'])
MACRO_RAS4aF_NAMES = np.array(['mRAS4Aa','mRAS4Ab','zRAS4Aa'])

MACRO_RAS_NAMES  = np.concatenate((MACRO_RAS4a_NAMES,  MACRO_RAS4b_NAMES))
MACRO_RASF_NAMES = np.concatenate((MACRO_RAS4aF_NAMES, MACRO_RAS4bF_NAMES))
MACRO_RAF_NAMES  = np.concatenate((MACRO_RAF4a_NAMES,  MACRO_RAF4b_NAMES))

MACRO_RAS_ALL_NAMES = np.concatenate((MACRO_RAS_NAMES, MACRO_RASF_NAMES))
MACRO_PROTEIN_NAMES = np.concatenate((MACRO_RAS_NAMES, MACRO_RAF_NAMES, MACRO_RASF_NAMES))

# ------------------------------------------------------------------------------
# COMPLEX NAMES
# ------------------------------------------------------------------------------
MACRO_RAS4b_COMPLEX_NAMES = np.array(['mRASa_mRAFa', 'mRASb_mRAFb', 'zRASa_zRAFa'])
MACRO_RAS4a_COMPLEX_NAMES = np.array(['mRAS4Aa_mRAF4Aa', 'mRAS4Ab_mRAF4Ab', 'zRAS4Aa_zRAF4Aa'])

MACRO_COMPLEX_NAMES = np.concatenate((MACRO_RAS4b_COMPLEX_NAMES, MACRO_RAS4a_COMPLEX_NAMES))

# ------------------------------------------------------------------------------
