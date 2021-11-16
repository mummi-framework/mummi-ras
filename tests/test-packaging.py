# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
#!/usr/bin/env python

# ------------------------------------------------------------------------------
# this is how we do the very basic tests we want
# ------------------------------------------------------------------------------
def test_module(module):
    import importlib
    print('\n> import ' + module)
    try:
        mod = importlib.import_module(module)
        print(mod)
        print(dir(mod))
        return True
    except ImportError as e:
        print(e)
        print(' ===> Error:: module {} not found'.format(module))
        return False


# ------------------------------------------------------------------------------
def test_packaging(mod, submods):
    f = test_module(mod)
    if not f:
        exit(1)

    failed = []
    for _submod in submods:
        f = test_module(f'{mod}.{_submod}')
        if not f:
            failed.append(f'{mod}.{_submod}')

    print('')
    if len(failed) == 0:
        print('successfully imported {}'.format(mod))
    else:
        print('failed to import modules:', failed)


# ------------------------------------------------------------------------------
if __name__ == '__main__':
    mname = 'mummi_ras'
    modules = ['datastructures', 'feedback', 'interfaces', 'ml', 'online',
               'parsers', 'transformations', 'workflow',
               'online.aa', 'online.cg']

    test_packaging(mname, modules)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
