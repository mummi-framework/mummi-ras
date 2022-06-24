## MuMMI RAS v1.0
#### Released: Jun 24, 2022

This guide describes the installation of MuMMI RAS.


#### Step 0: Assumptions and requirements

1. MuMMI dependencies can be broadly categorized into two.
  * Workflow dependencies (e.g., `python`, `flux`, `dynim`, `keras`, etc.)
  * Application dependencies (e.g., `ddcMD`, `gromacs`, `amber`, etc.)

2. Whereas workflow dependencies are easy to provide, the application
dependencies include simulation and analysis tools, not all of which are currently
public or free. Some of these are not part of MuMMI and we assume these are
accessible to the user.

3. For collaborators on `lassen` and `summit`, we share all these installations
when we are allowed to. On these machines, we have created a working software stack
using [`spack`](https://github.com/spack/spack). In this repository, we include scripts (`mummi-ras/setup/envs`)
that allow using our shared software stack.

3. Again on `lassen` and on `summit`, we have created a small set of test data,
which can be used to launch MuMMI at small scales. This (and the larger dataset)
will be made public through NCI website. Until then, we can make this data available
upon request.


#### Step 1. Download the code

```
$ git clone -b main https://github.com/mummi-framework/mummi-core.git
$ git clone -b main https://github.com/mummi-framework/mummi-ras.git
$ git clone -b binned_sampling https://github.com/LLNL/dynim.git
```

**Note:** `dynim` is a requirement of MuMMI right now that is to be installed
explicitly by the user.


#### Step 2. Copy and edit config file

```
$ mkdir ~/.mummi
$ cp mummi-ras/setup/config/config.mummi.sh ~/.mummi/
```

MuMMI uses this config file to set some key environment variables. Please
copy the template file and edit to set:
* `MUMMI_APP` and `MUMMI_CORE`: the paths where you clone the code
* `MUMMI_ROOT`: the path/workspace where you want to run the simulation
* `MUMMI_RESOURCES`: additional resource files (please reach out to get access)
* `MUMMI_VENV`: virtual environment name of your choice


#### Step 3: Setup the MuMMI environment

As of now, this is the trickiest step of all. We are working towards
making the entire environment configuration available through `spack`.

For now, we provide a setup script that can be used to create this environment.
Here, we rely on the custom setup on `lassen` and `summit` (see Step 0) created
using `spack`.

```
$ source mummi-ras/setup/setup.env.sh
```

For folks who do not have access to our shared software stack, we recommend
installing the following dependencies through `spack` or from source.
* `ddcMD`:  https://github.com/LLNL/ddcMD
* `AMBER`: https://ambermd.org
* `GROMACS`: https://manual.gromacs.org/current/download.html
  * Note that we have a patch for gromacs installation for customization. To be open-sourced soon.
* `gridsim2d`: to be released shortly


#### Step 4: Install MuMMI

Once you run the setup script in Step 3, a virtual environment will be created
for you. Next, we will install MuMMI in this virtual environment.
```
$ pushd dynim; pip install .; popd

$ ./mummi-ras/setup/install/install_mummi_editable.sh
```


#### Step 5: Test installation

If all went well, you should be able to import `mummi_core` and `mummi_ras`
packages. Note the paths below pointing to the code repository, instead of an
install prefix because we installed in editable mode.

```
$ pushd ~
~ /usr/workspace/mummiusr/mummi_framework

$ python3
Python 3.7.3 (default, May 13 2021, 18:04:11)
[GCC 8.3.1 20190311 (Red Hat 8.3.1-3)] on linux
Type "help", "copyright", "credits" or "license" for more information.

>>> import mummi_core, mummi_ras
>>> print (mummi_core.__file__)
/usr/WS1/mummiusr/mummi_framework/mummi-core/mummi_core/__init__.py
>>> print (mummi_ras.__file__)
/usr/WS1/mummiusr/mummi_framework/mummi-ras/mummi_ras/__init__.py
>>> exit()

$ popd
```

You can also run a test to see if all submodules are properly installed.
```
$ ./mummi-ras/tests/test-packaging.py
```
