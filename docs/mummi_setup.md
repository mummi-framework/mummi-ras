## Procedure for setting up MuMMI

##### Updated as of Sep 2022 (on `lassen`)

#### Summary of steps


1. Clone `mummi` repositories, and make your custom config file.
2. Load `mummi` environment (`mummi` and its dependencies).
2. Install `mummi` (again) in a virtual environment -- this is the copy you will use.
3. Test `mummi` installation, and enjoy!! :-)

#### Detailed steps

#### Step 1.

First, we need to download the `mummi` codebase.
```
git clone git@code.ornl.gov:admirral/mummi-core.git
git clone git@code.ornl.gov:admirral/mummi-ras.git
git clone git@code.ornl.gov:admirral/mummi_resources.git
```

Next, we create a config file in the user space. A template is provided with the `mummi-ras` repository. Please copy this file and make edits.
```
mkdir -p ~/.mummi
cp $MUMMI_APP/setup/config/config.mummi.sh ~/.mummi
```

Let's look at the template and make changes. Here, you need to provide the path to the code (where you cloned the three repositories). You also need to provide `MUMMI_ROOT` -- a location where you will run the simulation.

###### Note: It is important that you change the paths pointing to your directories, rather than shared ones (as is in the template).
```
cat ~/.mummi/config.mummi.sh

#!/bin/bash

# ------------------------------------------------------------------------------
# likely to be set for each run
# ------------------------------------------------------------------------------
# path where you want to run mummi (workspace)
export MUMMI_ROOT=/usr/workspace/mummiusr/mummi_demoroot_20220724

# ------------------------------------------------------------------------------
# likely to be set only once
# ------------------------------------------------------------------------------
# virtual environment name of your choice
# no need to change (unless you want to do a complete reinstall)
export MUMMI_VENV=mummi-ras-demo

# path where you clone mummi repositories
# no need to change (unless you move your repositories)
export MUMMI_RESOURCES=/usr/workspace/mummiusr/mummi_resources
export MUMMI_CORE=/usr/workspace/mummiusr/mummi_framework/mummi-core
export MUMMI_APP=/usr/workspace/mummiusr/mummi_framework/mummi-ras

# ------------------------------------------------------------------------------
```

#### Step 2.

Setting up `mummi` involves several steps. A single script can do it for you though.
```
source $MUMMI_APP/setup/setup.env.sh
```

Here is an overview of what this script does.
1. First, it sources your config file (`~/.mummi/config.mummi.sh`) to figure out the paths to code and root directories, and the name of the virtual env.
2. Next, it figures out the machine (`lassen` or `summit`) and loads some host-specific packages from `mummi-ras/setup/envs/env.lassen.sh` etc.
3. Then, it loads the required `spack` packages (location of `spack` is specified in the host file above).
4. After this, a virtual environment is created (if needed) and activated. This virtual env's name is given in your config file.

The first time you run through these steps, you should see errors to the effect that `mummi` was not found. This is okay. This leads to the next step!

#### Step 3.

With all other packages available, you're now ready to install `mummi` itself. There is a script for this as well.
```
$MUMMI_APP/setup/install/install_mummi_editable.sh
```

This script works entirely within your virtual environment and will install both `mummi-core` and `mummi-ras`. This installation is made in the editable mode, which means you can edit your code without needing to have to reinstall.

###### Note. Installing `mummi` yourself without this script may interfere with `spack` installations. So, please use this script.

#### Step 4.

Run the test for `mummi` installation.
```
$MUMMI_APP/tests/test-packaging.py
```


#### Next steps
You're ready to use `mummi`. Check out the instructions in `mummi_launch.md`.
