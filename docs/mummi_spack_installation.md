## Procedure for spack installation of MuMMI

##### Important
This documentation is deprecated, we are now using a more modern Spack
installation: see https://code.ornl.gov/admirral/mummi-spack .

##### Updated as of Sep 2022 (on `lassen`)


#### Step 0. Download `spack`

We are using a custom version of `spack` from the following address.
```
git clone --branch mummi_sv0.18.1-13e6f87ef6 git@code.ornl.gov:admirral/mummi-ras.git
```



#### Step 1. Use the system modules as needed

Certain modules, such as `cuda` and `mpi` are available as system packages, and
are optimized by the vendor. We will simply use those rather than try to build
ourselves.

##### Step 1a. Load the system modules.

These modules are available on `lassen` as of Sep 07, 2022.
```
module load gcc/8.3.1
module load cmake/3.23.1
module load cuda/10.2.89
module load spectrum-mpi/2020.08.19
module load fftw/3.3.9
```

##### Step 1b. Ask `spack` to look for system modules

```
spack compiler find --scope site
spack external find --scope site

## sometimes, it helps to specify the name of the modules
spack external find --scope site  cuda
spack external find --scope site  spectrum-mpi
spack external find --scope site  fftw
```

##### Step 1c. Configure `spack` to not build these modules.

With the command above, system packages are added to the file
`$SPACK_ROOT/etc/spack/packages.yaml`. Next, we want to edit this file, and add
a flag to ensure `spack` never tries to build these. The snippet below shows the
format of this file. Manually add `buildable=False` to these packages.

```
cuda:
    buildable: False
    externals:
       spec: cuda@10.2.89
       prefix: /usr/tce/packages/cuda/cuda-10.2.89

spectrum-mpi:
    buildable: False
    externals:
    - spec: spectrum-mpi@10.3.1.03rtm0%gcc@8.3.1
      prefix: /usr/tce/packages/spectrum-mpi/spectrum-mpi-2020.08.19-gcc-8.3.1
    - spec: spectrum-mpi@10.3.1.03rtm0%xl@16.1
      prefix: /usr/tcetmp

fftw:
     buildable: False
     externals:
     - spec: fftw@3.3.9
       prefix: /usr/tcetmp/packages/fftw/fftw-3.3.9
     - spec: fftw@3.3.8
       prefix: /usr/tcetmp/packages/fftw/fftw-3.3.8
```

##### Step 1d. "Install" system modules

If the above configuration is correct, the following command will simply link
the system modules in `spack`. This is needed only as a sanity check.

```
spack install cuda fftw spectrum-mpi
```


#### Step 2. Install MuMMI and all its dependencies

The branch of `spack` you downloaded above contains the changes needed to install
`mummi`.

##### Step 2a. Install `mummi-ras`
The following command will install most of the required packages. It
can take some time though.

```
spack install  mummi-ras  ^python@3.7.3  ^py-scipy@1.6.3  %gcc@8.3.1
```

Note here that we are using specific versions of `gcc` and `python`.
* As of Sep 2022, `python@3.7.3` is needed due to some changes in later versions
in how binary streams are handled. Let's stick to `3.7.3` for now.
* We are using most up-to-date `gcc` available.
* There are some issues in installing newer versions of `scipy`, so we are using
an older and tried out version.

##### Step 2b. Install additional simulation packages

Certain simulation packages are not installed as part of "`mummi-ras`", to lend
the flexibility of updating their versions without disturbing the rest as well as
to have a cleaner installation. You need to issue individual commands for the following. Note the "variants" and dependencies here.

First, we install `gromacs`
```
spack install  gromacs@2019.6.classic   +blas+lapack~cuda~mpi         build_type=Release
spack install  gromacs@2019.6.mummifix  +blas+lapack~cuda~mpi+double  build_type=Release
spack install  gromacs@2019.6.mummifix  +blas+lapack~cuda~mpi~double  build_type=Release
```


These are our custom packages. So we install them separately.

**Note:** Sep 2022, there are issues with `spack` installation of `gridsim2d`.
So we will install this separately (as in Step 3).
```
spack install  ddcmd@20210520  +cuda
spack install  gridsim2d@v2022-08-25-precampaign
```

If you need to add a new version, do `spack edit ddcmd` or `spack edit
gridsim2d`, and add a new line for the new version. As long as the build
procedure remains the same, adding the new version is all that is needed.


#### Step 3. Install non-spack packages

We install `amber` from source outside the `spack`. We simply need to specify the
path to executable in `$MUMMI_APP/setup/envs/env.lassen.sh` (or same file for
`summit`).
