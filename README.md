# pygetm

This repository contains a rewrite of the General Estuarine Transport Model (GETM).
It is mostly written in Python; only performance-critical sectons of the code are implemented in Fortran.

## Installing

You will need the [Anaconda Python distribution](https://www.anaconda.com/products/individual). On many systems that is already installed: try running `conda --version`.
If that fails, you may need to load an anaconda module first: try `module load anaconda` or `module load anaconda3`. If that still does not give you a working `conda` command,
you may want to install [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

Before using conda for the very first time, you will need to initialize its environment:

```
conda init bash
```

If you are using a different shell than bash, replace `bash` with the name of your shell (see `conda init -h` for supported ones), or use
`conda init --all`.

This needs to be done just once, as it modifies your `.bashrc` that is sourced every time you login.
After this, restart your shell by logging out and back in.

### Installation with conda (currently Linux/Windows only)

To install or update pygetm:

```
conda install pygetm -c bolding-bruggeman -c conda-forge
```

### Manual build and install

If you need a customized version of pygetm, for instance, built with specific compiler options, or with specific biogeochemical models that are not part of the standard [FABM](https://fabm.net) distribution, you can manually obtain the pygetm source code, build it, and then install it.

#### Linux/Mac

To obtain the repository with setups and scripts, set up your conda environment, and build and install pygetm:

```
git clone --recursive https://github.com/BoldingBruggeman/getm-rewrite.git
cd getm-rewrite
conda env create -f environment.yml
conda activate pygetm
source ./install
```

If you are using a different shell than bash, you may need to replace `source` in the last line  by `bash`. If you are installing on an HPC system that already has a Fortran compiler and MPI libraries that you would like to use, replace `environment.yml` with `environment-min.yml` in the above.

To customize the build step, you typically edit [the install script](https://github.com/BoldingBruggeman/getm-rewrite/blob/devel/install). For instance, you can:
* Specify a specific Fortran compiler to use by adding `-DCMAKE_Fortran_COMPILER=<EXE>` in the first call to cmake
* Customize compiler flags by adding `export FFLAGS=<FLAGS>` before the first call to cmake
* Configure FABM by adding options (e.g., `-DFABM_INSTITUTES`, `-DFABM_<INSTITUTE>_BASE`) to the first call to cmake

#### Windows

As on other platforms, you need [Anaconda](https://www.anaconda.com/products/individual) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html). In addition, you need to ensure that software to obtain and build Fortran code is available. Therefore, install:

* [Git for Windows](https://git-scm.com/download/win)
* [Visual Studio Community 2019](https://my.visualstudio.com/Downloads?q=visual%20studio%202019&wt.mc_id=o~msft~vscom~older-downloads)
* [Intel Fortran Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran)
* [Microsoft MPI](https://www.microsoft.com/en-us/download/details.aspx?id=100593) - you need both the runtime library and the Software Development Kit

Now obtain the repository with setups and scripts, set up your conda environment, and build and install pygetm:

```
git clone --recursive https://github.com/BoldingBruggeman/getm-rewrite.git
cd getm-rewrite
conda env create -f environment-min.yml
conda activate pygetm

mkdir build
cd build
cmake ..\python
cmake --build . --config Release
cmake --install .
```

#### Staying up to date

To update this repository including its submodules (GOTM, FABM, etc.), make sure you are in the getm-rewrite directory and execute:

```
git pull
git submodule update --init --recursive
conda env update -f <ENVIRONMENT_YML>
conda activate pygetm
source ./install
```

In the above, replace `<ENVIRONMENT_YML>` with the name of the environment file you used previously: `environment.yml` for stand-alone conda environments, or `environment-min.yml` for a setup that uses the local MPI implementation and Fortran compiler.

## Using the model

You should always activate the correct Python environment before you use the model with `conda activate pygetm`.
This needs to be done any time you start a new shell.

### Jupyter Notebooks

The best place to start with the model is the [`python/examples`](https://github.com/BoldingBruggeman/getm-rewrite/tree/devel/python/examples) directory with Jupyter Notebooks that demonstrate the functionality of the model:

```
cd python/examples
python -m jupyterlab
```

### Simulations

Some of the original GETM test cases have been ported to pygetm:

* [north_sea](https://github.com/BoldingBruggeman/getm-rewrite/blob/devel/python/examples/north_sea_legacy.py) - including [an extended version](https://github.com/BoldingBruggeman/getm-rewrite/blob/devel/python/examples/north_sea_legacy.py) that shows new pygetm features such as command-line configurability.
* [box_spherical](https://github.com/BoldingBruggeman/getm-rewrite/blob/devel/python/examples/box_spherical.py)
* [seamount](https://github.com/BoldingBruggeman/getm-rewrite/blob/devel/python/examples/seamount.py)

To run a simulation:

```
python <RUNSCRIPT.py> [OPTIONS]
```

To run in parallel:

```
mpiexec -n <NCPUS> python <RUNSCRIPT.py> [OPTIONS]
```

## Generating optimal subdomain division

A tool is included to calculated an optimal subdomain division provided only a valid bathymetry file _topo.nc_ and the number of processes the setup is to be run on. The tool searches for a solution with the smallest subdomain size that still covers the entire calculation domain. Output of the command (a python pickle file) is used directly as an argument when the domain object is created in the Python run-script.

The command to generate the subdomain division is:
```bash
pygetm-subdiv optimize --legacy --pickle subdiv_7.pickle topo.nc 7
```
The calculated layout can be shown - and plotted - via:
```bash
pygetm-subdiv show subdiv_7.pickle --plot
```
Resulting in in a layout as:

<img src="./images/northsea_subdiv_7.png" alt="northsea_subdiv_7" style="zoom:72%;" />

The example is for the standard North Sea case provided by legacy GETM.

To use a given subdomain division at runtime the call to creating the domain object should be similar to:

```python
domain = pygetm.legacy.domain_from_topo(os.path.join(getm_setups_dir, 'NorthSea/Topo/NS6nm.v01.nc'), nlev=30, z0_const=0.001, tiling='subdiv.pickle')
```

Note the *tiling*-argument.

## Contributing

How to contribute to the development:

  1. Make a [fork](https://github.com/BoldingBruggeman/getm-rewrite/projects/4) - upper right -  of the repository to you private GitHub account(\*) 
  2. Make a [pull request](https://docs.github.com/en/free-pro-team@latest/github/collaborating-with-issues-and-pull-requests/about-pull-requests)

Note that all communication in relation to development of GETM is done via GitHub - using [tickets](https://github.com/BoldingBruggeman/tickets/issues), [code](https://github.com/BoldingBruggeman/getm-rewrite), [issues](https://github.com/BoldingBruggeman/getm-rewrite/issues), [projects](https://github.com/BoldingBruggeman/getm-rewrite/projects).


(\*) If you use a service other than GitHub for your daily work - please have a look [here](https://stackoverflow.com/questions/37672694/can-i-submit-a-pull-request-from-gitlab-com-to-github)

https://yarchive.net/comp/linux/collective_work_copyright.html
