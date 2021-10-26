# GETM

This repository contains a rewrite of the General Estuarine Transport Model (GETM) to modern Fortran.

When the first stable version is released (quite some time away) this will be the official version of GETM.

# Installing

You will need the [Anaconda Python distribution](https://www.anaconda.com/products/individual). On many systems that is already installed: try running `conda --version`.
If that fails, you may need to load an anaconda module first: try `module load anaconda` or `module load anaconda3`. If that still does not give you a working `conda` command,
you may want to install [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

Now make sure your conda environment is initialized:

```
conda init bash
```

This needs to be done just once, as it modifies your `.bashrc` that is sourced every time you login.
After this, restart your shell by logging out and back in.

Now obtain the repository with setups and scripts:

```
git clone --recursive git@github.com:BoldingBruggeman/getm-rewrite.git
cd getm-rewrite
conda env create -f environment.yml
conda activate pygetm
source ./install
```

# Using the model

The best place to start is the [`python/examples`](https://github.com/BoldingBruggeman/getm-rewrite/tree/devel/python/examples) directory with Jupyter Notebooks that demonstrate the functionality of the model:

```
cd python/examples
python -m jupyterlab
```

# Contributing

How to contribute to the development:

  1. Make a [fork](https://github.com/BoldingBruggeman/getm-rewrite/projects/4) - upper right -  of the repository to you private GitHub account(\*) 
  2. Make a [pull request](https://docs.github.com/en/free-pro-team@latest/github/collaborating-with-issues-and-pull-requests/about-pull-requests)

Note that all communication in relation to development of GETM is done via GitHub - using [tickets](https://github.com/BoldingBruggeman/tickets/issues), [code](https://github.com/BoldingBruggeman/getm-rewrite), [issues](https://github.com/BoldingBruggeman/getm-rewrite/issues), [projects](https://github.com/BoldingBruggeman/getm-rewrite/projects).


(\*) If you use a service other than GitHub for your daily work - please have a look [here](https://stackoverflow.com/questions/37672694/can-i-submit-a-pull-request-from-gitlab-com-to-github)

https://yarchive.net/comp/linux/collective_work_copyright.html
