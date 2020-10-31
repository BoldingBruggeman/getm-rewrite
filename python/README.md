Here a Python driver for the GETM model is implemented.

To build the driver follow the steps below - be sure to have made a [fork](https://github.com/BoldingBruggeman/getm-rewrite) - upper right.

In the following it is assumed the code is installed in source/repos/GETM/getm. If not - the path given below must be changed accordingly.

First we need to check out the python development branch:
  1. cd source/repos/GETM/getm
  2. git checkout python

Then we prepare for building and do the build:

  1. cd source/repos/GETM
  2. mkdir -p build/getm\_python
  3. cd build/getm\_python
  4. cmake ../../repos/GETM/getm/python && make install

If the build was sucessful the advection test can be done:
  1. cd source/repos/GETM/getm/examples
  2. python test_parallel_advection.py -h
  3. mpiexec -n 4 python test_parallel_advection.py

If you reach here - all is working. If not please provide a ticket here:
https://github.com/BoldingBruggeman/tickets/issues
