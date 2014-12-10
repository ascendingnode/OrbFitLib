OrbFitLib
=========

A Python library for fitting orbits.

## Requirements

* Python2 2.7 or greater, or Python3 3.4 or greater
* [Cython 0.17 or greater](http://cython.org) Library to link the C++ and Python code
* [NumPy & SciPy](http://www.numpy.org) General purpose Python numerical & scientific libraries
* [emcee](http://dan.iel.fm/emcee/current) For Markov-Chain Monte Carlo PDF estimation
* [NAIF CSPICE N0065 or greater](http://naif.jpl.nasa.gov/naif/index.html) NASA JPL Ephemeris library; see below

## Installation

You can install the latest versions of the required Python libraries locally (just for you) with:

```bash
pip install --upgrade cython numpy scipy emcee --user
```

or globally (all users):

```bash
sudo pip install --upgrade cython numpy scipy emcee
```

Installing NumPy with pip will require some additional C libraries; [see here](http://www.scipy.org/install.html). On Ubuntu this can accomplished with:

```bash
sudo apt-get install libblas-dev liblapack-dev
```

### Download CSPICE 

To avoid having to edit the SimSpice.pyx file:

1. Download the appropriate version of CSPICE [from the NAIF website](http://naif.jpl.nasa.gov/naif/toolkit_C.html) into the PyNBody directory
2. Extract the tarball: `tar xzvf cspice.tar.Z`

## Build and Install

To build and install locally:

```bash
python setup.py build_ext install --user
```

To build and install globally:

```bash
sudo python setup.py build_ext install
```
