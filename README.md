[![Build Status](https://travis-ci.org/aopp-pred/rpe.svg?branch=master)](https://travis-ci.org/aopp-pred/rpe) [![DOI (paper)](https://img.shields.io/badge/DOI%20%28paper%29-10.5194%2Fgmd--10--2221--2017-blue.svg)](https://doi.org/10.5194/gmd-10-2221-2017) [![DOI](https://zenodo.org/badge/49064493.svg)](https://zenodo.org/badge/latestdoi/49064493)

# Reduced Precision Emulator

An emulator for reduced-precision floating-point calculations in Fortran.
Licensed under the Apache 2.0 license (http://www.apache.org/licenses/LICENSE-2.0.html).

* Documentation: http://rpe.readthedocs.io
* Mailing list:  https://groups.google.com/d/forum/aopp-pred-rpe


## Documentation

See the [online documentation](http://rpe.readthedocs.io) for a user guide,
developer guide and API reference. The source code for the documentation is
distributed along with rpe in the `doc/` directory.


## Citation

The emulator is described in a [GMD paper](https://doi.org/10.5194/gmd-10-2221-2017).
If you use the software for published research please cite this paper, as well as the
DOI for the software repository itself.


## Building the library

The library is built using GNU Make and a fortran compiler. The default action
when invoking `make` is to build the static library `lib/librpe.a` and the
Fortran module `modules/rp_emulator.mod`.

The code is tested with and set-up for building with GNU gfortran or Intel
Fortran (ifort) compilers. Which compiler to use is determined by the value
of the `F90` environment variable. If this variable is not set the Makefile
will use `gfortran`. The default build is a simple invocation of `make`:

    make

It is also possible to build a shared library `lib/librpe.so` using the
`shared` target of the Makefile:

    make shared

You can optionally specify a compiler on the command line:

    F90=ifort make

If you want to use a compiler other than `gfortran` or `ifort` you will
need to specify both the appropriate `F90` variable and the correct `FFLAGS`
for your compiler, bearing in mind the source code contains some lines over
132 characters in length. Compiler flags are also used to ensure the module
`rp_emulator.mod` is placed in the `modules/` directory in the source tree.

### Unified source

The source code for the emulator is split across several files, and makes use
of the C preprocessor to combine them during the build process. If you want to
generate a unified source file for ease of use you can use the `source` target
of the Makefile:

    make source

This will generate the file `src/processed/rp_emulator.f90` which can be
integrated into the source of other projects.

