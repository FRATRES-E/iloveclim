# Getting started

## General Philosophy

### Install methodology

The installation of iLOVECLIM depends on a few packages to function correctly. For simplicity, these packages have been ported within the [LIPaS](https://github.com/dmr-dj/LIPaS) package installer. What is left for the user to do is to ensure the presence of a C/C++ compiler, a FORTRAN compiler and a netCDF library that is compiled with those two (both the C and the FORTRAN part of the netCDF).

Definition of the compilers and netCDF libraries is done through a machine-specific configuration file that can be found in LIPaS itself. With modern netCDF versions, it is easy to generate using the appropriate commands: nc-config nf-config

### Compiler choice
The install process should rely on one version of a complete suite of netCDF library, C compiler and FORTRAN compiler. The LIPaS install may prompt you to use one or more of them. Whatever your choice, *be consistent*: if you start the installation with one, carry it throughout.

## Install process

At a glance, installation is three steps:

1. Install the [prerequisites](prerequisites.md) (compilers and libraries)
2. Obtain the source code [there](obtaining-the-code.md)
3. Configure and compile using the [instructions](installation-example.md)

## Next steps

Once the above steps are completed, you can proceed to [run a first run](first-run.md), to produce a minimal experiment experiment that can be checked against the reference run experiment.


