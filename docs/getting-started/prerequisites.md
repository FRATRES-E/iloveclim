# Prerequisites

## Base aspects 

The base model language of the model is FORTRAN. The base requisite are thus a working FORTRAN / C compiler (three are regurlarly tested within the iLOVECLIM group and should work more or less readily: [Intel Fortran](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html), [GNU Fortran](https://gcc.gnu.org/fortran) and [NVIDIA HPC](https://developer.nvidia.com/hpc-sdk) compiling chain. 

Additionally, the input/output of the model relies on the netCDF format ; a working [netCDF library](https://www.unidata.ucar.edu/software/netcdf) (compiled in a version compatible with the FORTRAN / C compiler above mentioned) is thus also required. 

Optionaly, some aspects of iLOVECLIM can be run un parallel using [OpenMP](https://www.openmp.org) ; this is not a requirement but should be considered to speedup the computations (especially in experiments where you have a lot of oceanic tracers). 

!!! tip "HPC machines"
    Before doing the LIPaS installation steps, you want to check out the machine-specific pages to see if your computer is supported by default. The latter being the case, that will speedup your installation process tremendously.

## LIPaS Package installation manager

First, simply retrieve the LIPaS installer:

```bash
git clone https://github.com/dmr-dj/LIPaS
```
Then edit the configuration in configs where each directory is a machine name (or part of a machine name following regex) and each file defines a compiler-specific configuration file. Look at those present for inspiration. Once appropriately edited, you can run LIPaS once to see if this was correct:

```bash
./lipas.sh
```
(you can add the -v option to get verbose output).

The ask LIPaS to install the packages you need:

```bash
./lipas.sh -p semver,ncio,UUID-fortran,lapack,face
```

!!! warning "Version remark"
	These instructions apply to the iLOVECLIM version v1.3.0 and following. The list of above package given to LIPaS requires "face" since commit [c6968ad](https://github.com/FRATRES-E/iloveclim/commit/c6968ad0f713c3f1fdf070a7a36acd939c132eeb)

