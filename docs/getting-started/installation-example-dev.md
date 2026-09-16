# Installation

Lorem ipsum dolor sit amet, consectetur adipiscing elit. This page walks
through getting iLOVECLIM compiled and ready to run. Read it top to bottom the
first time; later you can jump straight to the section you need.

!!! warning "Version"
    These instructions apply to **iLOVECLIM vX.X**. For earlier versions, see
    the archived docs. *(Admonition boxes like this are how you flag caveats.)*

## Overview

Sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. A one-
paragraph orientation: what the build produces, roughly how long it takes, and
what the reader should have ready before starting. Ut enim ad minim veniam,
quis nostrud exercitation ullamco.

At a glance, installation is three steps:

1. Install the prerequisites (compilers and libraries).
2. Obtain the source code.
3. Configure and compile.

## Prerequisites

Duis aute irure dolor in reprehenderit. Before compiling you need a Fortran
compiler, an MPI implementation, and the NetCDF libraries.

### Compilers

Excepteur sint occaecat cupidatat non proident. iLOVECLIM is built with a
Fortran compiler; the following are known to work:

- `gfortran` — lorem ipsum dolor sit amet
- `ifort` / `ifx` — consectetur adipiscing elit

### Libraries

Sunt in culpa qui officia deserunt mollit. You will need:

| Library | Minimum version | Notes                          |
| ------- | --------------- | ------------------------------ |
| NetCDF  | X.X             | Lorem ipsum dolor              |
| MPI     | X.X             | OpenMPI or MPICH, sit amet     |

!!! tip
    On managed HPC systems these are usually provided as environment modules —
    see [HPC machines](machines/index.md) rather than installing by hand.

## Obtaining the code

Lorem ipsum dolor sit amet. Clone the repository from GitHub:

```bash
git clone https://github.com/FRATRES-E/iloveclim.git
cd iloveclim
```

!!! note
    Consectetur adipiscing elit — note anything about submodules, branches, or
    access here.

## Configuration

Ut labore et dolore magna aliqua. Set the paths and options the build needs.
The exact steps differ by machine, so pick your environment below.

=== "Generic Linux"

    ```bash
    export NETCDF_ROOT=/usr/local
    ./configure --with-netcdf=$NETCDF_ROOT
    ```

    Duis aute irure dolor in reprehenderit in voluptate.

=== "IRENE (TGCC)"

    ```bash
    module load gcc netcdf-fortran openmpi
    ./configure --preset=irene
    ```

    Excepteur sint occaecat cupidatat non proident.

=== "Jean-Zay (IDRIS)"

    ```bash
    module load gcc netcdf-fortran openmpi
    ./configure --preset=jean-zay
    ```

    Sunt in culpa qui officia deserunt mollit anim.

## Compiling

Lorem ipsum dolor sit amet, consectetur. Once configured, build the model:

```bash
make -j 4
```

The compiled executable appears in `bin/`. Sed do eiusmod tempor incididunt.

### Verifying the build

Ut enim ad minim veniam. Confirm the executable was produced and runs:

```bash
./bin/iloveclim --version
```

You should see the version string printed. If instead you see an error, jump to
[Troubleshooting](#troubleshooting) below.

## Troubleshooting

Quis nostrud exercitation ullamco laboris. Common issues:

??? question "`configure` cannot find NetCDF"
    Lorem ipsum dolor sit amet — check that `NETCDF_ROOT` points at the
    directory containing `include/netcdf.mod`. *(This is a collapsible
    admonition — handy for a long FAQ that would otherwise clutter the page.)*

??? question "Linker errors mentioning MPI symbols"
    Consectetur adipiscing elit — ensure the MPI module is loaded and that the
    compiler wrapper (`mpif90`) is on your `PATH`.

## Next steps

Duis aute irure dolor. With the model compiled, continue to
[Your first run](first-run.md) to launch a minimal experiment.
