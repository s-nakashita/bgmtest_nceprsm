# BGM for NCEP RSM

This repository contains 

- a program for generating perturbed member by adding rescaling perturbations to the unperturbed member (`addprtb.f90`)

- a program for calculating perturbation total energy between two states (`calcte.f90`)

- a driver for addprtb (`run_addprtb.sh`)

- a driver for calcte (`run_calcte.sh`)

## Prerequisites

- CMake (>=3.16)

- Fortran90 compiler (tested with gfortran12)

- The date module (common/date_module.f90) requires the [NCEPLIBS-w3emc](https://github.com/NOAA-EMC/NCEPLIBS-w3emc), whose older version is included in this repository (lib/src/w3lib-1.7).
If you use this older version, you need C compiler (tested with gcc12).

## Compile

To compile the programs, edit and run `compile.sh`.

## Test

Run following commands to obtain the sample perturbed 10 members:

```default
## the first cycle
./run_addprtb.sh 1
## the second cycle
./run_addprtb.sh 2
```

The resulting data will be contained in `rsm2rsm27_bgm/20220618{00,12}`.

The sample data is contained in `rsm2rsm27_bgmtest/20220618{00,12}` to compare the results.

If you run a simulation from these perturbed state, you can calculate the perturbation energy evolution using `run_calcte.sh`.