#!/bin/bash

#   This script runs a single simulation of a TMD nanoribbon and copies the results into the
#   directory containing the simulation data.

#   ARGUMENTS
#   TMD : Transition metal dichalcogenide
#   NX  : Length of the ribbon
#   NY  : Width of the ribbon
#   BETA: Inverse Temperature
#   U   : On-site Interaction
#   MU  : Chemical Potential
#   Frequency of recomputing the Green's function set to 4 sweeps, Trotter Error set to 0.125

N=3*$NX*$NY

src_file=simulation_eq_time

cd examples/data
mkdir ${NX}x${NY}
cd ${NX}x${NY}
mkdir BETA${BETA}-U$U-MU${MU}
cd BETA${BETA}-U$U-MU${MU}

cd ../../../..
make clean

make nx=$NX ny=${NY} beta=${BETA} dt_inv=8 source=$src_file

time ./simulation ${TMD} $U ${MU} 1024 256 2

cp -r tmp/* examples/data/${NX}x${NY}/BETA${BETA}-U$U-MU${MU}

make clean
