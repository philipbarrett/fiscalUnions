#!/bin/bash
#PBS -j oe
#PBS -V
#PBD -v np=25
#PBS -l procs=25
#PBS -o mpi_test.out
#PBS -e mpi_test.err

cd $PBS_O_WORKDIR

julia -p 25 parEg.jl 
