#!/usr/bin/env bash
set -e  # do not continue if any command fails

HELP_STRING="Usage: ./run.sh <num_mpi_ranks>"
if [[ $# -ne 1 ]]; then
  echo -e $HELP_STRING
  exit 0
fi

num_mpi_ranks=$1
num_mesh_parts=$(($num_mpi_ranks * 2))

cd grid;
echo $num_mesh_parts
mpirun -np $num_mesh_parts ../../build3/gridGen/buildGrid
cd ..
mpirun -np $num_mpi_ranks ../build3/driver/tioga.exe
