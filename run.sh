#!/bin/bash
rm -r ./Results
rm -r ./DEM
rm -r ./CFD
rm ./error.log
mkdir CFD
mkdir DEM
mkdir Results
cd ./DEM
mkdir log
mkdir result
np=16
for((i=0;i<np;i++));do
  mkdir "proc"$i
done
cd ..
cd ./CFD
for((i=0;i<np;i++));do
  mkdir "proc"$i
done
mkdir tecplot
cd ..
make

mpirun -n $np ./ibm-dem -pc_type hypre -ksp_rTol 1e-10 #-mat_view #-ksp_monitor
