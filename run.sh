#!/bin/bash

export LC_ALL="C"

# set parameters
N_a="12"
N_b="12"
L="60"

MaxOcc=$N_a

U="1"

MaxBondDim="512"
dir="./data/"
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export ITENSOR_USE_OMP=1


#r_vec="$(seq 0.05 0.025 0.5)"
t_vec="1 0.95 0.9 0.85 0.8"
r_vec="0.475 0.05"


for t in $t_vec
do
 for r in $r_vec
 do
 nohup nice -19 ./droplet $N_a $N_b $L $MaxOcc $t $U $r $MaxBondDim $dir > "out_L."$L"_N."$N_a"_N."$N_b"_t."$t"_r."$r &
 done
done
