#!/bin/bash

src_dir=$HOME/Code/ratchet_resetting/src_nd_mfpt

## Parameters
ell=4.0
kappa=1.0
delta=0.2

x0=2.0
dt=0.0001
samples=$((10 ** 4))
seed=1


run_dir=$HOME/Code/ratchet_resetting/simulation_data_mfpt/ell${ell}_kappa${kappa}_delta${delta}
if [ ! -d "$run_dir" ]; then
    mkdir -p $run_dir
    echo "Creating directory : $run_dir"
else
    echo "Directory : $run_dir"
fi

cd $run_dir

if [ -e "inpar" ]; then
    rm 'inpar'
fi
touch 'inpar'

echo ${ell} >> 'inpar'
echo ${kappa} >> 'inpar'
echo ${delta} >> 'inpar'
echo ${x0} >> 'inpar'
echo ${dt} >> 'inpar'
echo ${samples} >> 'inpar'
echo ${seed} >> 'inpar'

time ${src_dir}/ratchet_reset_mfpt inpar

printf "\nSimulation completed."
