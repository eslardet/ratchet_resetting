#!/bin/bash

src_dir=$HOME/Code/ratchet_resetting/src_nd_stop_start

## Parameters
ell=4.0
kappa=1.0
delta=0.2

x0=0
dt=0.0001
dt_save=10.0
t=20.0
samples=$((10 ** 3))
seed=1


run_dir=$HOME/Code/ratchet_resetting/simulation_data_stop_start/ell${ell}_kappa${kappa}_delta${delta}
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
echo ${dt_save} >> 'inpar'
echo ${t} >> 'inpar'
echo ${samples} >> 'inpar'
echo ${seed} >> 'inpar'

time ${src_dir}/ratchet_reset_stop_start inpar

printf "\nSimulation completed."
