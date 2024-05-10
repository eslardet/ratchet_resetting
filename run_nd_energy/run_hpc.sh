#!/bin/bash

src_dir=$HOME/ratchet_resetting/src_nd_energy


## Parameters
ell=$1
kappa=$2
delta=$3

x0=$4
dt=$5
dt_save=$6
t=$7
samples=$8
seed=$9


run_dir=$HOME/ratchet_resetting/simulation_data/ell${ell}_kappa${kappa}_delta${delta}
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

time ${src_dir}/ratchet_reset_energy inpar

printf "\nSimulation completed."