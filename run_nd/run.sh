#!/bin/bash

src_dir=$HOME/Code/ratchet_resetting/src_nd

## Parameters
alpha=4.0
beta=1.0
gamma=0.2

x0=0
dt=0.001
dt_save=0.01
t=1.0
samples=$((10 ** 3))
seed=1


run_dir=$HOME/Code/ratchet_resetting/simulation_data/alpha${alpha}_beta${beta}_gamma${gamma}
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

echo ${alpha} >> 'inpar'
echo ${beta} >> 'inpar'
echo ${gamma} >> 'inpar'
echo ${x0} >> 'inpar'
echo ${dt} >> 'inpar'
echo ${dt_save} >> 'inpar'
echo ${t} >> 'inpar'
echo ${samples} >> 'inpar'
echo ${seed} >> 'inpar'

time ${src_dir}/ratchet_reset inpar

printf "\nSimulation completed."
