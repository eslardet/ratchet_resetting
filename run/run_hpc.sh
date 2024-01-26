#!/bin/bash

src_dir=$HOME/ratchet_resetting/src

## Parameters
L=$1
a=$2
h=$3
D=$4
r=$5

x0=$6
dt=$7
dt_save=$8
t=$9
samples=${10}
seed=${11}


run_dir=$HOME/ratchet_resetting/simulation_data/h${h}_a${a}/D${D}/r${r}
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

echo ${L} >> 'inpar'
echo ${a} >> 'inpar'
echo ${h} >> 'inpar'
echo ${D} >> 'inpar'
echo ${r} >> 'inpar'
echo ${x0} >> 'inpar'
echo ${dt} >> 'inpar'
echo ${dt_save} >> 'inpar'
echo ${t} >> 'inpar'
echo ${samples} >> 'inpar'
echo ${seed} >> 'inpar'

time ${src_dir}/ratchet_reset_hpc inpar

echo "Simulation completed."
