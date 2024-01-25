#!/bin/bash

src_dir=$HOME/Code/Ratchet_Resetting/src

## Parameters
L=1
a=0.2
h=1
D=0.1
r=0.2

x0=0
dt=0.001
t=10
samples=$((10 ** 3))
seed=1


run_dir=$HOME/Code/Ratchet_Resetting/simulation_data/h${h}_a${a}/D${D}/r${r}/t${t}
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
echo ${t} >> 'inpar'
echo ${samples} >> 'inpar'
echo ${seed} >> 'inpar'

time ${src_dir}/ratchet_reset inpar

echo "Simulation completed."
