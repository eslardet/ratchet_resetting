This repository contains code to run and analysis a stochastic resetting model with a ratchet potential.

Once compiled, the C++ code in `src` can be run using shell scripts in the `run` folder, or run on the HPC/ Cluster using the bash scripts. 
This program will generate position files at different time save points for phase `d` and phase `r` under the folder `simulation_data`.

The final position files can then be analysed using the python functions in the `analysis` folder, which will generate image files in the `plots` folder.
