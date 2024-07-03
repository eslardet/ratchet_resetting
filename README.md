This repository contains code to run and analyse a stochastic resetting model mediated by a ratchet potential for the following paper: 

*Roberts, Connor, Emir Sezik, and Eloise Lardet. "Ratchet-mediated resetting: Current, efficiency, and exact solution." arXiv preprint arXiv:2405.10698 (2024).*

Please include the above citation if you reuse this code.

Once compiled, the C++ code in `src` can be run using shell scripts in the `run` folder, or run on the HPC/ Cluster using the bash scripts. Ensure that you compile the source code locally on the machine where you will be running the code from.
This program will generate position files at different time save points for phase `d` and phase `r` under the folder `simulation_data`.

Note that there are various different versions of the source code for different models and measurements. Most of the source code is run in C++, but Python was also used for small tests to generate the schematic plots. The majority of the analysis was performed for the non-dimensionalised model, but the original dimensional model is kept for comparison.
- `src` is a fully dimensional model in C++
- `src_nd` is the non-dimensionalised model in C++
- `src_nd_energy` is the non-dimensionalised model in C++ for energy calculations
- `src_nd_mfpt` is the non-dimensionalised model in C++ for mean first-passage time calculations
- `src_python` is the fully dimensional model in Python
- `src_nd` is the non-dimensionalised model in Python

The resulting text files can then be analysed using the python functions in the `analysis` or `analysis_nd` folder, which will generate image files in a folder called `plots`.
