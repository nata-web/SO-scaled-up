# Scaling up the self-optimization model by means of on-the-fly computation of weights

This repository contains the code used for the paper <b>"Scaling up the self-optimization model by means of on-the-fly computation of weights"</b> ([IEEE](https://ieeexplore.ieee.org/document/10022074), [arXiv](https://arxiv.org/abs/2211.01698)). The code is provided in two forms:

* Jupyter Notebooks using the Julia and Python languages. 
* Julia and Python code with a call to Fortran routine.

Julia and Python jupyter notebooks are mostly intended for visualization purposes of the self-optimization model with and without the 'on-the-fly' computation of weights for small networks. To run the code for large networks (N>400) it is recommended to use the "SO_scaled_up" files (in Julia or Python) with the FORTRAN routine enabled. 

## Content:

### Jupyter notebooks (without FORTRAN)
* **SO_scaled_up_P.ipynb** - Python notebook. All included.
* **SO_scaled_up_J.ipynb** - Julia notebook. Main file used for setting all the parameters and running the simulation (imports par_init.jl and func.jl)

### Code with FORTRAN
* **SO_scaled_up.jl** - Same as SO_scaled_up_J.ipynb, but with additional option to call for a Fortran routine.
* **SO_scaled_up.py** - Same as SO_scaled_up.jl, but for Python

### Other
* **timingF.jl** - Benchmark the functions with and without learning (with FORTRAN)
* **func.jl** - all the functions (used by Julia files)
* **par_init.jl** - initiazation (used by Julia files)
* **par_init_timing.jl** - initiazation (used by timingF.jl)
* **SO_logPlot.jl** - loglog plots of the execution times.
* **times/** - Contains the execution times presented in the paper

## Test the code
To see that the code produces the same result with and without the 'on-the-fly' computation of weights, run the code for the same simulation seed when speed=off in Julia (speed=0 in Python) and then again when speed=true (speed=1).


## To run the code with FORTRAN:
Make sure your system has gfortran and f2py for Python or just gfortran for Julia. Run the following commands before the execution of the python / julia code to compile the FORTRAN file:

**For Python:**

`f2py3 --f90flags="-fdefault-integer-8 -O3" -m soFortranF -c SO_fort.F90`

Or alternative optimized for AMD CPUs:

`f2py3 --f90flags="-fdefault-integer-8 -O3 -march=znver1 -mtune=znver1 -mfma -mavx2 -m3dnow -fomit-frame-pointer" -m soFortranF -c SO_fort.F90`

**For Julia:**

`gfortran SO_fort.F90 -o SOfortF.so -shared -fPIC -fdefault-integer-8 -O3 -g`

Or:

`gfortran SO_fort.F90 -o SO_fortF.so -shared -fPIC -fdefault-integer-8 -O3 -g -march=znver1 -mtune=znver1 -mfma -mavx2 -m3dnow -fomit-frame-pointer`

## 

If you have any questions, feel free to open an issue or send me an email: natalya.weber (at) oist.jp.

If you use our code for your research, consider citing our paper:
```
 @inproceedings{Weber_Koch_Froese_2022, 
 title = {Scaling up the self-optimization model by means of on-the-fly computation of weights}, 
 DOI = {10.1109/SSCI51031.2022.10022074}, 
 abstractNote = {The Self-Optimization (SO) model is a useful computational model for investigating self-organization in “soft” Artificial life (ALife) as it has been shown to be general enough to model various complex adaptive systems. So far, existing work has been done on relatively small network sizes, precluding the investigation of novel phenomena that might emerge from the complexity arising from large numbers of nodes interacting in interconnected networks. This work introduces a novel implementation of the SO model that scales as $mathcalO(N^2)$ with respect to the number of nodes N, and demonstrates the applicability of the SO model to networks with system sizes several orders of magnitude higher than previously was investigated. Removing the prohibitive computational cost of the naive $mathcalO(N^3)$ algorithm, our on-the-fly computation paves the way for investigating substantially larger system sizes, allowing for more variety and complexity in future studies.}, 
 booktitle = {2022 IEEE Symposium Series on Computational Intelligence (SSCI)}, 
 author = {Weber, Natalya and Koch, Werner and Froese, Tom}, 
 year = {2022}, 
 month = {Dec}, 
 pages = {1276–1282} 
 }
```
