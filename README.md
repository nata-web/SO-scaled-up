# Large-scale networks for the self-optimization model enabled by on-the-fly computation of weights

This repository contains the code used for the paper [Large-scale networks for the self-optimization model enabled by on-the-fly computation of weights](http://url)(link TBA).

The code is provided in two forms:

* Jupyter Notebooks using the Julia and Python languages. 
* Julia and Python code with a call to Fortran routine.



## Content:

* **SO_scaled_up_J.ipynb** - jupyter notebook (Julia), used for setting all the parameters and running the simulation
* **SO_scaled_up_P.ipynb** - same as SO_scaled_up_J.ipynb, but for Python
* **SO_scaled_up.jl** - Same as SO_scaled_up_J.ipynb, but with additional option to call for a Fortran routine.
* **SO_scaled_up.py** - Same as SO_scaled_up.jl, but for Python
* **timingF.jl** - Benchmark the functions with and without learning (with FORTRAN)
* **par_init.jl** - initiazation (used by Julia files)
* **par_init_timing.jl** - initiazation (used by timingF.jl)
* **func.jl** - all the functions (used by Julia files)
* **SO_logPlot.jl** - loglog plots of the execution times.
* **times/** - Contains the execution times presented in the paper

## To run the code with FORTRAN:
Make sure your system has gfortran and f2py for Python or just gfortran for Julia. Run the following commands before the execution of the python / julia code to compile the FORTRAN file:

**For Python:**

f2py3 --f90flags="-fdefault-integer-8 -O3" -m SO_fortF -c SO_fort.F90 

Or alternative optimized for AMD CPUs:

f2py3 --f90flags="-fdefault-integer-8 -O3 -march=znver1 -mtune=znver1 -mfma -mavx2 -m3dnow -fomit-frame-pointer" -m SO_fortF -c SO_fort.F90 

**For Julia:**

gfortran SO_fort.F90 -o SOfortF.so -shared -fPIC -fdefault-integer-8 -O3 -g

Or:

gfortran SO_fort.F90 -o SO_fortF.so -shared -fPIC -fdefault-integer-8 -O3 -g -march=znver1 -mtune=znver1 -mfma -mavx2 -m3dnow -fomit-frame-pointer


If you have any questions, feel free to open an issue or send me an email: natalya.weber@oist.jp.
