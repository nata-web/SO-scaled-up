cd(@__DIR__)

"""Load libraries and additional files"""

using PyCall            # Used func: pyimport

# to avoid output from plots, uncomment the next two lines
matplotlib = pyimport("matplotlib")
matplotlib.use("Agg")

axes_grid1 = pyimport("mpl_toolkits.axes_grid1")

using PyPlot
using BSON
using LaTeXStrings
using Random            # Used func: rand
using LinearAlgebra     # Used func: Symmetric, diagind
using ProgressMeter     # Used func: @showprogress
using StatsBase         # Used func: fit
using SparseArrays      # Used func: sprand
using BenchmarkTools    # Used func: @btime
using DelimitedFiles

include("par_init.jl")
include("func.jl")

##.
""" Initialize """

mydict = Dict{String,Any}(
    "N"             => 	[10000]        , # network size
    "ts_multip"     =>  20             , # sets the size of steps
    "resets" 		=> 	2000		   ,
    "alpha" 		=>	1.00e-9        ,
    "d"		  		=>	0.1		       , # sparse w: density of non-zero values
    "k"             =>  400            , # modelar w: size of module
    "p"             =>  0.1            , # modelar w: intermodule connections set to {-p,p} with equal probability
    "w_seed"		=>	123456   	   , # seed for specific w_orig
);

prt = Param(mydict);
F = true                # call for fortran routine
speed = true            # speed up by not explicitly adding w=w+dw
modular = true
name() = if modular "modular" else "sparse" end

plots = true
output_files = true
output_bson = false
path_out = "./output/output_julia/"

runs = 1
# sim_seeds = rand(1:2^16, runs)
sim_seeds =[47835]

N = prt.N                     # pointer
ts = prt.steps                # pointer
alpha = prt.alpha             # pointer

var = Variables(prt, N, ts)

if modular
    folder_name = "Modular/" * "N_" * string(N) * "/k_" * string(prt.k) * "/"
else
    folder_name = "Sparse/" * "N_" * string(N) * "/d_" * string(prt.d) * "/"
end
if isdir(path_out * folder_name) else mkpath(path_out * folder_name) end

for run in 1:runs
    println("\nRUN ", run)
    println("\nStarting simulation for ", name(), " matrix, N=", N, ", α=", alpha, " , speed=", string(speed))
    var.sim_seed = sim_seeds[run]

    Random.seed!(prt.w_seed) # make sure you start with same weight matrix
    if modular
        var.w_or = w_modular(N, prt, var)
    else
        var.w_or = w_sparse(N, prt, var)
    end

    eta = round(Int, 1 / alpha)

    var.w = round.(Int, var.w_or .* eta)
    Random.seed!(var.sim_seed) # set seed for simulation

    # simulate without learning (Initial random disrtibution)
    println("\nSimulation without learning...")
    num = 1
    @time simulate(N, eta, ts, prt, var, speed, F, false, num) # without learning
    # simulate with learning
    println("\nSimulation with learning...")
    num += 1
    @time simulate(N, eta, ts, prt, var, speed, F, true, num) # with learning
    # simulate without learning (for testing)
    println("\nSimulation without learning...")
    num += 1
    @time simulate(N, eta, ts, prt, var, speed, F, false, num) # without learning

    if plots
        fig_text = "N_" * string(N) * "_eta_" * string(eta) * "_run_" * string(run)
        figure(fig_text, figsize = (14, 9))
        plot_all_6(prt, var, num)
        subplots_adjust(top = 0.9)

        if modular
            title_text =
                L"N = " * string(N) *", steps = " *string(ts) *", resets = " *string(prt.resets) *
                ", α = " *string(alpha) *", k = " *string(prt.k) *", p = " *string(prt.p) *
                L"$, seed_{w} = $" *string(prt.w_seed) *L"$, seed_{sim} = $" *string(var.sim_seed)
        else
            title_text =
                L"N = " *string(N) *", steps = " *string(ts) *", resets = " *string(prt.resets) *
                ", α = " *string(alpha) *", d = " *string(prt.d) *
                L"$, seed_{w} = $" *string(prt.w_seed) *L"$, seed_{sim} = $" *string(var.sim_seed)
        end
        suptitle(title_text, y = 0.99)
    end

    if output_files
        filename = "N_" *string(N) *"_eta_" *string(eta) *"_s_" *string(var.sim_seed)
        fig_name = filename * ".png"
        savefig(path_out * folder_name * fig_name)
        # fig_name = filename * ".svg"
        # savefig(path_out * folder_name * fig_name)
        if output_bson
            bson_name = filename * "_main.bson"
            type_name = name()
            w_or = var.w_or
            BSON.@save joinpath(path_out, folder_name, bson_name) params = N type_name alpha ts F speed prt w_or
        end
    end

end
