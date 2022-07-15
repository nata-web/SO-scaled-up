### Benchmark the routines (in FORTRAN)
# Make sure your system has gfortran and prepare the SO_fort.F90 file
# before running this code

cd(@__DIR__)

"""Load libraries and additional files"""

using LaTeXStrings
using Random                # Used func: rand
using LinearAlgebra         # Used func: Symmetric, diagind
using ProgressMeter         # Used func: @showprogress
using StatsBase             # Used func: fit
using SparseArrays          # Used func: sprand
using BenchmarkTools        # Used func: @btime
using DelimitedFiles

include("par_init_timing.jl")
include("func.jl")

##.
""" Initialize """

mydict = Dict{String,Any}(
    # "N_sizes"       => 	Int.(range(100,10000,length=100))  , # network size
    "N_sizes"       => 	Int.(range(100,2000,length=20))  , # network size
    "runs"          => 	1			   , # number of runs per network
    "ts_multip"     =>  10             , # sets the size of steps
    # "resets" 		=> 	1000		       ,
    "resets" 		=> 	10		       ,
    "alphas" 		=>	1.00e-7 .* ones(20),
    "d"		  		=>	0.1		       , # sparse w: density of non-zero values
    "k"             =>  10             , # modelar w: size of module
    "p"             =>  0.1            , # modelar w: intermodule connections set to {-p,p} with equal probability
    "w_seed"		=>	123456   	   , # seed for specific w_orig
);

prt = Param(mydict);
modular = true
name() = if modular "modular" else "sparse" end

path_out = "./output/output_julia/times/"
if isdir(path_out) else mkpath(path_out) end

samples = 5
sim_seed = rand(1:2^16)
time_arr = zeros(Float64, 1, 4)
all_times = zeros(Float64, 3, samples+1+2)

@show ARGS                                  # takes command line arguments and
my_argument = parse(Int32,ARGS[1])          # parses them into Julia
N, ts, alpha = prt.N_sizes[my_argument], prt.steps[my_argument], prt.alphas[my_argument]
eta = round(Int, 1 / alpha)
time_arr[1,1] = N
all_times[:,1] .= N

var = Variables(prt, N, ts)
var.sim_seed = sim_seed

s = "Benchmark results for N = " * string(N) * " for " * string(prt.resets) * " reset/s"
file_name = "times_N" * string(N) * ".txt"
open(path_out*file_name, "w") do io
    write(io, s,'\n')
end

println("\nRunning simulation for ", name(), " matrix, N=", N, ", α=", alpha, " for ", prt.resets, " reset/s")

Random.seed!(prt.w_seed) # make sure you start with same weight matrix
if modular w_modular(N, prt, var) else w_sparse(N, prt, var) end
var.w = round.(Int, var.w_or .* eta)
Random.seed!(var.sim_seed) # set seed for simulation

BenchmarkTools.DEFAULT_PARAMETERS.samples = samples
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 330000 # If the execution time exceeds this limit, only 1 sample will be analyzed

# Benchmark without learning
learnF = 0
speedF = 1
all_times[1,2] = learnF
all_times[1,3] = speedF

s = "\nLearning = " * string(learnF) * ", speed = " * string(speedF) * "\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n"
println(s)
open(path_out*file_name, "a") do io
    write(io, s,'\n')
end

ccall((:__sofort_MOD_doseed,"./SOfortF.so"),Cvoid,(Ref{Int32},),var.sim_seed) # set seed
energies = zeros(Float64, ts, prt.resets)

b = @benchmarkable ccall((:__sofort_MOD_run,"./SOfortF.so"),Cvoid,
(Ref{Int64},Ref{Float64},Ref{Float64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64}),
var.w, var.w_or, energies, Ref(ts), Ref(prt.resets), Ref(N), Ref(learnF), Ref(eta), Ref(speedF)) samples=samples  evals=1

p = run(b)
io = IOBuffer()
show(io, "text/plain", p)
s = String(take!(io))
println(s)

open(path_out*file_name, "a") do io
    write(io, s,'\n')
end

time_arr[1,2] = minimum(p).time / 10^9 # convert ns to s
sample_size = length(p.times)
all_times[1,4:(sample_size+3)] .= p.times

# Benchmark the direct routine with learning
learnF = 1
speedF = 0
all_times[2,2] = learnF
all_times[2,3] = speedF

s = "\nLearning = " * string(learnF) * ", speed = " * string(speedF) * "\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n"
println(s)
open(path_out*file_name, "a") do io
    write(io, s,'\n')
end

ccall((:__sofort_MOD_doseed,"./SOfortF.so"),Cvoid,(Ref{Int32},),var.sim_seed) # set seed
energies = zeros(Float64, ts, prt.resets)

b = @benchmarkable ccall((:__sofort_MOD_run,"./SOfortF.so"),Cvoid,
(Ref{Int64},Ref{Float64},Ref{Float64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64}),
var.w, var.w_or, energies, Ref(ts), Ref(prt.resets), Ref(N), Ref(learnF), Ref(eta), Ref(speedF)) samples=samples  evals=1

p = run(b)
io = IOBuffer()
show(io, "text/plain", p)
s = String(take!(io))
println(s)

open(path_out*file_name, "a") do io
    write(io, s,'\n')
end

time_arr[1,3] = minimum(p).time / 10^9 # convert ns to s
sample_size = length(p.times)
all_times[2,4:(sample_size+3)] .= p.times

# Benchmark the on-the-fly routine with learning
learnF = 1
speedF = 1
all_times[3,2] = learnF
all_times[3,3] = speedF

Random.seed!(prt.w_seed) # make sure you start with same weight matrix
if modular w_modular(N, prt, var) else w_sparse(N, prt, var) end
var.w = round.(Int, var.w_or .* eta)
Random.seed!(var.sim_seed) # set seed for simulation

s = "\nLearning = " * string(learnF) * ", speed = " * string(speedF) * "\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n"
println(s)
open(path_out*file_name, "a") do io
    write(io, s,'\n')
end

ccall((:__sofort_MOD_doseed,"./SOfortF.so"),Cvoid,(Ref{Int32},),var.sim_seed) # set seed
energies = zeros(Float64, ts, prt.resets)


b = @benchmarkable ccall((:__sofort_MOD_run,"./SOfortF.so"),Cvoid,
(Ref{Int64},Ref{Float64},Ref{Float64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64}),
var.w, var.w_or, energies, Ref(ts), Ref(prt.resets), Ref(N), Ref(learnF), Ref(eta), Ref(speedF)) samples=samples  evals=1

p = run(b)
io = IOBuffer()
show(io, "text/plain", p)
s = String(take!(io))
println(s)

open(path_out*file_name, "a") do io
    write(io, s,'\n')
end

time_arr[1,4] = minimum(p).time / 10^9 # convert ns to s
sample_size = length(p.times)
all_times[3,4:(sample_size+3)] .= p.times

file_name2 = "times_array.csv"
file_name3 = "all_times.csv"

open(path_out*file_name2, "a") do io
    writedlm(io, time_arr, ',')
end

open(path_out*file_name3, "a") do io
    writedlm(io, all_times, ',')
end
