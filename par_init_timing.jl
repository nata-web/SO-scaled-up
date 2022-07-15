# Parameters

struct Param
    N_sizes     	:: Vector{Int64}
    runs        	:: Int64
	ts_multip		:: Int64
    resets      	:: Int64
	steps			:: Vector{Int64}
    alphas      	:: Vector{Float64}
	d           	:: Float64
	k           	:: Int64
    p           	:: Float64
    w_seed     		:: Int32
    ranges      	:: Vector{UnitRange{Int64}}

    Param(input_dict :: Dict{String,Any}) = new(
        input_dict["N_sizes"],
		input_dict["runs"],
        input_dict["ts_multip"],
		input_dict["resets"],
        input_dict["N_sizes"] * input_dict["ts_multip"], # usually 10; steps should be 10 times larger than the network size
		input_dict["alphas"],
		input_dict["d"],
		input_dict["k"],
        input_dict["p"],
        input_dict["w_seed"],
        # [1:input_dict["resets"], input_dict["resets"]+1:input_dict["resets"]*2], # ranges
		[1:input_dict["resets"], input_dict["resets"]+1:input_dict["resets"]*2,
			input_dict["resets"]*2+1:input_dict["resets"]*3], # ranges
    )
end

# Variables

mutable struct Variables
	w_or			::	Array{Float64}
	w				::	Array{Int64}
	sim_seed 		:: 	Int32
	state			::	Array{Int8}
    E_L_history  	:: 	Array{Float64}
	E_NL_history  	:: 	Array{Float64}
	E_NL2_history  	:: 	Array{Float64}

    Variables(prt::Param, N::Int64, ts::Int64) = new(
		zeros(N, N), 							# w_or
		zeros(N, N), 							# w
		0, 										# sim_seed
		zeros(N),								# state
		zeros(ts, prt.resets),					# E_L_history
		zeros(ts, prt.resets),					# E_NL_history
		zeros(ts, prt.resets)					# E_NL2_history
	)

end
