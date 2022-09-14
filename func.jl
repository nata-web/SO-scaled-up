### Functions used in "SO_scaled_up_J.ipynb" and "SO_scaled_up.jl"

rf(n) = rand((-1,1), n)         # function that returns randomly -1 or 1 with equal probability
θ(x) = x >= 0 ? 1 : -1          # Heaviside threshold function

##.
"""Functions that generate different weight matrices """

function w_sparse(N::Int64, prt::Param, var::Variables)
    """Generates a random sparse symmetric matrix of w of size NxN with density k/N of nonzero connections"""
    A = sprand(N, N, prt.d, rf, Float64)
    var.w_or = Symmetric(A) # Symmetric view of the upper triangle of the matrix `A`
end

########################################################################

function w_modular(N::Int64, prt::Param, var::Variables)
    """Generates a modular consistent connectivity matrix of weights of size NxN where:
    - intramodule connections set to {-1,1} with equal probability
    - intermodule connections set to {-p,p} with equal probability
    - modules are of size k"""
    p = prt.p
    k = prt.k
    W = rand((-p,p),N,N)

    for i in 1:N
        for j in 1:N
            if floor((i-1)/k) == floor((j-1)/k)
                W[i,j] = rand((-1,1))
            end
        end
    end
    W = Matrix{Float64}(Symmetric(W))
end

##.
"""Simulation Functions"""

function Generate_state(N::Int, var::Variables)
    """Generate state of size Nx1 with {-1,1} randomly distributed"""
    var.state = rf(N)
    nothing
end

########################################################################

function Binary_update(N::Int64, var::Variables)
    """Eq. (1)"""
    idx = rand(1:length(var.state))
    oldState = var.state[idx]
    var.state[idx] = θ.(dot(var.w[idx,:], var.state)) # θ.(activation)

    return idx, oldState
end

########################################################################

function Binary_energy(N::Int64, step::Int64, reset::Int64, learning::Bool, idx::Int64, oldState::Int8, var::Variables, k::Int64)
    """ Energy computation for a state"""
    if step == 1
        """ Eq. (2) """
        E = - 0.5 * (transpose(var.state) * var.w_or * var.state)[1]
        if learning
            var.E_L_history[step, reset] = E
        else
            if(k==1)
                var.E_NL_history[step, reset] = E
            else
                var.E_NL2_history[step, reset] = E
            end
        end
    else
        """ Eq. (8) """
        dE = (oldState - var.state[idx]) *
            (dot(var.state, var.w_or[:,idx]) - var.state[idx]*var.w_or[idx,idx])
        if learning
            var.E_L_history[step, reset] = var.E_L_history[step-1, reset] + dE
        else
            if(k==1)
                var.E_NL_history[step, reset] = var.E_NL_history[step-1, reset] + dE
            else
                var.E_NL2_history[step, reset] = var.E_NL2_history[step-1, reset] + dE
            end
        end
    end
end

########################################################################
function learn(N::Int64, steps::Int64, reset::Int64, prt::Param, var::Variables, learning::Bool, num::Int64)
    """Run the dynamics with or without learning (Eq.3), the direct routine """
    Generate_state(N, var)  # Randomize initial discrete behaviours/states s_i={+-1}
    dw = zeros(Int8, N,N)

    for step in 1:steps
        idx, oldState = Binary_update(N, var)
        if learning
            """Implementation C in the Appendix"""
            if step == 1
                dw =  var.state * transpose(var.state)
            else
                if (var.state[idx] > 0)
                    dw[idx,:] = dw[:,idx] = var.state
                else
                    dw[idx,:] = dw[:,idx] = .- var.state
                end
            end
            var.w .+= dw # primarily memory bandwidth constraint

            # # Use upper triangle only
            #     if (var.state[idx] > 0)
            #         dw[1:idx,idx] = var.state[1:idx]
            #     else
            #         dw[1:idx,idx] = .- var.state[1:idx]
            #     end
            # end
            # for i in 1:N
            #     var.w[1:i,i] .+= dw[1:i,i]
            # end

        end
        Binary_energy(N, step, reset, learning, idx, oldState, var, num)
    end
end

function updateW(t::Int64, idx::Int64, idx2t::Vector{Int64}, t2idx::Vector{Int64}, t2state::Vector{Int64}, dw::Array{Int8}, var::Variables)
    w = var.w                           # pointer for shorter notation
    state = var.state                   # pointer

    w[:,idx] .+= dw[:,idx] .* (t-idx2t[idx])
    # new_dw = state[idx] .* t2state[idx2t[idx]+1:t-1]
    for i = idx2t[idx]+1:t-1
        new_dw = state[idx] * t2state[i]
        w[t2idx[i],idx] += (new_dw - dw[t2idx[i],idx]) * (t - i)
        dw[t2idx[i],idx] = new_dw
    end
end

function learnSpeed(N::Int64, steps::Int64, reset::Int64, prt::Param, var::Variables, learning::Bool, num::Int64)
    """Run the dynamics with or without learning (Algorithm 2), the on-the-fly routine """
    Generate_state(N, var)           # Randomize initial discrete behaviours/states s_i={+-1}
    w = var.w                        # pointer for shorter notation
    w_or = var.w_or                  # pointer
    state = var.state                # pointer

    idx2t = ones(Int64, N)           # Map between idx and time it was last changed ('zeros' in Python)
    t2idx = ones(Int64, steps)       # The idx for which at time t the state was changed ('zeros' in Python)
    t2state = ones(Int64, steps)     # History of all states after they were changed

    dw = zeros(Int8, (N,N))

    for t in 1:steps
        idx = rand(1:length(state))
        oldState = state[idx]

        if learning
            updateW(t, idx, idx2t, t2idx, t2state, dw, var)
            idx2t[idx] = t # at what time idx changed
            t2idx[t] = idx # what idx that was
        end

        state[idx] = θ.(dot(w[idx,:], state))

        if learning
            t2state[t] = state[idx] # save the state that was changed at time t

            if t==1
                dw =  state * transpose(state)
            else
                if state[idx] >= 0
                    dw[:,idx] = state
                else
                    dw[:,idx] = .- state
                end
            end
        end
        Binary_energy(N, t, reset, learning, idx, oldState, var, num)
    end

    if learning
        t = steps + 1
        for idx in 1:N
            updateW(t, idx, idx2t, t2idx, t2state, dw, var)
        end
    end
end

########################################################################

function learnJulia(N::Int64, steps::Int64, prt::Param, var::Variables, speed::Bool, learning::Bool, num::Int64)

    for reset in 1:prt.resets
        if speed
            learnSpeed(N, steps, reset, prt, var, learning, num)
        else
            learn(N, steps, reset, prt, var, learning, num)
        end
    end

    # # We update the upper triangle in learn(), so we need to copy it
    # # to lower triangle to make it symmetric again
    # if !speed
    #     for i in 1:N
    #         var.w[i+1:N,i] = var.w[i,i+1:N]
    #     end
    # end
end

########################################################################

# For Jupyter notebook
function simulate(N::Int64, eta::Int64, steps::Int64, prt::Param, var::Variables, speed::Bool, learning::Bool, num::Int64)
    learnJulia(N, steps, prt, var, speed, learning, num)
end

function simulate(N::Int64, eta::Int64, steps::Int64, prt::Param, var::Variables, speed::Bool, F::Bool, learning::Bool, num::Int64)

    if F # Call Fortran routines
        ccall((:__sofort_MOD_doseed,"./SOfortF.so"),Cvoid,(Ref{Int32},),var.sim_seed) # set seed

        energies = zeros(Float64, steps, prt.resets)

        if learning
            learnF = 1
        else
            learnF = 0
        end
        if speed
            speedF = 1
        else
            speedF = 0
        end

        ccall((:__sofort_MOD_run,"./SOfortF.so"),Cvoid,
        (Ref{Int64},Ref{Float64},Ref{Float64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64}),
        var.w, var.w_or, energies, Ref(steps), Ref(prt.resets), Ref(N), Ref(learnF), Ref(eta), Ref(speedF))

        if learning
            var.E_L_history = energies

        else
            if(num==1)
                var.E_NL_history = energies
            else
                var.E_NL2_history = energies
            end
        end
    else
        Random.seed!(var.sim_seed) # set seed for simulation
        learnJulia(N, steps, prt, var, speed, learning, num)
    end

end
##.
"""Plot Functions"""

function plot_all_6(prt::Param, var::Variables, num::Int64)

    subplot(231)
    ax1 = plt.gca()
    ax1.set_title("Initial W")
    im1 = ax1.matshow(var.w_or, aspect="auto")
    divider = axes_grid1.make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im1, cax=cax)

    subplot(232)
    ax2 = plt.gca()
    ax2.set_title("Final W")
    im2 = ax2.matshow(var.w, aspect="auto")
    divider = axes_grid1.make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im2, cax=cax)

    # E_range = Int(prt.resets/2):Int((prt.resets/2)/50):prt.resets
    E_range = (prt.resets-50):prt.resets

    subplot(233)
    ax3 = plt.gca()
    ax3.set_title("Energies without learning")
    steps = size(var.E_NL_history, 1)
    ax3.plot(1:steps, var.E_NL_history[:,E_range])
    ax3.set_xlabel(L"steps")
    ax3.set_ylabel(L"Energy")

    subplot(234)
    ax4 = plt.gca()
    scat_text = "Attractor states visited"
    ax4.set_title(scat_text)
    ax4.scatter(prt.ranges[1], var.E_NL_history[end,:], s=10, marker="o",c="red")
    ax4.scatter(prt.ranges[2], var.E_L_history[end,:], s=10, marker="o",c="lightblue")
    if(num==3)
        ax4.scatter(prt.ranges[3], var.E_NL2_history[end,:], s=10, marker="o",c="green")
    end
    ax4.set_xlabel(L"Relaxations")
    ax4.set_ylabel(L"Energy")

    subplot(235)
    ax5 = plt.gca()
    hist_text = "Histogram of attractor energies"
    ax5.set_title(hist_text)
    avg_fac = Int(prt.resets*num / 100) # num is the number of simulations
    if(num==3)
        energies = vcat(var.E_NL_history[end, :], var.E_L_history[end, :], var.E_NL2_history[end, :])
    else
        energies = vcat(var.E_NL_history[end, :], var.E_L_history[end, :])
    end
    eMin, eMax = minimum(energies), maximum(energies)

    avg_E = permutedims(reshape(energies, 100,avg_fac)) # since julia is column-major
    h = zeros(avg_fac,30)
    for i in 1:avg_fac
        h[i,:] = plt.hist(avg_E[i,:], bins=30, range=(eMin,eMax))[1];
    end
    ax5.clear()
    im5 = ax5.matshow(transpose(h), origin="lower", extent=[0,avg_fac,eMin,eMax],aspect="auto")
    divider = axes_grid1.make_axes_locatable(ax5)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im5, cax=cax)

    ax5.set_xlabel(L"Time window")
    ax5.set_ylabel(L"Energy")

    subplot(236)
    ax6 = plt.gca()
    ax6.set_title("Energies with learning")
    steps = size(var.E_L_history, 1)
    ax6.plot(1:steps, var.E_L_history[:,E_range])
    ax6.set_xlabel(L"steps")
    ax6.set_ylabel(L"Energy")

    PyPlot.tight_layout()
end
