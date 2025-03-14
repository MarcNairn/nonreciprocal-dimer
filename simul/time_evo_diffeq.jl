using DifferentialEquations
using PyPlot
using LaTeXStrings
using LinearAlgebra
using Statistics
using ProgressMeter
using PyCall
using FFTW
using DelimitedFiles

PyPlot.pygui(false)
println(PyCall.python)

# ENV["PYTHON"] = "C:\\Users\\admin\\Anaconda3\\envs\\data_vis\\python.exe"

cmap = pyimport("cmap")

lajolla = cmap.Colormap("crameri:lajolla").to_mpl();
balance = cmap.Colormap("cmocean:balance").to_mpl();


# Configure LaTeX rendering and fonts
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["mathpazo"], size=18)  # Base font size
PyPlot.matplotlib[:rc]("axes", titlesize=22)             # Axis title
PyPlot.matplotlib[:rc]("axes", labelsize=20)             # Axis labels
PyPlot.matplotlib[:rc]("xtick", labelsize=18)            # X-ticks
PyPlot.matplotlib[:rc]("ytick", labelsize=18)            # Y-ticks
PyPlot.matplotlib[:rc]("legend", fontsize=18)            # Legend
PyPlot.matplotlib[:rc]("figure", titlesize=24)           # Figure title
PyPlot.svg(true)
# LaTeX preamble packages
PyPlot.matplotlib[:rc]("text.latex", preamble="\\usepackage{amsmath}\\usepackage{amsfonts}\\usepackage{amssymb}\\usepackage{lmodern}")


""" 
#############################################################

                    USEFUL FUNCTIONS BELOW

#############################################################
"""


function compute_time_average(sol::ODESolution{}, variable_index::Int)
    """  
    ###################
    Compute time average of an ODE sol object using portion of data, requires indexed variables 
    ################### 
    """
    t = sol.t
    t_mid = 0.6*t[end]
    indices = findall(x -> x >= t_mid, t)
    vals = sol[variable_index, indices]
    return mean(vals), var(vals)
end


function compute_time_average(t::Vector{Float64}, series::Vector{Float64})
    """
    ################### 
    Compute time average of a generic series using portion of data, no indexing required
    ################### 
    """
    t_mid = 0.6 * maximum(t)
    indices = findall(x -> x >= t_mid, t)
    return mean(series[indices])
end

"""
################################
DYNAMICS 
################################
"""

function spin_cav_ode!(du, u, p, t)
    s1x, s2x, s1y, s2y, s1z, s2z, α_r, α_i = u
    Δ1, Δ2, Δ, g, φ, κ, η = p

    # Equations for spin 1
    du[1] = Δ1 * s1y - 2*g * s1z * α_i            # ds1x/dt
    du[3] = -Δ1 * s1x - 2*g * s1z * α_r            # ds1y/dt
    du[5] = 2g * (s1x * α_i + s1y * α_r)          # ds1z/dt

    # Equations for spin 2
    du[2] = Δ2 * s2y - 2g * s2z * (sin(φ)*α_r + cos(φ)*α_i)  # ds2x/dt
    du[4] = -Δ2 * s2x   - 2g * s2z * (cos(φ)*α_r - sin(φ)*α_i) # ds2y/dt
    du[6] = 2g * (s2x*(cos(φ)*α_i + sin(φ)*α_r) + s2y*(cos(φ)*α_r - sin(φ)*α_i)) # ds2z/dt

    # Cavity field equations
    du[7] = -Δ*α_i - (κ/2)*α_r - g/2*(s1y +  cos(φ)*(s2y + s2x)) + η  # dα_r/dt
    du[8] = Δ*α_r - (κ/2)*α_i - g/2*(s1x - sin(φ)*(s2y + s2x)) # dα_i/dt
    return nothing
end


function single_spin_cav_ode!(du, u, p, t)
    s1x, s1y, s1z, α_r, α_i = u
    Δ1, Δ2, Δ, g, φ, κ, η = p

    # Equations for spin 1
    du[1] = Δ1 * s1y - 2*g * s1z * α_i            # ds1x/dt
    du[2] = -Δ1 * s1x - 2*g * s1z * α_r            # ds1y/dt
    du[3] = 2g * (s1x * α_i + s1y * α_r)          # ds1z/dt

    # Cavity field equations
    du[4] = -Δ*α_i - (κ/2)*α_r - g/2*s1y + η  # dα_r/dt
    du[5] = Δ*α_r - (κ/2)*α_i - g/2*s1x  # dα_i/dt
    return nothing
end

function spinonly_ode!(du, u, p, t)
    s1x, s2x, s1y, s2y, s1z, s2z = u
    Δ1, Δ2, Δ, g, φ, κ, η = p

    # Derived parameters
    K = κ / 2
    Δ_plus = Δ * cos(φ) + K * sin(φ)
    Δ_minus = Δ * cos(φ) - K * sin(φ)
    K_plus = K * cos(φ) + Δ * sin(φ)
    K_minus = K * cos(φ) - Δ * sin(φ)
    𝒥 = g^2 / (Δ^2 + K^2)  # 𝒥 = g²/(Δ² + (κ/2)²)

    # Equations for s1x, s2x, s1y, s2y
    du[1] = Δ1 * s1y + s1z * 𝒥 * (Δ * s1y + K * s1x + K_plus * s2x + Δ_plus * s2y - 2 * K * η / g)
    du[2] = Δ2 * s2y + s2z * 𝒥 * (Δ * s2y + K * s2x + K_minus * s1x + Δ_plus * s1y - 2 * K_minus * η / g)
    du[3] = -Δ1 * s1x + s1z * 𝒥 * (-Δ * s1x + K * s1y - Δ_minus * s2x + K_plus * s2y + 2 * Δ * η / g)
    du[4] = -Δ2 * s2x + s2z * 𝒥 * (-Δ * s2x + K * s2y - Δ_plus * s1x + K_minus * s1y + 2 * Δ_plus * η / g)

    # Equations for s1z, s2z
    du[5] = -𝒥 * (
        -2 * η / g * (K * s1x - Δ * s1y) +
        K * (s1x^2 + s1y^2) +
        K_plus * s1x * s2x +
        Δ_plus * s1x * s2y -
        Δ_minus * s1y * s2x +
        K_plus * s1y * s2y
    )

    du[6] = -𝒥 * (
        -2 * η / g * (K_minus * s2x - Δ_plus * s2y) +
        K * (s2x^2 + s2y^2) +
        K_minus * s2x * s1x +
        Δ_plus * s2x * s1y -
        Δ_plus * s2y * s1x +
        K_minus * s2y * s1y
    )

    return nothing
end


function single_spinonly_ode!(du, u, p, t)
    s1x, s1y,s1z = u
    Δ1, Δ2, Δ, g, φ, κ, η = p #Keep with extra parameters just for simplicity (wont be used)

    # Derived parameters
    K = κ / 2

    𝒥 = g^2 / (Δ^2 + K^2)  # 𝒥 = g²/(Δ² + (κ/2)²)

    # Equations for s1x, s2x, s1y, s2y
    du[1] = Δ1 * s1y + s1z * 𝒥 * (Δ * s1y + K * s1x - 2 * K * η / g)
    du[2] = -Δ1 * s1x + s1z * 𝒥 * (-Δ * s1x + K * s1y + 2 * Δ * η / g)
    # Equations for s1z, s2z
    du[3] = -𝒥 * (-2 * η / g * (K * s1x - Δ * s1y) +
        K * (s1x^2 + s1y^2)
    )
    return nothing
end



""" 
##############################
MAGNETIZATION  
##############################
"""


function sz_phase_diagram(Δ1, Δ2, Δ, g, κ, φ_range, η_range; u0=[0,0,0,0,-1,-1], t_end=100000*g)
    # Parameter grid
    φ_vals = φ_range
    η_vals = η_range
    s1z_avg = zeros(length(η_vals), length(φ_vals))
    s2z_avg = similar(s1z_avg)
    s1z_vrc = zeros(length(η_vals), length(φ_vals))
    s2z_vrc = similar(s1z_vrc)
    # Progress meter
    prog = Progress(length(φ_vals)*length(η_vals), 1)

    # Main loop
    for (i, φ) in enumerate(φ_vals)
        for (j, η) in enumerate(η_vals)
            p = (Δ1, Δ2, Δ, g, φ, κ, η)
            prob = ODEProblem(spinonly_ode!, u0, (0.0, t_end), p)
            sol = solve(prob, Tsit5(), saveat=0.1, abstol=1e-8, reltol=1e-8)
            
            s1z_avg[j,i], s1z_vrc[j,i] = compute_time_average(sol, 5)
            s2z_avg[j,i], s2z_vrc[j,i] = compute_time_average(sol, 6)
            
            next!(prog)
        end
    end

    return φ_vals, η_vals, s1z_avg, s2z_avg, s1z_vrc, s2z_vrc
end



function plot_sz_phase_diagram(φ_vals, η_vals, s1z_avg, s2z_avg, s1z_vrc, s2z_vrc)
    # Create figure with subplots
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 8.5), dpi=300,
                            constrained_layout=true, sharex="col")

    # Calculate extent for imshow [xmin, xmax, ymin, ymax]
    extent = [φ_vals[1], φ_vals[end], η_vals[1]/g, η_vals[end]/g]

    # Top row plots
    ax1_top, ax2_top = axs[1, 1], axs[1, 2]

    # Plot s1z_avg
    im1 = ax1_top.imshow(s1z_avg, cmap=balance, vmin=-1, vmax=1,
                        extent=extent, aspect="auto", origin="lower")
    ax1_top.set_ylabel(L"\eta/g")
    # ax1_top.set_xlabel(L"\varphi")
    ax1_top.set_xticks([-pi/2, 0, pi/2])
    ax1_top.set_xticklabels([L"-\pi/2", L"0", L"\pi/2"])
    ax1_top.set_yticks([0.6, 0.8, 1.0, 1.2])
    ax1_top.set_ylim(0.5, 1.25)

    # Plot s2z_avg
    im2 = ax2_top.imshow(s2z_avg, cmap=balance, vmin=-1, vmax=1,
                        extent=extent, aspect="auto", origin="lower")
    # ax2_top.set_xlabel(L"\varphi")
    ax2_top.set_xticks([-pi/2, 0, pi/2])
    ax2_top.set_xticklabels([L"-\pi/2", L"0", L"\pi/2"])
    ax2_top.set_yticks([0.6, 0.8, 1.0, 1.2])
    ax2_top.set_yticklabels([])
    ax2_top.set_ylim(0.5, 1.25)

    # Add shared colorbar for top row
    cbar_top = fig.colorbar(im1, ax=axs[1, :], orientation="vertical", pad=0.01)
    cbar_top.ax.xaxis.set_ticks_position("top")
    cbar_top.ax.xaxis.set_label_position("top")
    cbar_top.set_label(L"\langle s_z\rangle", 
                labelpad=0,  
                rotation=0,
                horizontalalignment="center",
                verticalalignment="bottom")
    cbar_top.ax.yaxis.label.set_position((0.25, 1.0))
    # Bottom row plots
    ax1_bottom, ax2_bottom = axs[2, 1], axs[2, 2]

    # Plot s1z_vrc
    im3 = ax1_bottom.imshow(s1z_vrc, cmap=lajolla.reversed(), vmin=0,
                            extent=extent, aspect="auto", origin="lower")
    ax1_bottom.set_ylabel(L"\eta/g")
    ax1_bottom.set_xticks([-pi/2, 0, pi/2])
    ax1_bottom.set_xticklabels([L"-\pi/2", L"0", L"\pi/2"])
    ax1_bottom.set_yticks([0.6, 0.8, 1.0, 1.2])
    ax1_bottom.set_ylim(0.5, 1.25)

    # Plot s2z_vrc
    im4 = ax2_bottom.imshow(s2z_vrc, cmap=lajolla.reversed(), vmin=0,
                            extent=extent, aspect="auto", origin="lower")
    ax2_bottom.set_xlabel(L"\varphi")
    ax2_bottom.set_xticks([-pi/2, 0, pi/2])
    ax2_bottom.set_xticklabels([L"-\pi/2", L"0", L"\pi/2"])
    ax2_bottom.set_yticks([0.6, 0.8, 1.0, 1.2])
    ax2_bottom.set_yticklabels([])
    ax2_bottom.set_ylim(0.5, 1.25)

    # Add shared colorbar for bottom row
    cbar_bottom = fig.colorbar(im3, ax=axs[2, :], orientation="vertical", pad=0.01)
    cbar_bottom.ax.xaxis.set_ticks_position("top")
    cbar_bottom.ax.xaxis.set_label_position("top")
    cbar_bottom.set_label(L"\langle s_z^2\rangle", 
                labelpad=0,  
                rotation=0,
                horizontalalignment="center",
                verticalalignment="bottom")
    cbar_bottom.ax.yaxis.label.set_position((0.25, 1.0))

    return fig
end

""" 
######################################################
FOURIER SPECTRUM 
######################################################
"""

# Function to compute FFT spectrum
function compute_fft_spectrum(signal, Δt)
    N = length(signal)
    Fs = 1 / Δt
    signal_centered = signal .- mean(signal)  # Remove DC component
    fft_vals = fft(signal_centered)
    freqs = fftfreq(N, Fs)
    # Shift frequencies and power to center DC component
    freqs_shifted = fftshift(freqs)
    power_shifted = fftshift(abs2.(fft_vals)) ./ (N * Fs)  # Two-sided PSD
    return freqs_shifted, power_shifted
end


#######PHASE DIAGRAM OF THE SPECTRUM########
function compute_dominant_freq(signal, Δt)
    freqs, power = compute_fft_spectrum(signal, Δt)
    Fs = 1/Δt
    nyquist = Fs / 2
    
    # Exclude DC and frequencies above Nyquist
    valid = findall(x -> 0 < abs(x) ≤ nyquist, freqs)
    power_valid = power[valid]
    freqs_valid = freqs[valid]
    
    # Thresholding: Ignore peaks below 10% of max power
    threshold = 0.1 * maximum(power_valid)
    significant = power_valid .> threshold
    if !any(significant)
        return NaN  # No meaningful peak
    end
    
    _, idx = findmax(power_valid[significant])
    return abs(freqs_valid[significant][idx])
end


function frequency_phase_diagram(φ_range, η_range, g, κ; show_progress=true)
    s1z_peaks = zeros(length(η_range), length(φ_range))
    s2z_peaks = similar(s1z_peaks)
    
    prog = Progress(length(φ_range) * length(η_range), enabled=show_progress)
    
    for (i, φ) in enumerate(φ_range)
        for (j, η) in enumerate(η_range)
            # Solve ODE (adjust parameters as per your system's requirements)
            p = (0.0, 0.0, 0.0, g, φ, κ, η)
            prob = ODEProblem(spinonly_ode!, [0.0, 0.0, 0.0, 0.0, -1.0, -1.0], (0.0, 10000.0), p)
            sol = solve(prob, Tsit5(), saveat=0.1, abstol=1e-8, reltol=1e-8)
            
            Δt = sol.t[2] - sol.t[1]
            s1z_peaks[j, i] = compute_dominant_freq(sol[5, :], Δt)
            s2z_peaks[j, i] = compute_dominant_freq(sol[6, :], Δt)
            next!(prog)
        end
    end
    return s1z_peaks, s2z_peaks
end 


function plot_frequency_phase_diagram(s1z_peaks, s2z_peaks, φ_range, η_range)
    zmin = min(minimum(s1z_peaks), minimum(s2z_peaks))
    zmax = max(maximum(s1z_peaks), maximum(s2z_peaks))
    clims = (zmin, zmax)


    # Create the figure and subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.4, 4), dpi=300)

    # Plot the first heatmap
    im1 = ax1.imshow(s1z_peaks, extent=[minimum(φ_range), maximum(φ_range), minimum(η_range ./ g), maximum(η_range ./ g)],
                    origin="lower", aspect="auto", cmap=lajolla.reversed(), vmin=clims[1], vmax=clims[2])
    ax1.set_xlabel(L"\varphi")
    ax1.set_ylabel(L"\eta/g")
    ax1.set_title(L"s_1")
    ax1.set_xticks([-π/2, 0, π/2])
    ax1.set_xticklabels([L"-\pi/2", L"0", L"\pi/2"])
    ax1.set_yticks([0.5, 1, 1.5])
    ax1.set_yticklabels(["0.5", "1.0", "1.5"])

    # Plot the second heatmap
    im2 = ax2.imshow(s1z_peaks, extent=[minimum(φ_range), maximum(φ_range), minimum(η_range ./ g), maximum(η_range ./ g)],
                    origin="lower", aspect="auto", cmap=lajolla.reversed(), vmin=zmin, vmax=zmax)
    ax2.set_xlabel(L"\varphi")
    ax2.set_title(L"s_2")
    ax2.set_xticks([-π/2, 0, π/2])
    ax2.set_xticklabels([L"-\pi/2", L"0", L"\pi/2"])
    ax2.set_yticks([0.5, 1, 1.5])
    ax2.set_yticklabels(["0.5", "1.0", "1.5"])

    # Create an axis for the colorbar
    cax = fig.add_axes([0.87, 0.22, 0.02, 0.65])  # [left, bottom, width, height]

    # Add the colorbar to the dedicated axis
    cbar = fig.colorbar(im1, cax=cax, orientation="vertical")
    cbar.set_label(L"\omega^\mathrm{max}", 
                labelpad=-20,  
                rotation=0,
                horizontalalignment="center",
                verticalalignment="bottom")
    # Manually position the label at the top using axes coordinates
    cbar.ax.yaxis.label.set_position((0.5, 1.0))
    plt.tight_layout(rect=[0, 0, 0.9, 1])

    return fig
end

""" 
############################
RADIANCE WITNESS 
############################
"""

function compute_R(p::NamedTuple, adaga0::Int64, u0_two:: Vector{Float64}, u0_one:: Vector{Float64}, tspan::Float64, w_cavity::Bool)
    K = p.κ / 2
    denominator = im * p.Δ - K
    # For two spins

    if w_cavity==true
        u0_two = [u0_two..., sqrt(adaga0/2), sqrt(adaga0/2)]
        u0_one = [u0_one..., sqrt(adaga0/2), sqrt(adaga0/2)]
        prob_two = ODEProblem(spin_cav_ode!, u0_two, tspan, p)
        prob_one = ODEProblem(single_spin_cav_ode!, u0_one, tspan, p)
        sol_two = solve(prob_two, Tsit5(), saveat=tspan/100, maxiters=1e9, abstol=1e-9, reltol=1e-9)
        sol_one = solve(prob_one, Tsit5(), saveat=tspan/100, maxiters=1e9, abstol=1e-9, reltol=1e-9)

        # For two spins
        a_two = sol_two[7, :] .+ im * sol_two[8, :]
        adaga_two = abs2.(a_two)
        avg_two = compute_time_average(sol_two.t, adaga_two)

        # For one spin
        a_one = sol_one[4, :] .+ im * sol_one[5, :]
        adaga_one = abs2.(a_one)
        avg_one = compute_time_average(sol_one.t, adaga_one)
    else
        prob_two = ODEProblem(spinonly_ode!, u0_two, tspan, p)
        prob_one = ODEProblem(single_spinonly_ode!, u0_one, tspan, p)
        sol_two = solve(prob_two, Tsit5(), saveat=tspan/100,maxiters=1e9, abstol=1e-9, reltol=1e-9)
        sol_one = solve(prob_one, Tsit5(), saveat=tspan/100, maxiters=1e9, abstol=1e-9, reltol=1e-9)
        
        # For two spins
        s1x_two = sol_two[1, :]
        s2x_two = sol_two[2, :]
        s1y_two = sol_two[3, :]
        s2y_two = sol_two[4, :]
    
        a_two = @. ( ( g/2 * (im*(s1x_two + exp(-im*φ)*s2x_two) + (s1y_two + exp(-im*φ)*s2y_two)) - η ) / denominator )
        adaga_two = abs2.(a_two)
        avg_two = compute_time_average(sol_two.t, adaga_two)
    
        # For one spin    
        s1x_one = sol_one[1, :]
        s1y_one = sol_one[2, :]
    
        a_one = @. ( (g/2 * ((im*s1x_one + s1y_one)) - η ) / denominator )
        adaga_one = abs2.(a_one)
        avg_one = compute_time_average(sol_one.t, adaga_one)
    end
    # Compute R
    R = (avg_two - 2 * avg_one) / (2 * avg_one)
    return R
end


# function compute_R(p::NamedTuple, u0_two, u0_one, tspan, w_cavity::Bool)
#     K = p.κ / 2
#     denominator = im * p.Δ - K
    
#     if w_cavity
#         # Handle cavity case with extended initial conditions
#         avg_two = _solve_cavity_case(spin_cav_ode!, vcat(u0_two, [sqrt(0.5), sqrt(0.5)]), p, tspan)
#         avg_one = _solve_cavity_case(single_spin_cav_ode!, vcat(u0_one, [sqrt(0.5), sqrt(0.5)]), p, tspan)
#     else
#         # Precompute common terms for spin-only case 
#         g_2 = p.g / 2
#         φ_term = exp(-im * p.φ)
#         pre_g = g_2 / denominator
#         pre_η = p.η / denominator
        
#         # Solve spin-only cases
#         avg_two = _solve_spin_case(spinonly_ode!, u0_two, p, tspan, pre_g, pre_η, φ_term, true)
#         avg_one = _solve_spin_case(single_spinonly_ode!, u0_one, p, tspan, pre_g, pre_η, φ_term, false)
#     end
    
#     # Calculate final result
#     R = (avg_two - 2 * avg_one) / (2 * avg_one)
#     return R
# end

# # Helper function for cavity case calculations
# function _solve_cavity_case(ode_func, u0, p, tspan)
#     prob = ODEProblem(ode_func, u0, tspan, p)
#     sol = solve(prob, Tsit5(), saveat=0.1, abstol=1e-9, reltol=1e-9)

#     # Dynamically find cavity field indices (last two elements)
#     cavity_real_idx = length(u0) - 1
#     cavity_imag_idx = length(u0)
    
#     @views a = sol[cavity_real_idx, :] .+ im .* sol[cavity_imag_idx, :]
#     return compute_time_average(sol.t, abs2.(a))
# end

# # Helper function for spin-only case calculations
# function _solve_spin_case(ode_func, u0, p, tspan, pre_g, pre_η, φ_term, is_two_spin)
#     prob = ODEProblem(ode_func, u0, tspan, p)
#     sol = solve(prob, Tsit5(), saveat=0.1, abstol=1e-9, reltol=1e-9)
    
#     @views begin
#         if is_two_spin
#             s1x, s2x, s1y, s2y = sol[1,:], sol[2,:], sol[3,:], sol[4,:]
#             a = pre_g .* (im .* (s1x .+ φ_term .* s2x) .+ (s1y .+ φ_term .* s2y)) .- pre_η
#         else
#             s1x, s1y = sol[1,:], sol[2,:]
#             a = pre_g .* (im .* s1x .+ s1y) .- pre_η
#         end
#     end
    
#     return compute_time_average(sol.t, abs2.(a))
# end

"""
############################################
MAIN SCRIPT (NUMERICS)
############################################    
"""

# Parameters
Δ1 = 0
Δ2 = 0
Δ = 0
g = 0.5
κ = 1.0
φ_range = range(-π/2, π/2, 101)     # φ (x-axis)
η_range = range(0.5*g, 1.5*g, 101)  # η (y-axis)
# Generate phase diagram
φ_vals, η_vals, s1z_avg, s2z_avg, s1z_vrc, s2z_vrc = sz_phase_diagram(Δ1, Δ2, Δ, g, κ, φ_range, η_range)


#### Benchmarking single and cuopled spin, with and without cavity 


Δ1 = 0
Δ2 = 0
Δ = 0 
κ = 1.0
g = 0.5κ
φ = 0
η = g
t_end = 1000/κ
p = (Δ1=Δ1, Δ2=Δ2, Δ=Δ, g=g, φ=φ, κ=κ, η=η)
saveat = t_end/1000

u0c_coupled=[0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0, 0]
u0s_coupled=[0.0, 0.0, 0.0, 0.0, -1.0, -1.0]
u0c_single=[0.0, 0.0, -1.0, 0, 0]
u0s_single=[0.0, 0.0, -1.0]

function spin_cav_dyns(p::NamedTuple,t_end, saveat, u0c_coupled,u0s_coupled,u0c_single,u0s_single,
    maxiters=1e9, reltol=1e-9, abstol=1e-9)

    # Precompute terms used in post-processing
    K = p.κ / 2
    denominator = im * p.Δ - K
    pre_g = (p.g / 2) / denominator
    pre_η = p.η / denominator
    φ_term = exp(-im * p.φ)

    # Common solver options
    solver_opts = (maxiters=maxiters, reltol=reltol, abstol=abstol)

    # 1. Coupled Spins WITH Cavity
    prob_cav_coupled = ODEProblem(spin_cav_ode!, u0c_coupled, (0.0, t_end), p)
    sol_cav_coupled = solve(prob_cav_coupled, Tsit5(); saveat=saveat, solver_opts...)

    # 2. Coupled Spins WITHOUT Cavity
    prob_spin_coupled = ODEProblem(spinonly_ode!, u0s_coupled, (0.0, t_end), p)
    sol_spin_coupled = solve(prob_spin_coupled, Tsit5(); saveat=saveat, solver_opts...)

    # 3. Single Spin WITH Cavity
    prob_cav_single = ODEProblem(single_spin_cav_ode!, u0c_single, (0.0, t_end), p)
    sol_cav_single = solve(prob_cav_single, Tsit5(); saveat=saveat, solver_opts...)

    # 4. Single Spin WITHOUT Cavity
    prob_spin_single = ODEProblem(single_spinonly_ode!, u0s_single, (0.0, t_end), p)
    sol_spin_single = solve(prob_spin_single, Tsit5(); saveat=saveat, solver_opts...)

    # Calculate |a|² for all cases
    abs2_cav_coupled = abs2.(sol_cav_coupled[7,:] .+ im.*sol_cav_coupled[8,:])

    s1x, s2x = sol_spin_coupled[1,:], sol_spin_coupled[2,:]
    s1y, s2y = sol_spin_coupled[3,:], sol_spin_coupled[4,:]
    abs2_spin_coupled = abs2.(pre_g .* (im.*(s1x .+ φ_term.*s2x) .+ (s1y .+ φ_term.*s2y)) .- pre_η)

    abs2_cav_single = abs2.(sol_cav_single[4,:] .+ im.*sol_cav_single[5,:])

    s1x, s1y = sol_spin_single[1,:], sol_spin_single[2,:]
    abs2_spin_single = abs2.(pre_g .* (im.*s1x .+ s1y) .- pre_η)

    avg_twoc = compute_time_average(sol_spin_coupled.t, abs2_cav_coupled)
    avg_onec = compute_time_average(sol_spin_single.t, abs2_cav_single)
    R_cav = ( avg_twoc- 2 * avg_onec) / (2 * avg_onec)

    avg_twos = compute_time_average(sol_spin_coupled.t, abs2_spin_coupled)
    avg_ones = compute_time_average(sol_spin_single.t, abs2_spin_single)
    R_spin = ( avg_twos- 2 * avg_ones) / (2 * avg_ones)


    R_cav_str = string(round(R_cav, digits=3))
    R_spin_str = string(round(R_spin, digits=3))
    # Create plot
    fig = figure(figsize=(8,5))
    t_points = sol_cav_coupled.t

    plot(t_points, abs2_cav_coupled, lw=2.5, color="dodgerblue", label="Coupled spins (with cavity)")
    #plot(t_points, abs2_spin_coupled, lw=2.5, color="dodgerblue", linestyle="--", label="Coupled spins (spin-only)")
    plot(t_points, abs2_cav_single, lw=2.5, color="crimson", label="Single spin (with cavity)")
    #plot(t_points, abs2_spin_single, lw=2.5, color="crimson", linestyle="--", label="Single spin (spin-only)")

    anon1 = annotate(L"$\mathcal{R}_c = " * R_cav_str * L"$", 
    (0.5, 0.7), 
    xycoords="axes fraction",
    fontsize=18,
    color="black")

    anon2 = annotate(L"$\mathcal{R}_s = " * R_spin_str * L"$", 
    (0.5, 0.6), 
    xycoords="axes fraction",
    fontsize=18,
    color="black")

    xlabel(L"Time (1/$\kappa$)")
    ylabel(L"|\alpha|^2")
    legend()
    grid(true)
    tight_layout()

    return fig
end

fig = spin_cav_dyns(p, t_end, saveat, u0c_coupled,u0s_coupled,u0c_single,u0s_single)


compute_time_average(sol_cav_coupled.t, abs2_cav_coupled)
compute_time_average(sol_spin_coupled.t, abs2_spin_coupled)

avg_two = compute_time_average(sol_spin_coupled.t, abs2_cav_coupled)
avg_one = compute_time_average(sol_spin_single.t,  abs2_cav_single)

R = (avg_two - 2 * avg_one) / (2 * avg_one)


# Extract time series (using second half)
t = sol.t
signal1 = sol[5, :]
signal2 = sol[6, :]



# Compute spectra
Δt = t[2] - t[1]
freqs1, power1 = compute_fft_spectrum(signal1, Δt)
freqs2, power2 = compute_fft_spectrum(signal2, Δt)

fig = plot_time_and_freq_spectrum(t, signal1, freqs1, power1, signal2, freqs2, power2, g)

function plot_time_and_freq_spectrum(t, signal1, freqs1, power1, signal2, freqs2, power2, g)
    fig = figure(figsize=(12, 8))

    # First subplot (top-left)
    subplot(2, 2, 1)
    plot(t, signal1, lw=5, color="steelblue")
    xlabel(L"\text{Time}")
    ylabel(L"s_1^z")  # Using LaTeX math mode
    
    
    # Second subplot (top-right)
    subplot(2, 2, 2)
    plot(freqs1, normalize(power1), lw=5, color="steelblue")
    xlabel(L"\text{Frequency}")
    ylabel(L"S(\omega)")
    xlim(-g/5, g/5)
    
    # Third subplot (bottom-left)
    subplot(2, 2, 3)
    plot(t, signal2, lw=5, color="coral")
    xlabel(L"\text{Time}")
    ylabel(L"s_2^z")
    
    
    # Fourth subplot (bottom-right)
    subplot(2, 2, 4)
    plot(freqs2, normalize(power2), lw=5, color="coral")
    xlabel(L"\text{Frequency}")
    ylabel(L"S(\omega)")
    xlim(-g/5, g/5)
    
    tight_layout()
    return fig
end



# φ_range = range(-pi, pi, length=101)
# η_range = range(0.0*g, 1.5*g, length=101)  
# s1z_peaks, s2z_peaks = frequency_phase_diagram(φ_range, η_range,g, κ)
# savefig("spectrum.pdf")


Δ1 = 0
Δ2 = 0
Δ = 0 
κ = 1.0
g = 0.5κ
φ = 0
η = κ

u0_two = [0.0, 0.0, 0.0, 0.0, -1.0, -1.0]
u0_one = [0.0, 0.0, -1.0]
p = (Δ1=Δ1, Δ2=Δ2, Δ=Δ, g=g, φ=φ, κ=κ, η=η)

R = compute_R(p, 1, u0_two, u0_one, 100.0/κ, false)

φ_range = range(0, 2*pi, length=30)
η_range = range(0.1*g, 100*g, length=30)  
total_iterations = length(φ_range) * length(η_range)
prog = Progress(total_iterations, 1, "Computing R grid: ", 50)
R_grid = zeros(length(φ_range), length(η_range));
for (i, φ) in enumerate(φ_range)
    for (j, η) in enumerate(η_range)
        # Create new parameter set with current phi and eta
        current_p = merge(p,(φ=φ, η=η))
        # Compute R for this parameter set
        R = compute_R(current_p, 0, u0_two, u0_one, 5000.0/κ, true)
        R_grid[i, j] = R

        next!(prog)
    end
end


fig = figure(figsize=(8, 6));   
ax = gca();
extent = [minimum(φ_range)/pi, maximum(φ_range)/pi, minimum(η_range)/g, maximum(η_range)/g];
# Create the image plot with proper extent
img_plot = ax.imshow(R_grid', extent=extent,origin="lower",aspect="auto",cmap=lajolla.reversed());
# Add colorbar
cbar = colorbar(img_plot, ax=ax);
cbar.set_label(L"\mathcal{R}");

# Set labels and title
semilogy()
xlabel(L"\varphi")
ylabel(L"\eta/g")
tight_layout()
display(fig)

# open("R.dat", "w") do f
#     writedlm(f, R_grid)
# end

# open("R_wdetns.dat", "w") do f
#     writedlm(f, R_grid)
# end

# R_grid = readdlm("R.dat");
# R_grid = readdlm("R_wdetns.dat");



"""

####################################
PRECOMPUTED EVALUATION / INTEGRAL EVAL (IN PROGRESS)
####################################

""" 

# First, ensure trapezoidal integration is available
function trapz(t::AbstractVector, y::AbstractVector)
    integral = 0.0
    for i in 2:length(t)
        dt = t[i] - t[i-1]
        integral += (y[i] + y[i-1]) * dt / 2
    end
    return integral
end

# Modified helper functions to return time integrals instead of averages
function _solve_cavity_case_integral(ode_func, u0, p, tspan)
    prob = ODEProblem(ode_func, u0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=0.1, abstol=1e-9, reltol=1e-9)
    
    cavity_real_idx = length(u0) - 1
    cavity_imag_idx = length(u0)
    
    @views a = sol[cavity_real_idx, :] .+ im .* sol[cavity_imag_idx, :]
    return trapz(sol.t, abs2.(a))  # Return full time integral
end

function _solve_spin_case_integral(ode_func, u0, p, tspan, pre_g, pre_η, φ_term, is_two_spin)
    prob = ODEProblem(ode_func, u0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=0.1, abstol=1e-9, reltol=1e-9)
    
    @views begin
        if is_two_spin
            s1x, s2x, s1y, s2y = sol[1,:], sol[2,:], sol[3,:], sol[4,:]
            a = pre_g .* (im .* (s1x .+ φ_term .* s2x) .+ (s1y .+ φ_term .* s2y)) .- pre_η
        else
            s1x, s1y = sol[1,:], sol[2,:]
            a = pre_g .* (im .* s1x .+ s1y) .- pre_η
        end
    end
    
    return trapz(sol.t, abs2.(a))  # Return full time integral
end

# New functions to compute time integrals
function compute_integral_one(p::NamedTuple, u0_one, tspan, w_cavity::Bool)
    if w_cavity
        return _solve_cavity_case_integral(single_spin_cav_ode!, vcat(u0_one, [0.0, 0.0]), p, tspan)
    else
        K = p.κ / 2
        denominator = im * p.Δ - K
        g_2 = p.g / 2
        pre_g = g_2 / denominator
        pre_η = p.η / denominator
        return _solve_spin_case_integral(single_spinonly_ode!, u0_one, p, tspan, pre_g, pre_η, 0.0, false)
    end
end

function compute_integral_two(p::NamedTuple, u0_two, tspan, w_cavity::Bool)
    if w_cavity
        return _solve_cavity_case_integral(spin_cav_ode!, vcat(u0_two, [0.0, 0.0]), p, tspan)
    else
        K = p.κ / 2
        denominator = im * p.Δ - K
        g_2 = p.g / 2
        pre_g = g_2 / denominator
        pre_η = p.η / denominator
        φ_term = exp(-im * p.φ)
        return _solve_spin_case_integral(spinonly_ode!, u0_two, p, tspan, pre_g, pre_η, φ_term, true)
    end
end

# Main computation
function compute_NCC_integrated(integral_two::Float64, integral_one::Float64; epsilon=1e-12)
    numerator = integral_two - 2 * integral_one
    denominator = integral_two + 2 * integral_one + epsilon
    return numerator / denominator
end

# Parameter setup remains the same
Δ1 = 0
Δ2 = 0
Δ = 0 
g = 0.5*κ
κ = 1.0
φ = 0
η = 1

u0_two = [0.0, 0.0, 0.0, 0.0, -1.0, -1.0]
u0_one = [0.0, 0.0, -1.0]
p = (Δ1=Δ1, Δ2=Δ2, Δ=Δ, g=g, φ=φ, κ=κ, η=η)

φ_range = range(-pi, pi, length=101)
η_range = range(g, 10*g, length=101)
tspan = (0.0, 500.0/κ)  # Full time range for integration
w_cavity = true  # Set based on your configuration

# Precompute single-spin integrals for each η
integral_one_vec = zeros(length(η_range))
@showprogress "Precomputing single-spin integrals..." for (j, η) in enumerate(η_range)
    current_p = merge(p, (η=η, φ=0.0))  # φ irrelevant for single-spin case
    integral_one_vec[j] = compute_integral_one(current_p, u0_one, tspan, w_cavity)
end

# Compute full NCC integrated grid
NCC_integrated_grid = zeros(length(φ_range), length(η_range))
@showprogress "Computing integrated NCC..." for (i, φ) in enumerate(φ_range)
    for (j, η) in enumerate(η_range)
        current_p = merge(p, (φ=φ, η=η))
        integral_two = compute_integral_two(current_p, u0_two, tspan, w_cavity)
        NCC_integrated_grid[i, j] = compute_NCC_integrated(integral_two, integral_one_vec[j])
    end
end

# Plot the integrated NCC grid
fig = figure(figsize=(8, 6))
ax = gca()
extent = [minimum(φ_range), maximum(φ_range), minimum(η_range)/g, maximum(η_range)/g]
img_plot = ax.imshow(NCC_integrated_grid', extent=extent, origin="lower", aspect="auto", cmap=lajolla.reversed())
cbar = colorbar(img_plot, ax=ax)
cbar.set_label("NCC (integrated)")
xlabel(L"\varphi")
ylabel(L"\eta/g")
tight_layout()
display(fig)