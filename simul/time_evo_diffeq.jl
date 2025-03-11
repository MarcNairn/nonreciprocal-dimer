using DifferentialEquations
using PyPlot
using LaTeXStrings
using LinearAlgebra
using Statistics
using ProgressMeter
using PyCall
using FFTW


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

# function spin_ode!(du, u, p, t)
#     s1x, s2x, s1y, s2y, s1z, s2z, α_r, α_i = u
#     Δ1, Δ2, Δ, g, φ, κ, η = p

#     # Equations for spin 1
#     du[1] = Δ1 * s1y - g * s1z * α_i            # ds1x/dt
#     du[3] = -Δ1 * s1x - g * s1z * α_r            # ds1y/dt
#     du[5] = g * (s1x * α_i + s1y * α_r)          # ds1z/dt

#     # Equations for spin 2
#     du[2] = Δ2 * s2y - g * s2z * (sin(φ)*α_r + cos(φ)*α_i)  # ds2x/dt
#     du[4] = -Δ2 * s2x - g * s2z * (cos(φ)*α_r - sin(φ)*α_i) # ds2y/dt
#     du[6] = g * (s2x*(cos(φ)*α_i + sin(φ)*α_r) + s2y*(cos(φ)*α_r - sin(φ)*α_i)) # ds2z/dt

#     # Cavity field equations
#     du[7] = -Δ*α_i - (κ/2)*α_r - g*(sin(φ)*s2x - g*s1y - g*cos(φ)*s2y) + 2η  # dα_r/dt
#     du[8] = Δ*α_r - (κ/2)*α_i - g*s1x - g*cos(φ)*s2x + g*sin(φ)*s2y # dα_i/dt
#     return nothing
# end

""" 
#############################################################

                    USEFUL FUNCTIONS BELOW

#############################################################
"""


function compute_time_average(sol::Vector{Vector{Float64}}, variable_index::Int)
    """  
    ###################
    Compute time average of an ODE sol object using portion of data, requires indexed variables 
    ################### 
    """
    t = sol.t
    t_mid = 5*t[end]/10
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
    t_mid = 0.5 * maximum(t)
    indices = findall(x -> x >= t_mid, t)
    return mean(series[indices])
end

"""
################################
SPIN ONLY DYNAMICS 
################################
"""

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

function compute_R(p::NamedTuple, u0_two,u0_one, tspan)
    K = p.κ / 2
    denominator = im * p.Δ - K
    # For two spins

    prob_two = ODEProblem(spinonly_ode!, u0_two, tspan, p)
    sol_two = solve(prob_two, Tsit5(), saveat=0.1, abstol=1e-9, reltol=1e-9)

    s1x_two = sol_two[1, :]
    s2x_two = sol_two[2, :]
    s1y_two = sol_two[3, :]
    s2y_two = sol_two[4, :]

    a_two = @. ( ( g/2 * (im*(s1x_two + exp(-im*φ)*s2x_two) + (s1y_two + exp(-im*φ)*s2y_two)) - η ) / denominator )
    adaga_two = abs2.(a_two)
    avg_two = compute_time_average(sol_two.t, adaga_two)

    # For one spin
    prob_one = ODEProblem(single_spinonly_ode!, u0_one, tspan, p)
    sol_one = solve(prob_one, Tsit5(), saveat=0.1, abstol=1e-9, reltol=1e-9)

    s1x_one = sol_one[1, :]
    s1y_one = sol_one[2, :]

    a_one = @. ( (g/2 * ((im*s1x_one + s1y_one)) - η ) / denominator )
    adaga_one = abs2.(a_one)
    avg_one = compute_time_average(sol_one.t, adaga_one)

    # Compute R
    R = (avg_two - 2 * avg_one) / (2 * avg_one)
    return R
end


"""
############################################
MAIN SCRIPT (NUMERICS)
############################################    
"""

# Parameters
Δ1 = 0
Δ2 = 0
Δ = 0
g = 0.1
κ = 1.0
φ_range = range(-π/2, π/2, 101)     # φ (x-axis)
η_range = range(0.5*g, 1.5*g, 101)  # η (y-axis)


# Generate phase diagram
φ_vals, η_vals, s1z_avg, s2z_avg, s1z_vrc, s2z_vrc = sz_phase_diagram(Δ1, Δ2, Δ, g, κ, φ_range, η_range)


# Parameters for testing
Δ1 = 0
Δ2 = 0
Δ = 0
g = 0.1
κ = 1.0
φ = 0#  pi/8
η = 0.9*g
u0 = [0.0, 0.0, 0.0, 0.0, -1.0, -1.0]
t_end = 100000*g  # Reduced for testing

# Solve the system
p = (Δ1, Δ2, Δ, g, φ, κ, η)
prob = ODEProblem(spinonly_ode!, u0, (0.0, t_end), p)
sol = solve(prob, Tsit5(), saveat=0.1, abstol=1e-8, reltol=1e-8)

# Extract time series (using second half)
t = sol.t
signal1 = sol[5, :]
signal2 = sol[6, :]



# Compute spectra
Δt = t[2] - t[1]
freqs1, power1 = compute_fft_spectrum(signal1, Δt)
freqs2, power2 = compute_fft_spectrum(signal2, Δt)

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



φ_range = range(-pi, pi, length=101)
η_range = range(0.0*g, 1.5*g, length=101)  
s1z_peaks, s2z_peaks = frequency_phase_diagram(φ_range, η_range,g, κ)
savefig("spectrum.pdf")





Δ1 = 0.05
Δ2 = -0.05
g = 0.5

u0_two = [0.0, 0.0, 0.0, 0.0, -1.0, -1.0]
u0_one = [0.0, 0.0, -1.0]
p = (Δ1=Δ1, Δ2=Δ2, Δ=Δ, g=g, φ=φ, κ=κ, η=η)
R= compute_R(p, u0_two, u0_one, 1000000.0*p.g)


total_iterations = length(φ_range) * length(η_range)
prog = Progress(total_iterations, 1, "Computing R grid: ", 50)
R_grid = zeros(length(φ_range), length(η_range));
for (i, φ) in enumerate(φ_range)
    for (j, η) in enumerate(η_range)
        # Create new parameter set with current phi and eta
        current_p = merge(p,(φ=φ, η=η))
        # Compute R for this parameter set
        R = compute_R(current_p, u0_two, u0_one, 50000.0*current_p.g)
        R_grid[i, j] = R

        next!(prog)
    end
end

fig = figure(figsize=(8, 6));
ax = gca();
extent = [minimum(φ_range), maximum(φ_range), minimum(η_range)/g, maximum(η_range)/g];
# Create the image plot with proper extent
img_plot = ax.imshow(R_grid', extent=extent,origin="lower",aspect="auto",cmap=balance, vmin=-1, vmax=2);
# Add colorbar
cbar = colorbar(img_plot, ax=ax);
cbar.set_label(L"\mathcal{R}");

# Set labels and title
xlabel(L"\varphi")
ylabel(L"\eta/g")
tight_layout()
display(fig)