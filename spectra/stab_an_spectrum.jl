using LinearAlgebra
using Plots
using LaTeXStrings
using ColorSchemes


default(
    fontfamily = "Computer Modern", # TeX's default serif font
    tickfont = font(12, "Computer Modern"),  # Tick labels font
    guidefont = font(14, "Computer Modern"), # Axis labels font
    legendfont = font(12, "Computer Modern"), # Legend font
    titlefont = font(16, "Computer Modern")  # Title font
)

# # Define the parameters
# ω0 = 1.0
# κ= 10*ω0 
# Δ = 20*ω0
# g = 2*κ
# χR = 1/sqrt(Δ^2 + κ^2/4)

# θl = atan(2Δ,κ)
# Δ1 = ω0
# Δ2 = ω0
# η = 5*κ


# # Function to compute the matrix M for a given φ
# function compute_M(φ)
#     θl = atan(2Δ, κ)
#     ξp = 2 * g^2 * abs(χR) * cos(θl)
#     ξm = 2 * g^2 * abs(χR) * sin(θl)
    
#     χp(m) = 2 * g^2 * abs(χR) * cos(θl + (-1)^m * φ)
#     χm(m) = 2 * g^2 * abs(χR) * sin(θl + (-1)^m * φ)
    
#     χ1p = χp(1)
#     χ1m = χm(1)
#     χ2p = χp(2)
#     χ2m = χm(2)
    
#     M = [
#         ξp   χ1p   -Δ1 + ξm   χ1m   0       0
#         χ2p  ξp    χ2m       -Δ2 + ξp   0       0
#         ξm + Δ1  χ1m   -ξp      -χ1p   0       0
#         χ2m  ξm + Δ2  -χ2p     -ξp    0       0
#         η * χ1m / g  0   -η * χ1m / g  0   0   0
#         0   η * χ2m / g  0   -η * χ2m / g  0   0
#     ]
#     return M
# end

# # Range of φ values
# φ_values = range(00, stop=π/2+0.1, length=10000)

# # Compute eigenvalues for each φ
# eigenvalues = [eigen(compute_M(φ)).values for φ in φ_values]

# # Extract real parts of eigenvalues
# real_parts = hcat([real.(eig) for eig in eigenvalues]...)
# im_parts = hcat([imag.(eig) for eig in eigenvalues]...)

# eigenvec_part = hcat([eigen(compute_M(φ)).vectors for φ in φ_values]...)


# eigen_p = [real(eig[3]) for eig in eigenvalues]
# eigen_m = [real(eig[4]) for eig in eigenvalues]

# # Plot the real parts of the eigenvalues
# plot(φ_values, real_parts', xlabel=L"\varphi", ylabel=L"\Re(\lambda)", legend=false, linewidth=3, ylims=(-1,1), color=[:lightsalmon :lightsalmon :orange3 :orange3 :maroon :maroon], xticks=([-π/2, 0, π/2], [L"-\pi/2", L"0", L"\pi/2"]))

# plot(φ_values, eigen_p, xlabel="φ", ylabel=L"$\Re(\lambda)$", legend=false, linewidth=3, ylims=(-10, 10), color=:sienna , xticks=([π/4, π/2, 3π/4, π], ["π/4", "π/2", "3π/4", "π"]))
# plot!(φ_values, eigen_m, xlabel="φ", ylabel=L"$\Re(\lambda)$", legend=false, linewidth=3, ylims=(-10, 10), color=:maroon, xticks=([π/4, π/2, 3π/4, π], ["π/4", "π/2", "3π/4", "π"]))

# savefig("stab_an_spectrum.pdf")


# plot(φ_values/π, im_parts', xlabel=L"\varphi/\pi", ylabel=L"\lambda", linewidth=3, ylims=(-100, 100), xticks=(:auto), label=[L"\Im(\lambda_1^+)" L"\Im(\lambda_1^-)" L"\Im(\lambda_2^+)" L"\Im(\lambda_2^-)" L"\Im(\lambda_3^+)" L"\Im(\lambda_3^-)"], color=[:lightsalmon :lightsalmon :orange3 :orange3 :maroon :maroon], linestyle=:dot, xlims=(0,0.7))
# plot!(φ_values/π, real_parts', xlabel=L"\varphi/\pi", ylabel=L"\lambda", linewidth=3, ylims=(-100, 100), xticks=(:auto), legend=:topright, label=[L"\Re(\lambda_1^+)" L"\Re(\lambda_1^-)" L"\Re(\lambda_2^+)" L"\Re(\lambda_2^-)" L"\Re(\lambda_3^+)" L"\Re(\lambda_3^-)"], color=[:lightsalmon :lightsalmon :orange3 :orange3 :maroon :maroon], linestyle=:solid, alpha=0.7, xlims=(0,0.7))




# # Define the parameters
# Γ_T = 0.5
# κ = 1.0
# ω = 1.0
# ω_0 = 1.0
# g = 0.9

# #Define matrix


# function compute_M(σz_ns) 
#     M = 
#     [
#         -im*ω-κ/2 -im*g 
#         im*g*σz_ns -im*ω_0-Γ_T/2
#     ]
#     return M
# end 

# # Range of σ^z_{ns} values
# Γ_uprange = range(0, stop=1, length=100)
# σz_ns_values = (Γ_uprange .- (1 .- Γ_uprange)) ./ (Γ_uprange .+ (1 .- Γ_uprange)) 


# eigenvalues = [eigen(compute_M(σz)).values for σz in σz_ns_values]

# # Extract real parts of eigenvalues
# real_parts = hcat([real.(eig) for eig in eigenvalues]...)
# im_parts = hcat([imag.(eig) for eig in eigenvalues]...)

# # plot(Γ_uprange, real_parts', xlabel=L"\Gamma_\uparrow/\Gamma_T", ylabel=L"\lambda", linewidth=3, xticks=(:auto), label=[L"\Re(\lambda_+)" L"\Re(\lambda_-)"], color=[:navy :orange3], linestyle=:dot)


# # Plot the real parts of the eigenvalues
# plot(Γ_uprange, real_parts', xlabel=L"\Gamma_\uparrow/\Gamma_T", ylabel=L"\lambda", linewidth=3.5, label=[L"\Re(\lambda_+)" L"\Re(\lambda_-)"], color=[:navy :orange3], linestyle=:dot, ylims=(-2,0.5), xticks=([0, 0.5, 1], ["0", "0.5", "1"]))
# plot!(Γ_uprange, range(0,0,100),linewidth=3, color=:black, linestyle=:solid, label="", alpha=0.7)
# # Add annotation in the top right corner
# annotate!(0.55, 0.40, text(L"\kappa = 1.0 \omega_0", :right, 16, :black))

# plot!(Γ_uprange, im_parts', xlabel=L"\Gamma_\uparrow/\Gamma_T", ylabel=L"\lambda", linewidth=3, label=[L"\Im(\lambda_+)" L"\Im(\lambda_-)"], color=[:navy :orange3], linestyle=:solid, alpha=0.7)


# savefig("stab_spectrum_TC_k1.pdf")



""" 
Below we study the full system with two distinct spin species and species dependent coupling. 

"""

# Define the parameters
δ = 0.0
Δ = 1.0
Δ0 = 1.0 + δ
Δ1 = 1.0 - δ
g = 3.0
κ = 10.0
Γ_T = 1.0

function compute_M(σz, ϕ, g) 
    χR = 1/(im*Δ-κ/2)
    θ = atan(2Δ, κ)
    χ = g^2*abs(χR)*exp(im*θ)
    ξ₋ = g^2*abs(χR)*exp(im*θ)*exp(-im*ϕ)
    ξ₊ = g^2*abs(χR)*exp(im*θ)*exp(im*ϕ)

    M = 
    [
       -im*Δ0 - Γ_T/2 + χ*σz  ξ₊*σz
        ξ₋*σz -im*Δ1 - Γ_T/2 + χ*σz
    ]
    return M
end 

function count_eigenreals(eigenlist)
    no_eigenreals = []
    for entry in 1:size(eigenlist, 2)
        count = sum(eigenlist[:, entry] .> 1e-8)
        push!(no_eigenreals, count)
    end
    return no_eigenreals
end

# Range of σ^z_{ns} values
Γ_uprange = range(0, stop=1, length=100)
σz_ns_values = (Γ_uprange .- (1 .- Γ_uprange)) ./ (Γ_uprange .+ (1 .- Γ_uprange))
g_range = range(0, κ, length=10)

# Range of φ values
φ_values = range(-π/2, stop=π/2, length=10)

eigenvalues = [eigen(compute_M(σz, φ, g)).values for σz in σz_ns_values, φ in φ_values]

# Extract real parts of eigenvalues
real_parts = hcat([real.(eig) for eig in eigenvalues]...)

# Plot the number of eigenvalues with positive real part as a function of φ and Γ_up
positive_real_parts = [sum(real_parts .> 0, dims=1)']


plot(φ_values, real_parts', xlabel=L"\varphi", ylabel=L"\Re(\lambda)", legend=false, linewidth=3 , xticks=([-π/2, 0, π/2], [L"-\pi/2", L"0", L"\pi/2"]))
plot!(φ_values, im_parts', linewidth=3)



z = reshape(count_eigenreals(real_parts)', length(Γ_uprange),:)
heatmap(φ_values./π, Γ_uprange, z, xlabel=L"\varphi/\pi", ylabel=L"\Gamma_\uparrow/\Gamma_T", colormap = :blues)
annotate!(-0.5, 0.75, text(L"g = \Delta", :right, 16, :black))

heatmap(g_range, Γ_uprange, reshape(count_eigenreals(real_parts)', length(Γ_uprange), :), xlabel=L"g", ylabel=L"\Gamma_\uparrow/\Gamma_T", colormap = :blues)


""" 

Below we study the same system but moving to the Pauli basis instead, the stability matrix now is 4×4

"""


function compute_M_xy(σz, ϕ, g)
    χR = 1/(im*Δ-κ/2)
    θ = atan(2Δ, κ)
    χ = 2g^2*abs(χR)
    Γ₊ = -Γ_T/2 + χ*σz*cos(θ)
    Γ₋ = -Γ_T/2 - χ*σz*cos(θ)
    M = [
        Γ₋  -χ*σz*cos(θ+ϕ)  -Δ0-χ*σz*sin(θ)  -χ*σz*sin(θ+ϕ)
        -χ*σz*cos(θ-ϕ)  Γ₋  -χ*σz*sin(θ-ϕ)  -Δ1-χ*σz*sin(θ)
        Δ0-χ*σz*sin(θ)  χ*σz*sin(θ+ϕ)  Γ₊  χ*σz*cos(θ+ϕ)
        χ*σz*sin(θ-ϕ)  Δ1-χ*σz*sin(θ)  χ*σz*cos(θ-ϕ)  Γ₊
    ]
return M
end

"""
 To compare against the "bare" single species model 

"""
function compute_M_xy_singlespin(σz, g)
    χR = 1/(im*Δ-κ/2)
    θ = atan(2Δ, κ)
    χ = 2g^2*abs(χR)
    Γ₊ = -Γ_T/2 + χ*σz*cos(θ)
    Γ₋ = -Γ_T/2 - χ*σz*cos(θ)
    M = [
        Γ₋ 0  -Δ0-χ*σz*sin(θ)  0    
        0  0  0  0
        Δ0-χ*σz*sin(θ) 0  Γ₊ 0
        0  0  0  0
    ]
    return M 
end

ϕ = π/2
eigenvaluesxy = [eigen(compute_M_xy(σz, ϕ, g)).values for σz in σz_ns_values, g in g_range]
#eigenvaluesxy = [eigen(compute_M_xy(σz, ϕ, g)).values for σz in σz_ns_values, ϕ in φ_values]
real_partsxy = hcat([real.(eig) for eig in eigenvaluesxy]...)



z = reshape(count_eigenreals(real_partsxy)', length(Γ_uprange), length(φ_values))
z = reshape(count_eigenreals(real_partsxy), length(Γ_uprange), length(g_range))

heatmap(φ_values ./ π, Γ_uprange, z, xlabel=L"\varphi/\pi", ylabel=L"\Gamma_\uparrow/\Gamma_T", colormap = :blues)
heatmap(g_range, Γ_uprange, z, xlabel=L"g", ylabel=L"\Gamma_\uparrow/\Gamma_T", colormap = :blues)
annotate!(3, 0.75, text(L"\varphi = 0", :right, 16, :black))


eigenvalues_single = [eigen(compute_M_xy_singlespin(σz, g)).values for σz in σz_ns_values, g in g_range]
real_parts_single = hcat([real.(eig) for eig in eigenvalues_single]...)

zsingle = reshape(count_eigenreals(real_parts_single), length(Γ_uprange), length(g_range))
heatmap(g_range, Γ_uprange, zsingle, xlabel=L"g", ylabel=L"\Gamma_\uparrow/\Gamma_T", colormap = :blues)

"""
For a fixed g value:
"""
eigenvalues_g = [eigen(compute_M_xy(σz, 0, 5)).values for σz in σz_ns_values]
eigenvalues_gs = [eigen(compute_M_xy_singlespin(σz, 5)).values for σz in σz_ns_values]
# Extract real parts of eigenvalues
real_parts_g = hcat([real.(eig) for eig in eigenvalues_g]...)
real_parts_gs = hcat([real.(eig) for eig in eigenvalues_gs]...)

im_parts_g = hcat([imag.(eig) for eig in eigenvalues_g]...)
im_parts_gs = hcat([imag.(eig) for eig in eigenvalues_gs]...)
# Plot the real parts of the eigenvalues for g = 3
plot(Γ_uprange, real_parts_g', xlabel=L"\Gamma_\uparrow/\Gamma_T", ylabel=L"\Re(\lambda)", legend=false, linewidth=4, grid=false, color=:rebeccapurple, ylims=(-5,5))
plot!(Γ_uprange, im_parts_g', linewidth = 3.5, linestyle = :dashdot, color=:rebeccapurple)
plot!(Γ_uprange, real_parts_gs', linewidth=4, color=:seagreen)
plot!(Γ_uprange, im_parts_gs', linewidth=3.5, linestyle=:dashdot, color=:seagreen)



eigenvecs = [eigen(compute_M(σz,0.1, 1)).vectors for ϕ in σz_ns_values]


# Extract components of the eigenvectors
eigenvecs_real = hcat([real.(vec) for vec in eigenvecs]...)
eigenvecs_imag = hcat([imag.(vec) for vec in eigenvecs]...)

# Ensure the dimensions match
eigenvecs_real = eigenvecs_real[:, 1:length(φ_values)]
eigenvecs_imag = eigenvecs_imag[:, 1:length(φ_values)]

eigenvecs_real = eigenvecs_real[:, 1:length(σz_ns_values)]
eigenvecs_imag = eigenvecs_imag[:, 1:length(σz_ns_values)]


# Plot the real parts of the eigenvector components as a function of φ
plot(φ_values, eigenvecs_real', xlabel=L"\varphi", ylabel=L"\Re", legend=false, linewidth=4)
plot!(φ_values, eigenvecs_imag', xlabel=L"\varphi", ylabel=L"\Im)", legend=false, linewidth=4, linestyle=:dash)
plot(σz_ns_values, eigenvecs_real', xlabel=L"\langle\sigma^z\rangle", ylabel=L"\Re", legend=false, linewidth=4)


