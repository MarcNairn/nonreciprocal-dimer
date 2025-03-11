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

# Define parameters
g = 0.75
η = 0.5
κ = 1.0
Δ = 1.0  
χR = (im*Δ-κ/2)^(-1)
θ = atan(2Δ/κ)
χ = g^2*abs(χR)*exp(im*θ) + η    


ϕ_range = range(-π, π, length=1000)  # Values for ϕ
g_range = range(0, 2κ, length = 1000)
η_range = range(0, 2κ, length = 1000)
# Compute λ_±
λ_plus = [-im*Δ - χ + sqrt(abs2(χ) + η^2 + 2η*χ*cos(ϕ)) for ϕ in ϕ_range]
λ_minus = [-im*Δ - χ - sqrt(abs2(χ) + η^2 + 2η*χ*cos(ϕ)) for ϕ in ϕ_range]

# Extract real and imaginary parts
real_plus = real.(λ_plus)
imag_plus = imag.(λ_plus)
real_minus = real.(λ_minus)
imag_minus = imag.(λ_minus)

# Plot
plot(
    ϕ_range./π, real_plus, label=L"\Re(\lambda_+)", linewidth=3, color=:blue, xlabel = L"\varphi/\pi", ylabel = L"\lambda"
)
plot!(
    ϕ_range./π, imag_plus, label=L"\Im(\lambda_+)", linewidth=3, linestyle=:dash, color=:blue
)
plot!(
    ϕ_range./π, real_minus, label=L"\Re(\lambda_-)", linewidth=3, color=:red
)
plot!(
    ϕ_range./π, imag_minus, label=L"\Im(\lambda_-)", linewidth=3, linestyle=:dash, color=:red
)

""" 3D plotting below """
# Prepare matrices for 3D plot
X = repeat(ϕ_range, 1, length(g_range))  # ϕ as x-axis
Z = repeat(g_range', length(ϕ_range), 1)  # g as z-axis
real_plus_vals = zeros(length(ϕ_range), length(η_range))
imag_plus_vals = zeros(length(ϕ_range), length(η_range))
real_minus_vals = zeros(length(ϕ_range), length(η_range))
imag_minus_vals = zeros(length(ϕ_range), length(η_range))

# Compute λ_± for different η values
for (j, g) in enumerate(η_range)
    χR = (im*Δ - κ/2)^(-1)
    θ = atan(2Δ/κ)
    χ = g^2 * abs(χR) * exp(im*θ) + η

    for (i, ϕ) in enumerate(ϕ_range)
        λ_plus = -im*Δ - χ + sqrt(abs2(χ) + η^2 + 2η*χ*cos(ϕ))
        λ_minus = -im*Δ - χ - sqrt(abs2(χ) + η^2 + 2η*χ*cos(ϕ))

        real_plus_vals[i, j] = real(λ_plus)
        imag_plus_vals[i, j] = imag(λ_plus)
        real_minus_vals[i, j] = real(λ_minus)
        imag_minus_vals[i, j] = imag(λ_minus)
    end
end

# Plot surfaces
p1 = surface(X ./ π, Z, real_plus_vals, color=:magma, xlabel=L"\varphi/\pi", ylabel=L"\eta", zlabel=L"\Re(\lambda)", title=L"Real parts $\lambda_\pm$", colorbar=false)
p2 = surface(X ./ π, Z, imag_plus_vals, color=:magma, xlabel=L"\varphi/\pi", ylabel=L"\eta", zlabel=L"\Im(\lambda)", title=L"Imaginary part of $\lambda_+$")
p3 = surface(X ./ π, Z, real_minus_vals, color=:magma, xlabel=L"\varphi/\pi", ylabel=L"\eta", zlabel=L"\Re(\lambda)")
p4 = surface(X ./ π, Z, imag_minus_vals, color=:magma, xlabel=L"\varphi/\pi", ylabel=L"\eta", zlabel=L"\Im(\lambda)", title=L"Imaginary part of $\lambda_-$")

savefig("surfaces_eigenreals_sdrive.pdf")
plot(p1, p2, p3, p4, layout=(2, 2))

"""Eigenreal phase diagram"""

# Function to count the number of eigenvalues with positive real parts
function count_positive_real(λs)
    return sum(real(λ) > 1e-8 for λ in λs)  # Count eigenvalues with Re(λ) > 1e-8
end
g_range = range(0, 2κ, length=1000)  
ϕ_range = range(-π, π, length=1000)  
η = 0.5*κ

# Store number of eigenvalues with positive real parts
num_eigenvalues = zeros(length(ϕ_range), length(g_range))

# Compute eigenvalues and count positive real parts
for (k, g) in enumerate(g_range)
    χR = (im*Δ - κ/2)^(-1)
    θ = atan(2Δ/κ)
    χ = g^2 * abs(χR) * exp(im*θ) + η

    for (i, ϕ) in enumerate(ϕ_range)
        λ_plus = -im*Δ - χ + sqrt(abs2(χ) + η^2 + 2η*χ*cos(ϕ))
        λ_minus = -im*Δ - χ - sqrt(abs2(χ) + η^2 + 2η*χ*cos(ϕ))

        # Store count of positive real eigenvalues
        num_eigenvalues[i, k] = count_positive_real([λ_plus, λ_minus])
    end
end

# Create a heatmap for the phase diagram
heatmap(ϕ_range ./ π, g_range, num_eigenvalues', color=:magma,
        xlabel=L"\varphi/\pi", ylabel=L"g/\kappa",colorbar=:false)

savefig("eigenreals_sdrive.pdf")

