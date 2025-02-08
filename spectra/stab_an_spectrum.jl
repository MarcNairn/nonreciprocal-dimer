using LinearAlgebra
using Plots
using LaTeXStrings

default(
    fontfamily = "Computer Modern", # TeX's default serif font
    tickfont = font(12, "Computer Modern"),  # Tick labels font
    guidefont = font(14, "Computer Modern"), # Axis labels font
    legendfont = font(12, "Computer Modern"), # Legend font
    titlefont = font(16, "Computer Modern")  # Title font
)

# Define the parameters
ω0 = 1.0
κ= 10*ω0 
Δ = 20*ω0
g = 2*κ
χR = 1/sqrt(Δ^2 + κ^2/4)

θl = atan(2Δ,κ)
Δ1 = ω0
Δ2 = ω0
η = 5*κ


# Function to compute the matrix M for a given φ
function compute_M(φ)
    θl = atan(2Δ, κ)
    ξp = 2 * g^2 * abs(χR) * cos(θl)
    ξm = 2 * g^2 * abs(χR) * sin(θl)
    
    χp(m) = 2 * g^2 * abs(χR) * cos(θl + (-1)^m * φ)
    χm(m) = 2 * g^2 * abs(χR) * sin(θl + (-1)^m * φ)
    
    χ1p = χp(1)
    χ1m = χm(1)
    χ2p = χp(2)
    χ2m = χm(2)
    
    M = [
        ξp   χ1p   -Δ1 + ξm   χ1m   0       0
        χ2p  ξp    χ2m       -Δ2 + ξp   0       0
        ξm + Δ1  χ1m   -ξp      -χ1p   0       0
        χ2m  ξm + Δ2  -χ2p     -ξp    0       0
        η * χ1m / g  0   -η * χ1m / g  0   0   0
        0   η * χ2m / g  0   -η * χ2m / g  0   0
    ]
    return M
end

# Range of φ values
φ_values = range(00, stop=π/2+0.1, length=10000)

# Compute eigenvalues for each φ
eigenvalues = [eigen(compute_M(φ)).values for φ in φ_values]

# Extract real parts of eigenvalues
real_parts = hcat([real.(eig) for eig in eigenvalues]...)
im_parts = hcat([imag.(eig) for eig in eigenvalues]...)

eigenvec_part = hcat([eigen(compute_M(φ)).vectors for φ in φ_values]...)


eigen_p = [real(eig[3]) for eig in eigenvalues]
eigen_m = [real(eig[4]) for eig in eigenvalues]

# Plot the real parts of the eigenvalues
plot(φ_values, real_parts', xlabel=L"\varphi", ylabel=L"\Re(\lambda)", legend=false, linewidth=3, ylims=(-1,1), color=[:lightsalmon :lightsalmon :orange3 :orange3 :maroon :maroon], xticks=([-π/2, 0, π/2], [L"-\pi/2", L"0", L"\pi/2"]))

plot(φ_values, eigen_p, xlabel="φ", ylabel=L"$\Re(\lambda)$", legend=false, linewidth=3, ylims=(-10, 10), color=:sienna , xticks=([π/4, π/2, 3π/4, π], ["π/4", "π/2", "3π/4", "π"]))
plot!(φ_values, eigen_m, xlabel="φ", ylabel=L"$\Re(\lambda)$", legend=false, linewidth=3, ylims=(-10, 10), color=:maroon, xticks=([π/4, π/2, 3π/4, π], ["π/4", "π/2", "3π/4", "π"]))

savefig("stab_an_spectrum.pdf")


plot(φ_values/π, im_parts', xlabel=L"\varphi/\pi", ylabel=L"\lambda", linewidth=3, ylims=(-100, 100), xticks=(:auto), label=[L"\Im(\lambda_1^+)" L"\Im(\lambda_1^-)" L"\Im(\lambda_2^+)" L"\Im(\lambda_2^-)" L"\Im(\lambda_3^+)" L"\Im(\lambda_3^-)"], color=[:lightsalmon :lightsalmon :orange3 :orange3 :maroon :maroon], linestyle=:dot, xlims=(0,0.7))
plot!(φ_values/π, real_parts', xlabel=L"\varphi/\pi", ylabel=L"\lambda", linewidth=3, ylims=(-100, 100), xticks=(:auto), legend=:topright, label=[L"\Re(\lambda_1^+)" L"\Re(\lambda_1^-)" L"\Re(\lambda_2^+)" L"\Re(\lambda_2^-)" L"\Re(\lambda_3^+)" L"\Re(\lambda_3^-)"], color=[:lightsalmon :lightsalmon :orange3 :orange3 :maroon :maroon], linestyle=:solid, alpha=0.7, xlims=(0,0.7))




# Define the parameters
Γ_T = 0.5
κ = 1.0
ω = 1.0
ω_0 = 1.0
g = 0.9

#Define matrix


function compute_M(σz_ns) 
    M = 
    [
        -im*ω-κ/2 -im*g 
        im*g*σz_ns -im*ω_0-Γ_T/2
    ]
    return M
end 

# Range of σ^z_{ns} values
Γ_uprange = range(0, stop=1, length=100)
σz_ns_values = (Γ_uprange .- (1 .- Γ_uprange)) ./ (Γ_uprange .+ (1 .- Γ_uprange)) 


eigenvalues = [eigen(compute_M(σz)).values for σz in σz_ns_values]

# Extract real parts of eigenvalues
real_parts = hcat([real.(eig) for eig in eigenvalues]...)
im_parts = hcat([imag.(eig) for eig in eigenvalues]...)

# plot(Γ_uprange, real_parts', xlabel=L"\Gamma_\uparrow/\Gamma_T", ylabel=L"\lambda", linewidth=3, xticks=(:auto), label=[L"\Re(\lambda_+)" L"\Re(\lambda_-)"], color=[:navy :orange3], linestyle=:dot)


# Plot the real parts of the eigenvalues
plot(Γ_uprange, real_parts', xlabel=L"\Gamma_\uparrow/\Gamma_T", ylabel=L"\lambda", linewidth=3.5, label=[L"\Re(\lambda_+)" L"\Re(\lambda_-)"], color=[:navy :orange3], linestyle=:dot, ylims=(-2,0.5), xticks=([0, 0.5, 1], ["0", "0.5", "1"]))
plot!(Γ_uprange, range(0,0,100),linewidth=3, color=:black, linestyle=:solid, label="", alpha=0.7)
# Add annotation in the top right corner
annotate!(0.55, 0.40, text(L"\kappa = 1.0 \omega_0", :right, 16, :black))

plot!(Γ_uprange, im_parts', xlabel=L"\Gamma_\uparrow/\Gamma_T", ylabel=L"\lambda", linewidth=3, label=[L"\Im(\lambda_+)" L"\Im(\lambda_-)"], color=[:navy :orange3], linestyle=:solid, alpha=0.7)


savefig("stab_spectrum_TC_k1.pdf")



""" 
Below we study the full system with two distinct spin species and species dependent coupling. 

"""

# Define the parameters
δ = 0.0
Δ = 1.0
Δ0 = 1.0 + δ
Δ1 = 1.0 - δ
g = 2.0
κ = 10.0
Γ_T = 1.0

function compute_M(σz, ϕ) 
    χR = 1/(im*Δ-κ/2)
    θ = atan(2Δ, κ)
    χ = g^2*abs(χR)*exp(im*θ)
    ξ₋ = g^2*abs(χR)*exp(im*θ)*exp(-im*ϕ)
    ξ₊ = g^2*abs(χR)*exp(im*θ)*exp(im*ϕ)

    M = 
    [
       -im*Δ - Γ_T/2 + χ*σz  ξ₊*σz
        ξ₋*σz -im*Δ - Γ_T/2 + χ*σz
    ]
    return M
end 

# Range of σ^z_{ns} values
Γ_uprange = range(0, stop=1, length=1000)
σz_ns_values = (Γ_uprange .- (1 .- Γ_uprange)) ./ (Γ_uprange .+ (1 .- Γ_uprange))

# Range of φ values
φ_values = range(-π, stop=π, length=1000)

eigenvalues = [eigen(compute_M(σz, φ)).values for σz in σz_ns_values, φ in φ_values]

# Extract real parts of eigenvalues
real_parts = hcat([real.(eig) for eig in eigenvalues]...)

# Plot the number of eigenvalues with positive real part as a function of φ and Γ_up
positive_real_parts = [sum(real_parts .> 0, dims=1)']



eigenvalues_phi = [eigen(compute_M(0.15, φ)).values for φ in φ_values]


real_parts = hcat([real.(eig) for eig in eigenvalues_phi]...)
im_parts = hcat([imag.(eig) for eig in eigenvalues_phi]...)

plot(φ_values, real_parts', xlabel=L"\varphi", ylabel=L"\Re(\lambda)", legend=false, linewidth=3 , xticks=([-π/2, 0, π/2], [L"-\pi/2", L"0", L"\pi/2"]))
plot!(φ_values, im_parts', linewidth=3)

function count_eigenreals(eigenlist)
    no_eigenreals1 = []
    no_eigenreals2 = []
    for entry in 1:Int(length(eigenlist)/2)
        if eigenlist[1,entry] > 1e-8 
            push!(no_eigenreals1, 1)
        else
            push!(no_eigenreals1, 0)
        end
        if eigenlist[2,entry] > 1e-8
            push!(no_eigenreals2, 1)
        else   
            push!(no_eigenreals2, 0)
        end
    end
    no_eigenreals = no_eigenreals1 .+ no_eigenreals2
    return no_eigenreals
    end



heatmap(φ_values./π, Γ_uprange, reshape(count_eigenreals(real_parts)', length(Γ_uprange), :), xlabel=L"\varphi/\pi", ylabel=L"\Gamma_\uparrow/\Gamma_T", colormap = :blues)