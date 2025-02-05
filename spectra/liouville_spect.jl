
using QuantumOptics
using LinearAlgebra
using QuantumOptics
using Arpack 
using LaTeXStrings
using Plots

ω0 = 1.0
δ = 0.05*ω0;


κ= 20.5*ω0
g =1*ω0
Δ = 2*ω0
Δ1 = ω0 - δ
Δ2 = ω0 + δ
η = 1*ω0
Γ = 0 #2.5*ω0/100
ϕ = π/4

n_max = 5     # Maximum photon number in cavity
N=2
# Hilbert spaces
cavity_space = FockBasis(n_max)  # Cavity space
spin_space = [SpinBasis(1//2) for _ in 1:N]  # Each spin is a spin-1/2 system
full_space = tensor(cavity_space, spin_space...)  # Total Hilbert space

# Cavity operators
a = destroy(cavity_space)
adag = dagger(a)
a_op = embed(full_space, 1, a)  # Embed the cavity operator in the total space
adag_op = dagger(a_op)

# Spin operators
σz1 = embed(full_space, 2, sigmaz(SpinBasis(1//2)))  # σ(2,2,1) for spin species 1
σz2 = embed(full_space, 3, sigmaz(SpinBasis(1//2)))  # σ(2,2,2) for spin species 2

# Initialize the full Hamiltonian
H = -Δ * adag_op * a_op - Δ1 * σz1 - Δ2 * σz2

# Add coupling terms for all spins
for i in 1:N
    σ12 = embed(full_space, i + 1, sigmam(SpinBasis(1//2)))  # σ(1,2,i)
    σ21 = dagger(σ12)  # σ(2,1,i)
    phase_factor = exp(-im * ϕ * (i - 1))
    H += g * (phase_factor * adag_op * σ12 + conj(phase_factor) * a_op * σ21)
end

# Add the driving term
H += im * η * (adag_op - a_op)

# Define the jump operators (dissipation)
J = [sqrt(κ) * a_op]

# Build the Liouvillian
L = liouvillian(H, J)

# Compute eigenvalues of the Liouvillian
n_eigenvals = 1294
eigvals_L = eigs(L.data, nev=n_eigenvals, which=:LM)[1]

# Extract real and imaginary parts
real_parts = real(eigvals_L)
imag_parts = imag(eigvals_L)

# Plot Im(λ) vs Re(λ)
plot(real_parts, imag_parts, seriestype=:scatter, legend=false, xlabel=L"Re($\lambda$)", ylabel=L"Im($\lambda$)", title="Liouville spectrum", color=:black, markersize=3)

# Find smallest nonzero Re(λ)
sorted_real_parts = sort(real_parts)[end-20:end] # Sort by real part

sorted_imag_parts = imag_parts[sortperm(real_parts)][end-20:end] # Sort imaginary parts accordingly

plot(sorted_real_parts, sorted_imag_parts, seriestype=:scatter, legend=false, xlabel=L"Re($\lambda$)", ylabel=L"Im($\lambda$)", title="Liouville spectrum", color=:black, markersize=3)