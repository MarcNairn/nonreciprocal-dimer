using QuantumCumulants
using ModelingToolkit, OrdinaryDiffEq
using Plots
using LaTeXStrings
using PGFPlotsX
default(
    fontfamily = "Computer Modern", # TeX's default serif font
    tickfont = font(12, "Computer Modern"),  # Tick labels font
    guidefont = font(14, "Computer Modern"), # Axis labels font
    legendfont = font(12, "Computer Modern"), # Legend font
    titlefont = font(16, "Computer Modern")  # Title font
)

# Define parameters
κ, g, Δ, Δ1, Δ2, ϕ= cnumbers("κ g Δ Δ1 Δ2 ϕ")

N=2

# Define hilbert space
hf = FockSpace(:cavity)
ha = ⊗([NLevelSpace(Symbol(:atom,i),2) for i=1:N]...)
h = hf ⊗ ha

# Define the fundamental operators
a = Destroy(h,:a)
σ(i,j,k) = Transition(h,Symbol("σ_{$k}"),i,j,k+1)

# Define the Hamiltonian
function pickH(name::String)
    if name=="TC"
        H = -Δ * a' * a - Δ1 * σ(2,2,1) - Δ2*σ(2,2,2) + sum(g*(exp(-im*ϕ*(i-1))*a'*σ(1,2,i) + exp(im*ϕ*(i-1))*a*σ(2,1,i)) for i=1:N)
    elseif name=="Dicke"
        H = -Δ * a' * a - Δ1 * σ(2,2,1) - Δ2*σ(2,2,2) + sum(g*(exp(-im*ϕ*(i-1))*a'*σ(1,2,i) + exp(im*ϕ*(i-1))*a*σ(2,1,i) +exp(-im*ϕ*(i-1))*a'*σ(2,1,i) + exp(im*ϕ*(i-1))*a*σ(1,2,i)) for i=1:N)
    end
return H
end


H = pickH("Dicke")
# Collapse operators
J = [a]
rates = [κ]

# Define the operators of interest
ops = [a',a, σ(2,2,1), σ(2,2,2), σ(1,2,1), σ(1,2,2)]


# Derive equation for average photon number
eqs = meanfield(ops,H,J;rates=rates,order=1)
#eqs_expanded = cumulant_expansion(eqs,2)

me_comp = complete(eqs) #automatically complete the system

# Build an ODESystem out of the MeanfieldEquations
@named sys = ODESystem(me_comp)


# initial state
u0 = zeros(ComplexF64, length(eqs))
u0[1] = sqrt(10)
u0[2] = sqrt(10)
u0[3] = 0
u0[4] = 0

δ = 0.0025;

κn= 1.0
gn = 2.5*κn/20
Δn = 1.0*κn
Δ1n = κn/20 - δ
Δ2n = κn/20 + δ
ϕn = π/4

# list of parameters
ps = (Δ, g, κ, ϕ, Δ1, Δ2)
p0 = ps .=> (Δn, gn, κn, ϕn, Δ1n, Δ2n)
tend = 10000/κn

prob = ODEProblem(sys,u0,(0.0,tend),p0)# Solve the initial ODE problem
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# Extract the final state
final_state = sol.u[end]

# Update parameters with new ϕ value
p0_new = ps .=> (Δn, gn, κn, -1*ϕn, Δ1n, Δ2n)

# Solve the new ODE problem with the final state as the initial state
prob = ODEProblem(sys,final_state,(0.0,tend),p0_new)
sol_new = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# Concatenate time and solution arrays
t_combined = [sol.t; sol_new.t .+ tend]
sigmaz1_combined = [real.(2 * sol[σ(2, 2, 1)] .- 1); real.(2 * sol_new[σ(2, 2, 1)] .- 1)]
sigmaz2_combined = [real.(2 * sol[σ(2, 2, 2)] .- 1); real.(2 * sol_new[σ(2, 2, 2)] .- 1)]

# Plot
plot(t_combined, [sigmaz1_combined, sigmaz2_combined], xlabel = L"t\kappa", ylabel = L"\sigma^z", label = ["Atom 1" "Atom 2"])