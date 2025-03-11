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
κ, g, ω, ω0, η, Γu, Γd= cnumbers("κ g ω ω0 η Γ_u Γ_d")

N=1

# Define hilbert space
hf = FockSpace(:cavity)
ha = ⊗([NLevelSpace(Symbol(:atom,i),2) for i=1:N]...)
h = hf ⊗ ha

# Define the fundamental operators
a = Destroy(h,:a)
σ(i,j,k) = Transition(h,Symbol("σ_{$k}"),i,j,k+1)


H = ω * a' * a + ω0 * σ(2,2,1)  + g*a'*σ(1,2,1) + a*σ(2,1,1)  

J = [a, σ(2,1,1), σ(1,2,1)]
rates = [κ, Γu, Γd]


ops = [a, a', σ(1,2,1),σ(2,2,1)]

# Derive equation for average photon number
eqs = meanfield(ops,H,J;rates=rates,order=1)

me_comp = complete(eqs) 

@named sys = ODESystem(me_comp)

u0 = zeros(ComplexF64, length(me_comp))
u0[1] = sqrt(1)
u0[2] = sqrt(1)
u0[3] = sqrt(0)
u0[4] = (Γun-Γdn)/(Γun+Γdn)

ω0n = 1.0

κn= 1*ω0n
gn = 0.9*ω0n
ωn = 1.0*ω0n
Γun = 0.3*ω0n
Γdn = 1-Γun


# list of parameters
ps = (ω, ω0, κ, g, Γu, Γd)
p0 = ps .=> (ωn, ω0n, κn, gn, Γun, Γdn)
tend = 300/κn

prob = ODEProblem(sys, u0, (0.0, tend), p0)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

a_dag_a = [real.(sol[a'] .* sol[a])]
sz = [real.(sol[σ(2,2,1)])] 

plot(sol.t, a_dag_a, xlabel=L"Time ($\kappa$)", ylabel=L"\langle a^\dagger a \rangle", legend=false, linewidth=3)
plot(sol.t, sz, xlabel=L"Time ($\kappa$)", ylabel=L"\langle \sigma_z \rangle", legend=false, linewidth=3)