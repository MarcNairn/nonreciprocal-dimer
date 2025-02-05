import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# Define the symbols
omega, omega0, kappa, gamma, g, GammaT = sp.symbols('omega omega0 kappa gamma g GammaT', real=True)

# Define the matrix M_inc
M_inc = sp.Matrix([
    [-sp.I*omega - kappa/2, -sp.I*g], 
    [sp.I*g*gamma, -sp.I*omega0 - GammaT/2],
])

# Define numerical values for parameters
omega0_val = 1
kappa_val = 1
omega_val = 1
g_val = 0.9
Gamma_T_val = 1

# Define a range of gamma values
gammarange = np.linspace(-1, 1, 101)

# Initialize lists to store the real parts of the eigenvalues
real_eigenvalue1 = []
real_eigenvalue2 = []

# Loop over gamma values
for gamma_val in gammarange:
    # Substitute numerical values into the matrix
    M_num = M_inc.subs({omega0: omega0_val, kappa: kappa_val, omega: omega_val, g: g_val, GammaT: Gamma_T_val, gamma: gamma_val})
    
    # Compute the eigenvalues
    eigenvalues = M_num.eigenvals()
    
    # Extract the real parts of the eigenvalues and evaluate them numerically
    real_parts = [sp.re(eigen).evalf() for eigen in eigenvalues]
    
    # Append the real parts to the lists
    real_eigenvalue1.append(float(real_parts[0]))
    real_eigenvalue2.append(float(real_parts[1]))

# Convert the lists to numpy arrays for plotting
real_eigenvalue1 = np.array(real_eigenvalue1)
real_eigenvalue2 = np.array(real_eigenvalue2)

# Plot the real parts of the eigenvalues as a function of gammarange
plt.plot(gammarange, real_eigenvalue1, label=r'Re($\lambda_1$)', color='blue')
plt.plot(gammarange, real_eigenvalue2, label=r'Re($\lambda_2$)', color='red')
plt.xlabel(r'$\gamma$')
plt.ylabel(r'Re($\lambda$)')
plt.title(r'Real Part of Eigenvalues vs $\gamma$')
plt.legend()
plt.grid(True)
plt.show()