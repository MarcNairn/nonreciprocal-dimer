import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment

# Set up plot style
plt.rcParams.update({
    'text.usetex': True,             # Enable LaTeX for text rendering
    'font.family': 'serif',          # Use a serif font family
    'font.serif': ['Computer Modern Roman'],  # Set the serif font to Computer Modern
    'axes.titlesize': 18,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14
})

κ = 1
K= κ/2
η = 0.3
g = 0.3
J = g**2/K**2

def compute_eigenvalues(phi):
    γ = g/2
    # Compute spin variables
    if η < γ:
        sx = η*np.cos(phi)/γ
        sy = -η*np.sin(phi)/γ
        sz = -np.sqrt(1 - η**2/γ**2)
    else:
        sx = γ*np.cos(phi)/η - np.sin(phi)*np.sqrt(1 - γ**2/η**2)
        sy = -γ*np.sin(phi)/η - np.cos(phi)*np.sqrt(1 - γ**2/η**2)
        sz = 0

    # Matrix elements
    v1x = (1 + np.cos(phi))*sx + np.sin(phi)*sy - (2*η)/g
    w1x = np.cos(phi)*sx + np.sin(phi)*sy - (2*η)/g
    v1y = np.sin(phi)*sx + (1 + np.cos(phi))*sy
    w1y = -np.sin(phi)*sx + np.cos(phi)*sy

    v2x = (1 + np.cos(phi))*sx + np.sin(phi)*sy - (2*η*np.cos(phi))/g
    w2x = np.cos(phi)*sx + np.sin(phi)*sy - (2*η*np.cos(phi))/g
    v2y = -np.sin(phi)*sx + (1 + np.cos(phi))*sy + (2*η*np.sin(phi))/g
    w2y = -np.sin(phi)*sx + np.cos(phi)*sy

    ux = np.cos(phi)*sx - np.sin(phi)*sy
    uy = np.sin(phi)*sx + np.cos(phi)*sy

    # Build matrix
    M = J*K*np.array([
        [sz, 0, v1x, sz*np.cos(phi), sz*np.sin(phi), 0],
        [0, sz, v1y, -sz*np.sin(phi), sz*np.cos(phi), 0],
        [-w1x, -w1y, 0, -ux, -uy, 0],
        [sz*np.cos(phi), sz*np.sin(phi), 0, sz, 0, v2x],
        [-sz*np.sin(phi), sz*np.cos(phi), 0, 0, sz, v2y],
        [-ux, -uy, 0, -w2x, -w2y, 0]
    ])
    
    return np.linalg.eigvals(M)

# Generate phi values
phi_vals = np.linspace(-np.pi, np.pi, 300)

# Compute eigenvalues with tracking
def track_eigenvalues(eigvals_list):
    tracked = [eigvals_list[0]]
    
    for i in range(1, len(eigvals_list)):
        # Create cost matrix using absolute differences
        cost = np.abs(np.subtract.outer(tracked[i-1], eigvals_list[i]))
        
        # Use Hungarian algorithm for optimal assignment
        row_ind, col_ind = linear_sum_assignment(cost)
        
        # Reorder eigenvalues
        tracked.append(eigvals_list[i][col_ind])
        
    return np.array(tracked)

# Compute and track eigenvalues
raw_eigs = np.array([compute_eigenvalues(phi) for phi in phi_vals])
tracked_eigs = track_eigenvalues(raw_eigs)

# Plotting
color_palette = ['#E69F00', '#56B4E9', '#009E73', '#0072B2', '#D55E00', '#CC79A7']

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

for i in range(6):
    ax1.plot(phi_vals/np.pi, tracked_eigs[:, i].real, 
             color=color_palette[i], linewidth=3)
    ax2.plot(phi_vals/np.pi, tracked_eigs[:, i].imag,
             color=color_palette[i], linewidth=3)

ax1.set_ylabel(r'$\Re(\lambda)$')
ax2.set_ylabel(r'$\Im(\lambda)$')
ax2.set_xlabel(r'$\varphi/\pi$')

plt.tight_layout()
plt.show()

""" PHASE DIAGRAMS """

# # Parameter ranges
# eta_vals = np.linspace(0.01, 1.0, 100)
# g_vals = np.linspace(0.01, 1.0, 100)
# phi_vals = np.linspace(-np.pi/2, np.pi/2, 100)

# # Function to count unstable eigenvalues
# # Function to count unstable eigenvalues
# def count_unstable_eigenvalues(phi, eta, g):

#     gamma = g / 2
#     if eta < gamma:
#         sx = eta * np.cos(phi) / gamma
#         sy = -eta * np.sin(phi) / gamma
#         sz = -np.sqrt(1 - eta**2 / gamma**2)
#     else:
#         sx = gamma * np.cos(phi) / eta - np.sin(phi) * np.sqrt(1 - gamma**2 / eta**2)
#         sy = -gamma * np.sin(phi) / eta - np.cos(phi) * np.sqrt(1 - gamma**2 / eta**2)
#         sz = 0

#     v1x = (1 + np.cos(phi)) * sx + np.sin(phi) * sy - (2 * eta) / g
#     w1x = np.cos(phi) * sx + np.sin(phi) * sy - (2 * eta) / g
#     v1y = np.sin(phi) * sx + (1 + np.cos(phi)) * sy
#     w1y = -np.sin(phi) * sx + np.cos(phi) * sy

#     v2x = (1 + np.cos(phi)) * sx + np.sin(phi) * sy - (2 * eta * np.cos(phi)) / g
#     w2x = np.cos(phi) * sx + np.sin(phi) * sy - (2 * eta * np.cos(phi)) / g
#     v2y = -np.sin(phi) * sx + (1 + np.cos(phi)) * sy + (2 * eta * np.sin(phi)) / g
#     w2y = -np.sin(phi) * sx + np.cos(phi) * sy

#     ux = np.cos(phi) * sx - np.sin(phi) * sy
#     uy = np.sin(phi) * sx + np.cos(phi) * sy

#     M = J * K * np.array([
#         [sz, 0, v1x, sz * np.cos(phi), sz * np.sin(phi), 0],
#         [0, sz, v1y, -sz * np.sin(phi), sz * np.cos(phi), 0],
#         [-w1x, -w1y, 0, -ux, -uy, 0],
#         [sz * np.cos(phi), sz * np.sin(phi), 0, sz, 0, v2x],
#         [-sz * np.sin(phi), sz * np.cos(phi), 0, 0, sz, v2y],
#         [-ux, -uy, 0, -w2x, -w2y, 0]
#     ])
    
#     eigenvalues = np.linalg.eigvals(M)
#     return np.sum(eigenvalues.real > 0)  # Count positive real parts

# # Compute phase diagrams
# instability_eta_phi = np.zeros((len(eta_vals), len(phi_vals)))
# instability_g_phi = np.zeros((len(g_vals), len(phi_vals)))

# for i, eta in enumerate(eta_vals):
#     for j, phi in enumerate(phi_vals):
#         instability_eta_phi[i, j] = count_unstable_eigenvalues(phi, eta, g=0.5)

# for i, g in enumerate(g_vals):
#     for j, phi in enumerate(phi_vals):
#         instability_g_phi[i, j] = count_unstable_eigenvalues(phi, eta=0.4, g=g)

# # Define integer levels (here the expected number of unstable modes ranges from 0 to 6)
# levels = np.arange(-0.5, 7.5, 1)

# # Plot phase diagrams
# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# c1 = ax1.contourf(phi_vals, eta_vals, instability_eta_phi, levels=levels, cmap='viridis')
# cb1 = fig.colorbar(c1, ax=ax1, label='Number of Unstable Modes', ticks=np.arange(0, 7, 1))
# ax1.set_ylabel(r'$\eta$')
# ax1.set_title(r'Phase Diagram: $\eta$ vs $\varphi$')

# c2 = ax2.contourf(phi_vals, g_vals, instability_g_phi, levels=levels, cmap='plasma')
# cb2 = fig.colorbar(c2, ax=ax2, label='Number of Unstable Modes', ticks=np.arange(0, 7, 1))
# ax2.set_ylabel(r'$g$')
# ax2.set_xlabel(r'$\varphi$')
# ax2.set_title(r'Phase Diagram: $g$ vs $\varphi$')

# plt.tight_layout()
# plt.show()



def compute_eigenvalues(phi, η_val, g_val):
    γ = g_val / 2
    if η_val < γ:
        sx = η_val * np.cos(phi) / γ
        sy = -η_val * np.sin(phi) / γ
        sz = -np.sqrt(1 - η_val**2 / γ**2)
    elif η_val >= γ:
        sx = γ * np.cos(phi) / η_val - np.sin(phi) * np.sqrt(1 - γ**2 / η_val**2)
        sy = -γ * np.sin(phi) / η_val - np.cos(phi) * np.sqrt(1 - γ**2 / η_val**2)
        sz = 0

    v1x = (1 + np.cos(phi)) * sx + np.sin(phi) * sy - (2 * η_val) / g_val
    w1x = np.cos(phi) * sx + np.sin(phi) * sy - (2 * η_val) / g_val
    v1y = np.sin(phi) * sx + (1 + np.cos(phi)) * sy
    w1y = -np.sin(phi) * sx + np.cos(phi) * sy

    v2x = (1 + np.cos(phi)) * sx + np.sin(phi) * sy - (2 * η_val * np.cos(phi)) / g_val
    w2x = np.cos(phi) * sx + np.sin(phi) * sy - (2 * η_val * np.cos(phi)) / g_val
    v2y = -np.sin(phi) * sx + (1 + np.cos(phi)) * sy + (2 * η_val * np.sin(phi)) / g_val
    w2y = -np.sin(phi) * sx + np.cos(phi) * sy

    ux = np.cos(phi) * sx - np.sin(phi) * sy
    uy = np.sin(phi) * sx + np.cos(phi) * sy

    K = 0.5  # κ = 1, so K = 0.5
    J = (g_val**2) / (K**2)
    M = J * K * np.array([
        [sz, 0, v1x, sz * np.cos(phi), sz * np.sin(phi), 0],
        [0, sz, v1y, -sz * np.sin(phi), sz * np.cos(phi), 0],
        [-w1x, -w1y, 0, -ux, -uy, 0],
        [sz * np.cos(phi), sz * np.sin(phi), 0, sz, 0, v2x],
        [-sz * np.sin(phi), sz * np.cos(phi), 0, 0, sz, v2y],
        [-ux, -uy, 0, -w2x, -w2y, 0]
    ])
    
    return np.linalg.eigvals(M)

# Generate η vs φ phase diagram
g_fixed = 0.1
eta_vals = np.linspace(0.01, 0.3, 100)
phi_vals = np.linspace(-np.pi, np.pi, 100)
counts_eta = np.zeros((len(eta_vals), len(phi_vals)))

for i, eta in enumerate(eta_vals):
    for j, phi in enumerate(phi_vals):
        eigs = compute_eigenvalues(phi, eta, g_fixed)
        counts_eta[i, j] = np.sum(np.real(eigs) > 1e-12)

# Generate g vs φ phase diagram
eta_fixed = 0.3
g_vals = np.linspace(0.01, 1, 100)
counts_g = np.zeros((len(g_vals), len(phi_vals)))

for i, g in enumerate(g_vals):
    for j, phi in enumerate(phi_vals):
        eigs = compute_eigenvalues(phi, eta_fixed, g)
        counts_g[i, j] = np.sum(np.real(eigs) > 1e-12)

# Plotting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# η vs φ plot
c1 = ax1.pcolormesh(phi_vals, eta_vals, counts_eta, shading='auto', cmap='viridis')
ax1.set_xlabel(r'$\varphi$')
ax1.set_ylabel(r'$\eta$')
ax1.set_title(r'$\eta$ vs $\varphi$')
fig.colorbar(c1, ax=ax1, label='Unstable Modes')

# g vs φ plot
c2 = ax2.pcolormesh(phi_vals, g_vals, counts_g, shading='auto', cmap='viridis')
ax2.set_xlabel(r'$\varphi$')
ax2.set_ylabel(r'$g$')
ax2.set_title(r'$g$ vs $\varphi$')
fig.colorbar(c2, ax=ax2, label='Unstable Modes')

plt.tight_layout()
plt.show()