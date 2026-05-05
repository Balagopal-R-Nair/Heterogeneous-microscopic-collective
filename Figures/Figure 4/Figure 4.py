import numpy as np
import matplotlib.pyplot as plt




num_part = 60


# Loading required files obtained from running model ensemble for different parameters

# Average speeds
mean_sigma_02_If_00 = np.loadtxt(f'avg_speed_particle_{num_part}_part_sigma_0.2_If_0.0')
mean_sigma_02_If_004 = np.loadtxt(f'avg_speed_particle_{num_part}_part_sigma_0.2_If_0.04')
mean_sigma_04_If_00 = np.loadtxt(f'avg_speed_particle_{num_part}_part_sigma_0.4_If_0.0')
mean_sigma_04_If_004 = np.loadtxt(f'avg_speed_particle_{num_part}_part_sigma_0.4_If_0.04')
mean_sigma_08_If_00 = np.loadtxt(f'avg_speed_particle_{num_part}_part_sigma_0.8_If_0.0')
mean_sigma_08_If_004 = np.loadtxt(f'avg_speed_particle_{num_part}_part_sigma_0.8_If_0.04')

# Intrinsic speeds
v0_sigma_02_If_00 = np.loadtxt(f'intrinsic_speed_particle_{num_part}_part_sigma_0.2_If_0.0')
v0_sigma_02_If_004 = np.loadtxt(f'intrinsic_speed_particle_{num_part}_part_sigma_0.2_If_0.04')
v0_sigma_04_If_00 = np.loadtxt(f'intrinsic_speed_particle_{num_part}_part_sigma_0.4_If_0.0')
v0_sigma_04_If_004 = np.loadtxt(f'intrinsic_speed_particle_{num_part}_part_sigma_0.4_If_0.04')
v0_sigma_08_If_00 = np.loadtxt(f'intrinsic_speed_particle_{num_part}_part_sigma_0.8_If_0.0')
v0_sigma_08_If_004 = np.loadtxt(f'intrinsic_speed_particle_{num_part}_part_sigma_0.8_If_0.04')



# Latex parameters
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['font.size'] = 18
 


# define the datasets
data = {
    0.2: {
        0.0: (v0_sigma_02_If_00, mean_sigma_02_If_00),
        0.04: (v0_sigma_02_If_004, mean_sigma_02_If_004)
    },
    0.4: {
        0.0: (v0_sigma_04_If_00, mean_sigma_04_If_00),
        0.04: (v0_sigma_04_If_004, mean_sigma_04_If_004)
    },
    0.8: {
        0.0: (v0_sigma_08_If_00, mean_sigma_08_If_00),
        0.04: (v0_sigma_08_If_004, mean_sigma_08_If_004)
    }
}

sigmas = [0.2, 0.4, 0.8]
If_values = [0.0, 0.04]



# plotting

fig, axes = plt.subplots(3, 2, figsize=(10, 12))

for i, sigma in enumerate(sigmas):
    # --- Compute limits for each row ---
    all_v0 = np.concatenate([data[sigma][If][0] for If in If_values])
    all_mean = np.concatenate([data[sigma][If][1] for If in If_values])

    lim_max = max(all_v0.max(), all_mean.max())
    margin = 0.05 * lim_max

    # Offset slightly from zero only for bottom row
    if i == len(sigmas) - 1:
        lim_min = -0.05 * lim_max
    else:
        lim_min = 0.0
    lim_max += margin

    for j, If in enumerate(If_values):
        ax = axes[i, j]
        v0, mean_v = data[sigma][If]

        # Scatter points — blue with black edge
        ax.scatter(v0, mean_v, facecolors='#7f7fff', edgecolors='black',
                   s=40, linewidth=0.5, alpha=0.8)

        # Diagonal line
        ax.plot([lim_min, lim_max], [lim_min, lim_max],
                color='gray', linestyle='--', linewidth=1)

        # Apply limits (slightly below zero for bottom row)
        ax.set_xlim(lim_min, lim_max)
        ax.set_ylim(lim_min, lim_max)

        # Titles and labels
        if j == 0:
            ax.set_ylabel(fr'$\sigma = {sigma}$', rotation=0,
                          labelpad=40, fontsize=16, va='center')
        if i == 0:
            ax.set_title(fr'$I_f = {If}$', fontsize=16)

        # Visible ticks on all sides
        ax.tick_params(axis='both', which='both', direction='in',
                       top=True, right=True)

# --- Common labels ---
fig.text(0.6, 0.04, r'Intrinsic Speed', ha='center', fontsize=18)
fig.text(0.04, 0.5, r'Observed Average Speed',
         va='center', rotation='vertical', fontsize=18)

plt.tight_layout(rect=[0.05, 0.05, 1, 0.97])
plt.savefig('Figure 4.pdf', dpi = 600, transparent =True)
plt.show()
