import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


num_part = 60


# Loading required files obtained from running model ensemble for different parameters
########## sigma = 0 ##################
collision_angles_ensemble_60_If_00_sigma0 = pd.read_csv(f'trial1_anglecheck_collision_angle_beta_{num_part}_part_sigma_0.0_If_0.0',header = None)
collision_angles_ensemble_60_If_002_sigma0 = pd.read_csv(f'trial1_anglecheck_collision_angle_beta_{num_part}_part_sigma_0.0_If_0.02',header = None)
collision_angles_ensemble_60_If_004_sigma0 = pd.read_csv(f'trial1_anglecheck_collision_angle_beta_{num_part}_part_sigma_0.0_If_0.04',header = None)

collision_angles_ensemble_60_If_00_sigma0 = collision_angles_ensemble_60_If_00_sigma0.to_numpy()
collision_angles_ensemble_60_If_002_sigma0 = collision_angles_ensemble_60_If_002_sigma0.to_numpy()
collision_angles_ensemble_60_If_004_sigma0 = collision_angles_ensemble_60_If_004_sigma0.to_numpy()



# ########### sigma = 0.6 ##################
collision_angles_ensemble_60_If_00_sigma06 = pd.read_csv(f'trial1_anglecheck_collision_angle_beta_{num_part}_part_sigma_0.6_If_0.0',header = None)
collision_angles_ensemble_60_If_002_sigma06 = pd.read_csv(f'trial1_anglecheck_collision_angle_beta_{num_part}_part_sigma_0.6_If_0.02',header = None)
collision_angles_ensemble_60_If_004_sigma06 = pd.read_csv(f'trial1_anglecheck_collision_angle_beta_{num_part}_part_sigma_0.6_If_0.04',header = None)

collision_angles_ensemble_60_If_00_sigma06 = collision_angles_ensemble_60_If_00_sigma06.to_numpy()
collision_angles_ensemble_60_If_002_sigma06 = collision_angles_ensemble_60_If_002_sigma06.to_numpy()
collision_angles_ensemble_60_If_004_sigma06 = collision_angles_ensemble_60_If_004_sigma06.to_numpy()



plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['font.size'] = 14
 



  


angles = np.array(np.radians([0,90]))


#plotting
datasets_sigma0 = [collision_angles_ensemble_60_If_00_sigma0,
                   collision_angles_ensemble_60_If_002_sigma0,
                   collision_angles_ensemble_60_If_004_sigma0]

datasets_sigma06 = [collision_angles_ensemble_60_If_00_sigma06,
                    collision_angles_ensemble_60_If_002_sigma06,
                    collision_angles_ensemble_60_If_004_sigma06]

colors = ['tab:red', 'tab:green', 'tab:blue']
labels = [r'$I_f = 0.0$', '$I_f = 0.02$', '$I_f = 0.04$']

nbins = 36

fig, axes = plt.subplots(1, 2, subplot_kw={'projection': 'polar'}, figsize=(12, 6))

# sigma = 0.0
ax = axes[0]
for data, color, label in zip(datasets_sigma0, colors, labels):
    counts, bins = np.histogram(data, bins=nbins, range=(0, 2*np.pi), density=True)
    ax.bar(bins[:-1], counts, width=np.diff(bins),
           align='edge', color=color, alpha=0.4,
           edgecolor='black', linewidth=0.5, label=label)

ax.set_thetamin(0)
ax.set_thetamax(180)
ax.set_ylim(0, 0.6)
ax.set_xticks(angles)
ax.set_yticks([0.0,0.2,0.4,0.6])
ax.set_theta_zero_location('E')
ax.set_theta_direction(1)
ax.set_title(r'$\sigma = 0.0$', fontsize=40, pad=15)
ax.tick_params(labelsize=28, pad = 15)
ax.tick_params(axis = 'y', labelrotation = -60)

# sigma = 0.6
ax = axes[1]
for data, color, label in zip(datasets_sigma06, colors, labels):
    counts, bins = np.histogram(data, bins=nbins, range=(0, 2*np.pi), density=True)
    ax.bar(bins[:-1], counts, width=np.diff(bins),
           align='edge', color=color, alpha=0.4,
           edgecolor='black', linewidth=0.5, label=label)

ax.set_thetamin(0)
ax.set_thetamax(180)
ax.set_ylim(0, 0.6)
ax.set_xticks(angles)
ax.set_yticks([0.0,0.2,0.4,0.6])
ax.set_theta_zero_location('E')
ax.set_theta_direction(1)
ax.set_title(r'$\sigma = 0.6$', fontsize=40, pad=15)
ax.tick_params(labelsize=28,  pad = 15)
ax.tick_params(axis = 'y', labelrotation = -60)
ax.legend(loc='upper right', bbox_to_anchor=(1.2, 1.1), fontsize=20)



fig.text(0.5, 0.1, 'Density', ha='center', va='center', fontsize=40)
fig.subplots_adjust(wspace=3.0)

plt.tight_layout(rect=[0, 0.05, 1, 1])  # leave space for bottom label
# plt.savefig('Figure 6a.pdf', dpi =600, bbox_inches = 'tight')
plt.show()







