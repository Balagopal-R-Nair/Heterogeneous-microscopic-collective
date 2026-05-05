import numpy as np
import random as rn
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.mlab as mlab



num_part = 60
num_iter = 1000
realizations = 100
initialize = 100
mean = 1
If_pi = 0.01


# Loading required files obtained from running model ensemble for different parameters

######## OBSERVED SPEEDS##############
########## sigma = 0 ##################

observed_speed_60_If_00_sigma0 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.0_If_0.0',header = None)
observed_speed_60_If_001_sigma0 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.0_If_0.01',header = None)
observed_speed_60_If_002_sigma0 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.0_If_0.02',header = None)
observed_speed_60_If_003_sigma0 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.0_If_0.03',header = None)
observed_speed_60_If_004_sigma0 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.0_If_0.04',header = None)

observed_speed_60_If_00_sigma0 = observed_speed_60_If_00_sigma0.to_numpy()
observed_speed_60_If_001_sigma0 = observed_speed_60_If_001_sigma0.to_numpy()
observed_speed_60_If_002_sigma0 = observed_speed_60_If_002_sigma0.to_numpy()
observed_speed_60_If_003_sigma0 = observed_speed_60_If_003_sigma0.to_numpy()
observed_speed_60_If_004_sigma0 = observed_speed_60_If_004_sigma0.to_numpy()


# ############ sigma = 0.2 ##################
observed_speed_60_If_00_sigma02 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.2_If_0.0',header = None)
observed_speed_60_If_001_sigma02 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.2_If_0.01',header = None)
observed_speed_60_If_002_sigma02 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.2_If_0.02',header = None)
observed_speed_60_If_003_sigma02 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.2_If_0.03',header = None)
observed_speed_60_If_004_sigma02 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.2_If_0.04',header = None)

observed_speed_60_If_00_sigma02 = observed_speed_60_If_00_sigma02.to_numpy()
observed_speed_60_If_001_sigma02 = observed_speed_60_If_001_sigma02.to_numpy()
observed_speed_60_If_002_sigma02 = observed_speed_60_If_002_sigma02.to_numpy()
observed_speed_60_If_003_sigma02 = observed_speed_60_If_003_sigma02.to_numpy()
observed_speed_60_If_004_sigma02 = observed_speed_60_If_004_sigma02.to_numpy()


# # ########### sigma = 0.4 ##################
observed_speed_60_If_00_sigma04 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.4_If_0.0',header = None)
observed_speed_60_If_001_sigma04 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.4_If_0.01',header = None)
observed_speed_60_If_002_sigma04 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.4_If_0.02',header = None)
observed_speed_60_If_003_sigma04 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.4_If_0.03',header = None)
observed_speed_60_If_004_sigma04 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.4_If_0.04',header = None)

observed_speed_60_If_00_sigma04 = observed_speed_60_If_00_sigma04.to_numpy()
observed_speed_60_If_001_sigma04 = observed_speed_60_If_001_sigma04.to_numpy()
observed_speed_60_If_002_sigma04 = observed_speed_60_If_002_sigma04.to_numpy()
observed_speed_60_If_003_sigma04 = observed_speed_60_If_003_sigma04.to_numpy()
observed_speed_60_If_004_sigma04 = observed_speed_60_If_004_sigma04.to_numpy()



# # ########### sigma = 0.6 ##################
observed_speed_60_If_00_sigma06 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.6_If_0.0',header = None)
observed_speed_60_If_001_sigma06 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.6_If_0.01',header = None)
observed_speed_60_If_002_sigma06 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.6_If_0.02',header = None)
observed_speed_60_If_003_sigma06 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.6_If_0.03',header = None)
observed_speed_60_If_004_sigma06 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.6_If_0.04',header = None)

observed_speed_60_If_00_sigma06 = observed_speed_60_If_00_sigma06.to_numpy()
observed_speed_60_If_001_sigma06 = observed_speed_60_If_001_sigma06.to_numpy()
observed_speed_60_If_002_sigma06 = observed_speed_60_If_002_sigma06.to_numpy()
observed_speed_60_If_003_sigma06 = observed_speed_60_If_003_sigma06.to_numpy()
observed_speed_60_If_004_sigma06 = observed_speed_60_If_004_sigma06.to_numpy()



########### sigma = 0.8 ##################
observed_speed_60_If_00_sigma08 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.8_If_0.0',header = None)
observed_speed_60_If_001_sigma08 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.8_If_0.01',header = None)
observed_speed_60_If_002_sigma08 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.8_If_0.02',header = None)
observed_speed_60_If_003_sigma08 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.8_If_0.03',header = None)
observed_speed_60_If_004_sigma08 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.8_If_0.04',header = None)

observed_speed_60_If_00_sigma08 = observed_speed_60_If_00_sigma08.to_numpy()
observed_speed_60_If_001_sigma08 = observed_speed_60_If_001_sigma08.to_numpy()
observed_speed_60_If_002_sigma08 = observed_speed_60_If_002_sigma08.to_numpy()
observed_speed_60_If_003_sigma08 = observed_speed_60_If_003_sigma08.to_numpy()
observed_speed_60_If_004_sigma08 = observed_speed_60_If_004_sigma08.to_numpy()


########### INTRINSIC SPEEDS ##############
#  ########### sigma = 0.2 ##################
intrinsic_speed_60_If_00_sigma02 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.2_If_0.0',header = None)
intrinsic_speed_60_If_001_sigma02 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.2_If_0.01',header = None)
intrinsic_speed_60_If_002_sigma02 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.2_If_0.02',header = None)
intrinsic_speed_60_If_003_sigma02 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.2_If_0.03',header = None)
intrinsic_speed_60_If_004_sigma02 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.2_If_0.04',header = None)

intrinsic_speed_60_If_00_sigma02 = intrinsic_speed_60_If_00_sigma02.to_numpy()
intrinsic_speed_60_If_001_sigma02 = intrinsic_speed_60_If_001_sigma02.to_numpy()
intrinsic_speed_60_If_002_sigma02 = intrinsic_speed_60_If_002_sigma02.to_numpy()
intrinsic_speed_60_If_003_sigma02 = intrinsic_speed_60_If_003_sigma02.to_numpy()
intrinsic_speed_60_If_004_sigma02 = intrinsic_speed_60_If_004_sigma02.to_numpy()





 
############ sigma = 0.4 ################## 

intrinsic_speed_60_If_00_sigma04 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.4_If_0.0',header = None)
intrinsic_speed_60_If_001_sigma04 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.4_If_0.01',header = None)
intrinsic_speed_60_If_002_sigma04 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.4_If_0.02',header = None)
intrinsic_speed_60_If_003_sigma04 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.4_If_0.03',header = None)
intrinsic_speed_60_If_004_sigma04 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.4_If_0.04',header = None)

intrinsic_speed_60_If_00_sigma04 = intrinsic_speed_60_If_00_sigma04.to_numpy()
intrinsic_speed_60_If_001_sigma04 = intrinsic_speed_60_If_001_sigma04.to_numpy()
intrinsic_speed_60_If_002_sigma04 = intrinsic_speed_60_If_002_sigma04.to_numpy()
intrinsic_speed_60_If_003_sigma04 = intrinsic_speed_60_If_003_sigma04.to_numpy()
intrinsic_speed_60_If_004_sigma04 = intrinsic_speed_60_If_004_sigma04.to_numpy()





############ sigma = 0.6 ################## 
intrinsic_speed_60_If_00_sigma06 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.6_If_0.0',header = None)
intrinsic_speed_60_If_001_sigma06 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.6_If_0.01',header = None)
intrinsic_speed_60_If_002_sigma06 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.6_If_0.02',header = None)
intrinsic_speed_60_If_003_sigma06 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.6_If_0.03',header = None)
intrinsic_speed_60_If_004_sigma06 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.6_If_0.04',header = None)

intrinsic_speed_60_If_00_sigma06 = intrinsic_speed_60_If_00_sigma06.to_numpy()
intrinsic_speed_60_If_001_sigma06 = intrinsic_speed_60_If_001_sigma06.to_numpy()
intrinsic_speed_60_If_002_sigma06 = intrinsic_speed_60_If_002_sigma06.to_numpy()
intrinsic_speed_60_If_003_sigma06 = intrinsic_speed_60_If_003_sigma06.to_numpy()
intrinsic_speed_60_If_004_sigma06 = intrinsic_speed_60_If_004_sigma06.to_numpy()




 
############ sigma = 0.8 ##################

intrinsic_speed_60_If_00_sigma08 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.8_If_0.0',header = None)
intrinsic_speed_60_If_001_sigma08 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.8_If_0.01',header = None)
intrinsic_speed_60_If_002_sigma08 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.8_If_0.02',header = None)
intrinsic_speed_60_If_003_sigma08 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.8_If_0.03',header = None)
intrinsic_speed_60_If_004_sigma08 = pd.read_csv(f'trial1_intrinsic_speed_particle_{num_part}_part_sigma_0.8_If_0.04',header = None)

intrinsic_speed_60_If_00_sigma08 = intrinsic_speed_60_If_00_sigma08.to_numpy()
intrinsic_speed_60_If_001_sigma08 = intrinsic_speed_60_If_001_sigma08.to_numpy()
intrinsic_speed_60_If_002_sigma08 = intrinsic_speed_60_If_002_sigma08.to_numpy()
intrinsic_speed_60_If_003_sigma08 = intrinsic_speed_60_If_003_sigma08.to_numpy()
intrinsic_speed_60_If_004_sigma08 = intrinsic_speed_60_If_004_sigma08.to_numpy()




plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['font.size'] = 18
 

############################################################
############################################################

If_values = [0.00, 0.02, 0.04]
If_labels  = ["00", "002", "004"]

sigmas = [0.0,0.4,0.8]
sigma_labels = ["0","04","08"]

ylims = {
    0.0: (0, 2.6),
    0.4: (0, 1.5),
    0.8: (0, 1.2)
}

fig, axes = plt.subplots(3, 3, figsize=(25, 25), sharex=True)

for i, sigma in enumerate(sigmas):
    for j, If in enumerate(If_values):
        ax = axes[i, j]


        If_str = If_labels[j]
        sigma_str = sigma_labels[i]

        obs_name = f"observed_speed_60_If_{If_str}_sigma{sigma_str}"

        try:
            obs_var = globals()[obs_name]
        except KeyError:
            print(f"Missing variable: {obs_name}")
            continue

        ax.hist(obs_var, bins=100, fc=(0, 0, 1, 0.5),
                edgecolor='black', density=True, histtype='stepfilled')

        # Plot intrinsic if available
        if sigma > 0:
            intr_name = f"intrinsic_speed_60_If_{If_str}_sigma{sigma_str}"
            try:
                intr_var = globals()[intr_name]
                ax.hist(intr_var, bins=100, fc=(1, 0, 0, 0.5),
                        density=True, histtype='stepfilled')
            except KeyError:
                print(f" Missing variable: {intr_name}")

        ax.set_xlim(0, 8)
        if not (i == 0 and j == 0):
            ax.set_ylim(ylims[sigma])
        
        ax.tick_params(axis='x', labelsize=35, bottom=True, labelbottom=True)
        if i != len(sigmas) - 1:
            ax.tick_params(axis='x', labelbottom=False)
        if j == 0:
    # First plot in each row
            ax.tick_params(axis='y', labelsize=35, left=True, labelleft=True)
            
            ax.set_ylabel(fr"$\sigma = {sigma:.1f}$", fontsize=35, labelpad=25)
        elif i == 0 and j == 1:
    # Second plot in the first row only
            ax.tick_params(axis='y', labelsize=35, left=True, labelleft=True)
        else:
    # All other subplots: hide y-ticks
            ax.tick_params(axis='y', left=False, labelleft=False)

       

        # Row and column labels
        if j == 0:
            ax.set_ylabel(fr"$\sigma = {sigma:.1f}$", fontsize=40, labelpad=100,fontweight='bold',rotation = 'horizontal')
        if i == 0:
            ax.set_title(fr"$ I_f = {If:.2f}$", fontsize=40, pad=14,fontweight='bold')

# Global labels
fig.text(0.5, 0.04, "Speed", ha='center', fontsize=45)
fig.text(0.04, 0.5, "Density", va='center', rotation='vertical', fontsize=45)

plt.tight_layout(rect=[0.06, 0.06, 1, 1])
# plt.savefig("figure 2.png", transparent = True, dpi=600)
plt.show()

################################################################














