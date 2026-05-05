import numpy as np
import matplotlib.pyplot as plt



# Loading files

slow_speeds_sigma_04_If_00 = np.loadtxt('slow_particle_speeds_sigma_0.4_If_0.0')
slow_speeds_sigma_04_If_004 = np.loadtxt('slow_particle_speeds_sigma_0.4_If_0.04')

fast_speeds_sigma_04_If_00 = np.loadtxt('fast_particle_speeds_sigma_0.4_If_0.0')
fast_speeds_sigma_04_If_004 = np.loadtxt('fast_particle_speeds_sigma_0.4_If_0.04')



time = np.arange(0,900)
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['font.size'] = 18
 



time = np.arange(0, 9, 0.01)

fig, axes = plt.subplots(4, 2, figsize=(10, 16), gridspec_kw={'width_ratios': [1, 1.7]})
plt.subplots_adjust(hspace=0.4, wspace=0.05)

###################################################
# --- Row 1: If = 0, slow ---
axes[0, 1].plot(time, slow_speeds_sigma_04_If_00, linewidth=1)
axes[0, 1].axhline(y=0.5, color='r', linewidth=1, label='Intrinsic speed')
axes[0, 1].axhline(y=np.mean(slow_speeds_sigma_04_If_00), color='b', linewidth=1, label='Mean speed')
axes[0, 1].set_title(r"$\sigma = 0.4, I_f = 0.0$, Intrinsic speed = 0.5", fontsize=18)
axes[0, 1].set_yticks([])

# Rotated histogram (left column)
axes[0, 0].hist(slow_speeds_sigma_04_If_00, bins=100, fc=(0, 1, 0, 0.5),
                edgecolor='black', density=True, histtype='stepfilled', orientation='horizontal')
axes[0, 0].axhline(y=0.5, color='r', linewidth=1, label='Intrinsic speed')
axes[0, 0].axhline(y=np.mean(slow_speeds_sigma_04_If_00), color='b', linewidth=1,label='Mean speed')
axes[0, 0].set_xlim([0,5])
axes[0,0].set_ylabel('(a)', rotation = 'horizontal', labelpad = 20)
axes[0, 0].invert_xaxis()  # makes it "attached" to right plot visually
axes[0, 0].legend(fontsize=10, loc='upper left')

###################################################
# --- Row 2: If = 0, fast ---
axes[1, 1].plot(time, fast_speeds_sigma_04_If_00, linewidth=1)
axes[1, 1].axhline(y=3.0, color='r', linewidth=1)
axes[1, 1].axhline(y=np.mean(fast_speeds_sigma_04_If_00), color='b', linewidth=1)
axes[1, 1].set_title(r"$\sigma = 0.4, I_f = 0.0$, Intrinsic speed = 3.0", fontsize=18)
axes[1, 1].set_yticks([])

axes[1, 0].hist(fast_speeds_sigma_04_If_00, bins=100, fc=(0, 1, 0, 0.5),
                edgecolor='black', density=True, histtype='stepfilled', orientation='horizontal')
axes[1, 0].axhline(y=3.0, color='r', linewidth=1)
axes[1, 0].axhline(y=np.mean(fast_speeds_sigma_04_If_00), color='b', linewidth=1)
axes[1, 0].set_xlim([0,3])
axes[1, 0].set_ylabel('(b)',rotation = 'horizontal', labelpad = 20)
axes[1, 0].invert_xaxis()

###################################################
# --- Row 3: If = 0.04, slow ---
axes[2, 1].plot(time, slow_speeds_sigma_04_If_004, linewidth=1)
axes[2, 1].axhline(y=0.5, color='r', linewidth=1)
axes[2, 1].axhline(y=np.mean(slow_speeds_sigma_04_If_004), color='b', linewidth=1)
axes[2, 1].set_title(r"$\sigma = 0.4, I_f = 0.04$, Intrinsic speed = 0.5", fontsize=18)
axes[2, 1].set_yticks([])

axes[2, 0].hist(slow_speeds_sigma_04_If_004, bins=100, fc=(0, 1, 0, 0.5),
                edgecolor='black', density=True, histtype='stepfilled', orientation='horizontal')
axes[2, 0].axhline(y=0.5, color='r', linewidth=1)
axes[2, 0].axhline(y=np.mean(slow_speeds_sigma_04_If_004), color='b', linewidth=1)

axes[2, 0].set_ylabel('(c)', rotation = 'horizontal', labelpad = 20)
axes[2, 0].invert_xaxis()

###################################################
# --- Row 4: If = 0.04, fast ---
axes[3, 1].plot(time, fast_speeds_sigma_04_If_004, linewidth=1)
axes[3, 1].axhline(y=3.0, color='r', linewidth=1)
axes[3, 1].axhline(y=np.mean(fast_speeds_sigma_04_If_004), color='b', linewidth=1)
axes[3, 1].set_title(r"$\sigma = 0.4, I_f = 0.04$, Intrinsic speed = 3.0", fontsize=18)
axes[3, 1].set_yticks([])
axes[3, 1].set_xlabel("Time", fontsize=22, labelpad = 15)

axes[3, 0].hist(fast_speeds_sigma_04_If_004, bins=100, fc=(0, 1, 0, 0.5),
                edgecolor='black', density=True, histtype='stepfilled', orientation='horizontal')
axes[3, 0].axhline(y=3.0, color='r', linewidth=1)
axes[3, 0].axhline(y=np.mean(fast_speeds_sigma_04_If_004), color='b', linewidth=1)

axes[3, 0].set_ylabel('(d)',rotation = 'horizontal', labelpad = 20)
axes[3, 0].invert_xaxis()
axes[3, 0].set_xlabel("Density", fontsize=22, labelpad = 15)

###################################################
# --- Common labels ---
fig.text(0.06, 0.5, "Speed", va='center', rotation='vertical', fontsize=18)
plt.savefig('Figure 5.pdf', dpi = 600, transparent = True)
plt.show()



