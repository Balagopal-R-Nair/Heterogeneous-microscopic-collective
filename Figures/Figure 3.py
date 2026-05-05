import pandas as pd
import matplotlib.pyplot as plt




num_part = 60


## Loading required files obtained from running model ensemble for different parameters
########## sigma = 0 ##################
observed_speed_60_If_004_sigma0 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.0_If_0.04',header = None)

observed_speed_60_If_004_sigma0 = observed_speed_60_If_004_sigma0.to_numpy()


# ############ sigma = 0.2 ##################
observed_speed_60_If_003_sigma02 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.2_If_0.03',header = None)
observed_speed_60_If_004_sigma02 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.2_If_0.04',header = None)

observed_speed_60_If_003_sigma02 = observed_speed_60_If_003_sigma02.to_numpy()
observed_speed_60_If_004_sigma02 = observed_speed_60_If_004_sigma02.to_numpy()


# # ########### sigma = 0.4 ##################
observed_speed_60_If_002_sigma04 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.4_If_0.02',header = None)
observed_speed_60_If_004_sigma04 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.4_If_0.04',header = None)

observed_speed_60_If_002_sigma04 = observed_speed_60_If_002_sigma04.to_numpy()
observed_speed_60_If_004_sigma04 = observed_speed_60_If_004_sigma04.to_numpy()



# # ########### sigma = 0.6 ##################

observed_speed_60_If_004_sigma06 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.6_If_0.04',header = None)

observed_speed_60_If_004_sigma06 = observed_speed_60_If_004_sigma06.to_numpy()



########### sigma = 0.8 ##################
observed_speed_60_If_004_sigma08 = pd.read_csv(f'trial1_observed_speed_{num_part}_part_sigma_0.8_If_0.04',header = None)
observed_speed_60_If_004_sigma08 = observed_speed_60_If_004_sigma08.to_numpy()


 




plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['font.size'] = 18
 





################################################################


plt.hist(observed_speed_60_If_004_sigma0,bins=100,fc=(1, 0, 0, 0),edgecolor='red',density = True, histtype='stepfilled',label = r'$\sigma = 0.0$')
plt.hist(observed_speed_60_If_004_sigma02,bins=100,fc=(0, 1, 1, 0),edgecolor='green',density = True, histtype='stepfilled',label = r'$\sigma = 0.2$')
plt.hist(observed_speed_60_If_004_sigma04,bins=100,fc=(0, 1, 0, 0),edgecolor='blue',density = True, histtype='stepfilled',label = r'$\sigma = 0.4$')
plt.hist(observed_speed_60_If_004_sigma06,bins=100,fc=(1, 1, 0, 0),edgecolor='orange',density = True, histtype='stepfilled',label = r'$\sigma = 0.6$')
plt.hist(observed_speed_60_If_004_sigma08,bins=100,fc=(0, 0, 1, 0),edgecolor='black',density = True, histtype='stepfilled',label = r'$\sigma = 0.8$')
plt.xlim([0,8])
plt.ylim([0,1])
plt.title(r'$I_f=0.04$')
plt.ylabel('Density')
plt.xlabel('Speed')
plt.legend()
plt.tight_layout(pad=1.5)

# plt.savefig('If = 0.04_distinguishing.pdf',dpi=600, bbox_inches = 'tight')
plt.show()


plt.hist(observed_speed_60_If_002_sigma04,bins=100,fc=(0, 1, 1, 0),edgecolor='red',density = True, histtype='stepfilled',label = r'$\sigma = 0.4, I_f = 0.02$')
plt.hist(observed_speed_60_If_003_sigma02,bins=100,fc=(0, 0, 1, 0),edgecolor='green',density = True, histtype='stepfilled',label = r'$\sigma = 0.2, I_f = 0.03$')
plt.hist(observed_speed_60_If_004_sigma0,bins=100,fc=(1, 1, 0, 0),edgecolor='blue',density = True, histtype='stepfilled',label = r'$\sigma = 0.0, I_f = 0.04$')
plt.legend()
plt.title(r'Varying $I_f$')
plt.xlabel('Speed')
plt.xlim([0,8])
plt.ylim([0,1])
plt.tight_layout(pad=1.5)

# plt.savefig('varying If_distinguishing.pdf', dpi = 600, bbox_inches = 'tight')
plt.xlabel('Speed')
plt.show()
############################################################