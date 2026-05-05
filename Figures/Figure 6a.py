import numpy as np
import matplotlib.pyplot as plt



plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = "Helvetica"
plt.rcParams['font.size'] = 12



I_f = np.array([0,0.01,0.02,0.03,0.04])

# Obtained from running run_model_ensemble for each set of parameters and manually saving the average collision count that is printed
cont_coll_00 = np.array([986.92,870.63,819.27, 817.37,897.12])
cont_coll_02 = np.array([995.39,906.13,886.34,918.53,986.3])
cont_coll_04 = np.array([1037.69,981.13,987.84,1046.66,1154.33])
cont_coll_06 = np.array([1073.1,1057.81,1097.43,1220.68,1366.93])
cont_coll_08 = np.array([1131.56,1136.93,1210.19,1377.5,1523.55])


plt.plot(I_f,cont_coll_00,label = r'$\sigma = 0.0$')
plt.plot(I_f,cont_coll_02,label = r'$\sigma = 0.2$')
plt.plot(I_f,cont_coll_04,label = r'$\sigma = 0.4$')
plt.plot(I_f,cont_coll_06,label = r'$\sigma = 0.6$')
plt.plot(I_f,cont_coll_08,label = r'$\sigma = 0.8$')
plt.xlabel(r'Hydrodynamic strength $(I_f)$', fontsize =14, labelpad = 18)
plt.ylabel('Average number of collisions', fontsize = 14, labelpad = 18)
plt.legend(fontsize = 14)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.savefig('Figure 6a.pdf', dpi =600, bbox_inches = 'tight')
plt.show()










