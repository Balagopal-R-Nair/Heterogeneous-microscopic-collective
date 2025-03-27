import numpy as np

from model_functions import run_model_ensemble,run_model_once




####################################################

box_length = 3          # Length of the square arena
num_part = 30           # Number of agents in teh simulation
num_iter = 1000         # Number of time steps
realizations = 2        # Number of realizaions
initialize = 100        # Number of initialization time steps; no data is stored during these time steps
mean = 1                # Mean of the normal distribution from which intrinsic speeds are sampled
sigma = 0.4             # Standard deviation of the normal distribution from which intrinsic speeds are sampled (i.e., extent of heterogeneity)
eta = 0.1               # Size of uniform distribution from which the noise in orientation of agent is sampled
If = 0.02               # Strength of hydrodynamic field
spr_const = 50          # Spring constant of the collsion force 
rad_int = 0.1           # Radius of agent


# Stores corresponding quantities for a time step
position = np.empty([num_part,3]) # also stores intrinsic speeds
angles =  np.empty([num_part,])
velocity = np.empty([num_part,2])
speed = np.empty([num_part,1])
force = np.zeros([num_part,2])

####################################################


intrinsic_speeds, observed_speeds = run_model_once(If,sigma,mean,box_length,rad_int,eta,spr_const,num_part,num_iter,initialize,position,velocity,angles,force)


#intrinsic_speeds_ens, observed_speeds_ens = run_model_ensemble(If,sigma,mean,box_length,rad_int,eta,spr_const,num_part,num_iter,initialize,realizations,position,velocity,angles,force)
