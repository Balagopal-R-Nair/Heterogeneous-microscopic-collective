import numpy as np
import matplotlib.pyplot as plt



from model_functions import run_model_ensemble,run_model_once



####################################################

box_length = 3          # Length of the domain
num_part = 60           # Number of agents in the simulation
num_iter = 1000         # Number of time steps
realizations = 100      # Number of realizaions
initialize = 100        # Number of initialization time steps; no data is stored during these time steps
mean = 1                # Mean of the  log normal distribution from which intrinsic speeds are sampled
sigma = 0.2             # Standard deviation of the log normal distribution from which intrinsic speeds are sampled (i.e., extent of heterogeneity)
eta = 0.1               # Variance of normal distribution from which the noise in orientation of agent is sampled (mean zero)
If = 0.0               # Strength of hydrodynamic field
spr_const = 50          # Spring constant of the collsion force 
rad_int = 0.1           # Radius of agent


# Stores corresponding quantities for a time step
position = np.empty([num_part,3]) # also stores intrinsic speeds
angles =  np.empty([num_part,])
velocity = np.empty([num_part,2])
speed = np.empty([num_part,1])
force = np.zeros([num_part,2])

####################################################

print('Parameters \n')
print(f'No. of particles = {num_part}, No. of iterations = {num_iter}, No. of realizations = {realizations}, No. of initializations steps = {initialize} \n')
print(f'spring const = {spr_const}, mean = {mean}')
print(f'sigma = {sigma} , If = {If}')

# Running simulation once
intrinsic_speeds, observed_speeds = run_model_once(If,sigma,mean,box_length,rad_int,eta,spr_const,num_part,num_iter,initialize,position,velocity,angles,force)

# with open(f"slow_particle_speeds_sigma_{sigma}_If_{If}") as f:
#     print(observed_speeds)


# Running for multiple realizations
# intrinsic_speeds,observed_speeds,collision_angles = run_model_ensemble(If,sigma,mean,box_length,rad_int,eta,spr_const,num_part,num_iter,initialize,realizations,position,velocity,angles,force)



