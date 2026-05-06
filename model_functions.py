import numpy as np
import matplotlib as plt
from tqdm import tqdm



#Initialises positions of every agent at random positions within the box. 
#The position array also stores the intrinsic speeds of the agents in the 3rd column
def init_position(position,num_part,box_length,mean,sigma):
    
    position[:,0] = np.random.uniform(0, box_length, num_part)
    position[:,1] = np.random.uniform(0, box_length, num_part) 

    s = np.sqrt(np.log(1 + sigma**2 / mean**2))
    m = np.log(mean) - s**2 / 2 
    
    position[:,2] = np.random.lognormal(m,s,num_part)        # lognormal distribution for intrinsic speeds
    
    
#Initialises orientations of the agents
def init_angles(num_part,angles):
    angles[:] = np.random.uniform(-np.pi, np.pi, num_part)

 
#Takes distance between two agents (x and y) as input and outputs corresponding distance with periodic condition met   
def apply_periodic_boundary(dx, dy, box_length):
    dx = (dx - box_length * np.round(dx / box_length))
    dy = (dy - box_length * np.round(dy / box_length))

    return dx, dy  
     
# Calculates the translation and rotation velocity field acting on every agent    
def field(angles,position,num_part,box_length,If,R):
   d = 0.01
   
   # Creates NxN array with i,j containing the x distance and y distance between j and i
   dx_ij, dy_ij = position[:,np.newaxis,0] - position[np.newaxis,:,0], position[:,np.newaxis,1] - position[np.newaxis,:,1]
   dx_ij, dy_ij = apply_periodic_boundary(dx_ij, dy_ij, box_length)
    
   
   rho_ij = np.abs(dx_ij**2 + dy_ij**2)  # NxN array containing distance squared between i and j in i,j
   np.fill_diagonal(rho_ij,np.inf)       # At diagonals i=j which gives rho=0. We later use rho in the denominator, 
                                         # so to prevent it from blowing up we set it to inf 
   

   
   rho_ij[rho_ij <= 2*R*2*R] = 2*R*2*R           # To have the velocity field plateau for inter agent distance less than 2*R
   
   
   dot = dx_ij*np.cos(angles) + dy_ij*np.sin(angles)
   det = dx_ij*np.sin(angles) - dy_ij*np.cos(angles)

   theta_ji = np.arctan2(det, dot)       # The angle between heading of agent j and line joining i and j 


   cos = np.cos(theta_ji)
   sin = np.sin(theta_ji)
   
   # calculating the vectors in the equation
   ejrho = np.stack(([dx_ij/np.sqrt(rho_ij),dy_ij/np.sqrt(rho_ij)]), axis = -1)
     
   ejtheta = np.stack([ejrho[:,:,1], -ejrho[:,:,0]], axis=-1)
   
   term1 = (ejrho*cos[:,:,np.newaxis])
   term2 = (ejtheta)*sin[:,:,np.newaxis]
   
   v_0 = position[:,2]
   Uij = -If*(term1 + term2)/(rho_ij[:,:,np.newaxis]) # We take a minus because of the vectors are calculated in such a way that -ve needs to be taken. 
   Uij = v_0[:,np.newaxis,np.newaxis]*Uij             # The velocity field is linear to intrinsic speed of the particle causing it
   Uij = np.sum(Uij,axis=0)                           # Taking sum over all the agents
 


   # Calculating rotational field
   # To calculate del of u_ij, we take 4 points around the focal agent and calculate u_ij at those points 
   # to get del of u_ij at the focal agent via finite difference
   
   
     #              2
     #           1     4            
     #              3 
     
     
   dx1, dy1 = (dx_ij - d), dy_ij
   dx2, dy2 = dx_ij, (dy_ij + d)
   dx3, dy3 = dx_ij, (dy_ij - d)
   dx4, dy4 = (dx_ij + d), dy_ij
          
      
            
   rho1 = dx1**2 + dy1**2
   dot1 = (dx1*np.cos(angles) + dy1*np.sin(angles))
   det1 = (dx1*np.sin(angles) - dy1*np.cos(angles))
   theta_ji_1 = np.arctan2(det1,dot1)
   cos1 = np.cos(theta_ji_1)
   sin1 = np.sin(theta_ji_1)
            
   ejrho1 = np.stack(([dx1/np.sqrt(rho1),dy1/np.sqrt(rho1)]), axis = -1)  
   ejtheta1 = np.stack([ejrho1[:,:,1], -ejrho1[:,:,0]], axis=-1)
   term1_1 = (ejrho1*cos1[:,:,np.newaxis])
   term2_1 = (ejtheta1)*sin1[:,:,np.newaxis]
            
    ##
                        
   rho2 = dx2**2 + dy2**2
   dot2 = (dx2*np.cos(angles) + dy2*np.sin(angles))
   det2 = (dx2*np.sin(angles) - dy2*np.cos(angles))
   theta_ji_2 = np.arctan2(det2,dot2)
   cos2 = np.cos(theta_ji_2)
   sin2 = np.sin(theta_ji_2)
            
   ejrho2 = np.stack(([dx2/np.sqrt(rho2),dy2/np.sqrt(rho2)]), axis = -1)
   ejtheta2 = np.stack([ejrho2[:,:,1], -ejrho2[:,:,0]], axis=-1)
   term1_2 = (ejrho2*cos2[:,:,np.newaxis])
   term2_2 = (ejtheta2)*sin2[:,:,np.newaxis]        
            

    ##
            
   rho3 = dx3**2 + dy3**2
   dot3 = (dx3*np.cos(angles) + dy3*np.sin(angles))
   det3 = (dx3*np.sin(angles) - dy3*np.cos(angles))
   theta_ji_3 = np.arctan2(det3,dot3)
   cos3 = np.cos(theta_ji_3)
   sin3 = np.sin(theta_ji_3)

   ejrho3 = np.stack(([dx3/np.sqrt(rho3),dy3/np.sqrt(rho3)]), axis = -1)
   ejtheta3 = np.stack([ejrho3[:,:,1], -ejrho3[:,:,0]], axis=-1)
   term1_3 = (ejrho3*cos3[:,:,np.newaxis])
   term2_3 = (ejtheta3)*sin3[:,:,np.newaxis]


    ##        
            
   rho4 = dx4**2 + dy4**2
   dot4 = (dx4*np.cos(angles) + dy4*np.sin(angles))
   det4 = (dx4*np.sin(angles) - dy4*np.cos(angles))
   theta_ji_4 = np.arctan2(det4,dot4)
   cos4 = np.cos(theta_ji_4)
   sin4 = np.sin(theta_ji_4)

   ejrho4 = np.stack(([dx4/np.sqrt(rho4),dy4/np.sqrt(rho4)]), axis = -1)
   ejtheta4 = np.stack([ejrho4[:,:,1], -ejrho4[:,:,0]], axis=-1)
   term1_4 = (ejrho4*cos4[:,:,np.newaxis])
   term2_4 = (ejtheta4)*sin4[:,:,np.newaxis]

    ##    
            
   Uij1 = -If*(term1_1 + term2_1)/rho1[:,:,np.newaxis]
   Uij1 = v_0[:,np.newaxis,np.newaxis]*Uij1
   Uij2 = -If*(term1_2 + term2_2)/rho2[:,:,np.newaxis]
   Uij2 = v_0[:,np.newaxis,np.newaxis]*Uij2
   Uij3 = -If*(term1_3 + term2_3)/rho3[:,:,np.newaxis]
   Uij3 = v_0[:,np.newaxis,np.newaxis]*Uij3
   Uij4 = -If*(term1_4 + term2_4)/rho4[:,:,np.newaxis]
   Uij4 = v_0[:,np.newaxis,np.newaxis]*Uij4
              
            
   ei_par = np.stack([np.cos(theta_ji), np.sin(theta_ji)], axis=2)  

   ei_perp = np.stack([np.cos(theta_ji + np.pi / 2), np.sin(theta_ji + np.pi / 2)], axis=2) 
 
        
        
   del_uij = np.zeros((num_part, num_part, 2, 2))
   del_uij[:,:,0,0] = (Uij1[:,:,0] - Uij4[:,:,0])/(2*d)  
   del_uij[:,:,0,1] = (Uij2[:,:,0] - Uij3[:,:,0])/(2*d) 
   del_uij[:,:,1,0] = (Uij1[:,:,1] - Uij4[:,:,1])/(2*d) 
   del_uij[:,:,1,1] = (Uij2[:,:,1] - Uij3[:,:,1])/(2*d) 
        

   # Calculating omega
   deluij_eiperp = np.empty([num_part,num_part,2])
   
   
   deluij_eiperp[:,:,0] = del_uij[:,:,0,0]*ei_perp[:,:,0] + del_uij[:,:,0,1]*ei_perp[:,:,1]
   deluij_eiperp[:,:,1] = del_uij[:,:,1,0]*ei_perp[:,:,0] + del_uij[:,:,1,1]*ei_perp[:,:,1]
 
   omega = np.empty([num_part,num_part])

   omega[:,:] = ei_par[:,:,0]*deluij_eiperp[:,:,0] +  ei_par[:,:,1]*deluij_eiperp[:,:,1]
   np.fill_diagonal(omega, 0)
   omega = np.sum(omega,axis=0)
   

   return(Uij,omega)





# Spring collision force
def inter_force(force,position,num_part,box_length,rad_int,spr_const):
    k = spr_const
    for i in range(num_part):
        for j in range(num_part):
            if i == j:
                continue
            dx = (position[i,0] - position[j,0])
            dy = (position[i,1] - position[j,1])
            dx,dy = apply_periodic_boundary(dx, dy, box_length)
            n = np.arctan2(dy,dx) # Force acts along the line joining the two agents
            
            r = np.sqrt(dx**2 + dy**2)
            
            if (r <= 2*rad_int):                # force is only added to the focal agent if the other agent is within 2*rad_int 
                spr_F = k*(2*rad_int - r)       # magnitude of force 
                force[i,0]+= spr_F*np.cos(n)
                force[i,1]+= spr_F*np.sin(n)

          
            
# updates velocity of every agent
def velocity_update(velocity,angles,position,force,num_part,box_length,dt,If,rad_int,eta):
    R = rad_int
    Uij, omega = field(angles,position,num_part,box_length,If,R)

        
    random = np.random.normal(0,eta,num_part)
    angles += omega*dt + random                      # updating angular velocity with noise
    
    angles = (angles + np.pi) % (2 * np.pi) - np.pi  # resets angles of every agent to have value from -pi to pi
    
    velocity[:,0] = position[:,2]*np.cos(angles) + (force[:,0])  + Uij[:,0]
    velocity[:,1] = position[:,2]*np.sin(angles) + (force[:,1])  + Uij[:,1]
    
    

    

 # updates positions                
def position_update(velocity,angles,position,num_part,box_length,dt):

    l = box_length 
   
    position[:,0:2] += velocity*dt
    
    # Accounting for periodic boundaries
    position[:, 0] = np.mod(position[:, 0], l)
    position[:, 1] = np.mod(position[:, 1], l)
    


    

# Returns an NxN truth table for pairs of agents that are in contact with each other    
def collision_check(position,R,box_length):
    dx_ij, dy_ij = position[:,np.newaxis,0] - position[np.newaxis,:,0], position[:,np.newaxis,1] - position[np.newaxis,:,1]
    dx_ij, dy_ij = apply_periodic_boundary(dx_ij, dy_ij, box_length)
     

    rho_ij = np.abs(dx_ij**2 + dy_ij**2)  # NxN array containing distance squared between i and j in i,j
    np.fill_diagonal(rho_ij,np.inf)       # At diagonals i=j which gives rho=0. To prevent counting self such case we set it to inf
    
    collision_check = rho_ij <= (2*R)**2     
    
    
    return(collision_check)  


# runs the model once for <num_iter> time steps, initialized for <initialize> time steps 
# and stores speed, intrinsic speed, x,y,positions, x,y,velocities in arrays to be used for plotting/analysis
# also stores position data in .xyz file format to use for animating in OVITO
def run_model_once(If,sigma,mean,box_length,rad_int,eta,spr_const,num_part,num_iter,initialize,position,velocity,angles,force):
    dt = 0.01 # time step
    data = [] # stores data for animating in OVITO
    
    # Initializing arrays to store data in
    total_speed = np.empty([(num_iter - initialize),num_part])            # Stores speed of all the agents in all time steps
    total_intrinsic_speed = np.empty([(num_iter - initialize),num_part])  # Stores intrinsic speed of all the agents in all time steps
    
    # uncomment the data array you wish to store and uncomment corresponding ones below too
    # total_position_x = np.empty([(num_iter - initialize),num_part]) 
    # total_position_y = np.empty([(num_iter - initialize),num_part])
    # total_velocity_x = np.empty([(num_iter - initialize),num_part])
    # total_velocity_y = np.empty([(num_iter - initialize),num_part])
    # total_angle =  np.empty([(num_iter - initialize),num_part])

    
    # Initializing
    init_position(position,num_part,box_length,mean,sigma)     
    init_angles(num_part,angles)  
    
    
    # Arrays in which data is stored for each step to calculate for next step or to store in larger set
    speed = np.empty([num_part,1])          
    force = np.zeros([num_part,2])
    velocity = np.empty([num_part,2])
    
    


    ###########################
    
    velocity_update(velocity,angles,position,force,num_part,box_length,dt,If,rad_int,eta) # calculates velocity of each agent based on the initial position and orientation
    total_collisions = 0
     
    for i in tqdm(range(num_iter)): 
       force = np.zeros([num_part,2])
       inter_force(force,position,num_part,box_length,rad_int,spr_const)
                
       velocity_update(velocity,angles,position,force,num_part,box_length,dt,If,rad_int,eta)
       speed = np.sqrt(np.square(velocity[:,0]) + np.square(velocity[:,1]))   
        
                      
       coll_check_old = collision_check(position,rad_int,box_length)  # stores collisions that happened in the previous time step to compare with current
       position_update(velocity,angles,position,num_part,box_length,dt) 
   
       if i >= initialize:
           ii = i-initialize
           
           # This loop is for when we need to store the number of collisions that occur in the system
           if ii>0:
               coll_check_current = collision_check(position,rad_int,box_length)
               
               # checks whether whether the agents were in contact in the previous step, only if they were not are they considered a collision
               continuous_coll_check = (coll_check_current.astype(np.float32) - coll_check_old.astype(np.float32))
               continuous_coll_check[continuous_coll_check<0] = 0
               
               total_collisions = np.sum(continuous_coll_check) # can be printed later if required
       
            # storing data    
           total_speed[ii,:] = speed
           total_intrinsic_speed[ii,:] = position[:,2]  
           
           # total_position_x[ii,:] = position[:,0]  # uncomment the data array you wish to store
           # total_position_y[ii,:] = position[:,1] 
           # total_velocity_x[ii,:] = velocity[:,0]
           # total_velocity_y[ii,:] = velocity[:,1]
           # total_angle[ii,:] = angles
            
           # Storing data for animating in OVITO                                                                               
           snapshot = [f"{num_part}\n Properties=Particle Type:R:1pos:R:3:Velocity:R:3\n"]
           snapshot.extend([f"1 {position[0, 0]} {position[0, 1]} 0 {np.cos(angles[0])} {np.sin(angles[0])} 0  \n"])
           snapshot.extend([f"2 {position[1, 0]} {position[1, 1]} 0 {np.cos(angles[1])} {np.sin(angles[1])} 0  \n"])
           snapshot.extend([f"3 {position[j, 0]} {position[j, 1]} 0 {np.cos(angles[j])} {np.sin(angles[j])} 0  \n" for j in range(2,num_part)])
           data.append("".join(snapshot))
              
            
    with open(f"sim_{num_part}part_sigma_{sigma}_If_{If}.xyz", "w") as f:
        f.write("".join(data))   
    f.close()                        
    
    # flattens the data arrays to make them easier to print to files, or use in histograms. 
    # comment out and replace the return function with appropriate data arrays if you dont want to flatten
    flat_total_speed = (total_speed).flatten()
    flat_total_intrinsic_speed = (total_intrinsic_speed).flatten()
    # flat_total_position_x = total_position_x.flatten()   # uncomment the data array you wish to store
    # flat_total_position_y = total_position_y.flatten()   # also add the corresponding data array to the return
    # flat_total_velocity_x = total_velocity_x.flatten()
    # flat_total_velocity_y = total_velocity_y.flatten()
    # flat_total_angle =  total_angle.flatten()

    return(flat_total_intrinsic_speed,flat_total_speed) # returns data arrays, add arrays that you wish to use 
    

####################################################################################

#runs the model for <realizations> realizations, each running for <num_iter> time steps, initialized for <initialize> time steps 
#and stores speed, intrinsic speed, x,y,positions, x,y,velocities in arrays to be used for plotting/analysis
def run_model_ensemble(If,sigma,mean,box_length,rad_int,eta,spr_const,num_part,num_iter,initialize,realizations,position,velocity,angles,force):    
    
    dt = 0.01 # time step
    
    # Initializing arrays to store data in
    total_speed_for_avg = np.empty([realizations,num_part,(num_iter - initialize)])         # stores observed speeds of agents for entire ensemble that is later averaged for each realization, used in Figure 4
    total_intrinsic_speed_for_avg = np.empty([realizations,num_part])                       # stores intrinsic speeds of agents for each realization
    total_speed_ens = np.empty([realizations,(num_iter - initialize),num_part])             # stores observed speeds of entire ensemble
    total_intrinsic_speed_ens = np.empty([realizations,(num_iter - initialize),num_part])   # stores instrinsic speeds for entire ensemble
    total_angle_ens =  np.empty([realizations,(num_iter - initialize),num_part])            # stores heading angles of all the agents for the entire ensemble
    continuous_collision_counter_realization = np.empty([realizations])                     # number of collisions in a realization
    # total_position_x_ens = np.empty([realizations,(num_iter - initialize),num_part])
    # total_position_y_ens = np.empty([realizations,(num_iter - initialize),num_part])
    # total_velocity_x_ens = np.empty([realizations,(num_iter - initialize),num_part])
    # total_velocity_y_ens = np.empty([realizations,(num_iter - initialize),num_part])
    
    
    # Arrays in which data is stored each step to calculate for next step or store in larger array
    avg_speed = np.empty([num_iter - initialize])
    avg_speed_ens = np.empty([realizations])
    beta_angle_ensemble = []
    
    for I in tqdm(range(realizations)):
       
        
       init_position(position,num_part,box_length,mean,sigma)     
       init_angles(num_part,angles)
    
    
    
       velocity_update(velocity,angles,position,force,num_part,box_length,dt,If,rad_int,eta)
       
       
       #collision_counter = 0
       continuous_collision_counter = 0
       for i in (range(num_iter)): 
            
            force = np.zeros([num_part,2])
            inter_force(force,position,num_part,box_length,rad_int,spr_const)
            
            velocity_update(velocity,angles,position,force,num_part,box_length,dt,If,rad_int,eta)
            speed = np.sqrt(np.square(velocity[:,0]) + np.square(velocity[:,1]))
            
            coll_check_old = collision_check(position,rad_int,box_length)
            position_update(velocity,angles,position,num_part,box_length,dt)  
            
            
            v0 = position[:,2]
            total_intrinsic_speed_for_avg[I,:] = v0    
            
                        
            if i >= initialize:
                ii = i-initialize
                total_speed_for_avg[I,:,ii] = speed
                
                total_angle_ens[I,ii,:] = angles
                if ii>0:
                    coll_check_current = collision_check(position,rad_int,box_length)
                    
                    continuous_coll_check = (coll_check_current.astype(np.float32) - coll_check_old.astype(np.float32))
                    continuous_coll_check[continuous_coll_check<0] = 0
                    
                    if np.any(continuous_coll_check):
                        coll_pairs = np.where(continuous_coll_check == 1)
                      
                        beta_index1 = coll_pairs[0]
                        beta_index2 = coll_pairs[1]
                        
                        for k in range(len(beta_index1)):
                            # collision angle (beta) is defined as angle between agent heading and line joining 2 agents
                            dx = position[beta_index2[k],0] - position[beta_index1[k],0]
                            dy = position[beta_index2[k],1] - position[beta_index1[k],1]
                            dx, dy = apply_periodic_boundary(dx, dy, box_length)
                            
                            rho = np.sqrt(dx**2 + dy**2)
                            
                            dx = dx/rho
                            dy = dy/rho
                            
                            dot = dx*np.cos(total_angle_ens[I,ii-1,beta_index1[k]]) + dy*np.sin(total_angle_ens[I,ii-1,beta_index1[k]])
                            det = dx*np.sin(total_angle_ens[I,ii-1,beta_index1[k]]) - dy*np.cos(total_angle_ens[I,ii-1,beta_index1[k]])
                            
                            b = np.arctan2(det,dot)
                            beta_angle_ensemble.append(b) # storing collisions angle data
                        
                    # updating collision counter
                    continuous_collision_count = np.count_nonzero(continuous_coll_check)//2
                    continuous_collision_counter+= continuous_collision_count
                    

                
                # storing data 
                total_speed_ens[I,ii,:] = speed     
                avg_speed[ii] = np.average(speed)      
                total_intrinsic_speed_ens[I,ii,:] = position[:,2]
            
                # total_position_x_ens[I,ii,:] = position[:,0]  # uncomment the data array you wish to store
                # total_position_y_ens[I,ii,:] = position[:,1] 
                # total_velocity_x_ens[I,ii,:] = velocity[:,0]
                # total_velocity_y_ens[I,ii,:] = velocity[:,1]
                
                
            
            
            

       continuous_collision_counter_realization[I] = continuous_collision_counter
       avg_speed_ens[I] = np.average(avg_speed) 
    
    
    avg_speed_per_particle = np.average(total_speed_for_avg,axis = 2)  # Average speeds of each agent in every realization
    total_avg_speed_ens = np.average(avg_speed_ens)                    # Average of all speeds of all the agents across all realizations 
    beta_angle_ensemble = np.array(beta_angle_ensemble)                # converts the list of collision angles into numpy array for easier storing 
    
    
    # flattens the data arrays to make them easier to print to files, or use in histograms. 
    # comment out and replace the return array with appropriate data arrays if you dont want to flatten            
    flat_total_speed_ens = (total_speed_ens).flatten()
    flat_total_intrinsic_speed_ens = (total_intrinsic_speed_ens).flatten()
    flat_avg_speed_per_particle = (avg_speed_per_particle).flatten()
    flat_total_intrinsic_speed_for_avg = (total_intrinsic_speed_for_avg).flatten()
    flat_collision_angle_beta = (beta_angle_ensemble).flatten
    # flat_total_position_x_ens = total_position_x_ens.flatten()   # uncomment the data array you wish to store
    # flat_total_position_y_ens = total_position_y_ens.flatten()   # also add the corresponding data array to the return
    # flat_total_velocity_x_ens = total_velocity_x_ens.flatten()
    # flat_total_velocity_y_ens = total_velocity_y_ens.flatten()
    # flat_total_angle_ens =  total_angle_ens.flatten()
    
    # Printing avg no. of collisions per realization and avg speed of an agent
    avg_coll = np.average(continuous_collision_counter_realization)  # taking average of number of collisions in each realization
    print(f'avg no. of collisions: {avg_coll}')
    # print(f'avg speed of agents: {total_avg_speed_ens}')
    
    
    # Saving data to be used for plotting
    
    # Data for figure 2 
    np.savetxt(f"observed_speed_particle_60_part_sigma_{sigma}_If_{If}", flat_total_speed_ens) 
    np.savetxt(f"intrinsic_speed_particle_60_part_sigma_{sigma}_If_{If}", flat_total_intrinsic_speed_ens) 
    
    # Data for figure 4
    # np.savetxt(f"avg_speed_particle_60_part_sigma_{sigma}_If_{If}", flat_avg_speed_per_particle) 
    # np.savetxt(f"intrinsic_speed_particle_60_part_sigma_{sigma}_If_{If}", flat_total_intrinsic_speed_for_avg)
    
    # Data for figure 6b 
    # np.savetext(f"trial1_anglecheck_collision_angle_beta_60_part_sigma_{sigma}_If_{If}", flat_collision_angle_beta)
    return(flat_total_intrinsic_speed_ens,flat_total_speed_ens,flat_collision_angle_beta)           
               


              
