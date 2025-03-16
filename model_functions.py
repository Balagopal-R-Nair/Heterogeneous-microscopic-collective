import numpy as np
from tqdm import tqdm



#Initialises positions of every agent at random positions within the box. 
#The position array also stores the intrinsic speeds of the agents in the 3rd column
def init_position(position,num_part,box_length,mean,sigma):
    
    position[:,0] = np.random.uniform(0, box_length, num_part)
    position[:,1] = np.random.uniform(0, box_length, num_part) 
    
    
    for i in range(num_part):
      v_0 = np.random.normal(mean,sigma)
      while v_0 < 0: # To ensure that sampled values are positive
         v_0 = np.random.normal(mean,sigma)
      position[i,2] = v_0


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

    

   rho_ij = np.abs(dx_ij**2 + dy_ij**2)  # NxN array containing distance squared between i and j in i,j
   np.fill_diagonal(rho_ij,np.inf)       # At diagonals i=j which gives rho=0. We later use rho in the denominator, 
                                         # so to prevent it from blowing up we set it to inf 
   
   
   rho_ij[rho_ij <= 2*R] = 2*R           # To have the velocity field plateau for inter agent distance less than 2*R
   

   dot = dx_ij*np.cos(angles) + dy_ij*np.sin(angles)
   det = dx_ij*np.sin(angles) - dy_ij*np.cos(angles)

   theta_ji = np.arctan2(det, dot)       # The angle between heading of agent j and line joining i and j 


   cos = np.cos(theta_ji)
   sin = np.sin(theta_ji)
   
   # calculating the vectors in the equation
   ejrho = np.stack(([dx_ij/np.sqrt(rho_ij),dy_ij/np.sqrt(rho_ij)]), axis = -1)
   rot90 = np.array([[0,1],[-1,0]])
   ejtheta = np.empty([num_part,num_part,2])
   ejtheta[:,:,0] = rot90[0,0]*ejrho[:,:,0] + rot90[0,1]*ejrho[:,:,1]
   ejtheta[:,:,1] = rot90[1,0]*ejrho[:,:,0] + rot90[1,1]*ejrho[:,:,1]
    
   term1 = (ejrho*cos[:,:,np.newaxis])
   term2 = (ejtheta)*sin[:,:,np.newaxis]
   
   v_0 = position[:,2]
   Uij = -If*(term1 + term2)/(rho_ij[:,:,np.newaxis]) # We take a minus because of the vectors are calculated in such a way that -ve needs to be taken. 
   Uij = v_0[:,np.newaxis,np.newaxis]*Uij             # The velocity field is linear to intrinsic speed of the particle causing it
   Uij = np.sum(Uij,axis=0)                           # Taking sum over j
 


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
   ejtheta1 = np.empty([num_part,num_part,2])
   ejtheta1[:,:,0] = rot90[0,0]*ejrho1[:,:,0] + rot90[0,1]*ejrho1[:,:,1]
   ejtheta1[:,:,1] = rot90[1,0]*ejrho1[:,:,0] + rot90[1,1]*ejrho1[:,:,1]
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
   ejtheta2 = np.empty([num_part,num_part,2])
   ejtheta2[:,:,0] = rot90[0,0]*ejrho2[:,:,0] + rot90[0,1]*ejrho2[:,:,1]
   ejtheta2[:,:,1] = rot90[1,0]*ejrho2[:,:,0] + rot90[1,1]*ejrho2[:,:,1]
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
   ejtheta3 = np.empty([num_part,num_part,2])
   ejtheta3[:,:,0] = rot90[0,0]*ejrho3[:,:,0] + rot90[0,1]*ejrho3[:,:,1]
   ejtheta3[:,:,1] = rot90[1,0]*ejrho3[:,:,0] + rot90[1,1]*ejrho3[:,:,1]
   term1_3 = (ejrho3*cos3[:,:,np.newaxis])
   term2_3 = (ejtheta1)*sin3[:,:,np.newaxis]


    ##        
            
   rho4 = dx4**2 + dy4**2
   dot4 = (dx4*np.cos(angles) + dy4*np.sin(angles))
   det4 = (dx4*np.sin(angles) - dy4*np.cos(angles))
   theta_ji_4 = np.arctan2(det4,dot4)
   cos4 = np.cos(theta_ji_4)
   sin4 = np.sin(theta_ji_4)

   ejrho4 = np.stack(([dx4/np.sqrt(rho4),dy1/np.sqrt(rho4)]), axis = -1)
   ejtheta4 = np.empty([num_part,num_part,2])
   ejtheta4[:,:,0] = rot90[0,0]*ejrho4[:,:,0] + rot90[0,1]*ejrho4[:,:,1]
   ejtheta4[:,:,1] = rot90[1,0]*ejrho4[:,:,0] + rot90[1,1]*ejrho4[:,:,1]
   term1_4 = (ejrho4*cos4[:,:,np.newaxis])
   term2_4 = (ejtheta1)*sin4[:,:,np.newaxis]

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

        
    random = np.random.uniform(-eta/2,eta/2,num_part)
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
    
    

# runs the model once for <num_iter> time steps, initialized for <initialize> time steps 
# and stores speed, intrinsic speed, x,y,positions, x,y,velocities in arrays to be used for plotting/analysis
# also stores position data in .xyz file format to use for animating in OVITO
def run_model_once(If,sigma,mean,box_length,rad_int,eta,spr_const,num_part,num_iter,initialize,position,velocity,angles,force):
    dt = 0.01 # time step
    data = [] # stores data for animating in OVITO
    
    
    total_speed = np.empty([(num_iter - initialize),num_part])
    total_intrinsic_speed = np.empty([(num_iter - initialize),num_part])
    # total_position_x = np.empty([(num_iter - initialize),num_part])  # uncomment the data array you wish to store and uncomment corresponding ones below too
    # total_position_y = np.empty([(num_iter - initialize),num_part])
    # total_velocity_x = np.empty([(num_iter - initialize),num_part])
    # total_velocity_y = np.empty([(num_iter - initialize),num_part])
    # total_angle =  np.empty([(num_iter - initialize),num_part])

    init_position(position,num_part,box_length,mean,sigma)     
    init_angles(num_part,angles)
    
    
    velocity_update(velocity,angles,position,force,num_part,box_length,dt,If,rad_int,eta) # calculates velocity of each agent based on the initial position and orientation
    for i in tqdm(range(num_iter)): 
       force = np.zeros([num_part,2])
       inter_force(force,position,num_part,box_length,rad_int,spr_const)
            
       velocity_update(velocity,angles,position,force,num_part,box_length,dt,If,rad_int,eta)
       speed = np.sqrt(np.square(velocity[:,0]) + np.square(velocity[:,1]))
             
                        
                        
       position_update(velocity,angles,position,num_part,box_length,dt)            
                        
       if i >= initialize:
           ii = i-initialize
                
            # storing data    
           total_speed[ii,:] = speed
           total_intrinsic_speed[ii,:] = position[:,2]  
           # total_position_x[ii,:] = position[:,0]  # uncomment the data array you wish to store
           # total_position_y[ii,:] = position[:,1] 
           # total_velocity_x[ii,:] = velocity[:,0]
           # total_velocity_y[ii,:] = velocity[:,1]
           # total_angle[ii,:] = angles
            
           # Storing data for animating in OVITO                                                                    
           snapshot = [f"{num_part}\nParticles\n"]
           snapshot.extend([f"{j} {position[j, 0]} {position[j, 1]} 0  \n" for j in range(0,num_part)])
           data.append("".join(snapshot))
              
            
    with open(f"{num_part}part_sigma_{sigma}_If_{If}.xyz", "w") as f:
        f.write("".join(data))   
    f.close()                        
    
    # flattens the data arrays to make them easier to print to files, or use in histograms. 
    # comment out and replace the return array with appropriate data arrays if you dont want to flatten
    flat_total_speed = (total_speed).flatten()
    flat_total_intrinsic_speed = (total_intrinsic_speed).flatten()
    # flat_total_position_x = total_position_x.flatten()   # uncomment the data array you wish to store
    # flat_total_position_y = total_position_y.flatten()   # also add the corresponding data array to the return
    # flat_total_velocity_x = total_velocity_x.flatten()
    # flat_total_velocity_y = total_velocity_y.flatten()
    # flat_total_angle =  total_angle.flatten()


    return(flat_total_intrinsic_speed,flat_total_speed) # returns data arrays



#runs the model for <realizations> realizations, each running for <num_iter> time steps, initialized for <initialize> time steps 
#,nd stores speed, intrinsic speed, x,y,positions, x,y,velocities in arrays to be used for plotting/analysis
def run_model_ensemble(If,sigma,mean,box_length,rad_int,eta,spr_const,num_part,num_iter,initialize,realizations,position,velocity,angles,force):    
    
    dt = 0.01 # time step
    total_speed_ens = np.empty([realizations,(num_iter - initialize),num_part])
    total_intrinsic_speed_ens = np.empty([realizations,(num_iter - initialize),num_part])
    # total_position_x_ens = np.empty([realizations,(num_iter - initialize),num_part])
    # total_position_y_ens = np.empty([realizations,(num_iter - initialize),num_part])
    # total_velocity_x_ens = np.empty([realizations,(num_iter - initialize),num_part])
    # total_velocity_y_ens = np.empty([realizations,(num_iter - initialize),num_part])
    # total_angle_ens =  np.empty([realizations,(num_iter - initialize),num_part])
    
    for I in tqdm(range(realizations)):
       
        
       init_position(position,num_part,box_length,mean,sigma)     
       init_angles(num_part,angles)
    
    
    
       velocity_update(velocity,angles,position,force,num_part,box_length,dt,If,rad_int,eta)
    
       for i in (range(num_iter)): 
            force = np.zeros([num_part,2])
            inter_force(force,position,num_part,box_length,rad_int,spr_const)
            
            velocity_update(velocity,angles,position,force,num_part,box_length,dt,If,rad_int,eta)
            speed = np.sqrt(np.square(velocity[:,0]) + np.square(velocity[:,1]))
                  
            
                        
                        
            position_update(velocity,angles,position,num_part,box_length,dt)            
                        
            if i >= initialize:
                ii = i-initialize
                
                # storing data 
                total_speed_ens[I,ii,:] = speed     
                total_intrinsic_speed_ens[I,ii,:] = position[:,2]
                # total_position_x_ens[I,ii,:] = position[:,0]  # uncomment the data array you wish to store
                # total_position_y_ens[I,ii,:] = position[:,1] 
                # total_velocity_x_ens[I,ii,:] = velocity[:,0]
                # total_velocity_y_ens[I,ii,:] = velocity[:,1]
                # total_angle_ens[I,ii,:] = angles
                
    # flattens the data arrays to make them easier to print to files, or use in histograms. 
    # comment out and replace the return array with appropriate data arrays if you dont want to flatten            
    flat_total_speed_ens = (total_speed_ens).flatten()
    flat_total_intrinsic_speed_ens = (total_intrinsic_speed_ens).flatten()
    # flat_total_position_x_ens = total_position_x_ens.flatten()   # uncomment the data array you wish to store
    # flat_total_position_y_ens = total_position_y_ens.flatten()   # also add the corresponding data array to the return
    # flat_total_velocity_x_ens = total_velocity_x_ens.flatten()
    # flat_total_velocity_y_ens = total_velocity_y_ens.flatten()
    # flat_total_angle_ens =  total_angle_ens.flatten()

    return(flat_total_intrinsic_speed_ens,flat_total_speed_ens)           
               


              
