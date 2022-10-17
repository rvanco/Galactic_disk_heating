#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 17:22:05 2022

@author: renaudvanco
"""


#%%
###########################################################################################################################
            # PACKAGE :
###########################################################################################################################

import numpy as np
import matplotlib.pyplot as plt
import random
import imageio
import os

#%%
###########################################################################################################################
            # Usefull definition :
###########################################################################################################################

        # Limits of the graph

def lim_graph(list_X, list_Y, list_Z) :
    x_max = np.max(list_X) + 0.1*np.max(list_X)
    x_min = np.min(list_X) + 0.1*np.min(list_X)
    x_lim = max([x_max, abs(x_min)])
    
    y_max = np.max(list_Y) + 0.1*np.max(list_Y)
    y_min = np.min(list_Y) + 0.1*np.min(list_Y)
    y_lim = max([y_max, abs(y_min)])
    
    return x_lim, y_lim

        # Switch between cartesian and polar coordinate :

def cart2cyl(x, y, z):
    """
    Convert from cartesian coordinate to cylindrical coordinate
    """
    rho = np.sqrt(x**2 + y**2 + z**2)
    phi = np.arctan2(y, x) # for negative result of phi just add 2*np.pi
    if phi < 0 :
        phi += 2*np.pi 
    return rho, phi

def cyl2cart(rho, phi):
    """
    Convert from Cylindrical coordinate to Cartesian coordinate
    """
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return x, y

#%%
###########################################################################################################################
            # Constant :
###########################################################################################################################

        # Potential :
    
G = 6.67430*10**(-11)       # m^3.kg^−1.s^−2
M_pot = 10**(10) # Solar mass
a = 2.5
b = a/20

# Perturber :
M_pertub = M_pot*10**(-3)


#%%
###########################################################################################################################
            # Miyamoto-Nagai potential :
###########################################################################################################################

def Miyamoto_Nagai(x_lim, y_lim, a, b, z) :
    
    x = np.linspace(-x_lim, x_lim, 50)
    y = np.linspace(-y_lim, y_lim, 50)
    X_pot, Y_pot = np.meshgrid(x, y)
    
    R = np.sqrt(X_pot**2 + Y_pot**2)
    S = np.sqrt( R**2 + ( a**2 + np.sqrt( b**2 + z**2 ) )**2 )
    Z_pot = - (G*M_pot)/(S)
    
    return X_pot, Y_pot, Z_pot

        # 2D :

def R_pt_pt (R, theta_pt, Z, G, M, a, b) :
    S = np.sqrt( R**2 + ( a**2 + np.sqrt( b**2 + Z**2 ) )**2 )
    grad_phi_r = G*M * (R/S**3)
    R_pp = - grad_phi_r + R*(theta_pt**2)
    return R_pp

def theta_pt_pt (R, R_pt, theta_pt, G, M) :
    theta_pp = - 2*R_pt*theta_pt/R
    return theta_pp

def Z_pt_pt (R, Z, G, M, a, b) :
    S = np.sqrt( R**2 + ( a**2 + np.sqrt( b**2 + Z**2 ) )**2 )
    grad_phi_z = G*M * ((Z*(a + np.sqrt(b**2 + Z**2))) / ((S**3)*np.sqrt( b**2 + Z**2 )))
    Z_pp = - grad_phi_z
    
    return Z_pp

def ang_mom(R, theta_dot) :
    mom_anguaire = (R**2)*theta_dot
    return mom_anguaire


#%%
###########################################################################################################################
            # RK4 (2D) :
###########################################################################################################################

# Parameter of integration :

t_max = 500
dt = 0.1
nb_star = 100
r_min = 0.5
r_max = 16

# Integration :

def RK4_2D (r_min, r_max, t_max, dt, nb_star) :
    
    list_X = []
    list_Y = []
    list_Z = []
    angul_mom = []
    time = []
    
    for i in range(0, nb_star) :
        
        R = random.uniform(r_min, r_max)
        theta = random.uniform(0, 2*np.pi)
        Z = 0
        print("R = ", R, "theta = ", theta)#, "z = ", Z)
        V_R = 0
        V_theta = 0.1/(R+0.1)
        V_Z = 0.01/(R+0.1)
        t = 0
    
        list_X.append([])
        list_Y.append([])
        list_Z.append([])
        angul_mom.append([])
        list_V_R = []
        list_V_theta = []
        list_V_Z = []
        
        epsilon = 0.00001
        while t < t_max+epsilon :
            
            if i==0: #### We register the time####
                time.append(t)
            
            X, Y = cyl2cart(R, theta)
            
            
            list_X[i] += [X]
            list_Y[i] += [Y]
            list_Z[i].append([Z])
            list_V_R += [V_R]
            list_V_theta += [V_theta]
            list_V_Z += [V_Z]
            angul_mom[i].append(ang_mom(R, V_theta))
            
            
            k_R1 = V_R*dt
            k_theta1 = V_theta*dt
            k_Z1 = V_Z*dt
            #k_Z1 = 0
            k_VR1 = R_pt_pt(R, V_theta, Z, G, M_pot, a, b) * dt
            k_Vtheta1 = theta_pt_pt(R, V_R, V_theta, G, M_pot) * dt
            k_VZ1 = Z_pt_pt (R, Z, G, M_pot, a, b) * dt
            
            k_R2 = (V_R + (1/2)*k_VR1) * dt
            k_theta2 = (V_theta + (1/2)*k_Vtheta1) * dt
            k_Z2 = (V_Z + (1/2)*k_VZ1) * dt
            #k_Z2 = 0
            k_VR2 = R_pt_pt(R+(1/2)*k_R1, V_theta+(1/2)*k_Vtheta1, Z+(1/2)*k_Z1, G, M_pot, a, b) * dt
            k_Vtheta2 = theta_pt_pt(R+(1/2)*k_R1, V_R+(1/2)*k_VR1, V_theta+(1/2)*k_Vtheta1, G, M_pot) * dt
            k_VZ2 = Z_pt_pt (R+(1/2)*k_R1, Z+(1/2)*k_Z1, G, M_pot, a, b) * dt

            k_R3 = (V_R + (1/2)*k_VR2) * dt
            k_theta3 = (V_theta + (1/2)*k_Vtheta2) * dt
            k_Z3 = (V_Z + (1/2)*k_VZ2) * dt
            #k_Z3 = 0
            k_VR3 = R_pt_pt(R+(1/2)*k_R2, V_theta+(1/2)*k_Vtheta2, Z+(1/2)*k_Z2, G, M_pot, a, b) * dt
            k_Vtheta3 = theta_pt_pt(R+(1/2)*k_R2, V_R+(1/2)*k_VR2, V_theta+(1/2)*k_Vtheta2, G, M_pot) * dt
            k_VZ3 = Z_pt_pt (R+(1/2)*k_R2, Z+(1/2)*k_Z2, G, M_pot, a, b) * dt

            k_R4 = (V_R + k_VR3) * dt
            k_theta4 = (V_theta + k_Vtheta3) * dt
            k_Z4 = (V_Z + (1/2)*k_VZ3) * dt
            k_VR4 = R_pt_pt(R+k_R3, V_theta+k_Vtheta3, Z+(1/2)*k_Z3, G, M_pot, a, b) * dt
            k_Vtheta4 = theta_pt_pt(R+k_R3, V_R+k_VR3, V_theta+k_Vtheta3, G, M_pot) * dt
            k_VZ4 = Z_pt_pt (R, Z, G, M_pot, a, b) * dt

        
            R += 1/6 * (k_R1 + 2*k_R2 + 2*k_R3 + k_R4)
            theta += 1/6 * (k_theta1 + 2*k_theta2 + 2*k_theta3 + k_theta4)
            Z += 1/6 * (k_Z1 + 2*k_Z2 + 2*k_Z3 + k_Z4)
            V_R += 1/6 * (k_VR1 + 2*k_VR2 + 2*k_VR3 + k_VR4)
            V_theta += 1/6 * (k_Vtheta1 + 2*k_Vtheta2 + 2*k_Vtheta3 + k_Vtheta4)
            V_Z += 1/6 * (k_VZ1 + 2*k_VZ2 + 2*k_VZ3 + k_VZ4)


            t += dt
    
    return list_X, list_Y, list_Z, angul_mom, time

#%%
###########################################################################################################################
            # Perturbator :
###########################################################################################################################



#%%
###########################################################################################################################
            # Visualisation 2D :
###########################################################################################################################


# Choose the case :
list_X, list_Y, list_Z, angul_mom, time = RK4_2D(r_min, r_max, t_max, dt, nb_star)

        # Angular momentum :
#%%
plt.figure()
for i in range(0, nb_star) :
    plt.scatter(time, angul_mom[i], s=2) #, color = "black")

plt.show()

        # X/Y trajectories in a potential :

x_lim, y_lim = lim_graph(list_X, list_Y, list_Z)
z = 0
X_pot, Y_pot, Z_pot = Miyamoto_Nagai(x_lim, y_lim, a, b, z)

plt.figure()
plt.pcolormesh(X_pot, Y_pot, Z_pot, shading="gouraud")
plt.colorbar()  
for i in range(len(list_X)) :
    plt.scatter(list_X[i], list_Y[i], s = 4) #, color = "black")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()
#%%
        # X/Z trajectories in a potential :

plt.figure()
 
for i in range(len(list_X)) :
    plt.scatter(list_X[i], list_Z[i], s = 4, color = "grey", alpha = 0.5)

plt.xlim(-16.5, 16.5)
#plt.ylim(-0.02, 0.02)
plt.xlabel("X")
plt.ylabel("Z")
plt.show()

#%%

        # 2D GIF :

filenames = []

for j in range(len(list_X[0])): # go through each time step
    if j % 10 == 0: 
        plt.figure()
        
        plt.pcolormesh(X_pot, Y_pot, Z_pot, shading="gouraud")
        plt.colorbar()  
        
        for S in range(len(list_X)): # go through each star

            plt.scatter(list_X[S][j], list_Y[S][j], s = 10, color='red')

            plt.xlim(-x_lim, x_lim)
            plt.ylim(-y_lim, y_lim)
            
            plt.title('Orbites of '+str(nb_star)+' stars in a Miyamoto-Nagai potential')
									
        filename=f'/Users/renaudvanco/Documents/Etude/Master/S3/Numerical_simulation/Plot_traj/{j/10}.png'
        filenames.append(filename)
        plt.savefig(filename)
        plt.close()
        print('plt.scatter at t = '+str(j*dt), "over " + str(t_max), " done")

i = 0
with imageio.get_writer('/Users/renaudvanco/Documents/Etude/Master/S3/Numerical_simulation/traj.gif',mode='I') as writer: 
    for filename in filenames: 
        image=imageio.imread(filename)
        writer.append_data(image)
        print("file n° ", i," add to the gif")
        i += 1

for filename in set(filenames):
    os.remove(filename)
    print("plt.scatter deleted")


#%%





#%%

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

for i in range(len(list_X)) :
    ax.scatter(list_X[i], list_Y[i], list_Z[i], s=2)
plt.show()










































