#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 15:42:50 2022

@author: renaudvanco
"""


###########################################################################################################################
            # PACKAGE :
###########################################################################################################################

import numpy as np
import random
import matplotlib.pyplot as plt
import imageio
import os

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

def sph2cart(az, el, r):
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    return x, y, z

def distance_S_P(R_S, R_P, theta_S, theta_P, Z_S, Z_P):
    R_BH = np.sqrt(R_S**2 + R_P**2 - 2*R_S*R_P*np.cos(theta_S - theta_P) + (Z_S - Z_P)**2)
    return R_BH

def distance_S_P_R(R_S, R_P, theta_S, theta_P):
    R_BH_R = np.sqrt(R_S**2 + R_P**2 - 2*R_S*R_P*np.cos(theta_S - theta_P))
    return R_BH_R

###########################################################################################################################
            # Miyamoto-Nagai potential :
###########################################################################################################################

def Miyamoto_Nagai_XZ(G, M_pot, x_lim, z_lim, a, b) :
    
    x = np.linspace(-x_lim, x_lim, 50)
    y = np.linspace(-z_lim, z_lim, 50)
    R_pot, Z_pot = np.meshgrid(x, y)
    
    R = np.sqrt(R_pot**2)
    S = np.sqrt( R**2 + ( a**2 + np.sqrt( b**2 + Z_pot**2 ) )**2 )
    
    V_pot = - (G*M_pot)/(S) + abs(Z_pot)*(a+np.sqrt(Z_pot**2 +b**2))/(S**3*np.sqrt(Z_pot**2 +b**2))
    
    return R_pot, Z_pot, V_pot

def Miyamoto_Nagai_XY(G, M_pot, x_lim, y_lim, a, b) :
    
    x = np.linspace(-x_lim, x_lim, 50)
    y = np.linspace(-y_lim, y_lim, 50)
    X_pot, Y_pot = np.meshgrid(x, y)
    
    R = np.sqrt(X_pot**2 + Y_pot**2)
    S = np.sqrt( R**2 + ( a**2 + np.sqrt(b**2) )**2 )
    
    Z_pot = - (G*M_pot)/(S) 
    
    return X_pot, Y_pot, Z_pot

def accel_star(R_S, R_BH, R_BH_R, R_P, R_pt, theta_P, theta_S, theta_pt, Z_S, Z_P, G, M_pot, M_star, M_perturb, a, b) :
    S = np.sqrt( R_S**2 + ( a**2 + np.sqrt( b**2 + Z_S**2 ) )**2 )
    
    grad_phi_r = (G*M_pot) * (R_S/S**3)
    grad_phi_z = (G*M_pot) * ((Z_S*(a + np.sqrt(b**2 + Z_S**2))) / ((S**3)*np.sqrt( b**2 + Z_S**2 )))

    delta_Z = Z_P - Z_S
    R_R = np.sqrt(R_BH**2 - delta_Z**2)

    theta_S_BH = np.pi - np.arccos((R_S**2 + R_R**2 - R_P**2)/(2*R_S*R_R))

    if np.isnan(theta_S_BH) == True :
        if R_S >= R_P :
            theta_S_BH = np.pi
        else :
            theta_S_BH = 0
        
    theta_Z = np.arccos((delta_Z**2 + R_BH**2 - R_R**2)/(2*delta_Z*R_BH))


    Perturber_R = (G*M_perturb/(R_BH+0.1)**2) * np.sin(theta_Z) * np.cos(theta_S_BH)
    Perturber_Z = (G*M_perturb/(R_BH+0.1)**2) * np.cos(theta_Z)
    Perturber_theta = 0

    R_pp = - grad_phi_r + R_S*(theta_pt**2) + Perturber_R
    theta_pp = - 2*R_pt*theta_pt/R_S + Perturber_theta 
    Z_pp = - grad_phi_z + Perturber_Z
    
    return np.array([R_pp, theta_pp, Z_pp])

def accel_perturb (R, R_pt, theta_pt, Z, G, M_pot, M_perturb, a, b) :
    S = np.sqrt( R**2 + ( a**2 + np.sqrt( b**2 + Z**2 ) )**2 )
    
    grad_phi_r = (G*M_pot) * (R/S**3)
    grad_phi_z = (G*M_pot) * ((Z*(a + np.sqrt(b**2 + Z**2))) / ((S**3)*np.sqrt( b**2 + Z**2 )))

    
    R_pp = - grad_phi_r + R*(theta_pt**2)
    theta_pp = - 2*R_pt*theta_pt/(R+10**-10)
    Z_pp = - grad_phi_z
    
    return np.array([R_pp, theta_pp, Z_pp])

def ang_mom(R, theta_dot) :
    mom_anguaire = (R**2)*theta_dot
    return mom_anguaire

###########################################################################################################################
            # Setting list :
###########################################################################################################################

def Setting_list(nb_star_gal, nb_star_perturb) :
    list_X=[]
    list_Y=[]
    list_X_P = []
    list_Y_P = []
    list_X_SP = []
    list_Y_SP = []
    
    list_R=[]
    list_theta=[]
    list_Z=[]
    list_V_R=[]
    list_V_theta=[]
    list_V_Z=[]
    
    list_R_P = []
    list_theta_P = []
    list_Z_P = []
    list_V_R_P = []
    list_V_theta_P = []
    list_V_Z_P = []

    list_R_SP = []
    list_theta_SP = []
    list_Z_SP = []
    list_V_R_SP = []
    list_V_theta_SP = []
    list_V_Z_SP = []
    
    list_ang_mom_P = []
    list_ang_mom = []
    
    for i in range(0, nb_star_gal) :
        list_X.append([])
        list_Y.append([])
    
        list_R.append([])
        list_theta.append([])
        list_Z.append([])
        list_V_R.append([])
        list_V_theta.append([])
        list_V_Z.append([])
        list_ang_mom.append([])
        
    for i in range(0, nb_star_perturb) :
        list_X_SP.append([])
        list_Y_SP.append([])
    
        list_R_SP.append([])
        list_theta_SP.append([])
        list_Z_SP.append([])
        list_V_R_SP.append([])
        list_V_theta_SP.append([])
        list_V_Z_SP.append([])

        
    return list_X, list_Y, list_X_P, list_Y_P, list_X_SP, list_Y_SP, list_R, list_theta, list_Z, list_V_R, list_V_theta, list_V_Z, list_R_P, list_theta_P, list_Z_P, list_V_R_P, list_V_theta_P, list_V_Z_P, list_R_SP, list_theta_SP, list_Z_SP, list_V_R_SP, list_V_theta_SP, list_V_Z_SP, list_ang_mom, list_ang_mom_P    

###########################################################################################################################
            # Initial Position :
###########################################################################################################################


def position_stars_bulge(r_bulge, i, G, M, a, b) :
    R = random.uniform(0.5, r_bulge)
    theta = random.uniform(0, 2*np.pi)
    elevation = random.uniform(-np.pi/2, np.pi/2)
    X, Y, Z = sph2cart(theta, elevation, R)

    S = np.sqrt( R**2 + ( a**2 + np.sqrt( b**2 + Z**2 ) )**2 )

    V_R = 0
    V_theta = np.sqrt((G*M)/S**3) #0.1/(R+0.1)
    V_Z = 0
        
    print(f"R_{i} = ", np.round(R,2), "theta = ", np.round(theta,2), "Z = ", np.round(Z,2))
    R = np.sqrt(X**2 + Y**2)
        
    return R, theta, Z, V_R, V_theta, V_Z

def position_stars_disk(r_bulge, r_disk, i, G, M, a, b) :        
    R = random.uniform(r_bulge, r_disk)
    theta = random.uniform(0, 2*np.pi)
    Z = random.uniform(-10/R, 10/R)

    S = np.sqrt( R**2 + ( a**2 + np.sqrt( b**2 + Z**2 ) )**2 )

    print(f"R_{i} = ", np.round(R,2), "theta = ", np.round(theta,2), "Z = ", np.round(Z,2))
    V_R = 0
    V_theta = np.sqrt((G*M)/S**3) #0.14/R
    V_Z = 0
    
    return R, theta, Z, V_R, V_theta, V_Z

def position_stars_perturber(r_perturber, theta_perturber, Z_perturber, V_R_perturber, V_theta_perturber, V_Z_perturber, i, G, M_perturb, a, b) :        
    R = random.uniform(r_perturber-1, r_perturber+1)
    theta = random.uniform(theta_perturber-10, theta_perturber+10)
    Z = random.uniform(Z_perturber-1, Z_perturber+1)
    S = np.sqrt( R**2 + ( a**2 + np.sqrt( b**2 + Z**2 ) )**2 )
    print(f"R_{i} = ", np.round(R,2), "theta = ", np.round(theta,2), "Z = ", np.round(Z,2))
    V_R = V_R_perturber
    V_theta = V_theta_perturber + np.sqrt((G*M_perturb)/S**3)
    V_Z = V_Z_perturber
    
    return R, theta, Z, V_R, V_theta, V_Z

def position_perturber(r_disk, G, M, a, b) :
    R = 0
    theta = np.pi
    Z = 0
    
    S = np.sqrt( R**2 + ( a**2 + np.sqrt( b**2 + Z**2 ) )**2 )

    print("R = ", np.round(R,2), "theta = ", np.round(theta,2), "Z = ", np.round(Z,2))
    
    V_R = 0
    V_theta = 0 #np.sqrt((G*M)/S**3)
    V_Z = 0 #-np.sqrt((G*M)/S**3)
    
    return R, theta, Z, V_R, V_theta, V_Z

def position_strait_perturb() :
    X_P = 0
    Y_P = 15
    Z_P = 10
    V_X_P = 0
    V_Y_P = 0
    V_Z_P = -0.01
    return X_P, Y_P, Z_P, V_X_P, V_Y_P, V_Z_P

###########################################################################################################################
            # Integration :
###########################################################################################################################


def RK4_star(G, M_pot, M_star, M_perturb, a, b, dt, R_S, theta_S, Z_S, V_R, V_theta, V_Z, R_P, R_BH, R_BH_R, theta_P, Z_P) :
    
    k_R1 = V_R*dt
    k_theta1 = V_theta*dt
    k_Z1 = V_Z*dt
    k_VR1, k_Vtheta1, k_VZ1 = accel_star(R_S, R_BH, R_BH_R, R_P, V_R, theta_P, theta_S, V_theta, Z_S, Z_P, G, M_pot, M_star, M_perturb, a, b) * dt
    
    k_R2 = (V_R + (1/2)*k_VR1) * dt
    k_theta2 = (V_theta + (1/2)*k_Vtheta1) * dt
    k_Z2 = (V_Z + (1/2)*k_VZ1) * dt
    k_VR2, k_Vtheta2, k_VZ2 = accel_star(R_S+(1/2)*k_R1, R_BH, R_BH_R, R_P, V_R+(1/2)*k_VR1, theta_P, theta_S, V_theta+(1/2)*k_Vtheta1, Z_S+(1/2)*k_Z1, Z_P, G, M_pot, M_star, M_perturb, a, b) * dt
    
    k_R3 = (V_R + (1/2)*k_VR2) * dt
    k_theta3 = (V_theta + (1/2)*k_Vtheta2) * dt
    k_Z3 = (V_Z + (1/2)*k_VZ2) * dt
    k_VR3, k_Vtheta3, k_VZ3 = accel_star(R_S+(1/2)*k_R2, R_BH, R_BH_R, R_P, V_R+(1/2)*k_VR2, theta_P, theta_S, V_theta+(1/2)*k_Vtheta2, Z_S+(1/2)*k_Z2, Z_P, G, M_pot, M_star, M_perturb, a, b) * dt

    k_R4 = (V_R + k_VR3) * dt
    k_theta4 = (V_theta + k_Vtheta3) * dt
    k_Z4 = (V_Z + k_VZ3) * dt
    k_VR4, k_Vtheta4, k_VZ4 = accel_star(R_S+k_R3, R_BH, R_BH_R, R_P, V_R+k_VR3, theta_P, theta_S, V_theta+k_Vtheta3, Z_S+k_Z3, Z_P, G, M_pot, M_star, M_perturb, a, b) * dt


    R_S += 1/6 * (k_R1 + 2*k_R2 + 2*k_R3 + k_R4)
    theta_S += 1/6 * (k_theta1 + 2*k_theta2 + 2*k_theta3 + k_theta4)
    Z_S += 1/6 * (k_Z1 + 2*k_Z2 + 2*k_Z3 + k_Z4)
    V_R += 1/6 * (k_VR1 + 2*k_VR2 + 2*k_VR3 + k_VR4)
    V_theta += 1/6 * (k_Vtheta1 + 2*k_Vtheta2 + 2*k_Vtheta3 + k_Vtheta4)
    V_Z += 1/6 * (k_VZ1 + 2*k_VZ2 + 2*k_VZ3 + k_VZ4)
    return R_S, theta_S, Z_S, V_R, V_theta, V_Z

def RK4_perturber(G, M_pot, M_perturb, a, b, dt, R, theta, Z, V_R, V_theta, V_Z) :

    k_R1 = V_R*dt
    k_theta1 = V_theta*dt
    k_Z1 = V_Z*dt
    k_VR1, k_Vtheta1, k_VZ1 = accel_perturb(R, V_R, V_theta, Z, G, M_pot, M_perturb, a, b) * dt
    
    k_R2 = (V_R + (1/2)*k_VR1) * dt
    k_theta2 = (V_theta + (1/2)*k_Vtheta1) * dt
    k_Z2 = (V_Z + (1/2)*k_VZ1) * dt
    k_VR2, k_Vtheta2, k_VZ2 = accel_perturb(R+(1/2)*k_R1, V_R+(1/2)*k_VR1, V_theta+(1/2)*k_Vtheta1, Z+(1/2)*k_Z1, G, M_pot, M_perturb, a, b) * dt
    
    k_R3 = (V_R + (1/2)*k_VR2) * dt
    k_theta3 = (V_theta + (1/2)*k_Vtheta2) * dt
    k_Z3 = (V_Z + (1/2)*k_VZ2) * dt
    k_VR3, k_Vtheta3, k_VZ3 = accel_perturb(R+(1/2)*k_R2, V_R+(1/2)*k_VR2, V_theta+(1/2)*k_Vtheta2, Z+(1/2)*k_Z2, G, M_pot, M_perturb, a, b) * dt

    k_R4 = (V_R + k_VR3) * dt
    k_theta4 = (V_theta + k_Vtheta3) * dt
    k_Z4 = (V_Z + k_VZ3) * dt
    k_VR4, k_Vtheta4, k_VZ4 = accel_perturb(R+k_R3, V_R+k_VR3, V_theta+k_Vtheta3, Z+k_Z3, G, M_pot, M_perturb, a, b) * dt

        
    R += 1/6 * (k_R1 + 2*k_R2 + 2*k_R3 + k_R4)
    theta += 1/6 * (k_theta1 + 2*k_theta2 + 2*k_theta3 + k_theta4)
    Z += 1/6 * (k_Z1 + 2*k_Z2 + 2*k_Z3 + k_Z4)
    V_R += 1/6 * (k_VR1 + 2*k_VR2 + 2*k_VR3 + k_VR4)
    V_theta += 1/6 * (k_Vtheta1 + 2*k_Vtheta2 + 2*k_Vtheta3 + k_Vtheta4)
    V_Z += 1/6 * (k_VZ1 + 2*k_VZ2 + 2*k_VZ3 + k_VZ4)
    
    return R, theta, Z, V_R, V_theta, V_Z

def strait_perturb(X_P, Y_P, Z_P, V_X_P, V_Y_P, V_Z_P, dt) :
    X = X_P + V_X_P*dt
    Y = Y_P + V_Y_P*dt
    Z = Z_P + V_Z_P*dt

    return X, Y, Z

def rotating_perturb(R_P, theta_P, Z_P, V_theta_P, dt) :
    R = R_P
    theta = theta_P + V_theta_P*dt
    Z = Z_P

    return R, theta, Z

###########################################################################################################################
            # Visualisation 2D :
###########################################################################################################################



        # Angular momentum :

def graph_ang_mom(nb_star, time, list_ang_mom, list_ang_mom_P) :
    plt.figure()
    for i in range(0, nb_star) :
        plt.scatter(time, list_ang_mom[i]/np.max(list_ang_mom[i]), s=2, color="red")
    #plt.scatter(time, list_ang_mom_P/np.max(list_ang_mom_P), s=2, color="black")
    plt.title('Angular momentum of '+str(nb_star)+' stars \n in a Miyamoto-Nagai potential \n')
    plt.savefig("/home/rvancoellie/Bureau/project_num/Ang_mom_normalized.png")
    plt.show()

        # X/Y trajectories in a potential :

def XY_traj(nb_star, list_X, list_X_P, list_Y, list_Y_P, list_Z, G, M_pot, a, b, x_lim, y_lim) :

    X_pot, Y_pot, Z_pot = Miyamoto_Nagai_XY(G, M_pot, x_lim, y_lim, a, b)
    
    plt.figure()
    plt.pcolormesh(X_pot, Y_pot, Z_pot, shading="gouraud", zorder=-1)
    plt.colorbar()  
    for i in range(len(list_X)) :
        plt.scatter(list_X[i], list_Y[i], s = 2)
    
    plt.xlim(-x_lim, x_lim)
    plt.ylim(-y_lim, y_lim)
    
    plt.scatter(list_X_P, list_Y_P, s = 10, color="black")
    
    
    plt.xlabel("X, (kPc)")
    plt.ylabel("Y, (kPc)")
    
    plt.title('trajectories of '+str(nb_star)+' stars in a Miyamoto-Nagai potential')
    plt.savefig("/home/rvancoellie/Bureau/project_num/XY_traj.png")
    
    plt.show()

        # X/Z trajectories in a potential :

def XZ_traj(nb_star, list_X, list_Z, list_X_P, list_Z_P, G, M_pot, x_lim, z_lim, a, b) :
    
    X_pot, Y_pot, Z_pot = Miyamoto_Nagai_XZ(G, M_pot, x_lim, z_lim, a, b)
    
    plt.figure()

    plt.pcolormesh(X_pot, Y_pot, Z_pot, shading="gouraud", zorder=-1)
    plt.colorbar() 

    for i in range(len(list_X)) :
        plt.scatter(list_X[i], list_Z[i], s = 2)
        
    plt.scatter(list_X_P, list_Z_P, s = 10, color="black")
    
    plt.xlim(-x_lim, x_lim)
    plt.ylim(-z_lim, z_lim)
    
    plt.xlabel("X, (kPc)")
    plt.ylabel("Z, (kPc)")
    plt.title('Trajectories of '+str(nb_star)+' stars in a Miyamoto-Nagai potential')
    plt.savefig("/home/rvancoellie/Bureau/project_num/XZ_traj.png")
    plt.show()


            # XY GIF
def XY_GIF(nb_star, list_X, list_X_P, list_X_SP, list_Y, list_Y_P, list_Y_SP, list_Z, list_Z_P, list_Z_SP, G, M_pot, a, b, x_lim, y_lim, dt, t_max, time) :
    image = []
    
    X_pot, Y_pot, Z_pot = Miyamoto_Nagai_XY(G, M_pot, x_lim, y_lim, a, b)
    
    for j in range(len(list_X[0])): # go through each time step
        if j % 10 == 0: 
            plt.figure()
            
            red_nb = 0
            blue_nb = 0
            
            plt.pcolormesh(X_pot, Y_pot, Z_pot, shading="gouraud", zorder=-1)
            plt.colorbar()  
            
            for S in range(len(list_X)): # go through each star of the galaxy
                if abs(list_Z_P[j]) > 10 :
                    if abs(list_Z[S][j]) > abs((2/3)*list_Z_P[j]) :
                        plt.scatter(list_X[S][j], list_Y[S][j], s = 1, color='blue')
                        blue_nb += 1
                    else :
                        plt.scatter(list_X[S][j], list_Y[S][j], s = 1, color='red')
                        red_nb += 1
                else :
                    plt.scatter(list_X[S][j], list_Y[S][j], s = 1, color='red')
                    red_nb = nb_star
                    
            for S in range(len(list_X_SP)) : # go through each star of the perturber
                plt.scatter(list_X_SP[S][j], list_Y_SP[S][j], s = 1, color='violet')
                
            plt.xlim(-x_lim, x_lim)
            plt.ylim(-y_lim, y_lim)
            
            plt.scatter(list_X_P[j], list_Y_P[j], s = 20, color='black')

            plt.scatter(x_lim+1, y_lim+1, s = 1, color='red', label=f"in the disk ({red_nb})")
            plt.scatter(x_lim+1, y_lim+1, s = 1, color='blue', label=f"out of the disk ({blue_nb})")
            plt.scatter(x_lim+1, y_lim+1, s = 1, color='violet', label=f"perturber stars")
            
            plt.text(x_lim,y_lim, "time="+str(round(time[j]/1000, 4))+"Gy")

            plt.title('Orbites of '+str(nb_star)+' stars in a Miyamoto-Nagai potential \n')
            plt.xlabel("X, (kPc)")
            plt.ylabel("Y, (kPc)")		
            plt.legend(loc='upper right')
            
            filename=f'/home/rvancoellie/Bureau/project_num/Plot_traj/{j/10}.png'
            plt.savefig(filename)
            
            image.append(imageio.imread(filename))
            
            plt.close()
            print('XY scatter at t = '+str(j*dt), "over " + str(t_max), " done")
    

    exportname = "/home/rvancoellie/Bureau/project_num/XY.gif"
    kargs = { 'duration': 0.1 }
    imageio.mimsave(exportname, image, 'GIF', **kargs)
    print("The .gif is finish :)")



            # XZ GIF

def XZ_GIF(nb_star, list_X, list_X_P, list_X_SP, list_Y, list_Y_P, list_Y_SP, list_Z, list_Z_P, list_Z_SP, G, M_pot, a, b, x_lim, z_lim, dt, t_max, time) :
    image = []
    
    X_pot, Y_pot, Z_pot = Miyamoto_Nagai_XZ(G, M_pot, x_lim, z_lim, a, b)

    
    for j in range(len(list_X[0])): # go through each time step
        if j % 10 == 0: 
            plt.figure()
            
            red_nb = 0
            blue_nb = 0
            
            plt.pcolormesh(X_pot, Y_pot, Z_pot, shading="gouraud", zorder=-1)
            plt.colorbar()  
            
            for S in range(len(list_X)): # go through each star
                if abs(list_Z_P[j]) > 10 :
                    if abs(list_Z[S][j]) > abs((2/3)*list_Z_P[j]) :
                        plt.scatter(list_X[S][j], list_Z[S][j], s = 1, color='blue')
                        blue_nb += 1
                    else :
                        plt.scatter(list_X[S][j], list_Z[S][j], s = 1, color='red')
                        red_nb += 1
                else :
                    plt.scatter(list_X[S][j], list_Z[S][j], s = 1, color='red')
                    red_nb = nb_star
            
            for S in range(len(list_X_SP)) : # go through each star of the perturber
                plt.scatter(list_X_SP[S][j], list_Z_SP[S][j], s = 10, color='violet')
            
            plt.xlim(-x_lim, x_lim)
            plt.ylim(-z_lim, z_lim)
                
            plt.scatter(list_X_P[j], list_Z_P[j], s = 20, color='black', zorder = -1)

            plt.scatter(x_lim+1, z_lim+1, s = 1, color='red', label=f"in the disk ({red_nb})")
            plt.scatter(x_lim+1, z_lim+1, s = 1, color='blue', label=f"out of the disk ({blue_nb})")
            plt.scatter(x_lim+1, z_lim+1, s = 1, color='violet', label=f"perturber stars")
            
            plt.text(x_lim,z_lim, "time="+str(round(time[j]/1000, 4))+"Gy")
            plt.title('Orbites of '+str(nb_star)+' stars in a Miyamoto-Nagai potential \n')
            plt.xlabel("X, (kPc)")
            plt.ylabel("Z, (kPc)")
            plt.legend(loc='upper right')
    					
            filename=f'/home/rvancoellie/Bureau/project_num/Plot_traj/{j/10}.png'
            plt.savefig(filename)
            
            image.append(imageio.imread(filename))
            
            plt.close()
            print('XZ scatter at t = '+str(j*dt), "over " + str(t_max), " done")
    
    
    exportname = "/home/rvancoellie/Bureau/project_num/XZ.gif"
    kargs = { 'duration': 0.1 }
    imageio.mimsave(exportname, image, 'GIF', **kargs)
    
    print("The .gif is finish :)")


###########################################################################################################################
            # Visualisation 3D :
###########################################################################################################################


def XYZ_GIF (list_X, list_Y, list_Z, list_X_P, list_Y_P, list_Z_P, x_lim, y_lim, z_lim, nb_star, dt, t_max, time) :
    image = []
    
    for j in range(len(list_X[0])): # go through each time step
        if j % 10 == 0: 
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            
            for S in range(len(list_X)): # go through each star
                
                ax.scatter(list_X[S][j], list_Y[S][j], list_Z[S][j], s = 1, color='red')
                
                plt.xlim(-x_lim, x_lim)
                plt.ylim(-y_lim, y_lim)
                ax.set_zlim(-z_lim, z_lim)
                
            ax.scatter(list_X_P[j], list_Y_P[j], list_Z_P[j], s = 10, color='black')
            ax.scatter(0,0,0, s=30, color ='purple')
            
            plt.title('Orbites of '+str(nb_star)+' stars in a Miyamoto-Nagai potential \n')
            plt.xlabel("X, (kPc)")
            plt.ylabel("Y, (kPc)")
            
            ax.set_zlabel("Z, (kPc)")
            
            ax.text(x_lim/2, y_lim, z_lim, "time"+str(round(time[j]/1000, 4))+"Gy")
            
            filename=f'/home/rvancoellie/Bureau/project_num/Plot_traj/{j/10}.png'

            plt.savefig(filename)
            
            image.append(imageio.imread(filename))

            
            plt.close()
            print('XYZ scatter at t = '+str(j*dt), "over " + str(t_max), " done")
    
    exportname = "/home/rvancoellie/Bureau/project_num/XYZ.gif"
    kargs = { 'duration': 0.1 }
    imageio.mimsave(exportname, image, 'GIF', **kargs)

    print("The .gif is finish :)")

