#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 17:22:05 2022

@author: renaudvanco
"""


import Definition as Df
import numpy as np


        # Potential :
    
G = 4*10**(-12) # kPc^3.M_s^-1.By^-2,    6.67430*10**(-11) m^3.kg^−1.s^−2
M_pot = 4*10**(9) # Solar mass
a = 2.5   # kPc
b = a/20
# Perturber :
M_perturb = M_pot/100 #10**(9) # Solar mass
# Star :
M_star = 1 # Solar mass

t_max = 3000
dt = 1
epsilon = 0.00001
time = []

r_bulge = 3   # kPc
r_disk = 35   # kPc 

nb_star_bulge = 25
nb_star_disk = 2475
nb_star_gal = nb_star_bulge + nb_star_disk
nb_star_perturber = 0
nb_star = nb_star_gal + nb_star_perturber


# Setting up list :

list_X, list_Y, list_X_P, list_Y_P, list_X_SP, list_Y_SP, list_R, list_theta, list_Z, list_V_R, list_V_theta, list_V_Z, list_R_P, list_theta_P, list_Z_P, list_V_R_P, list_V_theta_P, list_V_Z_P, list_R_SP, list_theta_SP, list_Z_SP, list_V_R_SP, list_V_theta_SP, list_V_Z_SP, list_R_BH, list_R_BH_SP, list_ang_mom, list_ang_mom_P =  Df.Setting_list(nb_star_gal, nb_star_perturber)   


# Initial position
print("Perturbator generation :")
R_P, theta_P, Z_P, V_R_P, V_theta_P, V_Z_P = Df.position_perturber(r_disk, G, M_pot, a, b)

list_R_P.append(R_P)
list_theta_P.append(theta_P)
list_Z_P.append(Z_P)
list_V_R_P.append(V_R_P)
list_V_theta_P.append(V_theta_P)
list_V_Z_P.append(V_Z_P)
list_ang_mom_P.append(Df.ang_mom(R_P, V_theta_P))

print("Bulge generation : ")
for i in range(0, nb_star_bulge) :
    R, theta, Z, V_R, V_theta, V_Z = Df.position_stars_bulge(r_bulge, i, G, M_pot, a, b)
    R_BH = Df.distance_S_P(R, R_P, theta, theta_P, Z, Z_P)
    
    list_R[i].append(R)
    list_theta[i].append(theta)
    list_Z[i].append(Z)
    list_V_R[i].append(V_R)
    list_V_theta[i].append(V_theta)
    list_V_Z[i].append(V_Z)
    list_R_BH[i].append(R_BH)
#    list_ang_mom[i].append(Df.ang_mom(R, V_theta))

print("Disk generation : ")
for i in range(nb_star_bulge, nb_star_gal) :
    R, theta, Z, V_R, V_theta, V_Z = Df.position_stars_disk(r_bulge, r_disk, i, G, M_pot, a, b)
    R_BH = Df.distance_S_P(R, R_P, theta, theta_P, Z, Z_P)
    
    list_R[i].append(R)
    list_theta[i].append(theta)
    list_Z[i].append(Z)
    list_V_R[i].append(V_R)
    list_V_theta[i].append(V_theta)
    list_V_Z[i].append(V_Z)
    list_R_BH[i].append(R_BH)
#    list_ang_mom[i].append(Df.ang_mom(R, V_theta))

print("Perturber generation : ")
for i in range(0, nb_star_perturber) :
    R, theta, Z, V_R, V_theta, V_Z = Df.position_stars_perturber(R_P, theta_P, Z_P, V_R_P, V_theta_P, V_Z_P, i, G, M_perturb, a, b)
    R_BH = Df.distance_S_P(R, R_P, theta, theta_P, Z, Z_P)
    
    list_R_SP[i].append(R)
    list_theta_SP[i].append(theta)
    list_Z_SP[i].append(Z)
    list_V_R_SP[i].append(V_R)
    list_V_theta_SP[i].append(V_theta)
    list_V_Z_SP[i].append(V_Z)
    list_R_BH_SP[i].append(R_BH)
#    list_ang_mom[i].append(Df.ang_mom(R, V_theta))


			# Integration :

t = 0
j=0
time.append(t)

while t < t_max+epsilon :
    if round(j) % 500 == 0 : 
        print("integration for t = ", np.round(t, 2), "over t_max = ", t_max)
    
    # Motion of stars in the galaxies :
    for i in range(0, nb_star_gal) :
        R_S, theta_S, Z_S, V_R, V_theta, V_Z = list_R[i][-1], list_theta[i][-1], list_Z[i][-1], list_V_R[i][-1], list_V_theta[i][-1], list_V_Z[i][-1]
        R_BH_R = Df.distance_S_P_R(R_S, R_P, theta_S, theta_P)
        R_BH = Df.distance_S_P(R_S, R_P, theta_S, theta_P, Z_S, Z_P)

        R_S, theta_S, Z_S, V_R, V_theta, V_Z = Df.RK4_star(G, M_pot, M_star, M_perturb, a, b, dt, R_S, theta_S, Z_S, V_R, V_theta, V_Z, R_P, R_BH, R_BH_R, theta_P, Z_P)
        
        list_R[i].append(R_S)
        list_theta[i].append(theta_S)
        list_Z[i].append(Z_S)
        list_V_R[i].append(V_R)
        list_V_theta[i].append(V_theta)
        list_V_Z[i].append(V_Z)
        list_R_BH[i].append(R_BH)
        list_ang_mom[i].append(Df.ang_mom(R_S, V_theta))

    # Motion of stars starting with the perturber :
    for i in range(0, nb_star_perturber) :
        R_S, theta_S, Z_S, V_R, V_theta, V_Z = list_R_SP[i][-1], list_theta_SP[i][-1], list_Z_SP[i][-1], list_V_R_SP[i][-1], list_V_theta_SP[i][-1], list_V_Z_SP[i][-1]
        R_BH_R = Df.distance_S_P_R(R_S, R_P, theta_S, theta_P)
        R_BH = Df.distance_S_P(R_S, R_P, theta_S, theta_P, Z_S, Z_P)

        R_S, theta_S, Z_S, V_R, V_theta, V_Z = Df.RK4_star(G, M_pot, M_star, M_perturb, a, b, dt, R_S, theta_S, Z_S, V_R, V_theta, V_Z, R_P, R_BH, R_BH_R, theta_P, Z_P)
        
        list_R_SP[i].append(R_S)
        list_theta_SP[i].append(theta_S)
        list_Z_SP[i].append(Z_S)
        list_V_R_SP[i].append(V_R)
        list_V_theta_SP[i].append(V_theta)
        list_V_Z_SP[i].append(V_Z)
        list_R_BH_SP[i].append(R_BH)


    # Perturber motion :
    R_P, theta_P, Z_P, V_R_P, V_theta_P, V_Z_P = Df.RK4_perturber(G, M_pot, M_perturb, a, b, dt, R_P, theta_P, Z_P, V_R_P, V_theta_P, V_Z_P)

    list_R_P.append(R_P)
    list_theta_P.append(theta_P)
    list_Z_P.append(Z_P)
    
    
    list_V_R_P.append(V_R_P)
    list_V_theta_P.append(V_theta_P)
    list_V_Z_P.append(V_Z_P)
    #list_ang_mom_P.append(Df.ang_mom(R_P, V_theta_P))

    j += 1
    t += dt
    time.append(t)

    
for i in range(0, nb_star_gal):
    for j in range(0, len(list_R[i])):
        X, Y = Df.cyl2cart(list_R[i][j], list_theta[i][j])

        list_X[i].append(X)
        list_Y[i].append(Y)

for i in range(0, nb_star_perturber):
    for j in range(0, len(list_R_SP[i])):
        X, Y = Df.cyl2cart(list_R_SP[i][j], list_theta_SP[i][j])

        list_X_SP[i].append(X)
        list_Y_SP[i].append(Y)

for j in range(0, len(list_R_P)):
        X, Y = Df.cyl2cart(list_R_P[j], list_theta_P[j])

        list_X_P.append(X)
        list_Y_P.append(Y)



###########################################################################################################################
          # Visualisation 2D :
###########################################################################################################################
x_lim, y_lim, z_lim = 50, 50, 30

#Df.graph_ang_mom(nb_star, time, list_ang_mom, list_ang_mom_P)
#Df.XY_traj(nb_star, list_X, list_X_P, list_Y, list_Y_P, list_Z, G, M_pot, a, b, x_lim, y_lim)
#Df.XZ_traj(nb_star, list_X, list_Z, list_X_P, list_Z_P, G, M_pot, x_lim, z_lim, a, b)
#Df.radii_graph(list_R, time)

Df.XY_GIF(nb_star, list_X, list_X_P, list_X_SP, list_Y, list_Y_P, list_Y_SP, list_Z, list_Z_P, list_Z_SP, list_V_R, list_V_theta, list_V_Z, list_R_BH, G, M_pot, M_perturb, a, b, x_lim, y_lim, dt, t_max, time)
Df.XZ_GIF(nb_star, list_X, list_X_P, list_X_SP, list_Y, list_Y_P, list_Y_SP, list_Z, list_Z_P, list_Z_SP, list_V_R, list_V_theta, list_V_Z, list_R_BH, G, M_pot, M_perturb, a, b, x_lim, z_lim, dt, t_max, time)

###########################################################################################################################
            # Visualisation 3D :
###########################################################################################################################

#Df.XYZ_GIF(list_X, list_Y, list_Z, list_X_P, list_Y_P, list_Z_P, x_lim, y_lim, z_lim, nb_star, dt, t_max, time)


















