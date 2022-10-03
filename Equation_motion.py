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


#%%
    
###########################################################################################################################
            # Newton free fall in cartesian :
###########################################################################################################################

    # Constants :

V_0 = 10
X_0 = -5
Y_0 = 0
theta_deg = 60
theta_rad = theta_deg/360 * 2*np.pi
g = 9.81
t = 0
dt = 0.005 # (sec)
t_max = 2 # sec


    # Function :

def newton (t, g=g, theta_0=theta_rad, V_0=V_0, X_0=X_0, Y_0=Y_0) :
    """
    Implement the gravitationnal force on an object with no friction :

    Parameters
    ----------
    t : time, in second
    g : terrestrial gravitationnal constante
    theta_0 : angle at wich the point is thrown
    V_0 : Initial speed of the point
    X_0 : initial X position
    Y_0 : initial Y position

    Returns
    -------
    a_x : acceleration in x
    a_y : acceleration in y
    V_x : speed in x
    V_y : speed in y
    x : X position
    y : Y position
    """
    
    a_y = -g
    a_x = 0
    V_y = -g*t + V_0*np.sin(theta_0)*t
    V_x = V_0*np.cos(theta_0)
    y = -(1/2)*g*t**2 + V_0*np.sin(theta_0)*t + Y_0
    x = V_0*np.cos(theta_0)*t + X_0
    return a_x, a_y, V_x, V_y, x, y

    

list_a_x = []
list_a_y = []
list_V_x = []
list_V_y = []
list_x = []
list_y = []
T = []

    # Motion equation loop :

y = 1

while t < t_max+0.0001 and y > -0.1 :
    a_x, a_y, V_x, V_y, x, y = newton(t, g, theta_rad, V_0, X_0, Y_0)
    list_a_x += [a_x]
    list_a_y += [a_y]
    list_V_x += [V_x]
    list_V_y += [V_y]
    list_x += [x]
    list_y += [y]
    T += [t]
    t += dt


    # Graph :

fig, ax = plt.subplots()
plt.xlim(X_0-0.5, np.max(list_x)+1)
plt.ylim(Y_0-0.5, np.max(list_y)+1)
ax.grid()


plt.plot(list_x, list_y, label = f"t = {np.round(T[len(T)-1], 3)}sec")

plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.title(f"NEWTON : \n Free fall of a mass in a gravity field of g={g} \n V_0 = {V_0}m/s, theta = {theta_deg}°")
plt.legend()

plt.show()



#%%

###########################################################################################################################
            # Switch between cartesian and polar coordinate :
###########################################################################################################################


def cart2cyl(x, y):
    """
    Convert from cartesian coordinate to cylindrical coordinate
    """
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x) # for negative result of phi just add 2*np.pi
    return rho, phi

def cyl2cart(rho, phi):
    """
    Convert from Cylindrical coordinate to Cartesian coordinate
    """
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return x, y

list_x = []
list_y = []
rho = 1
for i in range(0, 360, 10) :
    phi = i/360 * 2*np.pi
    x, y = cyl2cart(rho, phi)
    list_x += [x]
    list_y+= [y]



plt.scatter(list_x,list_y)
plt.title("pol to cart")

list_rho = []
list_phi_rad = []
list_phi_deg = []

for i in range(0, len(list_x)) :
    rho, phi = cart2cyl(list_x[i], list_y[i])
    list_rho += [round(rho, 2)]
    list_phi_rad += [phi]

list_phi_rad = np.array(list_phi_rad)
list_phi_rad[list_phi_rad<0] += 2*np.pi
list_phi_deg = np.copy(list_phi_rad)/(2*np.pi) *360


for i in range(0, len(list_rho)) :
    print(f"rho = {list_rho[i]},        phi = {round(list_phi_rad[i], 3)}rad = {round(list_phi_deg[i], 3)}deg ")
    
    



#%%
###########################################################################################################################
            # Miyamoto-Nagai potential :
###########################################################################################################################


        # Constants :


G = 6.67430   #*10**(-11)       # m^3.kg^−1.s^−2
M = 10        #**(10) # Solar mass
a = 2.5
b = a/20
R_0 = 1
Z_0 = 0
V_0r = 0
V_0phi = 10
V_0z = 0
phi_0 = 0

        # Functions :

def motion_star (G, M, R, S, Z, V_0r=V_0r, V_0phi=V_0phi, V_0z=V_0z, R_0=R_0, phi_0=phi_0, Z_0=Z_0) :
    """
    Motion of an object in a Miyamoto-Nagai potential :
        Problème car j'ai integré ça comme en cartesien
    """
    a_r = - G*M * (R/S**3)
    a_phi = 0
    a_z = - G*M * ( Z * (a + np.sqrt(Z**2+b**2)) / (S**3 * np.sqrt(Z**2+b**2)) )
    V_r = a_r*t + V_0r
    V_phi = a_phi + V_0phi
    V_z = a_z*t + V_0z
    e_r = (1/2)*a_r*t**2 + V_0r*t + R_0
    e_phi = (1/2)*a_phi*t**2 + V_0phi*t + phi_0
    e_z = (1/2)*a_z*t**2 + V_0z*t + Z_0
    
    return a_r, a_phi, a_z, V_r, V_phi, V_z, e_r, e_phi, e_z

def S_compute (R, Z, a=a, b=b) :
    
    S = np.sqrt( R**2 + ( a**2 + np.sqrt( b**2 + Z**2 ) )**2 )
    
    return S

        # Motion equation loop :

list_r = []
list_phi = []
list_z = []
T = []

t = 0
dt = 0.0001 # (sec)
t_max = 0.02 # sec

R = R_0
Z = Z_0

while t < t_max+0.0001 :
    S = S_compute(R, Z)
    a_r, a_phi, a_z, V_r, V_phi, V_z, R, phi, Z = motion_star (G, M, S, R, Z)
    
    list_r += [R]
    list_phi += [phi]
    list_z += [Z]
    T += [t]
    t += dt


        # To cartesian :

list_x = []
list_y = []

for i in range(0, len(list_r)) :
    x, y = cyl2cart(list_r[i], list_phi[i])
    list_x += [x]
    list_y+= [y]



        # Graph

plt.figure()
plt.xlim(-1, 1)
plt.ylim(-1, 1)
plt.plot(list_x,list_y)
plt.show()














