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


G = 6.67430   #*10**(-11)       # m^3.kg^−1.s^−2
M = 10        #**(10) # Solar mass
a = 2.5
b = a/20
z = 0
epsilon = 0.00001

theta_zero = np.pi
V_R_zero = 0
V_theta_zero = 0.4
R_0 = 4

        # Runge kutta 4 :


def R_pt_pt (R, R_pt, theta, theta_pt, G, M, a, b, z=0) :
    S = np.sqrt( R**2 + ( a**2 + np.sqrt( b**2 + z**2 ) )**2 )
    grad_phi_r = G*M * (R/S**3)
    R_pp = - grad_phi_r + R*(theta_pt**2)
    return R_pp

def theta_pt_pt (R, R_pt, theta, theta_pt, G, M) :
    theta_pp = - 2*R_pt*theta_pt/R
    return theta_pp

t_max = 15
dt = 0.05
nb_star = 15
r_min = 1
r_max = 6

limit = r_max + 0.1*r_max


x = np.linspace(-limit, limit, 50)
y = np.linspace(-limit, limit, 50)
X_pot, Y_pot = np.meshgrid(x, y)

R = np.sqrt(X_pot**2 + Y_pot**2)
S = np.sqrt( R**2 + ( a**2 + np.sqrt( b**2 + z**2 ) )**2 )
Z_pot = - (G*M)/(S)


plt.figure()


plt.pcolormesh(X_pot, Y_pot, Z_pot, shading="gouraud")
plt.colorbar()

for i in range(0, nb_star-1) :
    
    R = random.uniform(r_min, r_max)
    theta = random.uniform(0, 2*np.pi)

    V_R = 0
    V_theta = 1/R
    t = 0

    
    list_X = []
    list_Y = []
    list_V_R = []
    list_V_theta = []
    
    while t < t_max+epsilon :
        
        X, Y = cyl2cart(R, theta)
        
        list_X += [X]
        list_Y += [Y]
        list_V_R += [V_R]
        list_V_theta += [V_theta]
        
        k_R1 = V_R*dt
        k_theta1 = V_theta*dt
        k_VR1 = R_pt_pt(R, V_R, theta, V_theta, G, M, a, b) * dt
        k_Vtheta1 = theta_pt_pt(R, V_R, theta, V_theta, G, M) * dt
        
        k_R2 = (V_R + (1/2)*k_R1) * dt
        k_theta2 = (V_theta + (1/2)*k_theta1) * dt
        k_VR2 = R_pt_pt(R+(1/2)*k_R1, V_R+(1/2)*k_VR1, theta+(1/2)*k_theta1, V_theta+(1/2)*k_Vtheta1, G, M, a, b) * dt
        k_Vtheta2 = theta_pt_pt(R+(1/2)*k_R1, V_R+(1/2)*k_VR1, theta+(1/2)*k_theta1, V_theta+(1/2)*k_Vtheta1, G, M) * dt
        
        k_R3 = (V_R + (1/2)*k_R2) * dt
        k_theta3 = (V_theta + (1/2)*k_theta2) * dt
        k_VR3 = R_pt_pt(R+(1/2)*k_R2, V_R+(1/2)*k_VR2, theta+(1/2)*k_theta2, V_theta+(1/2)*k_Vtheta2, G, M, a, b) * dt
        k_Vtheta3 = theta_pt_pt(R+(1/2)*k_R2, V_R+(1/2)*k_VR2, theta+(1/2)*k_theta2, V_theta+(1/2)*k_Vtheta2, G, M) * dt
        
        k_R4 = (V_R + k_R3) * dt
        k_theta4 = (V_theta + k_theta3) * dt
        k_VR4 = R_pt_pt(R+k_R3, V_R+k_VR3, theta+k_theta3, V_theta+k_Vtheta3, G, M, a, b) * dt
        k_Vtheta4 = theta_pt_pt(R+k_R3, V_R+k_VR3, theta+k_theta3, V_theta+k_Vtheta3, G, M) * dt
        
    
        R += 1/6 * (k_R1 + 2*k_R2 + 2*k_R3 + k_R4)
        theta += 1/6 * (k_theta1 + 2*k_theta2 + 2*k_theta3 + k_theta4)
        V_R += 1/6 * (k_VR1 + 2*k_VR2 + 2*k_VR3 + k_VR4)
        V_theta += 1/6 * (k_Vtheta1 + 2*k_Vtheta2 + 2*k_Vtheta3 + k_Vtheta4)
        
        
        
        t += dt

    plt.scatter(list_X, list_Y, s = 2)
    
    
plt.xlim(-limit,limit)
plt.ylim(-limit,limit)

plt.show()
#%%
















