#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 14:32:48 2022

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
            # Constants :
###########################################################################################################################

G = 6.67430   #*10**(-11)       # m^3.kg^−1.s^−2
M = 10        #**(10) # Solar mass
a_test = 16
b_test = 16
a = 2.5     # kpc
b = a/20    # kpc
phi_min = 0
phi_max = 2*np.pi
R_C = 10

#%%
###########################################################################################################################
            # Functions :
###########################################################################################################################

def S_compute (R, a=a_test, b=b_test, z=0) :
    
    S = np.sqrt( R**2 + ( a**2 + np.sqrt( b**2 + z**2 ) )**2 )
    
    return S

def Miyamoto_Nagai (G, M, S) :
    Phi = - (G*M/S)
    return Phi

def Grad_Potential_r (G, M, R, S) :
    
    r_potential = G*M * (R/(S**3))
    
    return r_potential

def Grad_Potential_z (G, M, R, S, a, b, z=0) :
    
    r_potential = G*M * ( z*(a+np.sqrt(z**2+b**2)) / (S**3*np.sqrt(z**2+b**2)) )
    
    return r_potential


def angular_speed (G, M, S) :
    
    Omega = np.sqrt(G*M) * S**(-3/2)
    
    return Omega

def Lz_compute (R, Omega) :
    
    Lz = R**2 * Omega
    
    return Lz

def V_phi_compute (K, Omega, R, R_c) :
    
    V_phi = (K**2)/(2*Omega) * (R-R_c)
    
    return V_phi

def Schwarzschild_Velocity_Distrib (V_R, V_phi, sigma_R, Omega, K, sigma_z=0, V_z=0, C=1) :
    
    gamma = 2*Omega/K
    SVD = C * np.exp(-(V_R**2 + gamma**2 * V_phi**2)/(2*sigma_R**2) - (V_z**2)/(2*sigma_z**2))
    
    return SVD 



#%%


def Frequency (G, M, R, S) :
    
    K = G*M * (R * d(S**(-3))/dR + 4*S**(-3))
    
    return K

def Star_position (r_min=0, r_max=a_test, phi_min=phi_min, phi_max=phi_max) :
    
    Rc = random.uniform(r_min, r_max)
    phi = random.uniform(phi_min, phi_max)
    
    return Rc, phi

#%%
###########################################################################################################################
            # Initial Condition :
###########################################################################################################################

#       1 - Define a pair (Rc , Lz) for a circular orbit ;
#       2 - Put the star at coordinates (Rc , φ, 0) where φ : [ 0, 2π ] is drawn from a uniform distribution
Rc, phi = Star_position(0, a_test, phi_min, phi_max)
print(f"Rc = {Rc}, phi = {phi}")

S = S_compute(Rc)
print(f"S = {S}")

Lz = Lz_compute(G, M, S)
print(f"Lz = {Lz}")


#%%
###########################################################################################################################
            # Potential visualisation :
###########################################################################################################################

x = np.arange(-100, 100, 1)
y = np.arange(-100, 100, 1)
list_R = np.zeros((len(x),len(y)))


#%%
for i in x :
    for j in y :
        R = np.sqrt(i**2 + j**2)
        list_R[i,j] = R


















