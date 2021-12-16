# -*- coding: utf-8 -*-
'''
Project A: An optical ray tracer
Aim: To generate collimated rays that are uniformly distributed.
Written by Son-Gyo Jung; Group F; Lab Group A9; 04.02.16
'''

from math import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def rtpairs(R, N):
    '''
    A generator that yields a sequence of rings of r, θ coordinate pairs.
        
    Parameters:
        R - a list of radii
        N - a list containing the number of angles θ for each value of radius.
    '''
    for r in range(len(R)):
    	rad=0.0
    	for t in range(N[r]):
    		rad += 2*(np.pi)/N[r]
    		yield R[r], rad
			

def rtuniform(n, rmax, m): 
    '''
    Generates a sequence of r, θ that are approximately uniformly distributed over a disk by employing the rtpairs(R, N) function. 
    
    Parameters:
        n - Number of rings (this includes the central ray).
        rmax - Parameter which determines the maximum radius. (Note: the array created starts with zero radius by default.)
        m - The step-size of number of rays in each subsequent ring.
    '''
    R = np.arange (0.0, rmax, rmax/n)
    N = np.arange(0, n*m, m) 
    N[0] = 1 #set the first digit as 1 by default

    return rtpairs(R, N)
    

#plt.close()	
fig=plt.figure()
ax=fig.add_subplot((111))#plot in 2-dimensional space


for r, t in rtuniform(n=6, rmax=6., m=6):
    ax.scatter (r*cos(t), r*sin(t))
    
ax.set_title("Input in the x-y Plane")
ax.set_xlabel("y (mm)")
ax.set_ylabel("x (mm)")    
plt.grid()
plt.show ()
