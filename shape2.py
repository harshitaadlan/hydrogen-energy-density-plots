#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 18:28:54 2024

@author: harshita
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 00:50:40 2020

@author: harshita
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import  numpy as np
z=1
a=5.29*10**(-11)
def dxy(r,theta,phi):
    rp =2/((81*(3.14)**(0.5))*(z/a)**(3/2))
    e = (z*r/a)**2*np.exp((-1)*z*r/a*3)
    ang =np.sin(theta)*np.cos(theta)*np.cos(phi)
    psi = rp*e*ang
    return psi
def HFunc(r,theta,phi,n,l,m):
    '''
    Hydrogen wavefunction // a_0 = 1

    INPUT
        r: Radial coordinate
        theta: Polar coordinate
        phi: Azimuthal coordinate
        n: Principle quantum number
        l: Angular momentum quantum number
        m: Magnetic quantum number

    OUTPUT
        Value of wavefunction
    '''

    coeff = np.sqrt((2.0/n)**3 * spe.factorial(n-l-1) /(2.0*n*spe.factorial(n+l)))
    laguerre = spe.assoc_laguerre(2.0*r/n,n-l-1,2*l+1)
    sphHarm = spe.sph_harm(m,l,phi,theta) # Note the different convention from doc

    return coeff * np.exp(-r/n) * (2.0*r/n)**l * laguerre * sphHarm


rl=[]
r_x=[]
r_y=[]
r_z=[] 
n=3
l=2
m=1
plot_range = int((4*n+4*l)*a)
fig = plt.figure()
ax = Axes3D(fig)
outMat = np.zeros((nl,nw,nh))
for x in range (-plot_range,plot_range):
    for y in range (-plot_range,plot_range):
        for z in range (-plot_range,plot_range):
            
            r = np.sqrt(x**2+y**2+z**2)
            rl.append(r)
            if (r == 0):
                    phi = 0
            else:
                    phi = np.arccos(z/r)  #this  is theta  

            if (x == 0):
                    theta = np.pi/2
            else:
                    theta = np.arctan(y/x)
            funcEval = HFunc(r,theta,phi,n,l,m)

            outMat[x,y,z] = np.real(funcEval*np.conj(funcEval)
print(funcEval)
#            if psi2==1:
#                x=r*np.sin(phi)*np.cos(theta)
#               y=r*np.sin(phi)*np.sin(theta)
#               z=r*np.cos(theta)
#               print(x,y,z)
#                plt.scatter(x, y, z)

plt.show()
              
#theta = np.r_[0:np.pi:50j]                  
#phi = np.r_[0:2*np.pi:50j]
#radius = np.r_[0:n:50j]
#x = np.array([dxy(rl,theta, phi) for r in radius])
#y = np.array([dxy(rl,theta, phi) for r in radius])
#z = np.array([dxy(rl,theta, phi) for r in radius])




#fig = plt.figure()
#ax = Axes3D(fig)