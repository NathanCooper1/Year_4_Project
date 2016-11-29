# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 14:48:11 2016

@author: Nathan
"""

from numpy import *
from matplotlib import pyplot as plt
from scipy.stats import norm
import matplotlib.mlab as mlab
import math
import project_functions
plt.close('all')
mu1 = 1500
sigma1 = 300
x = linspace(0, 3000, 10000)
g,ax=plt.subplots(1,3,sharex=True,sharey=True)
ax[0].plot(x,mlab.normpdf(x, mu1, sigma1),label='Prior distribution')
ax[0].set_title('Beam profile')
plt.show()
ax[0].set_ylabel('$Arbitrary$ $units$')
#ax[0].set_xlabel('Position')
ax[1].set_xlabel('Position')
#ax[2].set_xlabel('Position')
ax[0].set_xticks([1000, 2000, 3000])
x1=1000
x2=2000
ax[1].axvline(x=x1,ymin=0.,ymax=1,c="blue",linewidth=0.5)  
ax[1].axvline(x=x2,ymin=0.,ymax=1,c="blue",linewidth=0.5)  
ax[1].set_xlim((0,3000))
ax[1].set_title('Target')
#plt.figure()
y1=mlab.normpdf(x,x1, sigma1)+mlab.normpdf(x,x2, sigma1)
ax[2].plot(x,y1)
plt.xlim((0,3000))
ax[2].set_title('Convolved target')
