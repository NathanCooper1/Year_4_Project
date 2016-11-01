# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 15:34:55 2016

@author: Nathan
"""
'''
import matplotlib as plt
from numpy import *
from sympy import integrate,symbols
import numpy as np
from matplotlib.pyplot import *
import pylab as py
from astropy.io import fits

import pyfits
imgname = "/Users/Nathan/Desktop/practiseNGC4736.fits"
img = pyfits.getdata(imgname)

pixel_average = img.mean()
pixel_low = img.min()
pixel_high= img.max()

py.clf
py.imshow(img)

'''
import numpy as np
import scipy.ndimage
import matplotlib.pyplot as plt

#-- Generate some data...
x, y = np.mgrid[-5:5:0.1, -5:5:0.1]
z = np.sqrt(x**2 + y**2) + np.sin(x**2 + y**2)

#-- Extract the line...
# Make a line with "num" points...
x0, y0 = 5, 4.5 # These are in _pixel_ coordinates!!
x1, y1 = 60, 75
num = 1000
x, y = np.linspace(x0, x1, num), np.linspace(y0, y1, num)

# Extract the values along the line, using cubic interpolation
zi = scipy.ndimage.map_coordinates(z, np.vstack((x,y)))

#-- Plot...
fig, axes = plt.subplots(nrows=2)
axes[0].imshow(z)
axes[0].plot([x0, x1], [y0, y1], 'ro-')
axes[0].axis('image')

axes[1].plot(zi)

plt.show()
