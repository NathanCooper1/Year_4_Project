# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 15:34:55 2016

@author: Nathan
"""

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

