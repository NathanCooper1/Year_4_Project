# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 15:00:35 2016

@author: Nathan
"""
from numpy import *
from astropy.io import fits
from matplotlib import pyplot as plot
plot.ion()
import aplpy
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord

from matplotlib.cbook import MatplotlibDeprecationWarning
import warnings
warnings.simplefilter('ignore', MatplotlibDeprecationWarning)
plot.close('all')

#file='/Users/Nathan/Documents/UNI/Year 4 Project/Fits files/obsid_1342204057.refs["level2_5"].product.refs["extdPLW"].product.fits'

#data=fits.open(file)
data = loadtxt("/Users/Nathan/Documents/UNI/Year 4 Project/Fits files/fits titles.txt",dtype='string',delimiter=',')




#fig=aplpy.FITSFigure(file)
#
#fig.show_grayscale()
#fig.show_contour()
#fig.add_grid
#fig.add_colorbar
#x=SkyCoord('16h35m08.48s','-48d46m32.2s',unit=(u.hourangle, u.deg) )
#
#fig.recenter(x.ra.degree,x.dec.degree,0.06)


def smallplot(file,figure,subplot):
    
    fig=aplpy.FITSFigure(file,hdu=1,figure=figure,subplot=subplot)
    
    fig.show_grayscale(stretch='log')
    fig.show_contour()
    fig.add_grid
    fig.set_theme('publication')
    x=SkyCoord('16h35m08.48s','-48d46m32.2s',unit=(u.hourangle, u.deg) )    
    fig.recenter(x.ra.degree,x.dec.degree,0.06)
    return fig

f,a=plot.subplots(3,3)
f.clf()
for i in range(len(data)):
    smallplot(data[i],f,(3,3,(i+1)))
 
   
#f,a=plot.subplots(3,3)
#f=plot.figure(1)
#f.clf()
#
#fig1=smallplot(file,f,(3,3,1))
#smallplot(file,f,(3,3,3))
#smallplot(file,f,(3,3,5))
#
plot.show()
##a[0,0].smallplot('/Users/Nathan/Documents/UNI/Year 4 Project/Fits files/obsid_1342204057.refs["level2_5"].product.refs["hiresPLW"].product.fits')
