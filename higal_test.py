# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 15:00:35 2016

@author: Nathan
"""
from numpy import *
from astropy.io import fits
from matplotlib import pyplot as plt
plt.ion()
import aplpy
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord

from matplotlib.cbook import MatplotlibDeprecationWarning
import warnings
warnings.simplefilter('ignore', MatplotlibDeprecationWarning)
plt.close('all')

data = loadtxt("/Users/Nathan/Documents/UNI/Year 4 Project/Fits files/fits titles.txt",dtype='string',delimiter=',')


def smallplot(file,figure,subplot,Title):
    
    fig=aplpy.FITSFigure(file,hdu=1,figure=figure,subplot=subplot)
    
    fig.show_grayscale(stretch='log')
    fig.show_contour()
    fig.add_grid
    fig.set_title(Title)
    fig.tick_labels.set_font(size='x-small')
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.set_theme('publication')
    x=SkyCoord('16h35m08.48s','-48d46m32.2s',unit=(u.hourangle, u.deg) )    
    fig.recenter(x.ra.degree,x.dec.degree,0.06)
    return fig
'''
f,a=plot.subplots(3,3)
title=array(['extdPLW','extdPMW','extdPSW','hiresPLW','hiresPMW','hiresPSW'])
f.clf()
for i in range(len(data)):
    smallplot(data[i],f,(3,3,(i+1)),title[i])
 
plt.show()
'''

fig=aplpy.FITSFigure(data[0],hdu=1)   
fig.show_grayscale(stretch='log')
fig.show_contour()
fig.add_grid
fig.set_title('extdPLW')
fig.tick_labels.set_font(size='x-small')
fig.tick_labels.set_xformat('hh:mm:ss')
fig.set_theme('publication')


hdulist = fits.open(data[0])
scidata = hdulist[1].data

hdulist2 = fits.open(data[3])
scidata2 = hdulist[1].data

scidata = nan_to_num(scidata)
x=SkyCoord('16h35m14s','-48d47m29s',unit=(u.hourangle, u.deg) )    
#y=SkyCoord('16h34m58s','-48d45m29.2s',unit=(u.hourangle, u.deg) )  
w=wcs.WCS(data[0])

plt.axvline(x=(int(w.wcs_world2pix(x.ra.degree, x.dec.degree,0)[0])),ymin=0.25,ymax=0.75,c="blue",linewidth=0.5)

graph= scidata[int(w.wcs_world2pix(x.ra.degree, x.dec.degree,0)[0]),:]
graph2= scidata2[int(w.wcs_world2pix(x.ra.degree, x.dec.degree,0)[0]),:]

plt.figure()
plt.plot(graph)

plt.plot(graph2)

