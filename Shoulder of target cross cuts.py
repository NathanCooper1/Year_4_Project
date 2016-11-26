# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 12:42:20 2016

@author: Nathan
"""
###looking at the shoulder of the target

from numpy import *
from astropy.io import fits
from matplotlib import pyplot as plt
plt.ion()
import aplpy
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
#from project_functions.py import *

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
    
def verticalcut(file1,file2,x,figure,subplot):
    hdulist = fits.open(file1)
    scidata = hdulist[1].data
    hdulist2 = fits.open(file2)
    scidata2 = hdulist2[1].data
    coord=array([[x.ra.degree,x.dec.degree]])
    w=wcs.WCS(hdulist[1].header)
    w2=wcs.WCS(hdulist2[1].header)
    graph= scidata[150:250,int(w.wcs_world2pix(coord,1)[0,0])]
    graph2= scidata2[300:500,int(w2.wcs_world2pix(coord,1)[0,0])]
    graphx=linspace(1,100,100)*float(hdulist[1].header['CDELT1'])
    graph2x=linspace(1,200,200)*float(hdulist2[1].header['CDELT1'])
    g=figure.add_subplot(subplot)
    g.plot(graphx,graph,label='Native')
    g.plot(graph2x,graph2,label='Hi res')
    g.legend('loc=best')
    
def horizontalcut(file1,file2,x,figure,subplot):
    hdulist = fits.open(file1)
    scidata = hdulist[1].data
    hdulist2 = fits.open(file2)
    scidata2 = hdulist2[1].data
    coord=array([[x.ra.degree,x.dec.degree]])
    w=wcs.WCS(hdulist[1].header)
    w2=wcs.WCS(hdulist2[1].header)
    graph= scidata[int(w.wcs_world2pix(coord,1)[0,1]),:]
    graph2= scidata2[int(w2.wcs_world2pix(coord,1)[0,1]),:]
    graphx=linspace(1,100,100)*float(hdulist[1].header['CDELT1'])
    graph2x=linspace(1,200,200)*float(hdulist2[1].header['CDELT1'])
    g=figure.add_subplot(subplot)
    g.plot(graphx,graph,label='Native')
    g.plot(graph2x,graph2,label='Hi res')
    g.legend('loc=best')
    
x=SkyCoord('16h35m13.7s','-48d45m42s',unit=(u.hourangle, u.deg) )       

f=plt.figure()
smallplot(data[5],f,(1,1,1),('Target'))  

hdulist = fits.open(data[5])
w=wcs.WCS(hdulist[1].header)
coord=array([[x.ra.degree,x.dec.degree]])

plt.axvline(x=(int(w.wcs_world2pix(coord,1)[0,0])),ymin=0.,ymax=1,c="blue",linewidth=0.5)    
plt.axhline(y=(int(w.wcs_world2pix(coord,1)[0,1])),xmin=0.,xmax=1,c="red",linewidth=0.5)    
'''
g,ax=plt.subplots(1,3)
smallplot(data[5],g,(1,3,1),(''))
ax[0].axvline(x=(int(w.wcs_world2pix(coord,1)[0,0])),ymin=0.,ymax=1,c="blue",linewidth=0.5)
ax[0].axhline(y=(int(w.wcs_world2pix(coord,1)[0,1])),xmin=0.,xmax=1,c="red",linewidth=0.5)  
verticalcut(data[4],data[5],x,g,(1,3,2))
horizontalcut(data[4],data[5],x,g,(1,3,3))
'''