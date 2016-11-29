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
#from project_functions import *

from matplotlib.cbook import MatplotlibDeprecationWarning
import warnings
warnings.simplefilter('ignore', MatplotlibDeprecationWarning)
plt.close('all')

data = loadtxt("/Users/Nathan/Documents/UNI/Year 4 Project/Fits files/fits titles.txt",dtype='string',delimiter=',')

def fitsplot(file,figure,x,Title):    
    fig=aplpy.FITSFigure(file,hdu=1,figure=figure)
    fig.show_grayscale(stretch='log')
    fig.show_contour()
    fig.add_grid
    fig.set_title(Title)
    fig.tick_labels.set_font(size='x-small')
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.set_theme('publication')   
    fig.recenter(x.ra.degree,x.dec.degree,0.06)
    linearr=array([[x.ra.degree,x.ra.degree],[-48.83,-48.73]])
    fig.show_lines([linearr],linewidth=2,color='c')
    return fig

def smallplot(file,figure,x,subplot,Title):    
    fig=aplpy.FITSFigure(file,hdu=1,figure=figure,subplot=subplot)
    fig.show_grayscale(stretch='log')
    fig.show_contour()
    fig.add_grid
    fig.set_title(Title)
    fig.tick_labels.set_font(size='x-small')
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.set_theme('publication')   
    fig.recenter(x.ra.degree,x.dec.degree,0.06)
    linearr=array([[x.ra.degree,x.ra.degree],[-48.83,-48.73]])
    fig.show_lines([linearr],linewidth=2,color='c')
    return fig

def verticalcut(file1,file2,x,figure):
    hdulist = fits.open(file1)
    scidata = hdulist[1].data
    hdulist2 = fits.open(file2)
    scidata2 = hdulist2[1].data
    coord=array([[x.ra.degree,x.dec.degree]])
    w=wcs.WCS(hdulist[1].header)
    w2=wcs.WCS(hdulist2[1].header)
    graph= scidata[:,int(w.wcs_world2pix(coord,1)[0,0])]
    graph2= scidata2[:,int(w2.wcs_world2pix(coord,1)[0,0])]
    graphx=(arange(len(graph))*float(hdulist[1].header['CDELT1']))*60
    graph2x=(arange(len(graph2))*float(hdulist2[1].header['CDELT1']))*60
    plt.plot(graphx,graph,label='Native')
    plt.plot(graph2x,graph2,label='Hi res')
    plt.ylabel('$MJy/Sr$')
    plt.xlabel('Position/$Arcminutes$')
    plt.xlim((-49,-43))
    
    
def verticalcutsmall(file1,file2,x,figure,subplot):
    hdulist = fits.open(file1)
    scidata = hdulist[1].data
    hdulist2 = fits.open(file2)
    scidata2 = hdulist2[1].data
    coord=array([[x.ra.degree,x.dec.degree]])
    w=wcs.WCS(hdulist[1].header)
    w2=wcs.WCS(hdulist2[1].header)
    graph= scidata[:,int(w.wcs_world2pix(coord,1)[0,0])]
    graph2= scidata2[:,int(w2.wcs_world2pix(coord,1)[0,0])]
    graphx=(arange(len(graph))*float(hdulist[1].header['CDELT1']))*60
    graph2x=(arange(len(graph2))*float(hdulist2[1].header['CDELT1']))*60
    g=figure.add_subplot(3,3,subplot)
    g.plot(graphx,graph,label='Native')
    g.plot(graph2x,graph2,label='Hi res')
    g.set_ylabel('$MJy/Sr$')
    g.set_xlabel('Position/$Arcminutes$')
    g.set_xlim((-49,-43))
    


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
    g=figure.add_subplot(3,3,subplot)
    g.plot(graphx,graph,label='Native')
    g.plot(graph2x,graph2,label='Hi res')
    g.set_xlabel('Position')
    g.set_ylabel('$MJy/Sr$')
    


###Vertical cut x postion
x=SkyCoord('16h35m08s','-48d46m29s',unit=(u.hourangle, u.deg) )   

f,a=plt.subplots(3,3)
title=array(['extdPLW','extdPMW','extdPSW','hiresPLW','hiresPMW','hiresPSW'])
f.clf()

'''
for i in range(len(data)):
    #if i>0:
        #continue
    smallplot(data[i],f,x,(3,3,(i+1)),title[i])

for i in range(3):
    verticalcutsmall(data[i],data[i+3],x,f,(i+7))
    plt.legend()
f.tight_layout(pad=0)
'''
'''
#large scale image with line showing cut position
fig=aplpy.FITSFigure(data[4],hdu=1)   
fig.show_grayscale(stretch='log')
fig.show_contour()
fig.add_grid
fig.set_title('Target location')
fig.tick_labels.set_font(size='x-small')
fig.tick_labels.set_xformat('hh:mm:ss')
fig.set_theme('publication')
hdulist = fits.open(data[4])
scidata = hdulist[1].data
w=wcs.WCS(hdulist[1].header)
coord=array([[x.ra.degree,x.dec.degree]])
plt.axvline(x=(int(w.wcs_world2pix(coord,1)[0,0])),ymin=0.,ymax=1,c="blue",linewidth=0.5)
plt.axhline(y=(int(w.wcs_world2pix(coord,1)[0,1])),xmin=0.,xmax=1,c="blue",linewidth=0.5)
'''

h=plt.figure()
fitsplot(data[0],h,x,'Native')
h=plt.figure()
fitsplot(data[3],h,x,'Hi-res')
h=plt.figure()
verticalcut(data[0],data[3],x,h)
plt.title('Vertical cross cut')
plt.legend()


#linearr=array([[x.ra.degree,x.ra.degree],[-48.83,-48.73]])
#fig.show_lines([linearr],linewidth=2,color='c')



#

#
#
#hdulist2 = fits.open(data[3])
#scidata2 = hdulist2[1].data



##y=SkyCoord('16h34m58s','-48d45m29.2s',unit=(u.hourangle, u.deg) )  
#w=wcs.WCS(hdulist[1].header)
#whi=wcs.WCS(hdulist2[1].header)
#coord=array([[x.ra.degree,x.dec.degree]])
#pltx=w.wcs_world2pix(coord,1)
#
#plt.axvline(x=(int(w.wcs_world2pix(coord,1)[0,0])),ymin=0.,ymax=1,c="blue",linewidth=0.5)
#plt.plot(pltx[0,1],pltx[0,0],'mx',ms=20)
#graph= scidata[150:250,int(w.wcs_world2pix(coord,1)[0,0])]
#graph2= scidata2[300:500,int(whi.wcs_world2pix(coord,1)[0,0])]
#graphx=linspace(1,100,100)*float(hdulist[1].header['CDELT1'])
#graph2x=linspace(1,200,200)*float(hdulist2[1].header['CDELT1'])

##fig2,axes2=plt.subplots(3,3)
#g=f.add_subplot(3,3,7)
##g.title('Cross cut of target')
#g.plot(graphx,graph)
#g.plot(graph2x,graph2)
#
#h=f.add_subplot(3,3,8)
#
#
#plt.show()
