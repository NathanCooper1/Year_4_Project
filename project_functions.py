# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 12:43:57 2016

@author: Nathan
"""

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
    
def verticalcut(file1,file2,x,subplot):
    hdulist = fits.open(file1)
    scidata = hdulist[1].data
    hdulist2 = fits.open(file2)
    scidata2 = hdulist2[1].data
    coord=array([[x.ra.degree,x.dec.degree]])
    w=wcs.WCS(hdulist[1].header)
    w2=wcs.WCS(hdulist2[1].header)
    graph= scidata[150:250,int(w.wcs_world2pix(coord,1)[0,0])]
    graph2= scidata2[150:250,int(w2.wcs_world2pix(coord,1)[0,0])]
    g=f.add_subplot(3,3,subplot)
    g.plot(graph,label='native')
    g.plot(graph2,ls='--',label='Hi res')
