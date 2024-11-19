#This is the BOBcat skymap plotting function. This is the function that actually plots the list of queried candidates on top of the smoothed strain upper limit data.

#importing tools
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import colors
from pylab import cm
from matplotlib import patheffects
from matplotlib import text
%matplotlib inline
%config InlineBackend.figure_format = 'retina'

import astropy.units as u
import astropy.constants as const
from astropy.coordinates import SkyCoord

import healpy as hp

#creating function
def skymap_plot(skymap_smoothed, data, nside):
    """
    This is the BOBcat skymap plotting function. It will take inputs of PTA UL data, a dataframe of 
    candidates with coordinate data, and a pixel count parameter (nside) to generate complete skymaps. 
    """
    
    #SETUP SECTION
    plt.rc('xtick',**{'labelsize':16})
    plt.rc('ytick',**{'labelsize':16})
    plt.rc('axes',**{'labelsize':20,'titlesize':18})
    
    #mollview tick labels
    #tick_labels=['22h', '20h', '18h', '16h', '14h', '12h', '10h', '8h', '6h', '4h', '2h']
    #y_ticks = ['$-75^\circ$', '$-60^\circ$', '$-45^\circ$','$-30^\circ$','$-15^\circ$','$0^\circ$',
    #       '$15^\circ$','$30^\circ$','$45^\circ$','$60^\circ$','$75^\circ$',]
    
    #colorbar tick labels
    cbar_ticks = [2e-15, 3e-15, 4e-15, 5e-15, 6e-15, 7e-15, 8e-15, 9e-15, 1e-14]
    cbar_labels = ['$2 \times 10^{-15}$', '$3 \times 10^{-15}$', '$4 \times 10^{-15}$', ' ', '$6 \times 10^{-15}$', ' ',
               ' ' , ' ', '$10^{-14}$']
    
    
    #assigning number of pixels to plot
    npix = hp.nside2npix(nside)
    
    #creating mollview plot
    mv = hp.mollview(skymap_smoothed, title="", cbar = None , cmap = 'viridis_r', rot = 180, norm='log',
#                        min = min(skymap11[0:768,3]), max = max(skymap11[0:768,3]),   
#                        bgcolor='white', badcolor = 'white',
                    )
    #drawing graticule (lines of latitude and longitude)
    hp.graticule(15, 30)

    #plotting candidates as points on the map
    hp.visufunc.projscatter(data['Right Ascension (Decimal Degrees)'], data['Declination (Decimal Degrees)'], #data[2], data[3], 
                        lonlat=True, marker='*', color='w',edgecolors='r',s=200, zorder = 2
                           )
        
    #Making colorbar for background map
    fig = plt.gcf()
    ax = plt.gca()
    image = ax.get_images()[0]
    cmap = fig.colorbar(image, ax=ax, orientation = 'horizontal', label = 'GW Strain Upper Limit', 
                        ticks = cbar_ticks, pad = 0.05
                       )
        
    #plotting RA labels
    for i in range(2,24,2):
        text = hp.projtext( i*180/12+5, 2,  str(i)+'h', lonlat=True, coord='G', fontsize = 'large', fontweight = 100, zorder = 1)
        text.set_path_effects([patheffects.withStroke(linewidth=4, foreground='w')])
    #plotting dec labels
    for i in range(-75,0,15):
        text = hp.projtext( 360, i+7*i/75,  str(i)+r'$^{\circ}$', lonlat=True, coord='G', fontsize = 'large', fontweight = 100, zorder = 10)
        text.set_path_effects([patheffects.withStroke(linewidth=5, foreground='w')])
    for i in range(0,90,15):
        text = hp.projtext( 360, i-3*i/75,  str(i)+r'$^{\circ}$', lonlat=True, coord='G', fontsize = 'large', fontweight = 100, zorder = 10)
        text.set_path_effects([patheffects.withStroke(linewidth=4, foreground='w')])
        
    #saves plots as images
    #for ftype in ['pdf', 'png']:
    #    plt.savefig('skymaps/skymap_12p5.'+ftype, transparent = False, bbox_inches = 'tight', dpi = 300)
    
    return

#this plots the skymap using 'skymap_smoothed' and 'data', the outputs from the functions "skymap_ul()" and "skymap_data()"
