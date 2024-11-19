#This is the BOBcat skymap generator function. It acts as an all-in-one skymap generator, taking inputs of UL data file, candidaet query, and the nside parameter connecting the data file to pixel count.

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


#creating the function
def skymap_generator(ul_file, nside):
    """
    This is the BOBcat skymap plotter. It compiles queried data from the BOBcat database alongside 
    Strain Upper Limit data from several different PTA surveys to generate a skymap of candidate position
    relative to sky sensitivities. This version of the function calls upon two other functions, skymap_data()
    and skymap_ul(), in order to build the map. For inputs, please enter the name of the UL data file you
    wish to use followed by the nside parameter, which defines the number of pixels the map will contain.
    Please ensure that this nside value generates a number of pixels compatible with the chosen UL data file.
    """
    
    #calling other functions
    skymap_smoothed = skymap_ul(ul_file, nside) #calls the strain upper limit function
    data = skymap_data()                        #calls the query function  #this section is still WIP. As the query function expands to be able to handle queries, add a query parameter in the inputs for skymap_generator()
    
    #SETUP SECTION
    #plot label parameters
    plt.rc('xtick',**{'labelsize':16})
    plt.rc('ytick',**{'labelsize':16})
    plt.rc('axes',**{'labelsize':20,'titlesize':18})
    
    #colorbar tick labels
    cbar_ticks = [2e-15, 3e-15, 4e-15, 5e-15, 6e-15, 7e-15, 8e-15, 9e-15, 1e-14]
    cbar_labels = ['$2 \times 10^{-15}$', '$3 \times 10^{-15}$', '$4 \times 10^{-15}$', ' ', '$6 \times 10^{-15}$', ' ',
               ' ' , ' ', '$10^{-14}$']
    
    #PLOTTING SECTION
    #assigning number of pixels to plot
    npix = hp.nside2npix(nside)
    
    #creating mollview plot
    mv = hp.mollview(skymap_smoothed, title="", cbar = None , cmap = 'viridis_r', rot = 180, norm='log',
#                        min = min(skymap11[0:768,3]), max = max(skymap11[0:768,3]),   
                        bgcolor='white', badcolor = 'white',
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

#this should return a completed plot
