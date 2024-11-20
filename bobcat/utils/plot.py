import numpy as np
import healpy as hp
import pandas as pd
import os

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import colors
from matplotlib import patheffects
from matplotlib import text
%matplotlib inline
%config InlineBackend.figure_format = 'retina'

from pylab import cm
import astropy.units as u
import astropy.constants as const
from astropy.coordinates import SkyCoord



# Need this to reference specific data sources.
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
_install_dir = os.path.abspath(os.path.dirname(__file__))


JESSICA we need to add the upper limit h vs f plots here... this seems to only be skymap plots?


######!!!!! SARAH THINKS THIS FUNCTION NEEDS TO BE NOT PART OF THIS CODE.
def skymap_cands(cand_file):
    """.

    This is the BOBcat database query selection function. It is
    intended to parse the BOBcat database for candidates given a query
    of desired values (ex. Mass 1 > 1e8 solar masses, or 0.5 < q <
    0.6, etc.). This function will then create a new dataframe that
    can then be parsed through the skymap generator

    We have a test candidate file in data/updated_candidates.csv

    """
    #A more refined version of this function would call upon the actual BOBcat database and query it
    #The following code is temporary as a means of testing that these functions are compatible given a dataframe
    #of chosen candidates

    cand_file = os.path.join(_install_dir, "data/updated_candidates.csv")
    
    #reading csv of candidates
    df = pd.read_csv(cand_file)
    data = df[:]
    #querying of df? to select reduced list?
    return data
    
    
    #df = BOBcat csv
    #query appplied
    #data = df post-query
    #Eventually, this will return a proper dataframe of chosen candidates to be plotted in the skymap_plotter


def skymap_ul(ul_file, nside):
    """.
    
    Read upper-limit data from a PTA dataset and smooth it to nside
    pixelization.  Currently, it has only been tested with NANOgrav
    datasets. To choose a dataset, enter the name of the text file
    pertaining to that dataset's upper limit data (ex. For 11 year
    NANOgrav data in file called 11yr_skymap_v4.txt, enter the string:
    '11yr_skymap_v4.txt')

    Inputs:

       ul_file: Raw txt file of strain UL data for sky position
                JESSICA ADD MORE DETAIL HERE. HOW IS THE SKY DATA
                PACKED? IS A PARTICULAR CONVENTION EXPECTED? IS THIS
                HEALPIX?

       nside: pixel count matched to datapoints in the txt file. (Ex:
              nside=8 --> pixels=768) JESSICA ADD MORE DETAIL HERE. IS
              THIS nside THE SIZE OF THE INPUT FILE OR IS IT THE
              DESIRED SMOOTHING LEVEL?

       JESSICA SHOULD WE ADD A DESIRED SMOOTHING LEVEL AS AN INPUT (with default=0.1?)

    Outputs:

       Smoothed skymap JESSICA ADD FORMAT INFORMATION. IS THIS OUTPUT
                       A DICTIONARY? A LIST? A TUPLE?

    """
    
    #SKYMAP IMPORT
    skymap = np.loadtxt(ul_file, skiprows = 1)
    #calling UL column
    m=skymap[:,3]        #UL data is in fourth column
    
    #nside = 8 #corresponds to 768 pixels
    skymap_smoothed = hp.smoothing(m, 0.1)
    
    min_location = np.where(skymap_smoothed == skymap_smoothed.min())[0]
    theta_min_s, phi_min_s = hp.pix2ang(nside, min_location)
    max_location = np.where(skymap_smoothed == skymap_smoothed.max())[0]
    theta_max_s, phi_max_s = hp.pix2ang(nside, max_location)
    
    thetas = skymap[:,1] #theta data is in second column
    phis = skymap[:,2]   #phi data is in third column
    
    return skymap_smoothed

    


def skymap_plot(ul_file, cands, nside):
    """.
    
    Plots a list of candidates (cands) on top of the smoothed strain
    upper limit data (skymap_smoothed) with an imaging resolution
    nside.

    Inputs:
        skymap_smoothed: PTA Upper limit sky map data JESSICA DETAIL THE EXPECTED FORMAT AND FILE TYPE

        cands: a dataframe of candidates with coordinates JESSICA DETAIL THE EXPECTED FORMAT AND INPUT TYPE

        nside: a pixel count parameter (nside) to generate complete skymaps.

    """

    # Get smoothed skymap
    skymap_smoothed = skymap_ul(ul_file, nside) #calls the strain upper limit function

    
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
    hp.visufunc.projscatter(cands['Right Ascension (Decimal Degrees)'], cands['Declination (Decimal Degrees)'], #cands[2], cands[3], 
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




# !!!!!! SARAH THINKS THIS IS REDUNDANT.
def skymap_generator(ul_file, nside):
    """.

    This is the BOBcat skymap plotter. It acts as an all-in-one skymap
    generator, taking inputs of UL data file, candidate query, and the
    nside parameter connecting the data file to pixel count. It
    compiles queried candidates from the BOBcat database alongside Strain
    Upper Limit data from several different PTA surveys to generate a
    skymap of candidate position relative to sky sensitivities. This
    version of the function calls upon two other functions,
    skymap_cands() and skymap_ul(), in order to build the map.


    Inputs:

        ul_file: name of the UL data file you wish to use JESSICA ADD
                 IN SOME WORDS ABOUT WHAT THIS IS. ASCII? WHAT DOES IT
                 CONTAIN?  HOW IS THE DATA FORMATTED?

        nside: number of pixels the map will contain.  nside value
               must generate a number of pixels compatible with the
               chosen UL data file. JESSICA ADD IN SOME WORDS ABOUT
               HOW TO KNOW WHAT NUMBERS ARE COMPATIBLE.

    """
    
    #calling other functions
    skymap_smoothed = skymap_ul(ul_file, nside) #calls the strain upper limit function
    cands = skymap_cands()                        #calls the query function

    # This section is still WIP. As the query function expands to be
    # able to handle queries, add a query parameter in the inputs for
    # skymap_generator()
    
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
    hp.visufunc.projscatter(cands['Right Ascension (Decimal Degrees)'], cands['Declination (Decimal Degrees)'], #cands[2], cands[3], 
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
