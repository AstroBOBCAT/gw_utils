#This is the BOBcat Upper Limit smoothing and processing function. It takes raw UL data from PTA surveys and smooths it such it can be properly plotted.
#It takes two inputs: 1) the raw txt file of strain UL data for sky position, and 2) an "nside" parameter that matches pixel count to datapoints in the txt file. (Ex: nside=8 --> pixels=768)

#importing tools
import numpy as np
import healpy as hp

#function:
def skymap_ul(ul_file, nside):
    """
    This BOBcat function is intended to select the UL data from a chosen PTA dataset and then process it for use
    in the plotted map. Currently, it has only been tested with NANOgrav datasets. To choose a dataset, enter the 
    name of the text file pertaining to that dataset's upper limit data (ex. For 11 year NANOgrav data in file  
    called 11yr_skymap_v4.txt, enter the string: '11yr_skymap_v4.txt')
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

#from here, the "skymap_smoothed" can be plotted using the plotting functions
