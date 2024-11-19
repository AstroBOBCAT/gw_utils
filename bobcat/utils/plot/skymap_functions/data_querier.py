#This is the BOBcat data query function. It is still a WIP. Th goal is to have this function (in one form or another) interpret queries to BOBcat and apply them to the database so that queried candidates can be plotted in skymaps.

#importing tools
import pandas as pd

#function:
def skymap_data():
    """
    This is the BOBcat database query selection function. It is intended to parse the BOBcat database for 
    candidates given a query of desired values (ex. Mass 1 > 1e8 solar masses, or 0.5 < q < 0.6, etc.). This 
    function will then create a new dataframe that can then be parsed through the skymap generator
    """
    #A more refined version of this function would call upon the actual BOBcat database and query it
    #The following code is temporary as a means of testing that these functions are compatible given a dataframe
    #of chosen candidates
    
    #reading csv of candidates
    df = pd.read_csv('Test_Candidates_with_Coordinates.csv')
    data = df[:]
    #querying of df? to select reduced list?
    return data
    
    
    #df = BOBcat csv
    #query appplied
    #data = df post-query


#Eventually, this will return a proper dataframe of chosen candidates to be plotted in the skymap_plotter
