## NED_name.py contains the code for the function NED_name(). This function takes the ra and dec, in degrees,
## of an object and and quieries the NED database for the names associated with the objects that are within a 
## radius from the ra and dec passed to NED. Then the objects are sorted so that the object with the 
## smallest separation is selected. This is essentially picking out the object with ra and dec values
## closest to those passed to the function. The name of this object is returned by the function.


######
# Import the need libraries and modules for the function to work.
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.ipac.ned import Ned #needed to connect to the NED database
######


###############
def NED_name(ra, dec, radius_tol = 0.01):
    '''
    Use NED search library to determine name of nearest NED-listed galaxy
    within a narrow position search.

    Inputs:
        ra = J2000 right ascension, units = degrees
        dec = J2000 declination, units = degrees
        radius_tol = tolerance of the search radius used in the query to the NED database, units = degrees, default = 0.01
    
    Outputs:
        name = name of the object closest to the given ra and dec, units = NA
    '''

    # Check that the arugments are actually numerical, i.e. the ra and dec passed to the function
    # are in degrees.
    if isinstance(ra, (int, float)) and isinstance(dec, (int, float)):
        # Put coordinates into an array.
        coords = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')

        # Query the NED database for objects within radius_tol degrees, radially, 
        # of the chosen point (ra,dec), and creates a table.
        result_table = Ned.query_region(coords, radius=radius_tol*u.deg, equinox='J2000.0')

        # Find object closest to the center of the region and drops all other objects from the dataframe.
        result = result_table[result_table['Separation']== result_table['Separation'].min()]
        
        # Extract the name from dataframe and return it.
        name = result[0]['Object Name'] 
        return name
    
    # If the arguments are not numerical, i.e. of hmsdms, then raise an error.
    else:
        raise RuntimeError("Arguments given are incorrect. Make sure ra and dec are in degrees.")
###############


## NED_name_resolver.py contains code for the NED_name_resolver() funciton. This function allows the user
## to enter an name of an astronomical object and get back the first name associated with the object
## in the NED database. This is useful as there are many different naming schemes through astronomy
## and the results are that most objects have multiple names. The "NED name" that is returned from the
## function is just the first name in the list of names NED has associated with it. This is not a complete
## fail safe as there are some naming schemes or objects' names that aren't in the NED database even if the
## actual object is. 


######
# Import the need libraries and modules for the function to work.
from astroquery.ipac.ned import Ned #allows the connection to the NED database
from .coordinate_finder import coord_finder
from .coordinate_converter import coord_converter
from .NED_name import NED_name
######

###################
def NED_name_resolver(name):
    ''' 
    Take an astronomical object's name and return the associated NED name for the object.
     
    Inputs:
        name string of astronomical object
    Outputs:
        NED name (or the first name in the list of names) for the astronomical object
    '''

    # First we need to check that the name given is actually a string.
    # If it isn't a string we could change it into a string but it is probably better to raise an error to
    # make sure the user is aware of what is wanted for the function.
    if not isinstance(name, str):
        raise TypeError("Name must be a string.")

    # Try to find the name in the NED database.
    try:

        # Query the NED database for the object based on the name passed to the function.
        #result_table = Ned.query_object(name) #### This code seems to have weird things happening to it.
        #### Like it'll give me normal results for "3c66b" but completely different and werid results when I
        #### give it "3C66B". It also was except "sam", "abs", "asfsd", and other things that definitely
        #### aren't in the NED database.
        ra, dec = coord_finder(name)

        # Convert ra and dec from hmsdms to degrees.
        ra_deg, dec_deg = coord_converter(ra, dec)

        # Pull just the object's name out of the table (stripping the extra white spaces).
        #NED_name = str(result_table['Object Name']).replace('Object Name\n-----------\n', '').strip()
        res_name = NED_name(ra_deg, dec_deg)

        # Return the NED name
        return (res_name)
    # If the name is not in the NED database then raise an error. (There will most likely be another error raised 
    # by the NED query failing that says the same thing).
    except:
        raise SystemError("Name given is not in NED database.")
##################