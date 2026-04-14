from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.ipac.ned import Ned
import numpy as np
import requests

'''.

Wrappers that call the NASA Extragalactic Database.

'''

def coord_converter(ra, dec): 
    '''.

    Convert J2000 ra and dec with units of hmsdms into J2000 ra
    and dec with units of degrees. 

    Inputs:
        ra = J2000 right ascension, units = hms
        dec = J2000 declination, units = dms 

    Outputs:
        ra = J2000 right ascension, units = degrees
        dec = J2000 declination, units = degrees

    '''

    # Make sure the ra and dec passed to the function are strings,
    # i.e. they're most likely in hmsdms units.  If not raise an
    # error.
    if not isinstance(ra, str) or not isinstance(dec, str):
        raise TypeError("ra and dec must be strings.")
    
    # Try to convert the coordinates into decimal degrees.
    try:
        # Convert the coordinates into degrees.
        coords_arr = SkyCoord(ra, dec)

        # Split the ra and dec and return them.
        ra_deg, dec_deg = np.array([coords_arr.ra.degree, coords_arr.dec.degree]) 
        return ra_deg, dec_deg
    # If something goes wrong in the above code block it is most
    # likely because the ra and dec passed to the function weren't
    # entered exactly correctly. Therefore, raise an error.
    except:
        raise SystemError("Make sure ra and dec are entered correctly.")



def coord_finder(name):
    '''.
    
    Use astropy search library to determine J2000 ra and dec of an
    object. This function will take the name of an astronomical object and
    query the astropy database with the name. The name used in this
    function does NOT have to be the NED name of the object. The function
    will return the J2000 ra and dec of the object in units of hmsdms as
    long so the object is found within astropy.

    Inputs:
        name = string of the name of the object the coordinates are needed for

    Outputs:
        ra = J2000 right ascension, units = hms
        dec = J2000 declination, units = dms

    '''

    # First we need to check that the name given is actually a string.
    # If it isn't a string we could change it into a string but it is
    # probably better to raise an error to make sure the user is aware
    # of what is wanted for the function.
    if not isinstance(name, str):
        raise TypeError("Name must be a string.")

    # Try to find the name in the astropy databases.
    try:
        # Search for the coordinates of the object given the name. 
        coords = SkyCoord.from_name(name, parse = True).to_string('hmsdms').split()

        # Assign the ra and dec variables.
        ra = coords[0]
        dec = coords[1]

        # Return an array of the ra and dec.
        return (ra, dec)
    
    # If the name is not in the astropy databases then raise an error.
    # (There will most likely be another error raised by the query
    # failing that says the same thing).
    except:
        raise SystemError("Name of source is in shorthand or not in an established database.")



def name(ra, dec, radius_tol = 0.01):
    '''.
    
    Use NED search library to determine name of nearest NED-listed
    galaxy within a narrow position search. This function takes the ra
    and dec, in degrees, of an object and and quieries the NED
    database for the names associated with the objects that are within
    a radius from the ra and dec passed to NED. Then the objects
    are sorted so that the object with the smallest separation is
    selected. This is essentially picking out the object with ra and
    dec values closest to those passed to the function. The name of
    this object is returned by the function.

    Inputs:

        ra = J2000 right ascension, units = degrees

        dec = J2000 declination, units = degrees

        radius_tol = tolerance of the search radius used in the query
                     to the NED database, units = degrees, default =
                     0.01
    
    Outputs:

        name = name of the object closest to the given ra and dec

    '''

    # Check that the arugments are actually numerical, i.e. the ra and
    # dec passed to the function are in degrees.
    NEDInstance = Ned() # Establish an instance of NED to change its parameters

    # Change the timeout parameter of the NED instance to 150 seconds. Requires making a separate request.
    NEDInstance._session.request = lambda *args, **kwargs: requests.Session().request(
        *args,
        timeout=150,
        **kwargs
)
    if isinstance(ra, (int, float)) and isinstance(dec, (int, float)):
        # Put coordinates into an array.
        coords = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')

        # Query the NED database for objects within radius_tol
        # degrees, radially, of the chosen point (ra,dec), and creates
        # a table.
        result_table = NEDInstance.query_region(coords, radius=radius_tol*u.deg, equinox='J2000.0',timeout=150)

        # Find object closest to the center of the region and drops
        # all other objects from the dataframe.
        result = result_table[result_table['Separation']== result_table['Separation'].min()]
        
        # Extract the name from dataframe and return it.
        name = result[0]['Object Name'] 
        return name
    
    # If the arguments are not numerical, i.e. of hmsdms, then raise
    # an error.
    else:
        raise RuntimeError("Arguments given are incorrect. Make sure ra and dec are in degrees.")




def name_resolver(name):
    '''.
    
    Take an astronomical object's name and return the associated NED
    name for the object.  This function allows the user to enter an
    name of an astronomical object and get back the first name
    associated with the object in the NED database. This is useful as
    there are many different naming schemes through astronomy and the
    results are that most objects have multiple names. The "NED name"
    that is returned from the function is just the first name in the
    list of names NED has associated with it. This is not a complete
    fail safe as there are some naming schemes or objects' names that
    aren't in the NED database even if the actual object is.
     
    Inputs:
        name string of astronomical object
    Outputs:
        NED name (or the first name in the list of names) for the astronomical object

    '''

    # First we need to check that the name given is actually a string.
    # If it isn't a string we could change it into a string but it is
    # probably better to raise an error to make sure the user is aware
    # of what is wanted for the function.
    if not isinstance(name, str):
        raise TypeError("Name must be a string.")

    # Try to find the name in the NED database.
    try:

        # Query the NED database for the object based on the name
        #passed to the function. 

        #Ned.query_object(name)  
        #### Note, the query_object code seems to have weird things
        #### happening to it so we are not using it.  Like it'll give
        #### me normal results for "3c66b" but completely different
        #### and werid results when I give it "3C66B". It also was
        #### except "sam", "abs", "asfsd", and other things that
        #### definitely aren't in the NED database.
        ra, dec = coord_finder(name)

        # Convert ra and dec from hmsdms to degrees.
        ra_deg, dec_deg = coord_converter(ra, dec)

        # Pull just the object's name out of the table (stripping the extra white spaces).
        #name = str(result_table['Object Name']).replace('Object Name\n-----------\n', '').strip()
        res_name = name(ra_deg, dec_deg)

        # Return the NED name
        return (res_name)
    # If the name is not in the NED database then raise an error. (There will most likely be another error raised 
    # by the NED query failing that says the same thing).
    except:
        raise SystemError("Name given is not in NED database.")




def redshift(ra, dec, ned_name, radius_tol = 0.01):

    '''.
    
    Find the redshift associated with an J2000 ra and dec (in degrees)
    with some search radius tolerance from the NED database. This
    function takes the ra and dec (in the units of degrees) and
    quieries the NED database for the redshift values associated with
    the objects that are within a radius of 0.01 degrees from the ra
    and dec passed to NED. Then the objects are sorted so that the
    object with the smallest separation is selected. This is
    essentially picking out the object with ra and dec values closest
    to those passed to the function. The redshift of this object is
    returned by the function.

    Inputs:
        ra = J2000 right ascension, units = degrees

        dec = J2000 declination, units = degrees

        radius_tol = tolerance of the search radius used in the query
                     to the NED database, units = degrees, default = 0.01
    
    Outputs:

        z = redshift of the object closest to the given ra and dec

    '''

    NEDInstance = Ned()

    NEDInstance._session.request = lambda *args, **kwargs: requests.Session().request(
        *args,
        timeout=150,
        **kwargs
)

    # Check that the arugments are actually numerical, i.e. the ra and
    # dec passed to the function are in degrees.
    if isinstance(ra, (int, float)) and isinstance(dec, (int, float)):
        # Put coordinates into an array.
        coords = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')

        # Query the NED database for objects within radius_tol
        # degrees, radially, of the chosen point (ra,dec), and creates
        # a table.
        for attempts in range(3): # Try 3 times, in case of timeout or other issue that repetition will solve.
            try:
                result_table = NEDInstance.query_region(coords, radius=radius_tol*u.deg, equinox='J2000.0')
            except:
                print(f"NED query failed {attempts+1}/3, retrying...")
                if attempts == 2:
                    raise
                continue
        # Find object closest to the center of the region and drops
        # all other objects from the dataframe.
        result = result_table[result_table['Separation'] == result_table['Separation'].min()]
        # Extracts redshift value from dataframe and return it.
        z = result[0]['Redshift'] 
        z = float(z) # Convert to float. Will produce a warning if more than one redshift is found.
        if np.isnan(z): # If the redshift couldn't be found, try to find it by name.
            print("Coordinate search for redshift failed, attempting search by name.")
            z = name_redshift(ned_name, NEDInstance)
        return z
    
    # If the arguments are not numerical, i.e. of hmsdms, then raise
    # an error.
    else:
        raise RuntimeError("Arguments given are incorrect. Make sure ra and dec are in degrees.")

def name_redshift(ned_name, NEDInstance):

    '''
    Backup function for redshift. If the coordinate search for redshift fails, this function will be called.

    Inputs:
        ned_name = NED name of the object
        NEDInstance = NED Instance, from the redshift() function.
    
    Outputs:
        z = redshift of the object
    '''
    for attempt in range(3):
        try:
            result = NEDInstance.query_object(ned_name)
            z = result['Redshift'][0]
            z = float(z)
            if np.isnan(z):
                print("Name search for redshift returned masked/NaN.")
                raise
            return z
        except Exception as e:
            print(f"NED name query failed {attempt+1}/3: {e}")
    # This part of the code is only called if it fails 3 times.
    print("NED name query failed 3 times.")
    raise