## find_m1_m2.py holds the code for the find_m1_m2() function. This function takes in all
## six mass values used to describe a binary system. From any two out of the six values,
## m1 and m2 can be found. We prioritize m1 and m2 values because calculating any of the 
## other four mass values is easiest from m1 and m2. Therefore, this function will determine
## which two (minimum) mass values were passed to it and then calculate m1 and m2 from there.
## The elif statements are put in an order such that the easiest calculations are put at the
## beginning and should complete first. This means if more than two mass values are passed to 
## the function, the simplest calculations will occur to find m1 and m2 values.




##################
def find_m1_m2(m1 = None, m2 = None, Mtot = None, q = None, Mc = None, mu = None):
    '''
    Find the values of m1 and m2 for a binary system given at least two of the six common
    mass parameters, i.e. a minimum of two values from m1, m2, Mtot, q, Mc, and mu are 
    needed to calculate m1 and m2.

    Inputs:
        m1 = mass of first binary object, units = same as the other mass value given, default = None
        m2 = mass of second binary object, units = same as the other mass value given, default = None
        Mtot = total mass of the system, units = same as the other mass values given, default = None
        q = mass ratio of the system, units = NA, default = None
        Mc = chirp mass of the system, units = same as the other mass values given, default = None
        mu = reduced mass of the system, units = same as the other mass values given, default = None
    
    Outputs:
        m1 = mass of first binary object, units = same as the units of the masses passed to function
        m2 = mass of second binary object, units = same as the units of the masses passed to function
    '''
    
    # Check if m1 and m2 are actually passed to the function. If they are then we can skip the remaining 
    # code and just pass them back.
    if m1!=None and m2!=None:
        pass
    
    # If we have either m1 or m2 and Mtot, the other m1 or m2 can be found easily. Mtot = m1 + m2
    elif m1!=None and Mtot!=None:
        m2 = Mtot-m1
    elif m2!=None and Mtot!=None:
        m1 = Mtot-m2

    # If we have either m1 or m2 and q, the other m1 or m2 can be found easily. q = m2/m1
    elif m1!=None and q!=None:
        m2 = q*m1
    elif m2!=None and q!=None:
        m1 = m2/q

    # If we have both q and mu, we can easily find m1 and m2.
    elif q!=None and mu!=None:
        m1 = mu*(q+1)/q 
        m2 = mu*(q+1)

    # If we have either m1 or m2 and mu, the other m1 or m2 can be found easily. mu = (m1*m2) / (m1+m2)
    elif m1!=None and mu!=None:
        m2 = (mu*m1)/(m1-mu)
    elif m2!=None and mu!=None:
        m1 = (mu*m2)/(m2-mu)

    # If we have both Mtot and q, we can find m1 and m2 easily.
    elif Mtot!=None and q!=None:
        m1 = Mtot/(q+1)
        m2 = (q*Mtot)/(q+1)

    # If we have both q and Mc, we can rearrange their equation to easily find m1 and m2.
    elif q!=None and Mc!=None:
        m1 = Mc*((q+1)/q**3)**(1/5)
        m2 = Mc*((q+1)*q**2)**(1/5)

    # If we have both Mtot and mu, we can rearrange their equations to find m1 and m2.
    elif Mtot!=None and mu!=None:
        m1 = (1/2)*(Mtot + (Mtot**2 - 4*(Mtot*mu))**(1/2))
        m2 = (1/2)*(Mtot - (Mtot**2 - 4*(Mtot*mu))**(1/2))

    # If we have Mc and mu, we can rearrange their equations to find m1 and m2.
    elif Mc!=None and mu!=None:
        m1 = mu*(1 + ((Mc**(5/2)-mu**(3/2))/(mu**(3/2))))
        m2 = mu*(1 + (mu**(3/2)/(Mc**(5/2)-mu**(3/2))))

    # If we have both Mtot and Mc, we can rearrange their equations to find m1 and m2.
    # There is an additional limit of Mc<=(0.43527528165*Mtot) which will prevent the 
    # calculator from outputting complex numbers by ensuring the chirp mass is in 
    # agreement with the total mass first. If they do not agree, the calculator will cut to the error.
    elif Mtot!=None and Mc!=None and Mc<=(0.43527528165*Mtot): 
        m1 = (1/2)*(Mtot + (Mtot**2 - 4*(Mtot*(Mc**5))**(1/3))**(1/2))
        m2 = (1/2)*(Mtot - (Mtot**2 - 4*(Mtot*(Mc**5))**(1/3))**(1/2))

    # If we have Mc and Mtot but Mc doesn't follow the restriction given to if by the Mtot,
    # then raise an error.
    elif Mc!=None and Mtot!=None and Mc>=(0.43527528165*Mtot):
            raise RuntimeError("Chrip mass too large to be correct, please check") 

    # If we have either m1 or m2 and Mc, the other m1 or m2 can be found be rearranging the chirp mass equation.
    # However, the equation involves several radicals, so as values of Mc get further from the actual (given m1),
    # there will be an increasingly larger complex portion. However, this portion is typically small,
    # and so this ".real" simply deletes the complex part so calculations can continue as normal.
    # This is currently the best way we can get around this odd part of calculations.
    elif m1!=None and Mc!=None:
        m2 = (((2/3)**(1/3))*(Mc**5))/(((9*(Mc**5)*(m1**7))+(3*((27*(Mc**10)*(m1**14))-(4*(Mc**15)*(m1**9))))**(1/2))**(1/3)) + (((9*(Mc**5)*(m1**7))+(3*((27*(Mc**10)*(m1**14))-(4*(Mc**15)*(m1**9))))**(1/2))**(1/3))/(((18)**(1/3))*(m1**3))
        m2 = round(m2.real, 16)
    elif m2!=None and Mc!=None:
        m1 = (((2/3)**(1/3))*(Mc**5))/(((9*(Mc**5)*(m2**7))+(3*((27*(Mc**10)*(m2**14))-(4*(Mc**15)*(m2**9))))**(1/2))**(1/3)) + (((9*(Mc**5)*(m2**7))+(3*((27*(Mc**10)*(m2**14))-(4*(Mc**15)*(m2**9))))**(1/2))**(1/3))/(((18)**(1/3))*(m2**3))
        m1 = round(m1.real, 16)

    # If we don't have at least two mass parameters we need to raise error for the function.
    else:
        raise RuntimeError("Not enough mass parameters to calculate m1 and m2. Need at least two mass values.")

    # Before returning m1 and m2, we want to make sure q will be betweeon 0 and 1 so we check if we
    # need to switch them or not.
    if m1!=None and m2!=None and m2 > m1:
        m1, m2 = m2, m1
    else:
        pass

    # Finally return m1 and m2 values.
    return(m1,m2)
 #####################


 ## Mc_calculator.py holds the code for the function Mc_calc(). This function will calculate
## the chirp mass of a system given the two masses within the system. It uses the equation
## (((m1*m2)**3)/(m1+m2))**(1/5)


#############
def Mc_calc(m1,m2):
    '''
    Calculate the chirp mass of a system given two masses.

    Inputs:
        m1 = mass of first binary object, units = same as the other mass value given
        m2 = mass of second binary object, units = same as the other mass value given
    
    Outputs:
        Mc = chirp mass of the system, units = same as the units of the masses passed to function
    '''
    
    return (((m1*m2)**3)/(m1+m2))**(1/5)
##############


## Mtot_calculator.py holds the code for the function Mtot_calc(). This function will calculate
## the total mass of a system given the two masses within the system. It uses the equation
## m1+m2


#############
def Mtot_calc(m1,m2):
    '''
    Calculate the total mass of a system given two masses.

    Inputs:
        m1 = mass of first binary object, units = same as the other mass value given
        m2 = mass of second binary object, units = same as the other mass value given
    
    Outputs:
        Mtot = total mass of the system, units = same as the units of the masses passed to function
    '''
    
    return m1+m2
##############
  

  ## mu_calculator.py holds the code for the function mu_calc(). This function will calculate
## the reduced mass of a system given the two masses within the system. It uses the equation
## (m1*m2)/(m1+m2)


#############
def mu_calc(m1,m2):
    '''
    Calculate the reduced mass of a system given two masses.

    Inputs:
        m1 = mass of first binary object, units = same as the other mass value given
        m2 = mass of second binary object, units = same as the other mass value given
    
    Outputs:
        mu = reduced mass of the system, units = same as the units of the masses passed to function
    '''
    
    return (m1*m2)/(m1+m2) 
############ 

## q_calculator.py holds the code for the function q_calc(). This function will calculate
## the mass ratio of a system given the two masses within the system. It uses the equation
## m(smaller)/m(bigger) so that q is always between 0 and 1.


##################
def q_calc(m1,m2):
    '''
    Calculate the mass ratio of a system given two masses.

    Inputs:
        m1 = mass of first binary object, units = same as the other mass value given
        m2 = mass of second binary object, units = same as the other mass value given
    
    Outputs:
        q = mass ratio of the system, units = NA
    '''

    # Check to see if m1>m2, which is what we want to force so that q remains between 0 and 1.
    if m1 > m2:
        return m2/m1 
    # If m2>m1 then we just change the ratio so that q is still between 0 and 1.
    else:
        return m1/m2
###################




#This is the BOBcat SMBHB Candidate Mass Value Calculator. It will read provided inputs of mass quantities and can calculate unknown mass quantities.
#Currently, the BOBcat Mass Value Calculator is equipped to calculate Mass 1, Mass 2, Total Mass, Mass Ratio, Chirp Mass, and Reduced Mass.


#THE MASS CALCULATING FUNCTION:
#As of right now, this calculator only fills in empty values, and cannot be a judge of erroneous/incongruent inputs.
#If the inputs were found via different methods/models, try entering only the inputs that used a consistent model.
def mass_val_calc(m1 = None, m2 = None, Mtot = None, q = None, Mc = None, mu = None): #creating the function, reading the collection of inputs
    """
    This is the BOBcat SMBHB mass value calculator! It will calculate estimates of Mass 1, Mass 2, Total Mass, Mass Ratio,
    Chirp Mass, and Reduced Mass (in solar masses) given known mass parameters of the binary candidate (in solar masses). 
    To successfully calculate the unknown values, at least two known parameters are needed as inputs.

    Args: (all in solar masses)
    m1 = Larger Mass
    
    m2 = Smaller Mass
    Mtot = Total Mass
    q = Mass Ratio
    Mc = Chirp Mass
    mu = Reduced Mass
    These are also the outputs
    """
    #m1_new, m2_new = find_m1_m2(m1, m2, Mtot, q, Mc, mu)
    
    
    #Here, we calculate all other values using m1 and m2

    #Mtot_new = Mtot_calc(m1,m2) #solving for Total Mass using Mass 1 and Mass 2
   #q_new = q_calc(m1,m2) #solving for Mass Ratio using Mass 1 and Mass 2
    #Mc_new = Mc_calc(m1,m2) #solving for Chirp Mass using Mass 1 and Mass 2
    #mu_new = mu_calc(m1,m2) #solving for Reduced Mass using Mass 1 and Mass 2
        
    #mass_array, updated_masses = update_masses(m1, m2, Mtot, q, Mc, mu)
    
    none_counter = 0


    if m1 == None:
        none_counter+=1
    if m2 == None:
        none_counter+=1
    if Mtot == None:
        none_counter+=1
    if q == None:
        none_counter+=1
    if Mc == None:
        none_counter+=1
    if mu == None:
        none_counter+=1


    mass_array, updated_masses = update_masses(m1, m2, Mtot, q, Mc, mu)

    if none_counter != len(updated_masses):
        raise RuntimeError("Something is wrong. Please check that you entered the correct numbers for each argument.")
    else:
        return mass_array 


# #This is the BOBcat SMBHB Candidate Mass Value Calculator. It will read provided inputs of mass quantities and can calculate unknown mass quantities.
# #Currently, the BOBcat Mass Value Calculator is equipped to calculate Mass 1, Mass 2, Total Mass, Mass Ratio, Chirp Mass, and Reduced Mass.


# #THE MASS CALCULATING FUNCTION:
# #As of right now, this calculator only fills in empty values, and cannot be a judge of erroneous/incongruent inputs.
# #If the inputs were found via different methods/models, try entering only the inputs that used a consistent model.
# def mass_val_calc(m1,m2,Mtot,q,Mc,mu): #creating the function, reading the collection of inputs
#     """
#     This is the BOBcat SMBHB mass value calculator! It will calculate estimates of Mass 1, Mass 2, Total Mass, Mass Ratio,
#     Chirp Mass, and Reduced Mass (in solar masses) given known mass parameters of the binary candidate (in solar masses). 
#     To successfully calculate the unknown values, at least two known parameters are needed as inputs.

#     Args: (all in solar masses)
#     m1 = Larger Mass
#     # m2 = Smaller Mas
#     Mtot = Total Mass
#     q = Mass Ratio
#     Mc = Chirp Mass 
#     mu = Reduced Mass
#     These are also the outputs
#     """
#     #Setting all strings (strings replacing NaNs) to None
#     if type(m1)!=int and type(m1)!=float:
#         m1 = None
#     if type(m2)!=int and type(m2)!=float:
#         m2 = None
#     if type(Mtot)!=int and type(Mtot)!=float:
#         Mtot = None
#     if type(q)!=int and type(q)!=float:
#         q = None
#     if type(Mc)!=int and type(Mc)!=float:
#         Mc = None
#     if type(mu)!=int and type(mu)!=float:
#         mu = None
    
#     #Beginning calculations of m1 and m2:
#     if m1!=None and m2!=None: #m1 and m2 are our priority variables. These if/elif statements are used to find which parameters we are given. From there, we calculate m1, m2, or both.
#         pass #when given m1 and m2, we can jump straight to the calculation formulas.

#     #Below, the function will be choosing variables in a priority order, as in it will calculate any missing values using the first pair of inputs it detects
#     #given our series of elif statements. This priority order is guided by which parameters we believe to be most commonly reported and most easily measurable.
#     #The order is: mass 1, mass 2, total mass, mass ratio, chirp mass, and then reduced mass
    
#     elif m1!=None and Mtot!=None:
#         m2 = Mtot-m1
#     elif m2!=None and Mtot!=None:
#         m1 = Mtot-m2
#     elif m1!=None and q!=None:
#         m2 = q*m1
#     elif m2!=None and q!=None:
#         m1 = m2/q
#     elif m1!=None and Mc!=None:
#         m2 = (((2/3)**(1/3))*(Mc**5))/(((9*(Mc**5)*(m1**7))+(3*((27*(Mc**10)*(m1**14))-(4*(Mc**15)*(m1**9))))**(1/2))**(1/3)) + (((9*(Mc**5)*(m1**7))+(3*((27*(Mc**10)*(m1**14))-(4*(Mc**15)*(m1**9))))**(1/2))**(1/3))/(((18)**(1/3))*(m1**3))
#         m2 = round(m2.real, 16) #the above equation involves several radicals, so as values of Mc get further from the actual (given m1), there will be an increasingly larger complex portion. However, this portion is typically small, and so this ".real" simply deletes the complex part so calculations can continue as normal.
#     elif m2!=None and Mc!=None:
#         m1 = (((2/3)**(1/3))*(Mc**5))/(((9*(Mc**5)*(m2**7))+(3*((27*(Mc**10)*(m2**14))-(4*(Mc**15)*(m2**9))))**(1/2))**(1/3)) + (((9*(Mc**5)*(m2**7))+(3*((27*(Mc**10)*(m2**14))-(4*(Mc**15)*(m2**9))))**(1/2))**(1/3))/(((18)**(1/3))*(m2**3))
#         m1 = round(m1.real, 16)
#     elif m1!=None and mu!=None:
#         m2 = (mu*m1)/(m1-mu)
#     elif m2!=None and mu!=None:
#         m1 = (mu*m2)/(m2-mu)
#     elif Mtot!=None and q!=None:
#         m1 = Mtot/(q+1)
#         m2 = (q*Mtot)/(q+1)
#     elif Mtot!=None and Mc!=None and Mc<=(0.43527528165*Mtot): #this additional limit will prevent the calculator from outputting complex numbers by ensuring the chirp mass is in agreement with the total mass first. If they do not agree, the calculator will cut to "else"
#         m1 = (1/2)*(Mtot + (Mtot**2 - 4*(Mtot*(Mc**5))**(1/3))**(1/2))
#         m2 = (1/2)*(Mtot - (Mtot**2 - 4*(Mtot*(Mc**5))**(1/3))**(1/2))
#     elif Mtot!=None and mu!=None:
#         m1 = (1/2)*(Mtot + (Mtot**2 - 4*(Mtot*mu))**(1/2))
#         m2 = (1/2)*(Mtot - (Mtot**2 - 4*(Mtot*mu))**(1/2))
#     elif q!=None and Mc!=None:
#         m1 = Mc*((q+1)/q**3)**(1/5)
#         m2 = Mc*((q+1)*q**2)**(1/5)
#     elif q!=None and mu!=None:
#         m1 = mu*(q+1)/q 
#         m2 = mu*(q+1)
#     elif Mc!=None and mu!=None:
#         m1 = mu*(1 + ((Mc**(5/2)-mu**(3/2))/(mu**(3/2))))
#         m2 = mu*(1 + (mu**(3/2)/(Mc**(5/2)-mu**(3/2))))
#     else:
#         print('No Calculations Possible') #if we are given only one input or no inputs, the function will tell the user nothing could be calculated
#         if Mc!=None and Mtot!=None and Mc>=(0.43527528165*Mtot):
#             print('Chirp Mass too large') #flagging an error message if chirp mass is out of bounds relative to total mass and no calculations are possible because of it
#         array = np.array([m1, m2, Mtot, q, Mc, mu]) #compiles unchanged inputs into an array
#         return array #ends function here should there not be a sufficient number of inputs

#     if m2 > m1: #here, we want to ensure q remains between 0 and 1, so we need m1 to be larger than m2. If m2 is found to be larger than m1, we switch them before calculating.
#         m1_new = m2 #assigning old values to new variables
#         m2_new = m1
#         m1 = m1_new #reinserting new variables so they satisfy m1>m2 before performing calculations
#         m2 = m2_new
#     else: #if m1>m2, then nothing needs to be switched around
#         pass
    
#     #Here, we calculate all other values using m1 and m2
#     Mtot_new = m1+m2 #solving for Total Mass using Mass 1 and Mass 2
#     q_new = m2/m1 #solving for Mass Ratio using Mass 1 and Mass 2
#     Mc_new = (((m1*m2)**3)/(m1+m2))**(1/5) #solving for Chirp Mass using Mass 1 and Mass 2
#     mu_new = (m1*m2)/(m1+m2) #solving for Reduced Mass using Mass 1 and Mass 2
        
#     #with m1 and m2 known, here we check which inputs were blank, and fill them in with their corresponding values
#     #if they are known, we cross-compare inputs versus calculated values and only maintain inputs within reason (margin)
#         #- All those outside the margin are to be flagged for checking, or reassigned to a new model
#         #- As desired, change the multiplier in the nested if statement to the appropriate error bounds (default is 10%)
#     if Mtot==None:
#         Mtot = Mtot_new #assigning a calculated value if there is none inputted
#     else:
#         if abs(Mtot-Mtot_new)<=(0.1*Mtot_new):
#             Mtot = Mtot #maintaining inputted value
#         else:
#             Mtot_old = Mtot #reassign inputted (out of range) value to "old"
#             Mtot = Mtot_new #assign new calculation to output
#     if q==None:
#         q = q_new
#     else:
#         if q>1:
#             q_old = q
#             q = q_new #q_new will always be less than 1 because of the statements reassigning the larger mass to m1 for calculations
#         if abs(q-q_new)<=(0.1*q_new):
#             q = q #maintaining inputted value
#         else:
#             q_old = q #reassign old (out of range) value to "old"
#             q = q_new #assign new calculation to output
#     if Mc==None:
#         Mc = Mc_new
#     else:
#         if abs(Mc-Mc_new)<=(0.1*Mc_new):
#             Mc = Mc #maintaining inputted value
#         else:
#             Mc_old = Mc #reassign old (out of range) value to "old"
#             Mc = Mc_new #assign new calculation to output
#     if mu==None:
#         mu = mu_new
#     else:
#         if abs(mu-mu_new)<=(0.1*mu_new):
#             mu = mu #maintaining inputted value
#         else:
#             mu_old = mu #reassign old (out of range) value to "old"
#             mu = mu_new #assign new calculation to output    

#     mass_array = np.array([m1, m2, Mtot, q, Mc, mu]) #compiles unchanged inputs and newly calculated outputs into an array
#     return mass_array #returns the updated array with newly calculated values in place

# #With this array, the new values are available to be reinserted into the BOBcat database.






## update_m1.py has the code for the update_m1 function. This function will check to see if 
## the mass 1 value passed to the function is the same or within some tolerance to the 
## mass 1 value calculated. If there is not mass 1 value passed into the function, or the 
## tolerance is exceeded, then the calculated mass 1 value becomes the value used and 
## returned as the mass 1.


###############
def update_m1(m1 = None, m2 = None, Mtot = None, q = None, Mc = None, mu = None, tolerance = 1e5):
    '''
    Update the mass 1 value given to that calculated if the tolerance is not met.

    Inputs:
        m1 = mass of first binary object, units = same as the other mass value given, default = None
        m2 = mass of second binary object, units = same as the other mass value given, default = None
        Mtot = total mass of the system, units = same as the other mass values given, default = None
        q = mass ratio of the system, units = NA, default = None
        Mc = chirp mass of the system, units = same as the other mass values given, default = None
        mu = reduced mass of the system, units = same as the other mass values given, default = None
        tolarance = what the difference has to be less than or equal to in order to not replace
        the given value with that calculated, units = NA, default = 1e5
    
    Outputs:
        m1 = mass of first binary object, units = same as the units of the masses passed to function
    '''
    
    # Try to calculate m1 from whatever given mass values you have using the find_m1_m2 function
    try:
        m1_calced = find_m1_m2(m1, m2, Mtot, q, Mc, mu)[0]
    # If the find_m1_m2 function fails, then raise an error becuase none of the mass values can be calculated.
    except:
        raise ValueError("Not enough mass values available to calculate m1 and m2.")

    # If no m1 was passed to the function or the difference between the value passed and
    # the value calculated is greater than the tolerance, the calculated m1 value is the
    # the mass value used for m1.
    if m1 == None or abs(m1 - m1_calced) >= tolerance:
        m1 = m1_calced

    # If none of the above criteria are met then the m1 value passed to the function
    # stays as the m1 value.
    else:
        pass
    return m1
##############



## update_m2.py has the code for the update_m2 function. This function will check to see if 
## the mass 2 value passed to the function is the same or within some tolerance to the 
## mass 2 value calculated. If there is not mass 2 value passed into the function, or the 
## tolerance is exceeded, then the calculated mass 2 value becomes the value used and 
## returned as the mass 2.


###############
def update_m2(m1 = None, m2 = None, Mtot = None, q = None, Mc = None, mu = None, tolerance = 1e5):
    '''
    Update the mass 2 value given to that calculated if the tolerance is not met.

    Inputs:
        m1 = mass of first binary object, units = same as the other mass value given, default = None
        m2 = mass of second binary object, units = same as the other mass value given, default = None
        Mtot = total mass of the system, units = same as the other mass values given, default = None
        q = mass ratio of the system, units = NA, default = None
        Mc = chirp mass of the system, units = same as the other mass values given, default = None
        mu = reduced mass of the system, units = same as the other mass values given, default = None
        tolarance = what the difference has to be less than or equal to in order to not replace
        the given value with that calculated, units = NA, default = 1e5
    
    Outputs:
        m2 = mass of first binary object, units = same as the units of the masses passed to function
    '''
    
    # Try to calculate m2 from whatever given mass values you have using the find_m1_m2 function
    try:
        m2_calced = find_m1_m2(m1, m2, Mtot, q, Mc, mu)[1]
    # If the find_m1_m2 function fails, then raise an error becuase none of the mass values can be calculated.
    except:
        raise ValueError("Not enough mass values available to calculate m1 and m2.")

    # If no m2 was passed to the function or the difference between the value passed and
    # the value calculated is greater than the tolerance, the calculated m2 value is the
    # the mass value used for m2.
    if m2 == None or abs(m2 - m2_calced) >= tolerance:
        m2 = m2_calced

    # If none of the above criteria are met then the m2 value passed to the function
    # stays as the m2 value.
    else:
        pass
    return m2
##############


## update_Mc.py has the code for the update_Mc function. This function will check to see if 
## the chirp mass value passed to the function is the same or within some tolerance to the 
## chirp mass value calculated given the m1 and m2. If there is not chirp mass value passed
## into the function, or the tolerance is exceeded, then the calculated chirp mass value
## becomes the value used and returned as the chirp mass.

###############
def update_Mc(m1, m2, Mc = None, tolerance = 1e5):
    '''
    Update the chirp mass value given to that calculated if the tolerance is not met.

    Inputs:
        m1 = mass of first binary object, units = same as the other mass value given
        m2 = mass of second binary object, units = same as the other mass value given
        Mc = chirp mass of the system, units = same as the other mass values given, default = None
        tolarance = what the difference has to be less than or equal to in order to not replace
        the given value with that calculated, units = NA, default = 1e5
    
    Outputs:
        Mc = chirp mass of the system, units = same as the units of the masses passed to function
    '''
    
    # Calculate the chirp mass from the given m1 and m2
    Mc_calced = Mc_calc(m1, m2)

    # If no chirp mass was passed to the function or the difference between the value passed and
    # the value calculated is greater than the tolerance, the calculated chirp mass value is the
    # the mass value used for chirp mass.
    if Mc == None or abs(Mc - Mc_calced) >= tolerance:
        Mc = Mc_calced

    # If none of the above criteria are met then the chirp mass value passed to the function
    # stays as the chirp mass value.
    else:
        pass
    return Mc
##############

## update_Mtot.py has the code for the update_Mtot function. This function will check to see if 
## the total mass value passed to the function is the same or within some tolerance to the 
## total mass value calculated given the m1 and m2. If there is not total mass value passed
## into the function, or the tolerance is exceeded, then the calculated total mass value
## becomes the value used and returned as the total mass.


###############
def update_Mtot(m1, m2, Mtot = None, tolerance = 1e5):
    '''
    Update the total mass value given to that calculated if the tolerance is not met.

    Inputs:
        m1 = mass of first binary object, units = same as the other mass value given
        m2 = mass of second binary object, units = same as the other mass value given
        Mtot = total mass of the system, units = same as the other mass values given, default = None
        tolarance = what the difference has to be less than or equal to in order to not replace
        the given value with that calculated, units = NA, default = 1e5
    
    Outputs:
        Mtot = total mass of the system, units = same as the units of the masses passed to function
    '''
    
    # Calculate the total mass from the given m1 and m2
    Mtot_calced = Mtot_calc(m1, m2)

    # If no total mass was passed to the function or the difference between the value passed and
    # the value calculated is greater than the tolerance, the calculated total mass value is the
    # the mass value used for total mass.
    if Mtot == None or abs(Mtot - Mtot_calced) >= tolerance:
        Mtot = Mtot_calced

    # If none of the above criteria are met then the total mass value passed to the function
    # stays as the total mass value.
    else:
        pass
    return Mtot
##############


## update_mu.py has the code for the update_mu function. This function will check to see if 
## the reduced mass value passed to the function is the same or within some tolerance to the 
## reduced mass value calculated given the m1 and m2. If there is not reduced mass value passed
## into the function, or the tolerance is exceeded, then the calculated reduced mass value
## becomes the value used and returned as the reduced mass.

###############
def update_mu(m1, m2, mu = None, tolerance = 1e5):
    '''
    Update the reduced mass value given to that calculated if the tolerance is not met.

    Inputs:
        m1 = mass of first binary object, units = same as the other mass value given
        m2 = mass of second binary object, units = same as the other mass value given
        mu = reduced mass of the system, units = same as the other mass values given, default = None
        tolarance = what the difference has to be less than or equal to in order to not replace
        the given value with that calculated, units = NA, default = 1e5
    
    Outputs:
        mu = reduced mass of the system, units = same as the units of the masses passed to function
    '''
    
    # Calculate the reduced mass from the given m1 and m2
    mu_calced = mu_calc(m1, m2)

    # If no reduced mass was passed to the function or the difference between the value passed and
    # the value calculated is greater than the tolerance, the calculated reduced mass value is the
    # the mass value used for reduced mass.
    if mu == None or abs(mu - mu_calced) >= tolerance:
        mu = mu_calced

    # If none of the above criteria are met then the reduced mass value passed to the function
    # stays as the reduced mass value.
    else:
        pass

    return mu
##############


## update_q.py has the code for the update_q function. This function will check to see if 
## the mass ratio value passed to the function is the same or within some tolerance to the 
## mass ratio value calculated given the m1 and m2. If there is not a mass ratio value passed
## into the function, or the tolerance is exceeded, then the calculated mass ratio value
## becomes the value used and returned as the mass ratio.


###############
def update_q(m1, m2, q = None, tolerance = 0.01):
    '''
    Update the mass ratio value given to that calculated if the tolerance is not met.

    Inputs:
        m1 = mass of first binary object, units = same as the other mass value given
        m2 = mass of second binary object, units = same as the other mass value given
        q = mass ratio of the system, units = NA, default = None
        tolarance = what the difference has to be less than or equal to in order to not replace
        the given value with that calculated, units = NA, default = 0.01
    
    Outputs:
        q = mass ratio of the system, units = NA, default = None
    '''
    
    # Calculate the mass ratio from the given m1 and m2
    q_calced = q_calc(m1, m2)

    # If no mass ratio was passed to the function or the difference between the value passed and
    # the value calculated is greater than the tolerance, the calculated mass ratio value is the
    # the value used for mass ratio.
    if q == None or abs(q - q_calced) >= tolerance:
        q = q_calced

    # If none of the above criteria are met then the mass ratio value passed to the function
    # stays as the mass ratio value.
    else:
        pass

    return q
##############



## update_all_masses.py has the code for the function update_masses(). This function
## complies all of the update_* functions into one so that the user can fully update
## all six mass values easily in one command. There is one caveat to this function,
## the tolerance used in all of the update_* functions must be the default tolerance.
## The user cannot define the tolerances themselves for the individual update_* functions.


################
def update_masses(m1, m2, Mtot, q, Mc, mu):
    '''
    Update all six of the mass values used to fully describe a binary system to make sure the arugments
    passed to it are all compatibable with each other according to the equations. NOTE: This function 
    ONLY uses the default tolerances in all the individual update functions. There is no way to set
    user defined tolerances.

    Inputs:
        m1 = mass of first binary object, units = same as the other mass value given, default = None
        m2 = mass of second binary object, units = same as the other mass value given, default = None
        Mtot = total mass of the system, units = same as the other mass values given, default = None
        q = mass ratio of the system, units = NA, default = None
        Mc = chirp mass of the system, units = same as the other mass values given, default = None
        mu = reduced mass of the system, units = same as the other mass values given, default = None
    
    Outputs:
        m1 = mass of first binary object, units = same as the other mass value given
        m2 = mass of second binary object, units = same as the other mass value given
        Mtot = total mass of the system, units = same as the other mass values given
        q = mass ratio of the system, units = NA
        Mc = chirp mass of the system, units = same as the other mass values given
        mu = reduced mass of the system, units = same as the other mass values given
    
    '''
   
    # Update m1.
    m1_updated = update_m1(m1, m2, Mtot, q, Mc, mu)

    # Update m2.
    m2_updated = update_m2(m1, m2, Mtot, q, Mc, mu)

    # Update Mtot.
    Mtot_updated = update_Mtot(m1, m2, Mtot)

    # Update q.
    q_updated = update_q(m1, m2, q)

    # Update Mc.
    Mc_updated = update_Mc(m1, m2, Mc)

    # Update mu.
    mu_updated = update_mu(m1, m2, mu)

    # In order to know whether or not any of the mass values were updated we need to keep track
    # of the mass values that actually are updated from the value initially given to the function.
    # This is done by appending the names of the variables that were updated to an array that
    # will be returned along with the updated mass values.

    # First create the empty array.
    updated_masses = []

    # Next check if any of the mass values have been updated and if so append the mass name to the 
    # array.
    if m1_updated != m1:
        updated_masses.append("m1")
    if m2_updated != m2:
        updated_masses.append("m2")
    if Mtot_updated != Mtot:
        updated_masses.append("Mtot")
    if q_updated != q:
        updated_masses.append("q")
    if Mc_updated != Mc:
        updated_masses.append("Mc")
    if mu_updated != mu:
        updated_masses.append("mu")

    # Return all the updated mass values and the array with the names of the values that have been updated.
    return ([(m1_updated,m2_updated,Mtot_updated,q_updated,Mc_updated,mu_updated), updated_masses])
################