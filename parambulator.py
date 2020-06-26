# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 08:31:49 2020

@author: gsied
"""

def propellant_calc(r0x,r0y,r0z,v0x,v0y,v0z,Mo,Isp):
    """
    Parameters
    ----------
    r0x : Positive Number
        Initial position for i direction in km.
    r0y : Positive Number
        Initial position for j direction in km.
    r0z : Positive Number
        Initial position for k direction in km.
    v0x : Positive Number
        Initial velocity for i direction in km/s.
    v0y : Positive Number
        Initial velocity for j direction in km/s.
    v0z : Positive Number
        Initial velocity for k direction in km/s.
    Mp : Positive Number
        Maximum propellant mass in kg.
    Mo : Positive Number
        Mass of spacecraft in kg.
    Isp : Positive Number
        Specific impulse of engine in s.
    variable: String
        String of either m, t, or v to select desired variable to optimize

    Returns
    -------
    delta_total : array
        n arrays of length n giving results of specified variable to travel (m, t,or v)
        from satellite number = row number to satellite number = column number 

    """
    import numpy as np
    import math
    
    #read in initial conditions
    mu = 3.986e5 # Gm of earth in km^3/s^2
    r0x = r0x
    r0y = r0y
    r0z = r0z
    v0x = v0x
    v0y = v0y
    v0z = v0z
    size = (len(r0x),len(r0x))
    delta_total = np.zeros(size)
    
    #create matracies of initial conditions, each row is a satellite
    r0 = np.column_stack((r0x,r0y,r0z))
    v0 = np.column_stack((v0x,v0y,v0z))
    
    for point1 in range(len(r0)):
        for point2 in range((point1+1),len(r0)):
            
            #create empty vectors
            h0 = []
            r0_mag = []
            v0_mag = []
            orbit_params = []
            
            #create vectors of initial conditions for comparison
            r01 = r0[point1]
            v01 = v0[point1]
            r02 = r0[point2]
            v02 = v0[point2]
            
            #calculate angular momentum vectors
            h1 = np.cross(r01,v01)
            h0.append(h1)
            h2 = np.cross(r02,v02)
            h0.append(h2)
            
            #calculate magnitudes of vectors in km and km/s
            r0_mag.append(np.linalg.norm(r01))
            v0_mag.append(np.linalg.norm(v01))
            r0_mag.append(np.linalg.norm(r02))
            v0_mag.append(np.linalg.norm(v02))
            
            #calculate perigee, semi major axis, apogee (all in km),
            #and inclincation (radians) of orbits and pass them to orbit_params 
            for i in range(2): 
                h0_unit = h0[i]/np.linalg.norm(h0[i])
                x0 = r0_mag[i]*(v0_mag[i])**2/mu #dimensionless orbit variable
                beta0 = np.arccos(np.linalg.norm(h0[i])/(r0_mag[i]*v0_mag[i]))
                e = math.sqrt((x0-1)**2*(math.cos(beta0))**2\
                                + (math.sin(beta0))**2) #orbit eccentricity
                rp = np.linalg.norm(h0[i])**2/mu/(1+e) #periapse distance
                a = rp/(1-e) #semimajor axis distance
                ra = a*(1+e) #apoapse distance
                inc = math.acos(h0_unit[-1]) #orbit inclination
                if inc > math.pi/2:
                    inc = inc - math.pi # if inclincation is over 90 degrees use negative equivalent
                add = [rp,ra,inc,e,beta0]
                orbit_params.append(add)
            
            #Select starting orbit
            #initial assumtions: 1) no phasing needed
            #2)compare apogee values to reduce plane change deltav
            #3)do combined manuver from largest apogee to other perigee
            if orbit_params[0][1] > orbit_params[1][1]:
                rp = orbit_params[1][0]
                ra = orbit_params[0][1]
                a1 = (ra+orbit_params[0][0])/2
                a2 = (orbit_params[1][1]+rp)/2
                v2 = math.sqrt(mu*(2/rp-1/a2))
                v1 = math.sqrt(mu*(2/ra-1/a1))
                h = np.linalg.norm(h0[0])
                e = orbit_params[0][3]
                r = r0_mag[0]
                beta = orbit_params[0][4]
            else:
                rp = orbit_params[0][0]
                ra = orbit_params[1][1]
                a1 = (ra+orbit_params[1][0])/2
                a2 = (orbit_params[0][1]+rp)/2
                v1 = math.sqrt(mu*(2/ra-1/(a1)))
                v2 = math.sqrt(mu*(2/rp-1/(a2)))
                h = np.linalg.norm(h0[1])
                e = orbit_params[1][3]
                r = r0_mag[1]
                beta = orbit_params[1][4]
            
            #Calculate velocity of transfer orbit starting at apogee of transfer orbit
            a_t = (rp + ra)/2
            vt1 = math.sqrt(mu*(2/ra-1/a_t))
            delta_v1_in_plane = abs(v1-vt1) #km/s
            
            # if difference in inclination is greater than 5 degrees (value chosen at random)
            # assume a plane change is needed and include it in the calculations
            # use orbit with larger ra to minimize the needed delta v
            delta_i = abs(orbit_params[0][2]-orbit_params[1][2])
            if np.degrees(delta_i) > 5:
                delta_v_out_of_plane = h/ra*delta_i   #necessary inclination change    
            else:
                delta_v_out_of_plane = 0
            
            #calculate delta v needed for combined manuver (transfer orbits and plane change)
            delta_v_combined = math.sqrt((delta_v1_in_plane)**2\
                                         +(delta_v_out_of_plane)**2)
            
            #velocity at perigee of transfer orbit and required delta v to enter new orbit
            vt2 = vt1*ra/rp
            delta_v2 = abs(v2-vt2)
            
            #Total propellant needed for first transfer    
            delta_m1 = Mo*(1-math.exp(-delta_v_combined*1000/(9.81*Isp)))# same units as Mo
            
            #new mass after first burn
            Mo1 = Mo - delta_m1
           
            #Propellant needed to match second orbit
            delta_m2 = Mo1*(1-math.exp(-delta_v2*1000/(9.81*Isp)))# same units as Mo
            
            #new mass after burn 2
            Mo2 = Mo1 - delta_m2
            
            #Total time for travel from starting location of orbit 1 to apogee of orbit 1 +
            #transfer orbit time
            tau1 = 2*math.pi*math.sqrt(a1**3/mu) #period of orbit 1
            E1 = math.acos((a1-r)/(a1*e))
            
            #if beta is greater than 90 degrees it is equivalent to being negative
            #if beta is negative then the satellite is past apogee and affects
            #eccentric anamoly, E
            if beta > math.pi/2:
                E1 = 2*math.pi - E1
                
            tp1 = math.sqrt(a1**3/mu)*(E1 - e*math.sin(E1)) #time since periapse of orbit 1 in seconds
            time_to_a1 = tau1/2 - tp1 #time needed to reach apogee of orbit 1 in seconds
            
            tau_t = math.pi*math.sqrt(a_t**3/mu) #half period of transfer orbit
            delta_t = time_to_a1 + tau_t
            
            #assuming robot enters new orbit at perigee, do required phasing
            #to meet up with second satellite
            delta_v_phasing, phasing_TOF = phasing(a2, delta_t, E1, rp)
            
            #Propellant needed for phasing
            delta_m3 = Mo2*(1-math.exp(-delta_v_phasing*1000/(9.81*Isp)))# same units as Mo
            
            #populate return matrix based on desired variable choice
            delta_total[point1][point2] = delta_m1 + delta_m2 + delta_m3
            
            #graph is undirected so populate the columns with the rows
            delta_total[point2][point1] = delta_total[point1][point2]
    
    return delta_total

def phasing(a,delta_t, E, r):
    '''
    Parameters
    ----------
    a : Positive Float
        Semimajor axis of orbit to propagate
    delta_t : Positive Float
        Change of time in seconds to propogate orbit
    E : Float
        Eccentric anamoly of initial position to be moved forward in time
    r : Float
        Magnitude of position vector used to calculate the delta V needed to 
        enter phasing orbit

    Returns
    -------
    delta_v_phasing : Positive Float
        Necessary delta v in km/s for robot to phase within the orbit it is 
        currently in while taking the least amount of time
        
    TOF : Positive Float
        Time in seconds needed for phasing

    '''
    
    import math

    mu = 3.986e5 # Gm of earth in km^3/s^2
    
    omega = math.sqrt(mu/a**3)
    E2 = (E + omega*delta_t)%(2*math.pi)
    
    if E2 > math.pi:
        E2 = - E2
    a_c = (mu*((2*math.pi-E2)/(2*math.pi*omega))**2)**(1/3) #semimajor axis of phasing orbit
    t_phasing = (2*math.pi-E2)/omega #time of flight to achieve phasing in seconds
    
    #delta v needed to enter phasing orbit and return to original orbit
    delta_v_phasing = 2*abs(math.sqrt(mu*(2/r-1/a)) - math.sqrt(mu*(2/r-1/a_c))) #km/s
    
    return delta_v_phasing, t_phasing
