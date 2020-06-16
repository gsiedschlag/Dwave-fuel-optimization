# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 08:31:49 2020

@author: gsied
"""

def parambulator(r0x,r0y,r0z,v0x,v0y,v0z,Mp,Mo,Isp):
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

    Returns
    -------
    delta_m_total : list
        fuel used to travel between satellites

    """
    import numpy as np
    import math

    mu = 3.986e5 # Gm of earth in km^3/s^2
    r0x = r0x
    r0y = r0y
    r0z = r0z
    v0x = v0x
    v0y = v0y
    v0z = v0z
    size = (len(r0x),len(r0x))
    delta_m_total = np.zeros(size)
    
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
            
            #calculate periapse, semi major axis, and
            #apoapse of orbits in km
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
                add = [rp,ra,inc]
                orbit_params.append(add)
            #calculate delta v/m to move between orbits
            #initial assumtions: 1) no phasing needed
            #2)wait until apoapse of larger orbit or periapse of smaller orbit
            #3)travel to other satellites apoapse/periapse to minimize needed delta v
            if abs(orbit_params[0][0] - orbit_params[1][1]) <\
                abs(orbit_params[1][0] - orbit_params[0][1]):
                rp = orbit_params[0][0]
                ra = orbit_params[1][1]
                a1 = (rp+orbit_params[0][1])/2
                a2 = (orbit_params[1][0]+ra)/2
                v2 = math.sqrt(mu*(2/rp-1/a1))
                v1 = math.sqrt(mu*(2/ra-1/a2))
                h = np.linalg.norm(h0[0])
            else:
                rp = orbit_params[1][0]
                ra = orbit_params[0][1]
                a1 = (rp+orbit_params[1][1])/2
                a2 = (orbit_params[0][0]+ra)/2
                v1 = math.sqrt(mu*(2/rp-1/(a1)))
                v2 = math.sqrt(mu*(2/orbit_params[0][0]-1/(a2)))
                h = np.linalg.norm(h0[1])
                
            a_t = (rp + ra)/2
            vt1 = math.sqrt(mu*(2/rp-1/a_t))
            delta_v1 = abs(v1-vt1) #km/s
            delta_m1 = Mo*(1-math.exp(-delta_v1*1000/(9.81*Isp)))# same units as Mo
            Mo1 = Mo - delta_m1
            vt2 = vt1*rp/ra
            delta_v2_in_plane = abs(v2-vt2)
            # if difference in inclination is greater than 5 degrees
            # assume a plane change is needed and include it in the calculations
            # use orbit with larger ra to minimize the needed delta v
            delta_i = abs(orbit_params[0][2]-orbit_params[1][2])
            if np.degrees(delta_i) > 5:
                delta_v_out_of_plane = h/ra*delta_i   #stopping here on 6/15    
            else:
                delta_v_out_of_plane = 0
            delta_v2 = math.sqrt(delta_v2_in_plane**2 + delta_v_out_of_plane**2)
            delta_m2 = Mo1*(1-math.exp(-delta_v2*1000/(9.81*Isp)))# same units as Mo
            
            delta_m_total[point1][point2] = delta_m1 + delta_m2
            delta_m_total[point2][point1] = delta_m_total[point1][point2]
    return delta_m_total