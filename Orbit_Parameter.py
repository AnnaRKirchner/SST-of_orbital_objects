"""
Created on Mon Jan 30 17:39:59 2023
@author: Anna Kirchner

Description: This code contains a collection of methods for calculation of the orbital parameters 
            related to observational data.
            It includes a method for calculating the Earth's radius at a specific latitude of an observer.
            It computes weighted means for several seperate calculations of orbital parameters 
            for the same object.
            
            Uncertainty analyses (propagated) is performed for all data, 
            including determined orbital state vectors and fundamental vectors for orbit calculation: (h, n, e)

Methods in this file:

    1.get_earth_radius: Calculates the Earth's radius at a given observer location, considering the oblateness of the Earth
    
    2.weighted_mean: Computes the weighted mean of a set of values, internal and external uncertainty
    
    3.get_u_r_vector & 4.get_u_v_vector : Determines the uncertainty of orbital state vectors (r,v)

    5.get_hne_error: Computes uncertainties of fundamental vectors (h, n, e) used in orbit calculations
             
    6.get_orbital_params: Computes the Orbital parameters (for now: circular orbit (eccentricity==0), onedimensional),
                          calculated for each track middle point and weighted mean of all image tracks of the same object
"""

import numpy as np 
import math
import sympy as sp
import astropy.constants as const

"""
    Calculates Earth radius at given observer latitude
    
    Parameters:
        - obs (EarthLocation): Earth location coordinates (longitude, latitude, altitude) in degrees and meters. 
                               Default is (0, 0, 0), which represents the Greenwich Observatory.
        - u_coordinate: uncertainty of coordinate transformation
    
    Returns:
        - R_E, u_r_value: Earth radius and propagated uncertainty
"""
def get_earth_radius(obs, u_coordinate):
   
   u_lat=u_coordinate
   
   # equatorial Earth radius
   r_eq = 6378136.6 #[m]
   u_r_eq=0.1 #[m]
   # polar radius
   r_pole = 6356752.3142  #wert von geopy # nicht genügend RAM, um geopy zu installieren
   u_r_pol=0.1 #.[m]
   
   #actual radius at latitude
   R_E= np.sqrt(((r_eq**2 * np.cos(obs.lat.radian))**2 + 
                 (r_pole**2 * np.sin(obs.lat.radian))**2) / ((r_eq * np.cos(obs.lat.radian))**2 +
                                                             (r_pole * np.sin(obs.lat.radian))**2))
   
   'u(r)'
   #variables
   r_pol, lat, R_eq =sp.symbols('r_pol lat r_eq')   
   
   #function
   r_E= sp.sqrt(((R_eq**2 * sp.cos(lat))**2 + 
                 (r_pol**2 * sp.sin(lat))**2) / ((R_eq * sp.cos(lat))**2 +
                                                             (r_pol * sp.sin(lat))**2))
                                                
   #calc derivatives
   dr_dr_pol=sp.diff(r_E, r_pol)
   dr_dr_eq=sp.diff(r_E, R_eq)
   dr_dlat=sp.diff(r_E, lat)
   
   # Step 4: Substitute numerical values into the derivatives
   dr_dr_pol_val = dr_dr_pol.subs([(r_pol, r_pole), (lat, obs.lat.radian), (R_eq, r_eq)])
   dr_dr_eq_val=dr_dr_eq.subs([(r_pol, r_pole), (lat, obs.lat.radian), (R_eq, r_eq)])
   dr_dlat_val=dr_dlat.subs([(r_pol, r_pole), (lat, obs.lat.radian), (R_eq, r_eq)])

 
   u_r= sp.sqrt((dr_dr_pol_val*u_r_pol)**2+(dr_dr_eq_val*u_r_eq)**2+(dr_dlat_val*u_lat)**2)
   u_r_value=u_r.evalf()
   
   return  R_E, u_r_value

"""
    Calculates weighted mean, internal and external uncertainty
    
    Parameters:
        - values, uncertainties: Values with uncertainties to get weighted mean from
        
    Returns:
        - mean, u_internal, u_external
"""
def weighted_mean(values, uncertainties):
    
        
    # Calculate the weights as the inverse square of the uncertainties
    values= np.array(values)
    uncertainties=np.array(uncertainties)
    
    weights=1/ (uncertainties ** 2) 
        
    # Calculate the weighted mean
    mean = np.sum(values * weights) / np.sum(weights)
 
    # Calculate the uncertainties of weighted mean
    u_internal = np.sqrt(1/ np.sum(weights))

    u_external= np.sqrt(np.sum(weights*(values-mean)**2/((len(values)-1)*np.sum(weights))))
    
    return mean, u_internal, u_external


"""
    Calculates uncertainty of state vectors r in ECI coordinates
    
    Parameters:
        - r, u_r: orbit radius (height plus distance from Earth's radius) + uncertainty
        - Ra, Dec: Ra/Dec coordinates of position
        - u_coordinate:  uncertainty of coordinate transformation
    
    Returns:
        - u_r_vector: propagated uncertainty
"""
def get_u_r_vector(r, u_r, Ra, Dec, u_coordinate):
    r_s, a, d=sp.symbols('r, a, d')
    
    x= r_s * sp.cos(d)*sp.cos(a)
    y= r_s * sp.cos(d)*sp.sin(a)
    z= r_s * sp.sin(d)
    
    # Calculate the partial derivatives
    dx_dr = sp.diff(x, r_s)  
    dx_dra1 = sp.diff(x, d) 
    dx_dra0 = sp.diff(x, a)
    
    dy_dr = sp.diff(y, r_s)  
    dy_dra1 = sp.diff(y, d) 
    dy_dra0 = sp.diff(y, a)
    
    dz_dr=sp.diff(z, r_s)  
    dz_dra1= sp.diff(z, d) 
        
          
    #for r1 or r2 vector
    dx_dr_r = dx_dr.subs([(r_s, r), (d, Dec), (a, Ra)])
    dx_dra1_r = dx_dra1.subs([(r_s, r), (d, Dec), (a, Ra)])
    dx_dra0_r = dx_dra0.subs([(r_s, r), (d, Dec), (a, Ra)])
    
    dy_dr_r = dy_dr.subs([(r_s, r), (d, Dec), (a, Ra)])
    dy_dra1_r = dy_dra1.subs([(r_s, r), (d, Dec), (a, Ra)])
    dy_dra0_r = dy_dra0.subs([(r_s, r), (d, Dec), (a, Ra)])
    
    dz_dr_r = dz_dr.subs([(r_s, r), (d, Dec), (a, Ra)])
    dz_dra1_r = dz_dra1.subs([(r_s, r), (d, Dec), (a, Ra)])
    
    u_rx= sp.sqrt((dx_dr_r*u_r)**2+(dx_dra1_r*u_coordinate)**2+(dx_dra0_r*u_coordinate)**2)
    u_ry= sp.sqrt((dy_dr_r*u_r)**2+(dy_dra1_r*u_coordinate)**2+(dy_dra0_r*u_coordinate)**2)  
    u_rz= sp.sqrt((dz_dr_r*u_r)**2+(dz_dra1_r*u_coordinate)**2)
      
    u_rx= u_rx.evalf()
    u_ry= u_ry.evalf()
    u_rz= u_rz.evalf()
        
    u_r_vector=np.array([u_rx, u_ry, u_rz])
   
    return u_r_vector


"""
    Calculates uncertainty of state vector v in ECI coordinates
    
    Parameters:
        - displacement, u_disp: displacement between two positional vectors r1, r2 (ECI) and uncertainty
        - velocity, u_velocity: magnitude and uncertainty
    Returns:
        - u_v_vector: propagated uncertainty
"""
def get_u_v_vector(displacement, u_disp, velocity, u_velocity):
    
    #unit vector of direction times velocity
    #v_vec=displacement/disp_normed * velocity
    
    #components of displacement vector
    d0, d1, d2, v = sp.symbols('d0 d1 d2 v')
    
    v_vec=sp.Matrix([d0, d1, d2])/sp.sqrt(d0**2+d1**2+d2**2) *v
    

    dvec_dd0=sp.diff(v_vec, d0)
    dvec_dd1=sp.diff(v_vec, d1)
    dvec_dd2=sp.diff(v_vec, d2)
    dvec_dv=sp.diff(v_vec, v)
    
    substitutions = [(v, velocity),(d0, displacement[0]), (d1, displacement[1]), 
                     (d2, displacement[2])]
    
    dvec_dd0_val= dvec_dd0.subs(substitutions).evalf()
    dvec_dd1_val= dvec_dd1.subs(substitutions).evalf()
    dvec_dd2_val= dvec_dd2.subs(substitutions).evalf()
    dvec_dv_val= dvec_dv.subs(substitutions).evalf()

    squared_terms = np.square(np.multiply(dvec_dd0_val,u_disp[0])) + \
                    np.square(np.multiply(dvec_dd1_val,u_disp[1])) + \
                    np.square(np.multiply(dvec_dd2_val,u_disp[2])) + \
                    np.square(np.multiply(dvec_dv_val,u_velocity)) 
                    
    
    u_v_vec=np.sqrt(np.float64(squared_terms)).flatten()
    return u_v_vec

"""
    Computes  uncertainty of fundamental vectors h, n, e for orbit calculation'
    
    Parameters:
        - v_vector, u_v_vector: The velocity vector in ECI (Earth centered Inertial) + propagated uncertainty
        - r_vector, u_r_vector: The position vector in ECI  + propagated uncertainty
    
    Outputs:
        - h_vector, u_h: h_vector in ECI and uncertainty (propagated) of the h vector
        - n_vector, u_n: n_vector in ECI and uncertainty (propagated) of the h vector
        - e_vector, u_e: e_vector in ECI and uncertainty (propagated) of the h vector     
"""
def get_hne_error(v_vector,r_vector, u_v_vector, u_r_vector):
    
    global G
    global m_E
    global u_mu
    global mu
    
    north_vector= [0,0,1] #in ECI coordinates
    
    h_vector =np.cross(v_vector, r_vector)
    
    e_vector= (1/mu) * ((np.linalg.norm(v_vector)**2-(mu/np.linalg.norm(r_vector)))*r_vector-
                        (np.dot(v_vector, r_vector))*v_vector)

    n_vector= np.cross(h_vector, north_vector)
    
    'u(h)' 
    v1, v2, v3, r1, r2, r3= sp.symbols('v1 v2 v3 r1 r2 r3')
    
    h_vec=sp.Matrix([v2*r3-v3*r2, v3*r1-v1*r3, v1*r2-v2*r1])
    
    dhvec_dv1=sp.diff(h_vec, v1)
    dhvec_dv2=sp.diff(h_vec, v2)
    dhvec_dv3=sp.diff(h_vec, v3)
    dhvec_dr1=sp.diff(h_vec, r1)
    dhvec_dr2=sp.diff(h_vec, r2)
    dhvec_dr3=sp.diff(h_vec, r3)
    
    substitutions = [(v1, v_vector[0]),(v2, v_vector[1]), (v3, v_vector[2]), 
                     (r1, r_vector[0]), (r2, r_vector[1]), (r3, r_vector[2])]
    
    dhvec_dv1_val= dhvec_dv1.subs(substitutions).evalf()
    dhvec_dv2_val= dhvec_dv2.subs(substitutions).evalf()
    dhvec_dv3_val= dhvec_dv3.subs(substitutions).evalf()
    dhvec_dr1_val= dhvec_dr1.subs(substitutions).evalf()
    dhvec_dr2_val= dhvec_dr2.subs(substitutions).evalf()
    dhvec_dr3_val= dhvec_dr3.subs(substitutions).evalf()


    squared_terms = np.square(np.multiply(dhvec_dv1_val,u_v_vector[0])) + \
                    np.square(np.multiply(dhvec_dv2_val,u_v_vector[1])) + \
                    np.square(np.multiply(dhvec_dv3_val,u_v_vector[2])) + \
                    np.square(np.multiply(dhvec_dr1_val,u_r_vector[0])) + \
                    np.square(np.multiply(dhvec_dr2_val,u_r_vector[1])) + \
                    np.square(np.multiply(dhvec_dr3_val,u_r_vector[2]))
                    
    u_h= np.sqrt(np.float64(squared_terms)).flatten()
    
    'u(n)'
    h1, h2, h3, k1, k2, k3= sp.symbols('h1 h2 h3 k1 k2 k3')
    
    n=sp.Matrix([h2*k3-h3*k2, h3*k1-h1*k3, h1*k2-h2*k1])
    
    dn_dh1=sp.diff(n, h1)
    dn_dh2=sp.diff(n, h2)
    dn_dh3=sp.diff(n, h3)
    dn_dk1=sp.diff(n, k1)
    dn_dk2=sp.diff(n, k2)
    dn_dk3=sp.diff(n, k3)
    
    substitutions = [(h1, h_vector[0]),(h2, h_vector[1]), (h3, h_vector[2]), (k1, north_vector[0]), (k2, north_vector[1]), (k3, north_vector[2])]
    dn_dh1_val= dn_dh1.subs(substitutions).evalf()
    dn_dh2_val=dn_dh2.subs(substitutions).evalf()
    dn_dh3_val = dn_dh3.subs(substitutions).evalf()
    dn_dk1_val=dn_dk1.subs(substitutions).evalf()
    dn_dk2_val=dn_dk2.subs(substitutions).evalf()
    dn_dk3_val=dn_dk3.subs(substitutions).evalf()
    
    u_k=[0, 0, 0]
    squared_terms = np.square(np.multiply(dn_dh1_val,u_h[0])) + \
                    np.square(np.multiply(dn_dh2_val,u_h[1])) + \
                    np.square(np.multiply(dn_dh3_val,u_h[2])) + \
                    np.square(np.multiply(dn_dk1_val,u_k[0])) + \
                    np.square(np.multiply(dn_dk2_val,u_k[1])) + \
                    np.square(np.multiply(dn_dk3_val,u_k[2]))
                    
    u_n= np.sqrt(np.float64(squared_terms)).flatten()
   
     
    'u(e)'  
    m, v0, v1, v2, r0, r1, r2 = sp.symbols('mu v0 v1 v2 r0 r1 r2')
    e = (1/m) * (((v0**2+v1**2+v2**2) - (m/sp.sqrt(r0**2+r1**2+r2*+2)))*sp.Matrix([r0, r1, r2]).T - (v0*r0 + v1*r1 + v2*r2)*sp.Matrix([v0, v1, v2]).T)
    
    
    # Calculate the partial derivatives
    
    de_dm=sp.diff(e, m)
    
    # Berechne die partiellen Ableitungen
    de_dr = [sp.diff(e, var) for var in [r0, r1, r2]]
    de_dv= [sp.diff(e, var2) for var2 in [v0, v1, v2]]
    
    substitutions = [(m, G*m_E),
                    (r0, r_vector[0]),(r1, r_vector[1]), (r2, r_vector[2]), 
                    (v0, v_vector[0]), (v1, v_vector[1]), (v2, v_vector[2])]
    de_dmu_val = de_dm.subs(substitutions).evalf()
    de_dv0_val = de_dv[0].subs(substitutions).evalf()
    de_dr0_val = de_dr[0].subs(substitutions).evalf()
    de_dv1_val = de_dv[1].subs(substitutions).evalf()
    de_dr1_val = de_dr[1].subs(substitutions).evalf()
    de_dv2_val = de_dv[2].subs(substitutions).evalf()
    de_dr2_val = de_dr[2].subs(substitutions).evalf()

    # Calculate the squared terms
    squared_terms = np.square(np.multiply(de_dmu_val, u_mu)) + \
                    np.square(np.multiply(de_dv0_val, u_v_vector[0])) + \
                    np.square(np.multiply(de_dv1_val, u_v_vector[1])) + \
                    np.square(np.multiply(de_dv2_val, u_v_vector[2])) + \
                    np.square(np.multiply(de_dr0_val, u_r_vector[0])) + \
                    np.square(np.multiply(de_dr1_val, u_r_vector[1])) + \
                    np.square(np.multiply(de_dr2_val, u_r_vector[2]))
    
    # Calculate the square root element-wise
    u_e = np.sqrt(np.float64(squared_terms))
    u_e=u_e[0]
  
    return h_vector, u_h, n_vector, u_n, e_vector, u_e


"""
    Calculates orbital parameters for section from one track-middle point to another:
        
    Parameters:
        - observer (astropy.coordinates.EarthLocation): Earth location of the observer
        - sat_skycoord: Ra/Dec sky-coordinates list of a tracks start- and endpoints
        - sat_heights, u_height:  list of object heights above Earth's surface + uncertainty list
        - sat_velocity, u_vel: list of velocity magnitudes + uncertainty list
        - exposure_time: exposure time of one image
        - u_coordinate: uncertainty of coordinate transformation
    
    Returns:
        - mean_orbital: String with computed mean orbital parameters for one image series of the same object
        --> String contains:
                -orbital period T
                -mean motion n
                -inclination
                -RAAN
                -Argument of periapsis
                -true anomaly (not defined for e=0)
                -argument of latitude u (instead of true anomaly for circular orbit)
                -true longitude (for orbits with inclination approx. zero)
"""
def get_orbital_params(observer, sat_skycoord, sat_heights, u_height, sat_velocity, u_vel, exposure_time, u_coordinate):
    #observer location
    obs = observer
    
    #constants
    global G, u_G
    global m_E, u_m
    global r_E, u_r_E
    global mu, u_mu
    
    G=const.G.value  #± 0,000 85 [N*m2 / kg2]
    u_G= 6.7*10**(-15)
    
    m_E=const.M_earth.value   #[kg]
    u_m= 6*10**20 #[kg]
    
    mu= G*m_E
    u_mu=np.sqrt((m_E*u_G)**2+(G*u_m)**2)
    
    #specific radius at observer location (for Earth being a spheroid)  
    r_E, u_r_E= get_earth_radius(obs, u_coordinate)
    
    #eccentricity for circular orbit
    e=0 
    
    h = sat_heights # [m] from Orbit_Data.py
    v=sat_velocity  #[m/s] 
    
    #distance to Earth's center specified for observer location
    r = h + r_E 
    
    T=(2*np.pi*r)/v  #one period around earth
    
    n = 86400/ T  #mean motion
    
    #create lists to get the mean values for satellite and their uncertainties in follow up images
    (period,u_period, incl, u_incl, long_asc, u_long_asc,  arg_peri, u_arg_peri,true_anom, u_true_anom,
    arg_lat, u_arg_lat, tru_long,u_true_long, motion, u_mean_motion) =([] for i in range(16))

    
    #main loop gets orbital parameters for each image track
    k=0
    for i,point in enumerate(sat_skycoord):
        print('__________________________________________________________________')
        print(f'point: {point}')
        start=point[0]  
        end=point[1]
        ra1, dec1= start.ra.radian, start.dec.radian 
        ra2, dec2=end.ra.radian, end.dec.radian
        
    
        r_sat=r[i]
        u_r_sat=u_height[i]+u_r_E
        vel_sat= v[i]
        u_vel_sat= u_vel[i]
        
        period.append(T[i])
        u_T= np.sqrt(np.float64((2*np.pi/vel_sat*u_r_sat)**2+((2*np.pi*r_sat)/(vel_sat**2)*u_vel_sat)**2))
        u_period.append(u_T)
        
        motion.append(n[i])
        u_motion= 86400/(T[i])**2*u_T
        u_mean_motion.append(u_motion)
        
        #need this first in ECI coordinates
        r1_x= r_sat * np.cos(dec1)*np.cos(ra1)
        r1_y= r_sat * np.cos(dec1)*np.sin(ra1)
        r1_z= r_sat * np.sin(dec1)
        
        r2_x= r_sat * np.cos(dec2)*np.cos(ra2)
        r2_y= r_sat * np.cos(dec2)*np.sin(ra2)
        r2_z= r_sat * np.sin(dec2)
        
        r1_vector=np.array([r1_x, r1_y, r1_z])
        r2_vector=np.array([r2_x, r2_y, r2_z])
        
        #get r-vector uncertainty
        u_r1_vector=get_u_r_vector(r_sat, u_r_sat, ra1, dec1, u_coordinate)
        u_r2_vector=get_u_r_vector(r_sat, u_r_sat, ra2, dec2, u_coordinate)
     
        #calculate displacement vector, uncertainty and its norm
        displacement=np.abs(np.subtract(r2_vector,r1_vector))
        u_disp=np.abs(np.subtract(u_r2_vector,u_r1_vector))
        disp_normed=np.linalg.norm(displacement)
        
        
        #calc state vector v with two positional vectors r1, r2
        v_vector=displacement/disp_normed * vel_sat
        u_v_vector= get_u_v_vector(displacement, u_disp, vel_sat, u_vel_sat)
        
        #get fundamental vectors h,n,e
        #TODO:later also calc orbital parameter for endpoint at r2'
        h_vector, u_h_vector, n_vector, u_n_vector,e_vector, u_e_vector = get_hne_error(v_vector, r1_vector, u_v_vector, u_r1_vector)
        
        #print(f'r1, u_r1: {r1_vector}, {u_r1_vector}') 
        #print(f'r2, u_r2: {r2_vector}, {u_r2_vector}') 
        #print(f' v, u_v: {v_vector}, {u_v_vector}')
        #print(f'h, n, e: {h_vector} {n_vector} {e_vector}')
        
        inclination = np.arccos(h_vector[2]/np.linalg.norm(h_vector))
        inclination = math.degrees(inclination)
        incl.append(inclination)
        
        inclination_type =''
        
        if inclination <= 90:
                inclination_type = 'prograde'    
        elif inclination  >=90 and inclination <=180:
                inclination_type ='retrograde'
        else: 
            raise ValueError('inclination must be within 0-180 °')
            
        
        #RAAN, Omega
        if n_vector[1]>=0:
            longitude_of_ascending_node = math.degrees(np.arccos(n_vector[0]/np.linalg.norm(n_vector)))
        else:
            longitude_of_ascending_node = math.degrees(2*np.pi-np.arccos(n_vector[0]/np.linalg.norm(n_vector)))
        long_asc.append(longitude_of_ascending_node)
        
        #omega
        argument_of_periapsis =  math.degrees(np.arccos(np.dot(n_vector, e_vector)/(np.linalg.norm(n_vector)*np.linalg.norm(e_vector)))) 
        arg_peri.append(argument_of_periapsis)
        
        #true anomaly at epoch (undefined for circular orbit)
        true_anomaly = math.degrees(np.arccos(np.dot(r1_vector, e_vector)/(np.linalg.norm(r1_vector)*np.linalg.norm(e_vector))))
        true_anom.append(true_anomaly)
        
        #u
        argument_of_latitude = math.degrees(np.arccos(np.dot(n_vector, r1_vector)/(np.linalg.norm(n_vector)*np.linalg.norm(r1_vector))))
        arg_lat.append(argument_of_latitude)
        
        #if inclination is zero
        true_longitude = argument_of_latitude + longitude_of_ascending_node
        tru_long.append(true_longitude)
        
        
        
        ###uncertainty calculation for orbital parameters###
        'u(i)'    
        h0, h1, h2= sp.symbols('h0 h1 h2')
        
        i= sp.acos(h2/sp.sqrt(h0**2+h1**2+h2**2))
        
        substitutions=substitutions = [(h0, h_vector[0]),(h1, h_vector[1]), (h2, h_vector[2])]
        
        di_dh0=sp.diff(i, h0)
        di_dh1=sp.diff(i, h1)
        di_dh2=sp.diff(i, h2)
        
        di_dh0_val= di_dh0.subs(substitutions).evalf()
        di_dh1_val=di_dh1.subs(substitutions).evalf()
        di_dh2_val = di_dh2.subs(substitutions).evalf()
        
    
        u_i= math.degrees(np.sqrt(np.float64((di_dh0_val*u_h_vector[0])**2+(di_dh1_val*u_h_vector[1])**2+(di_dh2_val*u_h_vector[2])**2)))    
        u_incl.append(u_i)
        
        'u(Omega)'
        # is zero, because u_n is zero
        n0, n1, n2= sp.symbols('n0 n1 n2')
        Omega = sp.acos(n0/sp.sqrt(n0**2+n1**2+n2*+2)) #big Omega
        
        dOm_dn0=sp.diff(Omega, n0)
        dOm_dn1=sp.diff(Omega, n1)
        dOm_dn2=sp.diff(Omega, n2)
        
        substitutions= [(n0, n_vector[0]), (n1, n_vector[1]), (n2, n_vector[2])]
        
        dOm_dn0_val= dOm_dn0.subs(substitutions).evalf()
        dOm_dn1_val= dOm_dn1.subs(substitutions).evalf()
        dOm_dn2_val= dOm_dn2.subs(substitutions).evalf()
        
        u_longasc= math.degrees(np.sqrt(np.float64((dOm_dn0_val*u_n_vector[0])**2+(dOm_dn1_val*u_n_vector[1])**2+(dOm_dn2_val*u_n_vector[2])**2)))      
        u_long_asc.append(u_longasc)
        
        'u(omega)'
        n0, n1, n2, e0, e1, e2=sp.symbols('n0, n1, n2, e0, e1, e2')
        
        o = sp.acos((n0*e0 + n1*e1 + n2*e2) / (sp.sqrt(n0**2+n1**2+n2**2) * sp.sqrt(e0**2+e1**2+e2**2)))
        
        # Calculate the partial derivatives
        do_de0 = sp.diff(o, e0) 
        do_de1 = sp.diff(o, e1) 
        do_de2 = sp.diff(o, e2) 
        do_dn0 = sp.diff(o, n0) 
        do_dn1 = sp.diff(o, n1) 
        do_dn2 = sp.diff(o, n2) 
       
        substitutions=[(n0,n_vector[0]),(n1,n_vector[1]), (n2,n_vector[2]),(e0, e_vector[0]),(e1, e_vector[1]), (e2, e_vector[2])]
        # Step 4: Substitute numerical values into the derivatives
        do_dn0_val = do_dn0.subs(substitutions)
        do_dn1_val= do_dn1.subs(substitutions)
        do_dn2_val = do_dn2.subs(substitutions)
        do_de0_val= do_de0.subs(substitutions)
        do_de1_val = do_de1.subs(substitutions)
        do_de2_val= do_de2.subs(substitutions)
        
        
        u_argperi=sp.sqrt((do_de0_val*u_e_vector[0])**2+(do_de1_val*u_e_vector[1])**2+(do_de2_val*u_e_vector[2])**2+(do_dn0_val*u_n_vector[0])**2+(do_dn1_val*u_n_vector[1])**2+(do_dn2_val*u_n_vector[2])**2)  
        u_argperi=math.degrees(u_argperi) 
        u_arg_peri.append(u_argperi)
        
        'u(nu)'
        r0, r1, r2,e0, e1, e2=sp.symbols('r0 r1 r2 e0 e1 e2')
        
        nu=sp.acos((r0*e0 + r1*e1 + r2*e2) / (sp.sqrt(r0**2+r1**2+r2**2) * sp.sqrt(e0**2+e1**2+e2**2)))
    
        # Calculate the partial derivatives
        do_de0 = sp.diff(nu, e0) 
        do_de1 = sp.diff(nu, e1) 
        do_de2 = sp.diff(nu, e2) 
        do_dr0 = sp.diff(nu, r0) 
        do_dr1 = sp.diff(nu, r1) 
        do_dr2 = sp.diff(nu, r2) 
       
        substitutions=[(r0,r1_vector[0]),(r1,r1_vector[1]), (r2,r1_vector[2]),
                       (e0, e_vector[0]),(e1, e_vector[1]), (e2, e_vector[2])]
        # Step 4: Substitute numerical values into the derivatives
        do_dr0_val = do_dr0.subs(substitutions)
        do_dr1_val= do_dr1.subs(substitutions)
        do_dr2_val = do_dr2.subs(substitutions)
        do_de0_val= do_de0.subs(substitutions)
        do_de1_val = do_de1.subs(substitutions)
        do_de2_val= do_de2.subs(substitutions)
   
        u_trueanom_temp=sp.sqrt((do_de0_val*u_e_vector[0])**2+(do_de1_val*u_e_vector[1])**2+(do_de2_val*u_e_vector[2])**2
                                +(do_dr0_val*u_r1_vector[0])**2+(do_dr1_val*u_r1_vector[1])**2+(do_dr2_val*u_r1_vector[2])**2) 
        u_trueanom=math.degrees(u_trueanom_temp)
        u_true_anom.append(u_trueanom)
      
        
        'u(u)'
        r0, r1, r2,n0, n1, n2=sp.symbols('r0 r1 r2 n0 n1 n2')
        
        uu=sp.acos((r0*n0 + r1*n1 + r2*n2) / (sp.sqrt(r0**2+r1**2+r2**2) * sp.sqrt(n0**2+n1**2+n2**2)))
    
        # Calculate thn partial derivatives
        do_dn0 = sp.diff(uu, n0) 
        do_dn1 = sp.diff(uu, n1) 
        do_dn2 = sp.diff(uu, n2) 
        do_dr0 = sp.diff(uu, r0) 
        do_dr1 = sp.diff(uu, r1) 
        do_dr2 = sp.diff(uu, r2) 
       
        substitutions=[(r0,r1_vector[0]),(r1,r1_vector[1]), (r2,r1_vector[2]),
                       (n0, n_vector[0]),(n1, n_vector[1]), (n2, n_vector[2])]
        # Step 4: Substitute numerical values into the derivatives
        do_dr0_val = do_dr0.subs(substitutions)
        do_dr1_val= do_dr1.subs(substitutions)
        do_dr2_val = do_dr2.subs(substitutions)
        do_dn0_val= do_dn0.subs(substitutions)
        do_dn1_val = do_dn1.subs(substitutions)
        do_dn2_val= do_dn2.subs(substitutions)
        
        
        u_arglat_temp=sp.sqrt((do_dn0_val*u_n_vector[0])**2+(do_dn1_val*u_n_vector[1])**2+(do_dn2_val*u_n_vector[2])**2+
                              (do_dr0_val*u_r1_vector[0])**2+(do_dr1_val*u_r1_vector[1])**2+(do_dr2_val*u_r1_vector[2])**2)
        
        u_arglat=math.degrees(u_arglat_temp)
        u_arg_lat.append(u_arglat)
        
        'u(l)'
        u_truelong=u_arglat+u_longasc
        u_true_long.append(u_truelong)
        
        #print orbital parameters for single image track
        print(f'inclination [\u03B9]: {inclination} {u_i} °--> {inclination_type}')  
        print(f'longitude of ascending node [\u03A9]: {longitude_of_ascending_node} {u_longasc}°')
        print(f'argument of periapsis [\u03C9]: {argument_of_periapsis} {u_argperi} °')
        print(f'true anomaly [\u03BD]: {true_anomaly} {u_trueanom}°')
        print(f'argument of latitude [u]: {argument_of_latitude} {u_arglat}°')
        print(f'true_longitude [l]: {true_longitude} {u_truelong}°')

        k+=1
   
        ###end of main loop###
         
    #calc weighted mean orbital elements, their internal and external uncertainty  
    #check, if internal or external uncertainty is greater                     
    mean_incl, u_int_incl, u_ext_incl=weighted_mean(incl, u_incl)
    #print(weighted_mean(incl, u_incl))
    
    mean_long, u_int_long, u_ext_long=weighted_mean(long_asc, u_long_asc)
    #print(weighted_mean(long_asc, u_long_asc))
    
    mean_period, u_int_period, u_ext_period= weighted_mean(period, u_period) 
    #print(weighted_mean(period, u_period))
    
    mean_motion, u_int_motion, u_ext_motion=weighted_mean(motion, u_mean_motion)
    #print(weighted_mean(motion, u_mean_motion))
    
    mean_peri, u_int_peri, u_ext_peri=weighted_mean(arg_peri, u_arg_peri)
    #print(weighted_mean(arg_peri, u_arg_peri))
    
    mean_anom, u_int_anom, u_ext_anom=weighted_mean(true_anom, u_true_anom)
    #print(weighted_mean(true_anom, u_true_anom))
    
    mean_arg, u_int_arg, u_ext_arg=weighted_mean(arg_lat, u_arg_lat)
    #print(weighted_mean(arg_lat, u_arg_lat))
    
    mean_truelong, u_int_truelong, u_ext_truelong=weighted_mean(tru_long, u_true_long)
    #print(weighted_mean(tru_long, u_true_long))
    
    
    print('_______________________________________________________________')
    print('mean orbital parameters for detected satellite (all apparent points included):')  
    print(f'orbital period: {mean_period} (\u00B1 {u_ext_period}) s')
    print(f'mean motion: {mean_motion} (\u00B1 {u_ext_motion})rev/day')
    print(f'inclination [\u03B9]: {mean_incl} (\u00B1 {u_ext_incl})°--> {inclination_type}')  
    print(f'longitude of ascending node [\u03A9]: {mean_long} (\u00B1 {u_ext_long})°')
    print(f'argument of periapsis [\u03C9]: {mean_peri} (\u00B1 {u_ext_peri})°')
    print(f'true anomaly [\u03BD]: {mean_anom} (\u00B1 {u_ext_anom})°')
    print(f'argument of latitude [u]: {mean_arg} (\u00B1 {u_ext_arg})°')
    print(f'true_longitude [l]: {mean_truelong} (\u00B1 {u_ext_truelong})°')
    
    
    #set bigger uncertainty as mean uncertainty:
    
    mean_orbital=   ("mean orbital parameters for detected satellite (all apparent points included): \n "
                    f" orbital period: {mean_period} (\u00B1 {u_ext_period}) s \n "
                    f"mean motion: {mean_motion} (\u00B1 {u_ext_motion})rev/day \n"
                    f"inclination [$\iota$]: {mean_incl} (\u00B1 {u_ext_incl})°--> {inclination_type} \n"
                    f"longitude of ascending node [$\Omega$]: {mean_long} (\u00B1 {u_ext_long})° \n" 
                    f"argument of periapsis [$\omega$]: {mean_peri} (\u00B1 {u_ext_peri})° \n"
                    f"true anomaly [$\nu$]: {mean_anom} (\u00B1 {u_ext_anom})° \n"
                    f"argument of latitude [u]: {mean_arg} (\u00B1 {u_ext_arg})° \n"
                    f"true_longitude [l]: {mean_truelong} (\u00B1 {u_ext_truelong})° ")
    
    return mean_orbital


