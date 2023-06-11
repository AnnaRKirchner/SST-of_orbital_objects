"""
Created on Sat Jun  3 21:38:44 2023
@author: Anna Kirchner

Description: This script calculates and plots velocities and heights of satellites based on their start and endpoints.           
             The satellite lines are obtained either through 
             manual input or by using the satellite line detection API. 
             The results, including orbital parameters and velocity-heights plots, are saved in output files.
             
Settings:
    - input_file: Path to the input file containing the satellite positions in RA/Dec coordinates.
    - exposure_time: Exposure time in seconds for the satellite images.
    - satellite_lines: Path to the file containing the manually input satellite lines or 'API' to use the 
                         satellite line detection API.
    - output_directory: Path to the directory where the output files and plots will be saved.
    
    - api_key: API key for the satellite line detection service (required if `satellite_lines` is set to 'API').
    - Earth_location: Earth location coordinates (longitude, latitude, altitude) in degrees and meters. 
                        Default is (0, 0, 0), which represents the Greenwich Observatory.

Note: Make sure to have the required input files and packages installed before running the script.
"""
import numpy as np
import matplotlib.pyplot as plt
import os
from sympy.solvers import solve
from sympy import Symbol, I
import sympy as sp
from datetime import timedelta
from astropy import constants as const
from astropy import units as u
from astropy.coordinates import EarthLocation
from astropy.io import fits
from FinalLineDetection import satellite_LineDetection
from SkyCoord_Transform import img_track  
from OrbitParams_Copy import get_orbital_params, get_earth_radius, weighted_mean


#for LaTeX fond in plots
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "lmodern"
})


"""
    Calcualte velocity and height of image track by equalizing apparent velocity calculations
    (1. equation dependent on height and altitude, 2. equation from spherical distance of start-endpoint divided by exposure time)
    
    Parameters:
        - radec: A list containing skycoordinates with start and endpoints of image tracks from image series
        - Altaz: A list containing Alt, Az skycoordinates for start and endpoints of image tracks from image series
        - exposure_time: exposure of one single image
        - u_coordinate: uncertainty for coordinate transformation
    
    Returns:
        - v, u_v: velocity magnitude and propagated uncertainty for one image track
        - sat_height, u_h: height above Earth's surface and propagated uncertainty for one image track
"""
def get_velocity(radec, Altaz, exposure_time, u_coordinate): 
    G= const.G.value
 
    distance=[]
    u_dist=[]
    app_v=[]
    for idx in range(len(radec)):
        sat=radec[idx]
        sat_start=sat[0]
        sat_end=sat[1]
       
        a1, del1 = sat_start.ra.radian, sat_start.dec.radian
        a2, del2= sat_end.ra.radian, sat_end.dec.radian
        
        u_ra1=u_ra2=u_dec1=u_dec2=u_coordinate
    
        spherical_distance= np.arccos(np.cos(del1)*np.cos(del2)*np.cos(a2-a1)+ np.sin(del1)*np.sin(del2))
        
        # Calculate apparent velocity from distance and exposure time
        apparent_velocity = spherical_distance / exposure_time   #[arc/s]
        app_v.append(apparent_velocity) #one velocity for one imagetrack
        distance.append(spherical_distance)
        
        'u(dist): uncertainty of spherical distance'
        # Define the variables
        dec1, dec2, ra1, ra2 = sp.symbols('dec1 dec2 ra1 ra2 ')
        
        #function
        d= sp.acos(sp.cos(dec1)*sp.cos(dec2)*sp.cos(ra2-ra1)+ sp.sin(dec1)*sp.sin(dec2))
        
        # Calculate the partial derivatives
        dd_dra1 = sp.diff(d, ra1)  
        dd_ddec1 = sp.diff(d, dec1)  
        dd_dra2 = sp.diff(d, ra2)  
        dd_ddec2 = sp.diff(d, dec2)
        
        # Step 4: Substitute numerical values into the derivatives
        dd_dra1_value = dd_dra1.subs([(dec1, del1), (dec2, del2), (ra1, a1), (ra2, a2)])
        dd_ddec1_value = dd_ddec1.subs([(dec1, del1), (dec2, del2), (ra1, a1), (ra2, a2)])
        dd_dra2_value = dd_dra2.subs([(dec1, del1), (dec2, del2), (ra1, a1), (ra2, a2)]) 
        dd_ddec2_value = dd_ddec2.subs([(dec1, del1), (dec2, del2), (ra1, a1), (ra2, a2)])

        u_distance= sp.sqrt((dd_dra1_value*u_ra1)**2+(dd_ddec1_value*u_dec1)**2+(dd_dra2_value*u_ra2)**2+(dd_ddec2_value*u_dec2)**2)  
        u_dist.append(u_distance.evalf())
    
    
    #list with all detected altitudes (start and endpoints)
    altitudes=[]  
    # list with mean altitudes of track middle-point
    mean_altitudes=[]
    for sat in Altaz:
        start=sat[0]
        end=sat[1]
        
        Alt_s=start.alt.value  
        Alt_e=end.alt.value
    
        altitudes.append(Alt_s)
        altitudes.append(Alt_e)
        
        mean_alt= (Alt_s+Alt_e)/2
        mean_altitudes.append(mean_alt)
        
    #velocity for one image track
    v=[]
    # corresponding uncertainties
    u_v=[] 
    
    #height calculated for middle points
    sat_height=[]   
    # corresponding uncertainties
    u_h=[]
    
    
    for i in range(len(mean_altitudes)):
    #apparent velocity for middle point of track
    
        b=app_v[i] 
        dist=distance[i]
        u_dst=u_dist[i]
       
        
        h = Symbol('h')
        height_app=solve(np.sin(np.radians(mean_altitudes[i]))**4 *(G * m_E) / ((r_E+h)*h**2)-b**2, h) 
        
              
        'u(h)'
                        #define variables
        d, alt, gamma, m, r , exp=sp.symbols('d alt G m r exp')
       
        zaehler=sp.cbrt(sp.sqrt((27*(exp/d)**2*gamma*m*sp.sin(alt)**4 -2*r**3)**2 -4*r**6)+27*(exp/d)**2*gamma*m*sp.sin(alt)**4 -2*r**3)
        
        #check which solution of h is the positive one and choose the corresponding uncertainty calculation
        #h=zaehler/(3*np.cbrt(2))+ np.cbrt(2)*r**2/(3*zaehler) - r/3
        #h= -(0.13228 -0.22912*I)*zaehler - ((0.20999 +0.36371 *I)*r**2)/(zaehler) - r/3
        h=-(0.13228 -0.22912*I)*zaehler - ((0.20999 -0.36371 *I)*r**2)/(zaehler) - r/3   
        
        #calc derivatives
        dh_dgamma=sp.diff(h, gamma)
        dh_dm=sp.diff(h, m)
        dh_dalt=sp.diff(h, alt)
        dh_dr=sp.diff(h, r)
        #dh_dexp=sp.diff(h, exp) --> uncertainty is zero 
        dh_dd=sp.diff(h, d)
        
        # Step 4: Substitute numerical values into the derivatives
        dh_dG_val = dh_dgamma.subs([(d, dist), (alt, np.radians(mean_altitudes[i])), (gamma, G), (m, m_E), (r, r_E), (exp, exposure_time)])
        dh_dm_val=dh_dm.subs([(d, dist), (alt, np.radians(mean_altitudes[i])), (gamma, G), (m, m_E), (r, r_E), (exp, exposure_time)])
        dh_dalt_val=dh_dalt.subs([(d, dist), (alt, np.radians(mean_altitudes[i])), (gamma, G), (m, m_E), (r, r_E), (exp, exposure_time)])
        dh_dr_val=dh_dr.subs([(d, dist), (alt, np.radians(mean_altitudes[i])), (gamma, G), (m, m_E), (r, r_E), (exp, exposure_time)])
        #dh_dexp=sp.diff(h, exp) --> uncertainty is zero 
        dh_dd_val=dh_dd.subs([(d, dist), (alt, np.radians(mean_altitudes[i])), (gamma, G), (m, m_E), (r, r_E), (exp, exposure_time)])
        
        u_alt=u_coordinate
        u_d= u_dst
        u_height= sp.sqrt((dh_dG_val*u_G)**2+(dh_dm_val*u_m)**2+(dh_dalt_val*u_alt)**2+(dh_dr_val*u_r)**2+(dh_dd_val*u_d)**2)

        u_height_complex = u_height.evalf()
        u_height_real = sp.re(u_height_complex)    
        u_h.append(float(u_height_real))
        
        
        #get real component of  3 height solutions
        height_app= np.array(height_app).astype('complex128')
        height_real=np.real(height_app)
        
        for angular in height_real:
            #pass only positive solution
            if angular >0:
                
                velocity= np.sqrt((G * m_E) / (r_E + angular))
                
                v.append(velocity)
                sat_height.append(angular)
                
                'u(v)'
                #define variables for uncertainty calc ulation
                gamma, m, r, h=sp.symbols('G m r h')   
                
                #define function
                velocity= sp.sqrt((gamma * m) / (r + h))
                
                           #calc derivatives
                dv_dG=sp.diff(velocity, gamma)       
                dv_dm= sp.diff(velocity, m)
                dv_dr= sp.diff(velocity, r)
                dv_dh= sp.diff(velocity, h)
                
                # Step 4: Substitute numerical values into the derivatives
                dv_dG_val = dv_dG.subs([(gamma, G), (m, m_E), (r, r_E), (h, angular)])
                dv_dm_val = dv_dm.subs([(gamma, G), (m, m_E), (r, r_E), (h, angular)])
                dv_dr_val = dv_dr.subs([(gamma, G), (m, m_E), (r, r_E), (h, angular)])
                dv_dh_val = dv_dh.subs([(gamma, G), (m, m_E), (r, r_E), (h, angular)])

                
                u_velocity= sp.sqrt((dv_dG_val*u_G)**2+(dv_dm_val*float(u_m))**2+(dv_dr_val*u_r)**2+(dv_dh_val*u_height_real)**2) 
                u_v.append(float(u_velocity.evalf()))
 
    return v, sat_height , u_v, u_h

#///  Settings ///#
###input and output data###

#save textfile with orb.parameters and values in
output_folder=r'C:\Users\annar\Informatik\Standard_Frame\velocities\23_03_23\2'
#save plots in
output_plots=r'C:\Users\annar\Informatik\Standard_Frame\velocities\23_03_23\2' 

#if satellites still need to be detected with Line-Detection API
rectified_images= r'C:\Users\annar\Documents\BA\SkyImg\SORTED_SATELLITES_CanonG1X\Rectified\23_03_23\2'
# Dictionary with satellite start and enpoints as values with key: img, from "FinalLineDetection.py"
#satellites=satellite_LineDetection(rectified_images, output_folder) #calls line detection, opens gui windows and returns satellite lines dictionary

# or manual input of satellite lines, if linedetection has already been done and API can be skipped (time intensive)'
satellites={'IMG_2877_rectified.JPG': [(3418, 205), (3429, 362)], 'IMG_2878_rectified.JPG': [(3446, 608), (3453, 750)], 'IMG_2879_rectified.JPG': [(3472, 993), (3480, 1104)], 'IMG_2880_rectified.JPG': [(3495, 1296), (3502, 1389)], 'IMG_2881_rectified.JPG': [(3517, 1545), (3525, 1628)]} 

#for Andoya, jpgs or other image format
# contains metadata
distorted_folder = r'C:\Users\annar\Documents\BA\SkyImg\SORTED_SATELLITES_CanonG1X\Distorted\23_03_23\2'
#wcs file for one rectified image 
w = r'C:\Users\annar\Informatik\Standard_Frame\wcs_files\CanonG1X\23_03_23\2.fits'

###Location settings###
# observer location
#obs = EarthLocation(lat=69.27850*u.deg, lon=16.00930 * u.deg, height=380*u.m)  # ALOMAR, Andenes
#obs=EarthLocation(lat=-32.38072*u.deg, lon= 20.81078*u.deg, height=1761*u.m)   #South Africa (DLR Data) 
obs = EarthLocation(lat=69.306708*u.deg, lon=16.083471 * u.deg, height=14*u.m)  # DO veien, Andenes

# Eastern Daylight Time (winter),  UTC+2 in Andenes (summer)
utcoffset = timedelta(hours=1)
#//////////

wcs_errors=[]
#TODO: get code_error for more than only one wcs file
# can open 2  FITS files and calculate the mean error from it
hdul1 = fits.open(w)


#Initialize code error from the wcs transformation
code_error= float(hdul1[0].header['comment'][-13][12:])
wcs_errors.append(np.radians(code_error))
print("Code error (rad):", np.radians(code_error)) 

global scale, u_pix
scale= float(hdul1[0].header['comment'][-8][7:14])
print('Scale:', scale)

u_pix_arc= 30*scale  #rather high estimate
# Convert arcseconds to radians
u_pix = np.radians(u_pix_arc  / 3600)
print(f'u(pix) (rad): {u_pix}')


global u_wcs
u_wcs=np.mean(wcs_errors)
print(f'u(wcs_mean):{u_wcs}')

# u(transform)=u(coordinate)
global u_coordinate
u_coordinate=np.sqrt((u_pix)**2+ (u_wcs)**2)
print(f'u(px->radec) transformation (rad): {u_coordinate}')


#constants and uncertainties
global G, u_G, r_E, u_r, m_E, u_m, u_position, u_exp

G= const.G.value
u_G= 6.7*10**(-15) #[m^3/(kg*s^2)]

#get Earth radius and uncertainty for specific latitude of Earthlocation
r_E, u_r= get_earth_radius(obs, u_coordinate)

m_E=const.M_earth.value
u_m= 6*10**20 #[kg]

u_position=0 # u(latitude)  assumed to be 0, neglectable, but for more precision: get GPS error
u_exp=0 # u(exp), neglectable, but there is a time accuracy error for every camera

#get skycoordinates start and endpoint in horizontal system (az, alt)
radec, aa, exposure= img_track(satellites, distorted_folder, obs, utcoffset, w) 

#calculate velocity and height magnitude and their uncertainties
#h is only mean height from calculated middle point of the track
v, h, u_v, u_h=get_velocity(radec, aa, exposure, u_coordinate)
print(f'satellites velocity : {v} and u(v): {u_v}')
print(f'satellites height (above earth): {h} and u(h): {u_h}')


#get orbital parameters for each middle point
mean_orbital=get_orbital_params(obs, radec, h,u_h, v, u_v, exposure, u_coordinate) 


#///plots///#
#list with all detected altitudes (start and endpoints)
alt_list=[]   
mean_altitudes=[]
for sat in aa:
    start=sat[0]
    end=sat[1]
   
    Alt_s=start.alt.value  
    Alt_e=end.alt.value
    
    #alt_list.append(Alt_e)
    #alt_list.append(Alt_s)
    
    mean_alt= (Alt_s+Alt_e)/2
    mean_altitudes.append(mean_alt)
    
    #sorts AzAlt
    if aa[0][0].alt.value<aa[1][0].alt.value:
        if Alt_s<Alt_e:
           alt_list.append(Alt_s)
           alt_list.append(Alt_e)
        else:
            alt_list.append(Alt_e)
            alt_list.append(Alt_s)
    else:  
        if Alt_s<Alt_e:
           alt_list.append(Alt_e)
           alt_list.append(Alt_s)
        else:
            alt_list.append(Alt_s)
            alt_list.append(Alt_e)
                        
print(f'altitudes:{alt_list}') 

#save (satellite lines, code error, scale, height, velocity and their uncertainties, mean orbital parameters with uncertainties, altitudes) to textfile
#textfile= open(f"{output_folder}\\orbital_info.txt", "w+")
#textfile.write(f"satellites:{satellites} \n\n code error:{code_error} \n\n scale:{scale} \n\n velocities:{v}\n u(v):{u_v} \n\n heights:{h}\n u(h):{u_h} \n altitudes: {alt_list} \n\n "+ mean_orbital )
#textfile.close()


act_vel=[]
app_vel=[]
app_vel40=[]  
app_vel30=[]
app_vel_images=[] 
sat_heights=range(0, 20000000, 100)  #von 0 bis 20000 km 

for height in sat_heights:
    actual_velocity = np.sqrt((const.G.value * const.M_earth.value) / (r_E + height))
    act_vel.append(actual_velocity)
    
    apparent_vel= np.sin(np.radians(60))*actual_velocity  #in km/s
    apparent_vel40=np.sin(np.radians(40))*actual_velocity
    apparent_vel30=np.sin(np.radians(30))*actual_velocity
    app_vel40.append(apparent_vel40) 
    app_vel30.append(apparent_vel30)                  
    app_vel.append(apparent_vel)
    
temp=[]
for deg in alt_list:  
    if temp != []:
       app_vel_images.append(temp)
       temp.clear()
    for height in sat_heights:
        actual_velocity = np.sqrt((const.G.value * const.M_earth.value) / (const.R_earth.value + height))
        apparent_vel_images= np.sin(np.radians(deg))*actual_velocity  #with altitude angle from start and endpoints
        temp.append(apparent_vel_images)



#Plot 1: actual and apparent velocities (m/s) at different altitudes over satellite_heights
plt.figure(1, figsize=(10, 8))
keppler= plt.plot(sat_heights, act_vel, label=r"$v_{circ}$", color='green')       
plt.plot(sat_heights, app_vel, label='apparent velocity at 60° altitude', color="orange")
plt.plot(sat_heights, app_vel40, label='apparent velocity at 40° altitude', color="purple")
plt.plot(sat_heights, app_vel30, label='apparent velocity at 30° altitude', color="blue" )
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.grid(True, color='lightgray')
plt.xlabel(r'Height above Earth''s surface $h$(m)', fontsize=20, labelpad=10)
plt.ylabel(r'Orbit velocity $v_{circ} \quad (m\cdot s^{-1})$', fontsize=20, labelpad=15)  #write in latex
plt.ticklabel_format(style='sci', scilimits=(3,3), axis='x', useOffset=(False))
plt.ticklabel_format(style='sci', scilimits=(3,3), axis='y', useOffset=(False))
plt.legend(fontsize=16)
#plt.show()


#create different lists with apparent radial velocities
app_vel_60=[]
app_vel_30=[]
app_vel_40=[]
app_vel_90=[]

heights= range(2, 4000000, 100) #bis 4000 km
for height in heights:
        apparent_velocity_90= np.sin(np.radians(90))**2 *actual_velocity/ (height)
        apparent_velocity_60= np.sin(np.radians(60))**2 *actual_velocity/ (height)
        apparent_velocity_45= np.sin(np.radians(40))**2 *actual_velocity/ (height) 
        apparent_velocity_30= np.sin(np.radians(30))**2 *actual_velocity/ (height)
        app_vel_90.append(apparent_velocity_90*180/np.pi)
        app_vel_60.append(apparent_velocity_60*180/np.pi)
        app_vel_40.append(apparent_velocity_45*180/np.pi)
        app_vel_30.append(apparent_velocity_30*180/np.pi)
 
        
#Plot 2: Angular/Radial velocity for different altitudes
plt.figure(2, figsize=(10, 8))
plt.xlabel('Height above Earth\'s surface $h$(m)',fontsize=16, labelpad=10)
plt.ylabel(r'Angular velocity $\omega_{circ} \quad (^{\circ} \cdot s^{-1})$', fontsize=16, labelpad=15) 
plt.plot(heights, app_vel_90, color='green', label = "90°")
plt.plot(heights, app_vel_60, color='orange', label = "60°")
plt.plot(heights, app_vel_40, color='purple', label = "40°")
plt.plot(heights, app_vel_30, color='blue', label = "30°")
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.ticklabel_format(style='sci', scilimits=(3,3), axis='x', useOffset=(False))
plt.grid(True, color='lightgray')
plt.ylim([0, 1.25])
plt.xlim([0, 4000000])
plt.legend(fontsize=12)
#plt.show()


x= range(1, len(alt_list)+1)
alt_error= np.degrees(u_coordinate)  

#Plot 3: Altitude plot- start and enpoints altitude in degree over time t (t not scaled)
plt.figure(3, figsize=(10, 8))
plt.xlabel('Satellites start-and endpoints in chronological order', fontsize=50, labelpad=10)
plt.ylabel(' Altitude (°)', fontsize=50, labelpad=15)  
plt.scatter( x, alt_list, color='black', label='Altitude', s=100)
plt.errorbar(x, alt_list,  yerr=alt_error, fmt='none', capsize=6, color= 'red', label='u(Alt.)' , elinewidth=2)
#plt.ylim([50, 90])
#plt.xlim([0, 17])
plt.yticks(fontsize=40)
plt.xticks(fontsize=40)
ax = plt.gca()
ax.yaxis.get_offset_text().set_fontsize(40)
ax.xaxis.get_offset_text().set_fontsize(40)

figure = plt.gcf() # get current figure
figure.set_size_inches(19, 15)
plt.grid(True, color='lightgray')
plt.legend(fontsize=50)
#plt.show()

#save
if os.path.isfile(f'{output_plots}\\altitudes.png'):
    os.remove(f'{output_plots}\\altitudes.png')
plt.savefig(f'{output_plots}\\altitudes.png', dpi=150)


#Plot 4: v over h for all image tracks 
plt.figure(4, figsize=(10, 8))
plt.xlabel('Height above Earth\'s surface h(m)' , fontsize=50, labelpad=10)
plt.ylabel(r' Orbit velocity $v_{circ} \quad (m \cdot s^{-1})$',fontsize=50, labelpad=15)
plt.plot(sat_heights, act_vel, label=r"$v_{circ}$", color='green') 
plt.scatter(h, v, color='blue', label='detected satellite', s=100)
plt.errorbar(h, v,  yerr=u_v, xerr= u_h, fmt='none', capsize=6, color= 'red', label='deviation',  elinewidth=2)
plt.yticks(fontsize=40)
plt.xticks(fontsize=40)
plt.ticklabel_format(style='sci', scilimits=(3,3), axis='x', useOffset=(True))
plt.ticklabel_format(style='sci', scilimits=(3,3), axis='y', useOffset=(True))
ax = plt.gca()
ax.yaxis.get_offset_text().set_fontsize(40)
ax.xaxis.get_offset_text().set_fontsize(40)

figure = plt.gcf() # get current figure
figure.set_size_inches(19, 15)
plt.xlim([200000, 1400000])
plt.ylim([6500, 8000])    
plt.legend(loc='lower left',fontsize=50)
plt.grid(True, color='lightgray')
plt.show()


# save in output folder
if os.path.isfile(f'{output_plots}\\v_detected.png'):
    os.remove(f'{output_plots}\\v_detected.png')
plt.savefig(f'{output_plots}\\v_detected.png', dpi=150)




#get weighted mean of v and h
v_mean, u_int_v, u_ext_v=weighted_mean(v, u_v)
h_mean, u_int_h, u_ext_h=weighted_mean(h, u_h)

#checks, if internal or external uncertainties are bigger, sets bigger uncertainty and label
if u_int_v>u_ext_v:
    u_mean_v=u_int_v
    
    if u_int_h>u_ext_h:
        u_mean_h=u_int_h
        int_or_ext= 'internal uncertainty $u_{int}(\overline{v}_{sat})$ and $u_{int}(\overline{h}_{sat})$'

    else:
        u_mean_h=u_ext_h
        int_or_ext= 'internal uncertainty $u_{int}(\overline{v}_{sat})$ and $u_{ext}(\overline{h}_{sat})$'
else:
    u_mean_v=u_ext_v
    
    if u_int_h>u_ext_h:
        u_mean_h=u_int_h
        int_or_ext='external uncertainty $u_{ext}(\overline{v}_{sat})$ and $u_{int}(\overline{h}_{sat})$'

    else:
        u_mean_h=u_ext_h
        int_or_ext= 'internal uncertainty $u_{int}(\overline{v}_{sat})$ and $u_{ext}(\overline{h}_{sat})$'
    

#Plot 5: Weighted means v over h 
plt.figure(5, figsize=(10, 8))
plt.xlabel('Height above Earth\'s surface h(m)' , fontsize=50, labelpad=10)
plt.ylabel(r' Orbit velocity $v_{circ} \quad (m \cdot s^{-1})$',fontsize=50, labelpad=15)
plt.plot(sat_heights, act_vel, label=r"$v_{circ}$", color='green') 

plt.errorbar(h_mean, v_mean, color='red', xerr=float(u_mean_h), yerr=float(u_mean_v), fmt='none', capsize=6, label=int_or_ext , elinewidth=2)
plt.scatter(h_mean, v_mean, color='purple', label=' weighted mean $\overline{v}_{sat}$', s=100)
# Add a box annotation
plt.annotate(f'{h_mean/1000:.0f}({u_mean_h/1000:.0f}) $\cdot 10^3$, {v_mean/1000:.3f}({u_mean_v:2.0f}) $\cdot 10^3$', 
             xy=(h_mean, v_mean), xycoords='data', xytext=(0, 60),
             textcoords='offset points', bbox=dict(boxstyle='round,pad=0.4', fc='white', lw=0.5), size=40)

plt.ticklabel_format(style='sci', scilimits=(3,3), axis='x', useOffset=(True))
plt.ticklabel_format(style='sci', scilimits=(3,3), axis='y', useOffset=(True))
ax = plt.gca()
ax.yaxis.get_offset_text().set_fontsize(50)
ax.xaxis.get_offset_text().set_fontsize(50)

figure = plt.gcf() # get current figure
figure.set_size_inches(19, 15)
plt.xlim([200000, 1000000])
plt.ylim([6500, 8000])
plt.grid(True, color='lightgray')
plt.legend(fontsize=40)
plt.yticks(fontsize=40)
plt.xticks(fontsize=40)
plt.show()

 
if os.path.isfile(f'{output_plots}\\v_mean.png'):
    os.remove(f'{output_plots}\\v_mean.png')
plt.savefig(f'{output_plots}\\v_mean.png', dpi=150)







