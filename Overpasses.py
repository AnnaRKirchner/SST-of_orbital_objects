"""
Created on Tue May 23 21:45:19 2023
@author: Anna Kirchner

Description:
    This script plots detected satellite tracks and overpasses of satellites from a collection of over 4000 TLEs.
    The visualization is presented together over a specified observer location in a polar projection with the zenith at the center.
    The script utilizes multiprocessing for improved performance.

Settings:
    - Set dictionary with satellite start and endpoints with images as keys. This can be retrieved from a text file created by Height_Velocity.py.
      Example: satellites = {'IMG_2918_rectified.JPG': [(2808, 1238), (2895, 1368)],
                            'IMG_2919_rectified.JPG': [(2480, 767), (2600, 939)]}

    - Set timespan (t) for overpass.
    - Set path to distorted images (containing metadata): distorted_folder.
    - Set WCS filepath: wcs_filename.
    - Set path to TLE collection: tle_textfile.
    - Set destination path for saving the plots: dst_path.

    - Set observing location.
      Example: obs = EarthLocation(lat=69.27850*u.deg, lon=16.00930*u.deg, height=380*u.m)

    - Set UTC time with utcoffset.
"""

from multiprocessing import Pool
import os
import matplotlib.pyplot as plt
import numpy as np
from skyfield.api import wgs84, load, EarthSatellite
from SkyCoord_Transform import img_track
from astropy import units as u
from astropy.coordinates.earth import EarthLocation
import time
from datetime import timedelta

"""
    Generate a plot for a single TLE file element.
    
    Parameters:
        - tle_file_element: A list containing TLE data for a single satellite.
    
    No Returns
"""
def get_one_plot(tle_file_element):
    
    element=tle_file_element
    
    TLE = "".join(element[1:3])
    name= element[0]
    norad_nr=element[1][2:8]
    L1, L2 = TLE.splitlines()
    
    sat = EarthSatellite(L1, L2)
    do = wgs84.latlon(69.27850, 16.00930, elevation_m=14)
    difference = sat - do
    topocentric = difference.at(t)
    alt, az, distance = topocentric.altaz()
    
    # Check if the satellite is above the horizon at the given datetime
    if any(ele > 0 for ele in alt.degrees):
            
        polarplot(az.degrees, alt.degrees, name, t, norad_nr)    

"""
    Generate a polar plot for satellite overpasses.
    
    Parameters:
        - az: Azimuth values for the satellite.
        - alt: Altitude values for the satellite.
        - name: Name of the satellite.
        - timestamp: Timestamp for the overpass.
        - norad_nr: NORAD number of the satellite.
    
    No Returns
"""
def polarplot(az, alt, name, timestamp, norad_nr):  
    #Polar plot settings
    plt.style.use('ggplot')
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_ylim([0, 90])
    ax.set_yticks(range(0, 100, 10))
    ax.set_yticklabels(['90', '', '', '60', '', '', '30', '', '', '0'], fontsize=30)
    ax.grid(True, linewidth=0.5, linestyle='-', color='gray', alpha=0.5)
    ax.tick_params(axis='x', labelsize=30)
    ax.scatter(0, 0, c='green', s=70, alpha=1, label='Andenes')
   
    # to plot all objects from TLE file together
    # for row_idx in range(len(az)):
    #   Az=az[row_idx]
    #   Alt=alt[row_idx] 
    
    Az=az
    Alt=alt
    
    first_sat = True
    for idx in range(len(Az)):
        if Alt[idx]>0:
            try:
                ax.scatter(Az[idx]*np.pi/180, (90-Alt[idx]), c='blue', s=15, alpha=1, label= None if not first_sat else f'{name}\u200B{norad_nr}')
                first_sat= False
            except: 
                print('could not read satellite name')
                print(name)
                ax.scatter(Az[idx]*np.pi/180, (90-Alt[idx]), c='blue', s=15, alpha=1, label= None if not first_sat else 'Satellite  ?')
                first_sat= False
                
    
    # Plot satellite start and end positions from image track
    radec, skycoord_altaz, exposure = img_track(satellites, distorted_folder, obs, utcoffset, wcs_filename) 
    
    first_img_sat=True
    for sat in skycoord_altaz:
        start=sat[0]
        
        Az_s=start.az.value
        Alt_s=start.alt.value
        
        
        end=sat[1]
        Az_e=end.az.value
        Alt_e=end.alt.value
        
        #scatter satellites start and end position
        ax.scatter(Az_s*np.pi/180, (90-Alt_s), c='skyblue', s=10, alpha=0.5 if not first_sat else 1, label= None )
        
        ax.scatter(Az_e*np.pi/180, (90-Alt_e), c='darkblue', s=10, alpha=0.5 if not first_sat else 1, label= None )
        
        first_img_sat=False   
        
        
        #line from start to endpoint
        x_values = [Az_s*np.pi/180,Az_e*np.pi/180]
        y_values = [ 90-Alt_s, 90-Alt_e]
        plt.plot(x_values, y_values, linestyle="-", linewidth=2, color='red')
            
    ax.legend(loc=(1.1,0.6), fontsize=20)
    
    #formats timestamp
    timestamp_min= min(timestamp.utc_datetime()).strftime("%Y-%m-%d %H:%M:%S")
    date=timestamp_min.split(' ')[0]
    time_min=timestamp_min.split(' ')[1]
    timestamp_max= max(timestamp.utc_datetime()).strftime("%Y-%m-%d %H:%M:%S")
    time_max=timestamp_max.split(' ')[1]
    
    #set title and figure size       
    plt.title(f'Satellite overpasses in Andenes on {date} from {time_min} to {time_max} ', fontsize=25, fontweight='bold', y=1.05)
    
    fig.set_size_inches(30, 20)
    
    
    #save the figure to destination folder
    destination = os.path.join(dst_path, f"overpass_{0}.png")
    
    i=0 
    while(os.path.isfile(destination)):
        i+=1
        destination= os.path.join(dst_path, f"overpass_{i}.png")
    
    plt.ioff()
    plt.savefig(destination)
    plt.close()
 
 

#///   Settings   ///#

##input and output data##
# Dictionary with satellite start and enpoints with key: img, from textfile created by Height_Velocity.py
satellites = {'IMG_2975.JPG': [(2162, 2005), (2214, 2035)], 
              'IMG_2976.JPG': [(1908, 1869), (1962, 1894)], 
              'IMG_2977.JPG': [(1570, 1691), (1681, 1747)],
              'IMG_2978.JPG': [(1282, 1548), (1373, 1597)],
              'IMG_2979.JPG': [(974, 1395), (1067, 1438)], 
              'IMG_2980.JPG': [(464, 1164), (580, 1204)], 'IMG_2981.JPG': [(77, 984), (213, 1054)]} 


ts = load.timescale()
seconds = np.arange(0, 180, 0.5) 
t = ts.utc(2023, 3, 26, 23, 14, seconds)

# Path to distorted images (contains metadata)
distorted_folder = r'C:\Users\annar\Documents\BA\SkyImg\SORTED_SATELLITES_CanonG1X\Distorted\27_03_23\5'  # contains metadata

# WCS filepath for one rectified image
wcs_filename = r'C:\Users\annar\Informatik\Standard_Frame\wcs_files\CanonG1X\27_03_23\5.fits' #for one rectified image

## Location settings ##

# Observer location #

#obs = EarthLocation(lat=69.27850*u.deg, lon=16.00930 * u.deg, height=380*u.m)   # ALOMAR, Andenes
obs = EarthLocation(lat=69.306708*u.deg, lon=16.083471 * u.deg, height=14*u.m)  # DO veien, Andenes
#ob3= EarthLocation(lat = 69.287063*u.deg, lon=16.020171*u.deg, height= 200*u.m)  #Lighthouse, Andenes

# Eastern Daylight Time (winter),  UTC+2 in Andenes (summer)
utcoffset = timedelta(hours=1)

#Tle data
tle_textfile = r'C:\Users\annar\Informatik\Standard_Frame\MatchTest\TLE_26_03.txt'

#Destination folder to save plots 
dst_path= r'C:\Users\annar\Informatik\Standard_Frame\MatchTest\27_03_23\5\final'
#///////////////////////////////////////////////////////////////////////////////////////////////////

if __name__=="__main__":
    # Set up the TLE file and load satellite data
    with open(tle_textfile) as f:
        lines = f.readlines()
    
    satellite_names = lines[0::3]
    tle_files = []
    
    i = 0
    while len(lines)>=i+3:
        tle_files.append(lines[i:i+3])
        i = i + 3

    
    #get_one_plot(tle_files)
    t1=time.time()
    p=Pool()
    result =p.map(get_one_plot, tle_files)
    p.close()
    p.join()
    print(time.time()-t1)  #runtime
