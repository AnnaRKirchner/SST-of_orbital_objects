"""
Created on Tue May  2 12:53:16 2023
@author: Anna Kirchner

Description: This script only contains one function and creates skycoordinates from detected lines in a series of images.
             It creates two lists, one with RaDec coordinates of start and endpoints, one with AltAz coordinates 
             corresponding to the time metadata of the image

Methods in this file:
1. img_track: Calculates the start and endpoint of the satellite lines given, converts them to sky coordinates with horizontal coordinates.

"""
from PIL import Image
import os
from PIL.ExifTags import TAGS
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import AltAz
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs.utils import pixel_to_skycoord
from astropy.time import Time
from dateutil import parser as ps 
from datetime import datetime


      
"""
    Calculates the start and endpoint of the satellite lines given, converts them to sky coordinates with horizontal coordinates.
    
    Parameters:
        - satellites (dict): Dictionary containing satellite names as keys and their corresponding start and endpoint coordinates as values.
        - distorted_folder (str): Path to the folder containing distorted images.
        - obs (astropy.coordinates.EarthLocation): Earth location of the observer.
        - utcoffset (datetime.timedelta): Offset of the observer's local time from UTC.
        - wcs_filename (str): Filename of the World Coordinate System (WCS) file.
    
    Returns:
        - skycoord (list): List of sky coordinates in the ICRS frame representing the start and endpoint of the satellite lines (RaDec coordinates).
        - skycord_altaz (list): List of satellite start and end locations in AltAz coordinates.
        - exposure (str): Exposure time of the images.
"""
def img_track(satellites, distorted_folder, obs, utcoffset, wcs_filename):  
    
    #This internal method gets image metadata
    def get_exif(fn):
        ret = {}
        i = Image.open(fn)
        info = i._getexif()
        
        for tag, value in info.items():
            decoded = TAGS.get(tag, tag)
            ret[decoded] = value
        return ret
    
    #This internal method returns key for any value of satellites dictionary
    def get_key(val):
        for key, value in satellites.items():
            if val == value:
                return key

        return "key doesn't exist"

    # List of start and endpoint in sky coordinates should have the same size as the number of images
    skycoord=[]
    skycord_altaz=[]  #list of satellite start and endcoordinate in alt az
    #list of altaz transforms for given location
    altaz = []
    
    # Iterate over each image file in the distorted folder
    for fn in os.listdir(distorted_folder):
        fn = os.path.join(distorted_folder, fn)
        
        if fn.endswith(".JPG") or fn.endswith(".jpg") or fn.endswith(".png") or fn.endswith(".JPEG"):
            metadata = get_exif(fn)
            exposure= metadata['ExposureTime']
            datetime_str = metadata['DateTime']
            
            # Convert the datetime string to a datetime object
            datetime_obj = datetime.strptime(datetime_str, '%Y:%m:%d %H:%M:%S')
            date_time = datetime_obj.strftime('%Y-%m-%dT%H:%M:%S')
            obs_time= datetime.strptime(date_time, '%Y-%m-%dT%H:%M:%S') - utcoffset #converted to utc
            
            #print(f'utc:{obs_time}')

            # azimuth, altitude transform for every image and its specific datetime
            aa = AltAz(location=obs, obstime=obs_time)
            # add transform to the list
            altaz.append(aa)
    
    
    # TODO: get corresponding wcs for every image from astrometry with API
   
    fn = get_pkg_data_filename(wcs_filename)
    hdu = fits.open(wcs_filename)[0]
    w = WCS(hdu.header)  # sets WCS for transformations
    
    '''
    #To check, whether the wcs transformation is valid
    px=(0,1)
    skycrd= pixel_to_skycoord(px[0], px[1], w)
    print(skycrd)
    
    pixel2= skycoord_to_pixel(skycrd, w)
    valid= np.max(np.abs(np.subtract(px, pixel2))) < 1e-6
    print(valid)
    '''
    
    for key in satellites:
    
        start = satellites[key][0]
        end = satellites[key][1]
       # print(f'start (px):{start}, end (px): {end}')
    
        # convert to skycoord:
        start_radec = pixel_to_skycoord(start[0], start[1], w)
        end_radec = pixel_to_skycoord(end[0], end[1], w)
        
        #add to list of skycoordinates (RaDec)
        new_satellite= [start_radec,end_radec]
        skycoord.append(new_satellite)
    
    
    #convert to ICRS frame
    for i, sat in enumerate(skycoord):
        for k, point in enumerate(sat):
            skycoord[i][k]=point.transform_to('icrs')
            
            
    # convert each satellite start and endpoint to AltAz 
    for idx in range(len(skycoord)):
        sat=skycoord[idx]
        sat_start=sat[0]
        sat_end=sat[1]
        
        start_aa = sat_start.transform_to(altaz[idx])
        end_aa = sat_end.transform_to(altaz[idx])
        new_aa=[start_aa, end_aa]
        skycord_altaz.append(new_aa)

    
    return skycoord, skycord_altaz, exposure


