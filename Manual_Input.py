"""
Created on Tue Mar 14 13:08:38 2023
@author: Anna Kirchner

Description: Class of Manual Input Stars, used in Calibration_Test.py
            This script gets px_positions and radec coordinates from a csv file,
            which is present in a textfile.
            Code has been tested with different file-formats, which are commented out below

Input: Currently: textfile with path to csv file in line 8 (beginning from 0)
"""

from os.path import exists
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd



class Input:
    
   def __init__(self, textfile):
       self.data= self.read_textfile(textfile)[0]
       self.ra= self.get_positions(self.data)[1][0]
       self.dec=self.get_positions(self.data)[1][1]
       self.px_positions=self.get_positions(self.data)[0]
       self.source=self.read_textfile(textfile)[1]
       self.SkyCoord=SkyCoord(ra= self.ra*u.degree, dec= self.dec*u.degree, frame='icrs') #SkyCoord containing list of several radec Coord
        
       
        #read and split up textfile 
   def read_textfile(self, textfile):
        data= open(textfile, 'r').readlines()
        data = [element.strip("\n").replace(" ", "") for element in data]
        data = list(filter(None, data))
        
     
        source = data[0]
        
        '''        
        pos_data = data[1].split(",")
        if len(pos_data) == 2:
            position = (float(pos_data[0]), float(pos_data[1]))

        date_data = data[2].split(",")
        if len(date_data) == 3:
            current_date = (int(date_data[0]), int(date_data[1]), int(date_data[2]))

        time_data = data[3].split(",")
        if len(time_data) == 3:
            current_time = (int(time_data[0]), int(time_data[1]), int(time_data[2]))

        arcsec = float(data[4])

        #zenith_ap_data = data[5].split(",")
        #if len(zenith_ap_data) == 2:
        #    zenith_ap = (int(zenith_ap_data[0]), int(zenith_ap_data[1]))

        #north_offset_ap = float(data[6])

        height = float(data[7])
        '''
        return data, source


   def get_positions(self,data):
        px_positions=[]
       # eq_positions=[]
        
        #for Andoya.txtformat: eq position given in hourangle
        '''
        star_px_cord_temp=[]
        star_eq_cord_temp=[]
        for index in range(8, len(data)):
            star = data[index].split("/")
            img_cord = star[0].split(",")
            img_cord = (int(img_cord[0]), int(img_cord[1]))
            star_px_cord_temp.append(img_cord)
        
            eq_cord = star[1].split(",")
            
            eq_cord_temp = ( int(eq_cord[0]), int(eq_cord[1]), int(eq_cord[2]), int(eq_cord[3]), int(eq_cord[4]), int(eq_cord[5]))
            star_eq_cord_temp.append(eq_cord_temp)
        
        
        if len(star_eq_cord_temp) > 4:
            px_positions = star_px_cord_temp
            eq_positions = star_eq_cord_temp
        else:
            print("Enter more Stars")
                    
            
        #transform the 6 eq coordinates to ra dec
        ra=[]
        dec=[]
           
        for cord in eq_positions:
            ra_temp=cord[0:3]
            dec_temp=cord[3:6]
            #translate ra hourangle to degreee
            ra.append((ra_temp[0]+ra_temp[1]/60+ra_temp[2]/3600)*15)
            dec.append(dec_temp[0] + (dec_temp[1] / 60) + (dec_temp[2] / 3600))
            
        radec=[ra, dec]
        '''
        
        #for Andoya2.txt format: eq position already given in ra dec (from fits file)
        '''
        ra=[]
        dec=[]
        for index in range(8, len(data)):
            star = data[index].split("/")
            img_cord = star[0].split("," )
        
            img_cord = (float(img_cord[0]), float(img_cord[1]))
            px_positions.append(img_cord)
     
            eq_cord = star[1].split(",")
            
            ra.append(float(eq_cord[0]))
            dec.append( float(eq_cord[1]))
    
        radec=[ra, dec]
        '''
        
        #radec from ods file
        base_path = data[8]
        
        df = pd.read_csv(base_path , usecols=['x', 'y', 'ra', 'dec'])
        
        ra=[]
        dec=[]
        for row in range(0,len(df)):
            img_cord= (float(df.loc[row, "x"]), float(df.loc[row, "y"]))
            px_positions.append(img_cord)
            ra.append(float(df.loc[row, "ra"]))
            dec.append(float(df.loc[row, "dec"]))
            
        radec=[ra, dec]
        
        
        return px_positions, radec


