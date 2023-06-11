"""
Created on Sat Mar 11 16:40:21 2023
@author: Anna Kirchner

Description:
    Test of Calibration with Astrometry, but without the need of Checkerboard Images.
    Script produces wrong distortion parameters so far.
    The script uses the astropy library to perform transformations and calculations with a World Coordinate System (WCS).
    Includes functions to plot and visualize the distortion parameters.

Settings:
    - textfile: txt file containing the path to the image file and other relevant information (see example)
    - wcsfile: path to the corresponding wcs file for the image
    - rescale_window: set aspect ratio of rescaled window. Default= 2000
    - starcatalog (ra, dec, mag): set right ascension, declination and magnitude from HYG star database
    - magnitude_threshold: set magnitude_threshold for star detection
    - arcsec_per_px: (TODO) should be withdrawn from wcs file
    
Methods in this file: 
1.  undistort_point: This function takes in distortion parameters, distorted points, and tangential distortion parameters and corrects the distorted points using the provided distortion parameters.
2.  ResizeWithAspectRatio: This function resizes an image while maintaining the original aspect ratio.
3.  resize_pixel: This function translates a pixel position after the image has been resized.
4.  get_tangential_params: This function calculates the tangential distortion parameters.
5.  plot_distortion: This function computes radial and tangential distortion between sets of points.
  
"""

import cv2
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs.utils import skycoord_to_pixel, pixel_to_skycoord
from astropy.coordinates import match_coordinates_sky, concatenate
from scipy.optimize import curve_fit
from sympy import symbols, solve
from statistics import fmean
from Manual_Input import Input

"""
    Gets list of distorted points (px coordintes in image) and corrects them with given distortion parameters, 
    which fits to the original img size(camera dependent)

    Parameters:
        - radial_dist:
        - distorted_points:
        - tangential_dist_param:
    
    Returns:
        - corrected: 
"""
def undistort_point(radial_dist, distorted_points, tangential_dist_param):  # 
    
    image_center=(int(w.wcs.crpix[0]),  int(w.wcs.crpix[1]))
    k_1, k_2, k_3= radial_dist
    p_1, p_2=tangential_dist_param
    
    corrected=[]
    for star in distorted_points:
        star2center_vector = np.subtract(star, image_center)
        distance2center=np.linalg.norm(star2center_vector)
        x=star[0]
        y=star[1]
        tangential_distortion= [p_1*(distance2center**2+2*x**2)+ 2*p_2*x*y, p_2*(distance2center**2+2*y**2) + 2*p_1*x*y]
        radial_distortion= k_1*distance2center+(k_2 * distance2center**3) +  (k_3 * distance2center ** 5)
       
        #https://www.tangramvision.com/blog/camera-modeling-exploring-distortion-and-distortion-models-part-ii
        if distance2center!= 0:
            ratio = radial_distortion / distance2center   
        else:
            ratio=0
            
        radial_distortion_vector = star2center_vector +(star2center_vector* ratio)
        corrected_point = image_center+(star2center_vector+ (star2center_vector*ratio)+tangential_distortion) #wiki distortion calc
        
        
        corrected_point = (int(corrected_point[0]), int(corrected_point[1]))
        corrected.append(corrected_point)
          
    return corrected

"""
This function resizes an image with the same ratio as the original image
"""
def ResizeWithAspectRatio(image, rescale_window=None, height=None, inter=cv2.INTER_AREA):
    dim = None
    (h, w) = image.shape[:2]

    if rescale_window is None and height is None:
        return image
    if rescale_window is None:
        r = height / float(h)
        dim = (int(w * r), height)
    else:
        r = rescale_window / float(w)
        dim = (rescale_window, int(h * r))

    return cv2.resize(image, dim, interpolation=inter)

"""
This function translates a pixel position after the image has been resized

Parameters: 
    -og_image: original image with original shape
    -resized_image: resized image with resized shape

Returns:
    -resized: resized pixel position
"""
def resize_pixel(og_image, resized_image, px_point):
    
    #bigger= og_image.shape[0]>resized_image.shape[0]
    
    x=px_point[0]
    y=px_point[1]
    ratio_height= resized_image.shape[0]/og_image.shape[0]
    ratio_width= resized_image.shape[1]/og_image.shape[1]

    xS=x*ratio_width
    yS= y*ratio_height
  
        
    resized=(int(xS), int(yS))  
    return resized


    
"""
This function calculates the tangential distortion parameters

Parameters:
    - r:distanceB
    - position= (x= pointsB[i][0],y= pointsB[i][1])
    - t: tangential_dist_vec=(x_d= tangential_dist_vec[i][0], y_d=tangential_dist_vec[1])
Returns: 
    - mean_params: mean tangential parameters [p1, p2]
"""
def get_tangential_params(r,position,t):
    p1=[]
    p2=[]
    for i in range(len(r)):
        x=position[i][0]
        y=position[i][1]
        x_d=t[i][0]
        y_d=t[i][1]
        p_1, p_2= symbols('p_1 p_2')
        eq1= p_2*(r[i]**+2*y**2)+2*p_1*x*y-y_d
        eq2=p_1*(r[i]**2*x**2)+2*p_2*x*y-x_d
        params= solve((eq1,eq2), (p_1, p_2))
        p1.append(params[p_1])
        p2.append(params[p_2])
        
        
    p1=fmean(p1)
    p2=fmean(p2)
    
    mean_params=[p1, p2]
    return mean_params

"""
   Computes radial distortion between points/arrays of points A,B
           
   Parameters: 
       - pointsA: correct position
       - pointsB: stars to be rectified/undistorted
       - img: distorted image
       - wcs
    Returns: 
        - param_radial, param_cov : radial distortion parameters and covariance
        - param_tangential (param_cov_tang): tangential distortion parameters and covariance
"""
def plot_distortion(img, pointsA, pointsB, wcs):
    
    # calculates only radial distortion
    if type(pointsB[1]) == SkyCoord:
        
       # img_center= pixel_to_skycoord(int(img.shape[1] / 2), int(img.shape[0] / 2), w)
        img_center=pixel_to_skycoord(int(w.wcs.crpix[0]),  int(w.wcs.crpix[1]), w) #reference point of wcs
        
        pointsA=[] #radec points
        pointsB=[] #radec points
        for i in range(len(pointsA)):
             vecA=(pointsA[i].ra.degree, pointsA[i].dec.degree)
             vecB=(pointsB[i].ra.degree, pointsB[i].dec.degree)
             pointsA.append(vecA)
             pointsB.append(vecB)
                  
    else:
       # img_center= (int(img.shape[1] / 2), int(img.shape[0] / 2))
       img_center=(int(w.wcs.crpix[0]),  int(w.wcs.crpix[1])) 
       
    #y-axis: offset equals distortion of point A to B
    offset=[]
    for i in range(len(pointsA)):
        dist_vec = tuple(np.subtract(pointsA[i],pointsB[i]))
        distortion=np.linalg.norm(dist_vec)
        offset.append(distortion) #distance (norm of vector) from point A to B (undistorted to distorted)
    
   
    #x-axis: distance of distorted px position B to img center
    distanceB=[]
   
    for i in range(len(pointsB)):
            dist_vecB = tuple(np.subtract(img_center,pointsB[i]))
            dist=np.linalg.norm(dist_vecB)
            distanceB.append(dist) #distance (norm of vector) from ref point to B (undistorted to distorted)
        
    

    #uncertainties for distance: estimated for uncertainty of pixel-moments 
    #-->u_distance very neglectable, cause k-mean algorithm used in astrometry is very certain method to calc the moments
    #img position und transformation with wcs
    #4 pixel in each direction is probably alread too high estimate
    u_A=(4,4)  
    u_B=(4,4)
    u_dist=tuple(np.subtract(u_A, u_B))
    #should be same for distance and distortion value
    u=np.linalg.norm(u_dist) 
    
    #///   Regression analysis: distortion plot   ///#
    #tested different regression functions, this one fitted the Brown-Conrady model best   
    def regression(r, k_1, k_2, k_3):
        y=k_1*r+(k_2 * r**3) +  (k_3 * r ** 5)#+(k_3 * r ** 6)
        return y
    
     # curve_fit() function takes the test-function
     # x-data and y-data as argument and returns
     # the coefficients a,b,c,d, etc... in param and
     # the estimated covariance of param in param_cov
    param_radial, param_cov = curve_fit(regression, distanceB, offset) 
  
    k_1, k_2, k_3=param_radial
    
    r_line = np.arange(min(distanceB), max(distanceB), 1)
    y_line = regression(r_line, k_1, k_2, k_3)

    # tangential distortion : distortion in  x and y direction
    x_axis_temp=[]
    y_axis_temp=[]
    for point in pointsA:
        x_axis_temp.append(point[0])
        y_axis_temp.append(point[1])
    
    x_offset=[]
    for i in range(len(pointsA)):
        dist_vec=tuple(np.subtract(pointsA[i], pointsB[i]))
        
        x_offset.append(abs(dist_vec[0]))
        
    distortion_vec=[]
    y_offset=[]
    for i in range(len(pointsA)):
        dist_vec=tuple(np.subtract(pointsA[i], pointsB[i]))
        distortion_vec.append(dist_vec)
        y_offset.append(abs(dist_vec[1])) 
        
    
     #plot y-distortion  
    plt.figure(1)
    plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.6)
    plt.subplot(221)
    plt.scatter(distanceB, y_offset, c="red", s=10, label="$\Delta$ y")
    #plt.plot(r_line, y2_line, '--', color ='blue', label ="nonlinear regression with approximation function: $f(x)= y+p_1*(2*y**2)+2*p_2*y$ ")
    plt.title(label="Distance to reference point (px)")
    plt.xlabel('Image y_axis(px)')
    plt.ylabel('Distortion of y-values(px)')
    plt.legend()
    
    #plot x-distortion
    plt.subplot(223)
    plt.scatter(distanceB, x_offset, c="green", s=10, label="$\Delta$ x")
    #plt.plot(r_line, y3_line, '--', color ='blue', label ="nonlinear regression with approximation function: $f(x)= x+p_3*(2*x**2)+2*p_4*x$ ")
    plt.title(label="Distance to reference point (px)")  # before: Distortion in x-direction
    plt.xlabel('Image x_axis(px)')
    plt.ylabel('Distortion of x-values(px)')
    plt.legend()
    
    #radial distortion
    plt.subplot(222)
    plt.scatter(distanceB, offset, c="red", s=10, label="distorted image points")
    plt.plot(r_line, y_line, '--', color ='blue', label ="nonlinear regression with approximation function: $f(x)=k_1\cdot x+k_2\cdot x^3+k_3\cdot x^5$ ")
    plt.errorbar(distanceB, offset, xerr=u, yerr=u, fmt='.k')
    plt.title(label="Distortion Function")
    plt.xlabel('Distance to reference projection point (px)')
    plt.ylabel('Distortion (px)')
    plt.legend()
    plt.show()
    
      
    #list of radial dist vectors in same order as points
    radial_dist_vec=[]  
    
    ###Get tangential_parameter###
    for star in pointsB:
        star2center_vector = np.subtract(star, img_center)
        distance2center=np.linalg.norm(star2center_vector)

        radial_distortion=k_1*distance2center +  (k_2 * distance2center ** 3)+(k_3 * distance2center ** 5) 
        
        if distance2center!= 0:
            ratio = radial_distortion / distance2center
        else:
            ratio=0
            
        radial_distortion_vector = star2center_vector +(star2center_vector* ratio)
        radial_dist_vec.append(radial_distortion_vector)
    
        
    tangential_dist_vec=[]
    for i in range(len(distortion_vec)):
        tangential_distortion_vec=distortion_vec[i]-radial_dist_vec[i]
        tangential_dist_vec.append(tangential_distortion_vec)
        
    
    param_tangential= get_tangential_params(distanceB, pointsB, tangential_dist_vec)
    
    #TODO: Uncertainty Measurement of radial parameters
    
    return param_radial, param_cov , param_tangential #, param_cov_tang



#///   Settings   ///#

textfile= r'C:\Users\annar\Informatik\Standard_Frame\Rest\Test4.txt'
wcsfile= r'C:\Users\annar\Informatik\Standard_Frame\Rest\wcs_4.fits'

#get imagepath and position data from textfile 
manualInput=Input(textfile)   # contains textfile data
print(manualInput.source)
im = cv2.imread(manualInput.source)
rescale_window=2000 #width


#infos from astropy calibration
filename = get_pkg_data_filename(wcsfile)  
hdu = fits.open(filename)[0]  
w = WCS(hdu.header)  #sets WCS for transformations

#equatorial lat and lon  #distance=x*u.kpc: assigns distance to origin frame: zenith
center_RaDec=  SkyCoord(ra=74.460*u.degree, dec= 48.231*u.degree, frame='icrs', distance=0*u.kpc) 
#zenith_skyframe= center_RaDec.skyoffset_frame(center_RaDec)
arcsec_per_px= 51.3 



#///   Instantiate Star Catalog   ///#
# open database
ra = open(r'C:\Users\annar\Informatik\Niels\Final_Codes\Final\star_database\hygfull (2) (1).txt', 'r')
dec = open(r'C:\Users\annar\Informatik\Niels\Final_Codes\Final\star_database\hygfull (1) (1).txt', 'r')
mag = open(r'C:\Users\annar\Informatik\Niels\Final_Codes\Final\star_database\hygfull (2) (2).txt', 'r')

magnitude_threshold=6

right_ascension = ra.readlines()
declination = dec.readlines()

#create instances of CatalogStar from database
catalog_stars=[]  #list of SkyCoord objecta
for i, line in enumerate(mag):

            mag_float = float(line)
            if mag_float < magnitude_threshold:
                
                ra_temp = float(right_ascension[i]) * 15   #conversion from hourangle to degree
                dec_temp = float(declination[i])

                CatalogStar= SkyCoord(ra= ra_temp*u.degree, dec= dec_temp*u.degree, frame='icrs')
                #new_CatalogStar = CatalogStar(ra_temp, dec_temp, mag_float)  #, zenith_skyframe
                catalog_stars.append(CatalogStar)
                

#///   Image Processing   ///#

imS = ResizeWithAspectRatio(im, rescale_window)
imS_frame=np.array([(0,0), (0, imS.shape[1]), (imS.shape[1], imS.shape[0]), (imS.shape[1], 0)])

img_gray = cv2.cvtColor(imS, cv2.COLOR_BGR2GRAY)
cv2.imshow('gray', img_gray)
cv2.waitKey()


# Gaussian Blurring
# Again, you can change the kernel size
gausBlur = cv2.GaussianBlur(img_gray, (3,3),0) 


# Bilateral Filtering- better for star background, bad for line detection
bilFilter = cv2.bilateralFilter(img_gray,5,50,60)
#cv2.imshow('Bilateral Filtering', bilFilter)
#cv2.waitKey(0)
#cv2.destroyAllWindows()


kernel1 = np.array([[-1, -1, -1], [-1, 8, -1], [-1, -1, -1]]) 
sharp_bil= cv2.filter2D(src=bilFilter, ddepth=-1, kernel=kernel1)
sharp_Gauss= cv2.filter2D(src=gausBlur, ddepth=-1, kernel=kernel1)
#cv2.imshow('BilFilter Sharpened', sharp_bil)
#cv2.imshow('GaussBlur Sharpened', sharp_Gauss)
#cv2.waitKey(0)


#///   Star Detection   ///#

th, bw = cv2.threshold(sharp_bil, 20, 255,  cv2.THRESH_BINARY)

stars, hierarchy = cv2.findContours(bw,cv2.RETR_LIST,cv2.CHAIN_APPROX_SIMPLE)

cv2.imshow('binary', bw)
cv2.waitKey()

#display wcs reference projection point
wcs_refPoint= (int(w.wcs.crpix[0]),  int(w.wcs.crpix[1])) #central pixel point of wcs projection
res_wcs_refPoint=resize_pixel(im, imS, wcs_refPoint)
cv2.circle(imS, res_wcs_refPoint, 3, (255, 255, 255), thickness=5)

#list of detected stars as tuples, consisting of centerpoint and radius, blue
detected_stars= []  
for star in stars:
   # calculate moments for each contour
   M = cv2.moments(star)
 
   if M["m00"] != 0:
     cX = int(M["m10"] / M["m00"])
     cY = int(M["m01"] / M["m00"])
     #bigger radius
     cv2.circle(imS, (cX, cY), 3, (255, 0, 0), thickness=2)  
     detected_stars.append(((cX,cY), 3))
   else:
      cX=star[0][0][0]
      cY=star[0][0][1]
      # smaller stars, just one pixel in binary
      cv2.circle(imS, (cX,cY), 2, (255, 0, 0), thickness=1)
      detected_stars.append(((cX,cY), 2))

print(f'{len(detected_stars)} detected stars')    


#create circles around stars from manual input SkyCoord
manual_stars_eq=manualInput.SkyCoord #SkyCoord  objects, green
manual_points1=[] # NOT resized px coord from  manual eq (radec) inputs  #green ( where they should be )
Res_manual_points1=[]
for star in manual_stars_eq:
    xp, yp =skycoord_to_pixel(star,w, mode='wcs')
    center= (float(xp), float(yp))
    res_center=resize_pixel(im, imS, center) #just to display
    
    cv2.circle(imS, res_center, 3, (0, 255, 0), thickness=2)
    manual_points1.append(center)
    Res_manual_points1.append(res_center)
    

#manual Input pixel positions
manual_stars_px=manualInput.px_positions
manual_points2=[] # px coord from  manual pixel estimation  #red (slighty wrong)
Res_manual_points2=[] #resized
for point in manual_stars_px:
    res_center=resize_pixel(im, imS, point) #resized to display
    Res_manual_points2.append(res_center)

    manual_points2.append(point)
    cv2.circle(imS, res_center, 3, (0, 0, 255), thickness=2)

# draw distortion lines for visualisation of offset between catalog stars and rectified detected stars
for i in range(len(Res_manual_points1)):
    cv2.line(imS, Res_manual_points1[i], Res_manual_points2[i], (255, 255, 255), thickness=2)
    
###distortion between manual_points1 and manual_points2, not resized points for "real distortion"
#distortion after transformation of manual radec coord to pixel (uncertainty also comes from transformation)

# get radial and tangential distortion parameters and covariance in not resized image 
radial_dist_param, param_covariance, param_tangential= plot_distortion(im, manual_points1, manual_points2, w) 
print(f"Regression function coefficients: radial:{radial_dist_param}")  #, tangential: {tangential_dist_param}
print(f"Covariance of coefficients: {param_covariance}")
print(f'p1, p2: {param_tangential}')


# creates circles from catalog_stars image positions on imS 
# resize and distort them to match with the detected stars

database_stars_distorted=[] # distorted and resized pixel coordinates of catalog
for star in catalog_stars:
    if star.dec>= 0:   #only stars on northern hemisphere
        xp, yp = skycoord_to_pixel(star,w, mode='wcs')
        
        #filter stars from catalog, that are  visible in im (NaN, NaN)
        if pd.notna(xp) and pd.notna(yp):
            center= (int(xp), int(yp))
            res_center=resize_pixel(im, imS, center)
            if cv2.pointPolygonTest(imS_frame, res_center, False)==True:
                database_stars_distorted.append(res_center)
                cv2.circle(imS, res_center, 6, (64, 145, 6), thickness=1)
    
print(f' database stars in image:{len(database_stars_distorted)}')


#calc eq coord of detected stars-instantiate SkyCoord object by WCS transformation
#create list of detected stars, resized to im, cause they are detected in imS 

DetStars=[] #SkyCoord objects of detected stars
im_detected_stars=[]
for star in detected_stars:
    xp=star[0][0]  #x coordinate
    yp=star[0][1]  #y coordinate
    px=(xp, yp)
    rev_px= resize_pixel(imS, im, px)  
    im_detected_stars.append(rev_px)
    DetStars.append(pixel_to_skycoord(rev_px[0], rev_px[1], w))  
 
#sort detected_stars and database stars by x,y coordinates:
detected_centers=[]
for star in detected_stars:
    center=star[0]
    detected_centers.append(center)


#test if catalog stars and detected stars can be matched. 
#problem: more detected stars than catalog stars in the image and different order
#TODO: implement filter and comparing system
'''
catalog=concatenate(catalog_stars)
detected=concatenate(DetStars)
max_sep = 5.0 * u.degree
idx, d2d, d3d = detected.match_to_catalog_3d(catalog)
sep_constraint = d2d < max_sep
d_matches = detected[sep_constraint]
catalog_matches = catalog[idx[sep_constraint]]
print(f'number of  matches in detected star list = {len(d_matches)}')
print(f'number of catalog matches = {len(catalog_matches)}')
print(catalog_matches[0])
print(d_matches[0])


for i in range(len(d_matches)):   
    px_det = skycoord_to_pixel(d_matches[i], w, mode='wcs')
    px_cat = skycoord_to_pixel(catalog_matches[i], w, mode='wcs')
    px_det=(int(px_det[0]), int(px_det[1]))
    px_cat=(int(px_cat[0]), int(px_cat[1]))
    res_det=resize_pixel(im, imS, px_det)
    res_cat=resize_pixel(im, imS, px_cat)
    lines= cv2.line(imS, res_det, res_cat, (255, 255, 255), thickness=2)# draw distortion lines for visualisation
print(len(lines))

 # display the image with lines between matched stars
cv2.imshow("Detected Stars", imS)
cv2.waitKey(0)
'''


