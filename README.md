# **Satellite Surveillance and Tracking**: 
**Line Detection in Images and Computation of Orbital Data for satellites and space debris in Andenes, Norway**

## Table of Contents
* Introduction & Functionalities
* Supplementary Material
* Technological Setup
* Workflow
* Status


## Introduction & Functionalities 
The provided codes have been created in the context of a scientific Bachelor thesis within the field of Space Situational Awareness (SSA) and Satellite Surveillance and Tracking (SST). Night-sky images have been taken with a customary camera in a fixed setup with exposure times around 6s. Thus satellites or other orbiting objects can be seen as straight lines in the images with relatively fixed star-positions in the background.
The images have been taken from Andenes in Norway with a timestamp in Central European Time in winter (UTC+1).
The work includes a slightly modified OpenCV camera calibration file (OpenCV_Cal.py) for batch calibration of an imagestock and a code to calibrate a camera without additional calibration images (Automated_Calibration.py). The latter needs to be modified and can not calibrate a camera sufficiently yet. 
The repository also includes a code enabeling the calculation of the height, speed and the orbital parameters of a satellite or space-debris fragment (Orbit_Data.py). 
Further the overpasses of the captured orbiting object can be calculated for a specific location on Earth (Overpasses.py) and be compared to satellites and debris from a database of TLE files.

## Supplementary Material
 In the folder "Plots_and_Params" supplementary material is stored, which has been worked with in this study. It contains velocity over height plots and altitude plots showing the change in altitude of detected objects, and their calculated orbital parameters.
The folder "Img_Processing" contains example images to the example provided in the thesis, showing the image processing steps that have been preceeding the calculations.
The "wcs_files" folder contains the world-coordinate-system fits files generated with the Astrometry engine (https://nova.astrometry.net/) for image sequences of the same satellite captured on the 23rd and 27th of March, 2023 at different times.
"Calibration" contains checkerboard-images used for the calibration of the camera Canon G1X as well as processed calibration images with detected corners for display.
With he files in "Auto_Cal" the Automated Calibration code can be run exemplarily.

## Technological Setup
* Python 3.9
* additional packages for astrometry:
  * astropy
  * skyfield
* (IDE: Spyder)
* Any digital camera can be used as long as the lens properties and therefore the distortion characteristics do not change
(not usable for fish-eye-lenses). With a different lens, the calibration method has to be adjusted.

## Workflow
In this section, the order in which the files need to be run to get all mentioned results is given. 

0. Automated_Calibration.py:
  Uses class of ManualInputs in Manual_Input.py to calculate the distortion parameters of a camera as the offset between detected           background stars and reprojected stars from the HYG database.
  To run the code, one has to use the files in the folder "Auto_Cal" or create equally structured files.
  * ***Example of use***:
      To test the automated calibration code for the Canon EOS 1300D, download Andoya2.txt, Astrometry_Test.csv, wcs.fits and Test.jpg from "Auto_Cal"  and set the right paths to the files in the code. Change also the path to Test.jpg and Astrometry_Test.csv in the downloaded  textfile (line 0 and line 8). 
      For usage for another camera, the text- and csv file need to be created with the same structure.
  
1. OpenCV_Cal.py:
   Open CV batch calibration. Checkerboard images are needed to get accurate results. The more calibration images of the checkerboard-pattern, the better is the reprojection error (see https://docs.opencv.org/4.x/dc/dbb/tutorial_py_calibration.html).
Here, a pattern with 9 x 6 squares has been used ( can be found here: https://github.com/opencv/opencv/blob/4.x/doc/pattern.png)
    
2. Line_Detection.py:
  The line detection algorithm can be run seperately from Orbit_Data.py, if one just wants to detect the lines in an image.
  Works by selecting a cropped part of the image. Detection method parameters are then manually adjustable within a GUI until the desired   outcome is reached, the line is detected as accurate as possible.

3. Orbit_Data.py
   Computes magnitude of velocity and height above surface of the detected object. Line Detection of start- and endpoints can be accessed    from this code or has to be inserted manually after running Line_Detection.py in a seperate step.
    This script makes use of Orbit_Parameter to calcualte the orbital parameters and do the corresponding uncertainty calculations;
    it uses Sky_Coord_Transform.py to make coordinate transformations
    Make sure to adjust the paths and EarthLocation mentioned in the settings in the code-header and to install all necessary astrometry      packages.

5. Overpasses.py
   Works seperately from Orbit_Data.py, but one needs to know the satellites start- and endpoints.
   So first run Line_Detection.py.
   Tries to match over 4000 TLE files to the detected overpass in a polar plot with the observer location as origin.
   Make sure to look into the settings and especially do not forget to change the timespan for the overpass, and the UTC time according      to your observer location. Also SkyCoord_Transform.py needs to be downloaded in the same working directory for coordinate         transformations.

   
Code-specific settings are further described in the code-headers.

## Status
* Automated Calibration:
  Code needs to be modified as the distortion parameters do not coincide with the corresponding distortion model.
  
* Line Detection:
  Right now it performs on resized image, should change that to normal cv2 format for better quality and bigger images.
  Additionally the algorithm should be able to make an error-estimation for the pixel-deviation for each image seperately.
  This will work a lot better using uncompressed RAW images.

* Overpasses: comparison to more than 4000 TLE files is necessary

  The limitations in accuracy of the resulting velocities and heights as well as the orbital parameters need to be tested by processing RAW images.
  
