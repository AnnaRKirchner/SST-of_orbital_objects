Satellite Surveillance and Tracking: 
Line Detection in Images and Computation of Orbital Data for satellites and space debris in Andenes, Norway

##Table of Contents
*[Introduction](#introduction)
*[Technological Setup](#tech-setup)
*[Functionalities](#functionalities)
*[Example of use](#example-of-use)
*[Status](#status)


##Introduction
The provided codes have been created in the context of a scientific Bachelor thesis regarding space situational awareness (SSA) and satellite surveillance and tracking (SST). Night-sky images have been taken with a customary camera in a fixed setup with exposure times around 6s. Thus satellites or other orbiting objects can be seen as straight lines in the images with relatively fixed star-positions in the background.
The images have been taken from Andenes in Norway with a timestamp in Central European Time in winter (UTC+1).
The work includes a slightly modified OpenCV camera calibration file (OpenCV_Cal.py) and a code to calibrate a camera without additional calibration images (Automated_Calibration.py). The latter needs to be modified and can not calibrate a camera sufficiently yet. 
It also enables the user to calculate the height, speed and the orbital parameters of a satellite or space-debris fragment (Orbit_Data.py). 
Further the overpasses of the captured orbiting object can be calculated for a specific location on Earth (Overpasses.py) and be compared to satellites and debris from a database of TLE files.

In the folder "Plots_and_Params" supplementary material is stored, which has been worked with in the study. It contains velocity-height plots and altitude plots showing the change in altitude of detected objects, and their calculated orbital parameters.
The folder "Img_Processing" contains example images to the example provided in the thesis, showing the image processing steps that have been preceeding the calculations.
"wcs_files" contains the world-coordinate-system fits files generated with the Astrometry engine (https://nova.astrometry.net/) for image sequences of the same satellite captured on the 23rd and 27th of March, 2023 at different times.
"Calibration" contains checkerboard-images used for the calibration of the camera Canon G1X as well as processed calibration images with detected corners for display.

##Technological Setup
*Python 3.9
*Anaconda:Spyder as an IDE in this work
*Any digital camera can be used as long as the lens properties and therefore the distortion characteristics does not change heavily. 
(not usable for fish-eye-lenses). With a different lens, ,the calibration method has to be adjusted.

##Functionalities
In which order the files need to be run to get all functionalitites is described in this section:
1. OpenCV_Cal.py
(2. Line_Detection.py)
3. Orbit_Data.py
4. Overpasses.py

##Example of use
##Status
Automated Calibration needs to be modified. The distortion parameters do not coincide with the corresponding distortion model. 
