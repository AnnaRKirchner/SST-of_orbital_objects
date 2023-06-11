"""
Created on Wed Apr  5 17:33:36 2023
Open CV : Checkerboard Calibration

Adjustments for Rectification of input folder: 
    @author: Jessica McKenna & Anna Kirchner
    
Settings:
    -checkerboard_dir: path to folder with checkerboard images for calibration
    
    -input_dir:  path to folder with distorted images
     e.g.:f'C:\\Users\\annar\\Documents\\BA\\SkyImg\\SORTED_SATELLITES_CanonG1X\\Distorted\\Alomar_08_03_23\\{i}'
    
    -output_dir: path to output folder
     e.g.: f'C:\\Users\\annar\\Documents\\BA\\SkyImg\\SORTED_SATELLITES_CanonG1X\\Rectified\\08_03_23\\{i}'
    
"""

import numpy as np
import cv2 as cv
import glob
import os
import matplotlib.pyplot as plt

class ImageRectifier:
    def __init__(self, input_directory, output_directory, ret, mtx, dist, rvecs, tvecs):
        self.input_directory = input_directory
        self.output_directory = output_directory
        #self.calibrator = calibrator
        self.ret = ret
        self.mtx = mtx
        self.dist = dist
        self.rvecs = rvecs
        self.tvecs = tvecs
    
    def rectify_image(self, img_path, output_directory):
        # Read the image
        img = cv.imread(img_path)
        # Rectify the image
        h, w = img.shape[:2]
        newcameramtx, roi = cv.getOptimalNewCameraMatrix(self.mtx, self.dist, (w,h), 1, (w,h))
        dst = cv.undistort(img, self.mtx, self.dist, None, newcameramtx)
        x, y, w, h = roi
        dst = dst[y:y+h, x:x+w]
        mapx, mapy = cv.initUndistortRectifyMap(self.mtx, self.dist, None, newcameramtx, (w,h), 5)
        dst = cv.remap(img, mapx, mapy, cv.INTER_LINEAR)
        # crop the image
        x, y, w, h = roi
        dst = dst[y:y+h, x:x+w]
        
        filename, ext = os.path.splitext(os.path.basename(img_path))
        new_file_name = f"{filename}_rectified{ext}"
        # Saving the rectified image to a new directory
        cv.imwrite(os.path.join(self.output_directory, new_file_name), dst)
        return dst
    
    
# Define a function to stop the program after a certain amount of time
def stop_program():
    print("Time's up!")
    # Raise an exception to stop the program
    raise SystemExit

    

# termination criteria
criteria = (cv.TERM_CRITERIA_EPS + cv.TERM_CRITERIA_MAX_ITER, 30, 0.001)
#objp=(0,0,0), (1,0,0), (2,0,0) ....,(6,5,0)
objp = np.zeros((9*6,3), np.float32)
objp[:,:2] = np.mgrid[0:9,0:6].T.reshape(-1,2)
objp = objp * 0.0215
# Arrays to store object points and image points from all the images.
objpoints = [] # 3d point in real world space
imgpoints = [] # 2d points in image plane.
#fname=r'C:\Users\annar\Informatik\Standard_Frame\Calibration\IMG_8458.jpg'
#for fname in images:
# get the path/directory
checkerboard_dir = r'C:\Users\annar\Documents\BA\Calibration_Images\Try'
for file in os.listdir(checkerboard_dir):
    if (file.endswith(".jpg") or file.endswith(".JPG")):
        img_path = os.path.join(checkerboard_dir, file)
        img = cv.imread(img_path)
        gray = cv.cvtColor(img, cv.COLOR_BGR2GRAY)
        
        # Find the chess board corners
        
        ret, corners = cv.findChessboardCorners(gray, (9,6), None)
        # If found, add object points, image points (after refining them)
        if ret == True:
            objpoints.append(objp)
            corners2 = cv.cornerSubPix(gray,corners, (11,11), (-1,-1), criteria)
            imgpoints.append(corners2)
            # Draw and display the corners
            cv.drawChessboardCorners(img, (9,6), corners2, ret)
            res= cv.resize(img, (int(img.shape[1]*0.5), int(img.shape[0]*0.5)), interpolation = cv.INTER_LINEAR)
            plt.imshow(res)
            plt.waitforbuttonpress()
            plt.close()
            cv.imwrite('img_{}'.format(file), img)
            cv.waitKey(0)
        else:
            print('cant detect checkerboard corners')
            
      

ret, mtx, dist, rvecs, tvecs = cv.calibrateCamera(objpoints, imgpoints, gray.shape[::-1], None, None)
print(f'k1, k2, p1, p2, k3: {dist}, matrix:{mtx}')

#Single image to undistort
#img = cv.imread(r'C:\Users\annar\Documents\BA\SkyImg\SORTED_SATELLITES_CanonG1X\Distorted\23_03_23\9\IMG_2938.JPG')
#print(img.shape)
#h,  w = img.shape[:2]
#newcameramtx, roi = cv.getOptimalNewCameraMatrix(mtx, dist, (w,h), 1, (w,h))
#print(f'newcameramtx: {newcameramtx}')

#Batch rectification: Set the input and output directories for distorted and rectified images
for i in range(1,6):
    print(i)
    input_dir=f'C:\\Users\\annar\\Documents\\BA\\SkyImg\\SORTED_SATELLITES_CanonG1X\\Distorted\\Alomar_08_03_23\\{i}'
    output_dir = f'C:\\Users\\annar\\Documents\\BA\\SkyImg\\SORTED_SATELLITES_CanonG1X\\Rectified\\08_03_23\\{i}'
    
    # Initialize the ImageRectifier class with the required arguments
    rectifier = ImageRectifier(input_directory=input_dir, 
                               output_directory=output_dir, ret =ret, mtx=mtx, dist=dist, rvecs=rvecs, tvecs=tvecs)
    
    # Loop through all the images in the input directory and rectify them, saving them to the outputs
    images = glob.glob(os.path.join(input_dir, "*.jpg")) #f.endswith('.jpg') or f.endswith(".JPG") or f.endswith(".png")
    for img_path in images:
        image_rectifier = ImageRectifier(input_dir, output_dir, ret, mtx, dist, rvecs, tvecs)
        image_rectifier.rectify_image(img_path, output_dir)
