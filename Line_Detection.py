"""
@author: Anna Kirchner
Line Detection code using an API

1)select cropped picture manually, lets Houghline perform line detection
 in the cropped image and draw a line at the coordinates in the original image
2)find drawn line by second image processing and HoughLine detection / get correct coordinates from before for original scale
3)save image with drawn line in destination folder

right now it performs on resized image, should change that to normal cv2 format for better quality and bigger images
TODO: 0) make it work through a whole folder with still showing images
1) make error estimation for pixel deviation from actual start and end point and due to width of line due to compression

"""

import cv2
import tkinter as tk
import os
import numpy as np
import matplotlib.pyplot as plt

def satellite_LineDetection(rectified_images, output_folder):
    
    def show_gui_Canny():
       
       # Define your function
       def changeCanny():
           # Retrieve the input parameters
           lower = float(l_thresh.get())
           upper = float(u_thresh.get())
           apertureSize=float(aperture.get())
           
    
           # Redefine the function with the new parameters globally, so no need to return or call it again in main code!
           global edges
           edges = cv2.Canny(sharp_Gauss, lower, upper, apertureSize) 
       
           print(f" cv2.Canny(sharp_Gauss, {lower}, {upper}, {apertureSize} )")
           
           root.destroy()
           root.quit()   #goes back to perform Canny
           
    
       def presume():
           
           #jump out of the while loop in "perform Canny" function
           global marker1
           
           marker1=False
           
           #destroy gui window
           root.destroy()
           root.quit()
           
    
    
       # Create the GUI window
       root = tk.Tk()
       root.title("Parameter Editor")
    
       # Create a label to display the function
       function_label = tk.Label(root, text=" edges= cv2.Canny(sharp_Gauss,")
       function_label.grid(row=0, column=0)
    
       # Define the input parameter widgets
       l_thresh = tk.Entry(root)
       l_thresh.grid(row=0, column=1)
       l_thresh.insert(0, 'lower threshold')
    
    
       u_thresh = tk.Entry(root)
       u_thresh.grid(row=0, column=2)
       u_thresh.insert(0, 'upper threshold')
    
       aperture =tk.Entry(root)
       aperture.grid(row=0, column=3)
       aperture.insert(0, 'aperture size')
    
       function_close = tk.Label(root, text=')')
       function_close.grid(row=0, column=4)
    
    
       redefine_button = tk.Button(root, text="Redefine", command=changeCanny)
       redefine_button.grid(row=1, column=0)
    
       presume_button = tk.Button(root, text="Presume with code", command=presume)
       presume_button.grid(row=1, column=1)
    
       # Start the main loop
       root.mainloop() 
       
     
       
    def show_gui_HoughLinesP():
    
       # Define your function
       def changeHoughLinesP():
           # Retrieve the input parameters
           thresh = int(threshold.get())
           r = int(rho.get())
           t= np.pi/int(theta.get())
           minline=int(linelength.get())
           gap= int(maxGap.get())
           
    
           # Redefine the function with the new parameters
           global lines
           lines= cv2.HoughLinesP(edges, threshold=thresh, rho=r, theta=t, minLineLength=minline, maxLineGap=gap) #add image as first parameter
    
           print(f"cv2.HoughLinesP(edges, {thresh}, {r}, {t},{ minline}, {gap}) ")
           
           root.destroy()
           root.quit()
    
       def save():
           
           #go back to performHoughLines function and get out of loop
           global marker2
           marker2=False
           
           #destroy window
           root.destroy()
           root.quit()
           
  
    
       # Create the GUI window
       root = tk.Tk()
       root.title("HoughLinesP Parameter Editor")
    
       # Create a label to display the function
       function_label = tk.Label(root, text="lines = cv2.HoughLinesP(edges")
       function_label.grid(row=0, column=0)
    
       # Define the input parameter widgets
       threshold = tk.Entry(root)
       threshold.grid(row=0, column=1)
       threshold.insert(0, 'threshold')
    
    
       rho = tk.Entry(root)
       rho.grid(row=0, column=2)
       rho.insert(1, 'rho')
    
       numpybit = tk.Label(root, text='np.pi/')
       numpybit.grid(row=0, column=3)
       
       theta =tk.Entry(root)
       theta.grid(row=0, column=4)
       theta.insert(180, 'theta')
       
       
       linelength =tk.Entry(root)
       linelength.grid(row=0, column=5)
       linelength.insert(50, 'minLineLength')
       
       maxGap=tk.Entry(root)
       maxGap.grid(row=0, column=6)
       maxGap.insert(10,'maxLineGap')
       
    
       function_close = tk.Label(root, text=')')
       function_close.grid(row=0, column=7)
    
    
       redefine_button = tk.Button(root, text="Redefine", command=changeHoughLinesP)
       redefine_button.grid(row=1, column=0)
    
       presume_button = tk.Button(root, text="Save Image", command=save)
       presume_button.grid(row=1, column=1)
    
       # Start the main loop
       root.mainloop() 
       
    
    
    def performCanny():
        global edges
        edges= cv2.Canny(sharp_Gauss,10, 30, 3) #initial Canny Edge Detection
            
        print("edges= cv2.Canny(sharp_Gauss,10, 30, 3) ")
        
        global marker1
        marker1= True
        while marker1==True:
              plt.imshow(edges)
              plt.waitforbuttonpress()
              plt.close()
              show_gui_Canny()   #changes marker to false, if presume is pressed
        
        #edges should be defined globally with new parameters and marker is set to False or marker is set to false directly while keeping old edges version 
        #return
        
    def performHoughLinesP():
        # Here, we're detecting lines with a minimum length of 50 pixels and a maximum gap between segments of 10 pixels
        global lines
        lines = cv2.HoughLinesP(edges, threshold=50, rho=1 ,theta=np.pi/180, minLineLength=20, maxLineGap=10)
        
        print("cv2.HoughLinesP(edges, threshold=50, rho=1 ,theta=np.pi/180, minLineLength=20, maxLineGap=10)")
        
        print(lines)
        
        global marker2
        marker2= True
        while marker2==True:
            
            img_copy = imS.copy()   #change to img for bigger screens
            if lines is not None and len(lines) > 0:
                # Draw the lines on the original image
                
                #check, if lines are similar and combine them, if slope is the same
                for line in lines:
                    x1, y1, x2, y2 = line[0]
                    # Add the offset of the ROI to the coordinates of the line
                    x1 += int(r[0])
                    y1 += int(r[1])
                    x2 += int(r[0])
                    y2 += int(r[1])
                    cv2.line(img_copy, (x1, y1), (x2, y2), (0, 0, 255), 2)
            
            #show image
            #plt.imshow(img_copy)
            #plt.waitforbuttonpress()
            #plt.close()
            cv2.imshow('houghLines', img_copy)
            print(lines)
            
            #redefine parameters, until presume is pressed (presume sets marker2==False)
            show_gui_HoughLinesP()  
            
        #when out of loop print finallines of final Image imS
        satpoints=[] #list of start and endpoint of the satellite
        if lines!=[]:
            for line in lines:
                x1, y1, x2, y2 = line[0]
                # Add the offset of the ROI to the coordinates of the line
                x1 += int(r[0])
                y1 += int(r[1])
                x2 += int(r[0])
                y2 += int(r[1])
                start=(x1,y1)
                end=(x2,y2)
                satpoints.append(start)
                satpoints.append(end)
                cv2.line(img, (x1, y1), (x2, y2), (0, 0, 255), 2)  #imS
            
        print(len(lines))
        #when final parameters are set , "save"  needs to be pressed to jump out of the while loop
        return img, satpoints  #imS
        
    folder= rectified_images
    output_dir=output_folder
    satellites={} #dictionary for img as key, start and endpoints of satlines as values
    
    
    for imagefile in os.listdir(folder):
        file=os.path.join(folder, imagefile)       
    
        # Load the image
        img = cv2.imread(file)          # Read image
        
        #store original image dimesions
        original_height, original_width = img.shape[:2]        
        imS = cv2.resize(img, (960, 540))        # Resize image
        
        resized_height, resized_width =imS.shape[:2]
        
        # Select ROI
        r = cv2.selectROI(imS)    #imS
        
        # Close the selection window
        cv2.destroyAllWindows()
        
        # Crop the image based on the ROI
        cropped_img = imS[int(r[1]):int(r[1]+r[3]), int(r[0]):int(r[0]+r[2])]    #imS
        
        # Show the cropped image
        plt.imshow(cropped_img)
        
        #Process cropped image to detect Line of interest
        gray_img = cv2.cvtColor(cropped_img, cv2.COLOR_BGR2GRAY)
        
        blurred = cv2.GaussianBlur(gray_img, (1, 1), 0) 
        
        #create kernel for noise reduction
        kernel_size = 3
        sigma = 2
        kernel = cv2.getGaussianKernel(kernel_size, sigma)
        
        print(kernel)
        
        #bilateral noise filter
        #bilFilter = cv2.bilateralFilter(gray,5,50,60)c
        
        sharp_Gauss= cv2.filter2D(src=blurred, ddepth=-1, kernel=kernel)
        
        plt.figure(1)
        plt.imshow(sharp_Gauss)
        plt.waitforbuttonpress()
        plt.close()
        
        
        performCanny() #makes user able to redefine edges (image)
        
        
        # Perform HoughLinesP line detection on the cropped image, returns image with drawn satellite lines and lines
        finalimage, finallines= performHoughLinesP()
        
        
        
        #2)convert to radec/coordinates for polar plot
        
        # Save the image with the final lines to the destination folder
        #TODO: define if else statement to check algorithm
        filename=os.path.basename(file)
        
        #finallines: satellite start and enpoints ...
        #1)save them in dictionary with the according image
        satellites[filename]= finallines
        print(satellites[filename])
        
        start_points_x=[]
        start_points_y=[]
        end_points_x=[]
        end_points_y=[]
        if len(satellites[filename])>1:
            for i in range(len(satellites[filename])-1):
                if i%2==0:
                    start_points_x.append(satellites[filename][i][0])
                    end_points_x.append(satellites[filename][i+1][0])
                    start_points_y.append(satellites[filename][i][1])
                    end_points_y.append(satellites[filename][i+1][1])
    
                
                
            
        start_x_mean=np.mean(start_points_x)
        end_x_mean=np.mean(end_points_x)
        start_y_mean=np.mean(start_points_y)
        end_y_mean=np.mean(end_points_y)
        
        start_x_original = int(start_x_mean * (original_width / resized_width))
        start_y_original = int(start_y_mean * (original_height / resized_height))
        end_x_original = int(end_x_mean * (original_width / resized_width))
        end_y_original = int(end_y_mean * (original_height / resized_height))

        
        start=(start_x_original, start_y_original)
        end=(end_x_original, end_y_original)
        
        satellites[filename]=[start, end]
        
        dst_path=os.path.join(output_dir, filename)
        cv2.imwrite(dst_path, finalimage)
        print("saved")   
    
  
    print(satellites)
    return satellites
