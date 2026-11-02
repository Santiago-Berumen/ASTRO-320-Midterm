"""
For Astro 320

You break it you buy it

IF YOU MAKE CHANGES, SAVE YOUR OWN COPY

THIS VERSION IS FOR TAKING THE STANDARD DEVIATIONS OF A SMALL SQUARE REGION AS OPPOSED TO THE WHOLE BIAS FRAME.

Run this in a unix environment with python 3.6 and the below packages installed by entering 'python A320_bf_stdev_ROI.py'
"""

import numpy as np
import os
from astropy.io import fits



filepath = input(r'Enter the filepath (include the final /): ')
pixel_x = input('Enter the x value for the top left pixel of the region you want to search: ')
pixel_y = input('Enter the y value for the top left pixel of the region you want to search: ')
search_size = input('Enter how much you want to search from that point (ie 5 for a 5x5 pixel region, or 25 for a 25x25 region, etc...): ')

#MAKE SURE YOUR REGION FITS IN THE FRAME

pixel_x = int(pixel_x) #Turns the strings entered above into ints
pixel_y = int(pixel_y)
search_size = int(search_size)

files_in_dir = os.listdir(filepath) #this puts every file in the directory given by filepath into a list

bias_in_dir = []
for c in range(len(files_in_dir)): #this keeps only the files with 'bias', 'Bias', and the number '0' in it
    if 'bf' in files_in_dir[c] and ('.fit' in files_in_dir[c] or '.fits' in files_in_dir[c]):
        if '0' in files_in_dir[c]:
            bias_in_dir.append(files_in_dir[c])
"""
If for some reason you saved your bias frames without the word bias in the filename, change line 24 to the following:

if 'yourweirdfilename' in files_in_dir[c]:

MAKE SURE NOTHING ELSE IN THE FOLDER HAS 'yourweirdfilename' IN IT ANYWHERE, OTHERWISE IT WILL BE ADDED TO THE LIST AND MAY
SCREW UP YOUR DATA. This can be avoided by making a folder with only your bias frames in it.

If for some reason you used a different name for each bias frame, you get no sympathy. If you numbered your filenames in a format different from 00nn, then remove ' and ('.fit' in files_in_dir[c] or '.fits' in files_in_dir[c])' from line 24 and make sure that no other files in your directory have the words 'bias' or 'Bias' in it.
"""

filelocs = []
for a in range(len(bias_in_dir)): #this combines the filepath and file names into single strings and stores them in a list
    filelocs.append(filepath + bias_in_dir[a])
    
biasframes = []
for b in range(len(filelocs)): #this opens the fits files and gets the DN values. Don't touch this part
    hdul = fits.open(filelocs[b])
    biasframe = np.array(hdul[0].data)
    biasframes.append(biasframe)

biasframes = np.array(biasframes) #this is now a 3D list, with indices [z][y][x], where z is the 'page number' i.e. each frame

stdev_frames = np.zeros(biasframes[0].shape) #creates an empty list the size of a single frame to store the standard deviations

stdev_frames = [] #this gets the standard deviations of each pixel in your ROI
for d in range(search_size):
    for e in range(search_size):
        x = pixel_x + e
        y = pixel_y + d
        stdev_frames.append(np.std(biasframes[:,y,x]))


print(np.mean(stdev_frames)) #this gives the average standard deviation of the bias frames, aka the read noise
