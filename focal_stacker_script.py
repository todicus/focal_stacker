import os
import socket		# for getting system name
import Tkinter
from tkFileDialog import askopenfilename
from tkFileDialog import askdirectory
import struct
import sys
#import time
#import smtplib
from focal_stacker import single_stacker
from focal_stacker import folder_stacker
from focal_stacker import date_stacker
from focal_stacker import saveWarp

#workingPath = ''
workingPath = '/Users/logangrosenick/Documents/Data/20100810/Raw_20x'

#savePath = ''
savePath = '/Users/logangrosenick/Documents/Data/20100810/Raw_20x/stacks'

# warp can be specified, either as a .tif from which to calculate the .warp, or the .warp itself. If left blank (''), the warp will be calculated from the first tif in each folder.
warpFile = '/Users/logangrosenick/Documents/Data/20100810/Raw_20x/1.06_20x/1.06-minimum, bgsub.warp'
# uncomment next line to use a previously-generated warp file:
#warpFile = '/Volumes/Nautilus/20100714_gcamp3_fish_&_slice/Raw_20X_zfish/AVG_3.01.warp'

autoDetectObj = 0	# not used yet

# optical parameters
magnification = 20.
NA = 0.5

# light field image input parameters
pixLens = 9
linearWarp = 0
usePrevWarp = 0
startTP = 0

# stack output parameters
numSlices = 25
umPerSlice = 13
emailInterval = 10000

batchMode = 3	# what type of batch to perform: single (0), folder (1), date (3)

# OS parameters...
onMac 			= 1			# using mac (1) or linux (0) (my mac and my linux, in particular)
							# can I do this with a python function?
# set paths...
if onMac:
	# mac paths...
	pathInit = os.path.expanduser("~")
	#pathInit = '~/Dropbox/'
else:
	# linux paths...
	pathInit = '/media/Neuronaut/20081021/Raw LF/' #pathLF+'/'


# get path to folder full of raw LF images...
if workingPath=='':
	# ask for bout folder location with GUI...
	print 'Please choose input folder...'
	root = Tkinter.Tk()
	root.withdraw()
	workingPath = askdirectory(initialdir=pathInit)	# title='choose bout folder')

# get path to save folder for tiff stacks...
if savePath=='':
	# ask for bout folder location with GUI...
	print 'Please choose save folder...'
	root = Tkinter.Tk()
	root.withdraw()
	savePath = askdirectory(initialdir=pathInit)	# title='choose bout folder')

#if warpFile ext is .tif, calculate a .warp file
if os.path.splitext(warpFile)[1] == '.tif':
	print('calculating warp parameters from image: '+warpFile)
	warpFilePath = os.path.split(warpFile[:-len(os.path.splitext(warpFile)[1])]+'.warp')[1]
	saveWarp(warpFile, pixLens, warpFilePath, linearWarp)
	warpFile = warpFilePath
elif os.path.splitext(warpFile)[1] == '.warp':
	print('using previously-calculated warp parameters from file: '+warpFile)

print('Processing the sequence in folder: ' + workingPath)
print('Saving the sequence to folder: ' + savePath)

#---------------------------------------------------------------------------
if batchMode == 0:
	print('Processing raw LF image: ')
	single_stacker(magnification, NA, workingPath, pixLens, warpFile, numSlices, umPerSlice, savePath)
	
elif batchMode == 1:
	folder_stacker(workingPath, startTP, warpFile, linearWarp, emailInterval, magnification, NA, pixLens, numSlices, umPerSlice, savePath)
	
elif batchMode == 2:	#*** needs work
	# ask for folder of folders location with GUI...
	root = Tkinter.Tk()
	root.withdraw()
	pathFolder = askdirectory(initialdir=pathInit)	# title='choose folder folder')
	print('Processing folder: ' + pathFolder)
	folder_stacker(pathFolder, startTP, usePrevWarp, emailInterval, magnification, NA, pixLens, numSlices, umPerSlice)

elif batchMode == 3:
	pathDate = workingPath
	date_stacker(pathDate, warpFile, linearWarp, emailInterval, magnification, NA, pixLens, numSlices, umPerSlice, savePath)
