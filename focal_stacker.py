# Designed to make batch commands independent of GUI.  Not currently object-oriented, which could 
# cut down on the number of variables passed from function to function; could create a focal_stacker
# class, which contained various parameters.

import os
import sys
import time
import struct
#import Tkinter
#from tkFileDialog import askopenfilename
#from tkFileDialog import askdirectory
import smtplib
import socket        # for getting system name
#from focal_stack_GUI import setObjective

toNorm = 0     	# whether to use ImageStack's normalization procedure (default 0 now)
				#  this will avoid saturation in focal stacks built from really bright
				#  light field images.
numCores = '2' 	# processor cores, as a string; MAKE INPUT VARIABLE
#linearWarp		= 0			# to use lfrectify (1) or lfrectify2 (0)

# hard-coded paths which must be specified:
pathHome = os.path.expanduser("~")+'/'     # default to user's home dir
pathTemp = pathHome+'Documents/Code/tmp/'  # holds working files; ***could use python's mkdtemp?

# all files are in the executing folder, which are thus independent of the individual machine:
#*** could allow path selection in a GUI menu?
pathFocalProducer = 'produce_focal.py'
pathMacro = 'Build_Stack_Cmd.txt'
pathImageJ = 'ij.jar'
if os.uname()[0]=='Darwin': # osx paths
    pathImageStack = './ImageStack'
else:    # linux paths
    pathImageStack = './ImageStack_linux'


#-----------------------------------------------------------------------------
def single_stacker(magnification, NA, curTimePtPath, pixLens, warpFilePath, linearWarp, numSlices, umPerSlice, pathFocal):
    curTimePt = os.path.split(curTimePtPath)[1]    # get name of raw data

    # get start time...
    timePtStartTime = time.time()

    # to normalize at focal stack creation, which can avoid clipping in really bright LF images; not used currently
    if toNorm:
        normCMD = '-normalize'
        print('normalizing each stack to the brightest value in the stack')
    else:
        normCMD = ''

    # ----- image processing begins here -----

    # apply rectification...
    applyWarp(curTimePtPath, warpFilePath, pixLens, linearWarp)

    # get number of lenslets in x & y...
    frames, width, height, channels = struct.unpack('<IIII',open(pathTemp+'rectified.tmp','rb').read(16))
    print('number of lenslets x: '+str(width)+', y: '+str(height))

    # get focal stack temp file...
    print('converting '+curTimePt+' to a focal stack...')
    
    # send the following line to produce_focal.py...
    # *** should probably be in the objective parameter data file, which I should put here in focal_stacker (not GUI)
    arrayFocalLength = 0.
    if magnification == 40. or magnification == 20.:
        arrayFocalLength = 2500.
    elif magnification >= 60.:
        arrayFocalLength = 3750.
    
    focalCmd = 'python \''+pathFocalProducer+'\' '+numCores+' '+str(width)+' '+str(height)+' '+str(numSlices)+' '+str(umPerSlice)+' '+str(pixLens)+' '+str(magnification)+' '+str(NA)+' \''+curTimePtPath+'\' \''+pathTemp+'\' '+str(arrayFocalLength)
    #print('--- focal command: '+focalCmd)
    os.system('echo '+focalCmd+' >> \''+pathTemp+'log.txt\'')
    os.system(focalCmd+' >> \''+pathTemp+'log.txt\'')

    stackCmd = '\''+pathImageStack+'\' -load \''+pathTemp+'focal_stack.tmp\' '+normCMD+' -saveframes '+pathTemp+'%03d.tif' # can put '-normalize' before saveFrames
    #print('--- stacking command: '+stackCmd)
    os.system('echo '+stackCmd+' >> \''+pathTemp+'log.txt\'')
    os.system(stackCmd+' >> \''+pathTemp+'log.txt\'')

    # convert slices to stack...
    saveName = curTimePt[:-len(os.path.splitext(curTimePt)[1])]+'_stk('+str(numSlices)+'x'+str(umPerSlice)+')'
    print('converting slices to a single stack called '+saveName+'...')
    imagejCmd = 'java -jar '+pathImageJ+' -batch '+pathMacro+' \''+pathTemp+'000.tif,'+saveName+','+pathFocal+'\''
    #print(imagejCmd)
    os.system('echo '+imagejCmd+' >> \''+pathTemp+'log.txt\'')
    os.system(imagejCmd+' >> \''+pathTemp+'log.txt\'')
    
    # delete previous slices...
    for slice in os.listdir(pathTemp):
        if os.path.splitext(slice)[1] in ['.tif', '.tmp']:    # remove tifs & tmps
            os.remove(pathTemp + '/' + slice)

    # get timing...
    timePtEndTime = time.time()
    timePtDiff = timePtEndTime - timePtStartTime
    print('time point '+curTimePt+' took '+'%0.2f' % timePtDiff+' seconds')
    print('------------------------------------------------------------\n')


#-----------------------------------------------------------------------------
#*** placeholder for automatic pixels/lenslet calculator function.
#def getPixPerLens():
    

#-----------------------------------------------------------------------------
# Calculate the warp file used to rectify a raw LF image.  Returns 1 if successful.
# Saves the warp file to the savePath (should have extension .warp)
def saveWarp(warpImage, pixLens, savePath, linearWarp):
    pixLens = str(pixLens)
    
    # make warp file name & path...
    warpFilePath = savePath
    
    # generate rectify command...
    if linearWarp:
        rectifyCmd = '\''+pathImageStack +'\' -load \''+warpImage+'\' -gamma 0.5 -lfrectify '+pixLens+' '+pixLens+' \''+warpFilePath+'\''
    else:
        rectifyCmd = '\''+pathImageStack +'\' -load \''+warpImage+'\' -gamma 0.5 -lfrectify2 '+pixLens+' '+pixLens+' \''+warpFilePath+'\''
    
#    print(rectifyCmd)
    print('------ Rectification Output -------------------------')
    print('saving warp file as '+warpFilePath)
    os.system('echo '+rectifyCmd+' >> \''+pathTemp+'log.txt\'')
    os.system(rectifyCmd)#+' >> \''+pathTemp+'log.txt\'')
    print('-----------------------------------------------------')
    
    # PROBLEM: doesn't warn if rectify cannot be completed due to insufficient info. Check if warp file was created, or try to get feedback from LFrectify2?
    if not os.path.exists(warpFilePath):    # doesn't work if warp file already existed, and you expected to overwrite
        print('warp could not be created, either due to insufficient image intensity or permission problems')
        return 0
    else:
        return 1

#-----------------------------------------------------------------------------
def applyWarp(curTimePtPath, warpFilePath, pixLens, linearWarp):
    print('applying rectifiying warp...')
    pixLens = str(pixLens)
    
    if linearWarp:
    	applyRectCmd = '\''+pathImageStack +'\' -load \''+curTimePtPath+'\' -lfrectify \''+warpFilePath+'\''+' -unriffle '+pixLens+' '+pixLens+' -frametiles '+pixLens+' '+pixLens+' -save \''+pathTemp+'rectified.tmp\''
    else:
    	applyRectCmd = '\''+pathImageStack +'\' -load \''+curTimePtPath+'\' -lfrectify2 \''+warpFilePath+'\''+' -unriffle '+pixLens+' '+pixLens+' -frametiles '+pixLens+' '+pixLens+' -save \''+pathTemp+'rectified.tmp\''
    
#    print(applyRectCmd)
    os.system('echo '+applyRectCmd+' >> \''+pathTemp+'log.txt\'')
    os.system(applyRectCmd+' >> \''+pathTemp+'log.txt\'')

#-----------------------------------------------------------------------------
# Processes all raw LF images in a folder.
def folder_stacker(folderPath, startTP, warpFilePath, linearWarp, emailInterval, magnification, NA, pixLens, numSlices, umPerSlice, savePath):
    print('Generating '+str(numSlices)+' slices with '+str(umPerSlice)+' um per slice at every timepoint')
    
    # find date folder...
    curFolder = os.path.basename(folderPath)
    
    # check last folder name for parameters (could use filesize for bin mode)...
    #*** could make a separate method for parsing obj/bin out of file properties
    folderName = os.path.split(folderPath)[1]
    if folderName.find('1x1') > 0:
        print('this folder is not binned (1x1), switching to 17 px/lens')
        pixLens = 17
    
    # get start time...
    boutStartTime = time.time()    
    
    # get list of raw lightfield files (time points)...
    timePtList = sorted(os.listdir(folderPath))
    
    # setup folder to hold stacks...
    if not os.path.isdir(savePath):
        print('creating focal stack directory: ' + savePath)
        os.makedirs(savePath)
    else:
        print('save path already exists...\n')
    
    #if no warp file has been passed in (''), get the warp from the 1st or 2nd time point...
    if warpFilePath == '':
    	if len(timePtList) > 1:
            rectTiff = timePtList[1] 
        else:
            rectTiff = timePtList[0]
        print('rectifying on image: '+rectTiff)
        
    	warpFilePath = folderPath+'/'+rectTiff[:-len(os.path.splitext(rectTiff)[1])]+'.warp'
        print('calculating warp \n saving as: '+warpFilePath)
        
        # get the warp file...
        warp_flag = saveWarp(folderPath+'/'+rectTiff, pixLens, warpFilePath, linearWarp)
        if warp_flag == 0:
            print 'bad warp file'
            exit
    
    # loop on time points in bout folder...
    numTPs = 0;
    for curTimePt in timePtList:
        curTimePtPath = folderPath + '/' + curTimePt
        
        # skip non-tiff files, or those timepoints before startTP
        if not os.path.splitext(curTimePt)[1] == '.tif':    # skip non-tiffs
            continue
        elif startTP > 0:    # don't worry about skipping unless a number higher than 0 has been specified
            curTimePtNum = int(curTimePt[curTimePt.find('.tif')-4:curTimePt.find('.tif')])   # TROUBLE if images not numbered
            if curTimePtNum < startTP:
                numTPs = numTPs+1
                continue    # skip time points before startTP
        #print(curTimePt)
        
        # increment number of time points...
        numTPs = numTPs+1
        
        # get focal stack for this time point...
        single_stacker(magnification
          , NA
          , curTimePtPath
          , pixLens
          , warpFilePath
          , linearWarp
          , numSlices
          , umPerSlice
          , savePath)
        
        # notification email...
#        if numTPs%emailInterval==0:
#            session = smtplib.SMTP('smtp.stanford.edu')
#            smtpresult = session.sendmail('tanders@stanford.edu', 'todicus@gmail.com', '[python] '+curTimePt+' in trial '+curFolder+' [Batch_FocalStacker on '+socket.gethostname()+']')
    
    # delete temp files... (could do this after every stack, takes too much time?)
    for file in os.listdir(pathTemp):
        if os.path.splitext(file) != '.txt':    # don't remove log file
            os.remove(pathTemp + '/' + file)
    
    # get timing...
    boutEndTime = time.time()
    boutTimeDiff = boutEndTime - boutStartTime
    print('Folder '+curFolder+' took '+'%0.1f' % (boutTimeDiff/3600)+' hours')
    
    # record the fact that this folder was stacked
    stackLog = open('stackLog.txt','a')
    stackLog.write(folderPath+'\n')
    stackLog.close()
    
    print('------------------------------------------------------------\n')

#-----------------------------------------------------------------------------
# find all the raw files which have yet to be processed, given the files in the
# stack folder. 
def findNonStackedRaws(rawFolder, stackFolder):
    rawFiles = os.listdir(rawFolder)
    stackFiles = os.listdir(stackFolder)
    


#-----------------------------------------------------------------------------
# process a folder of folders containing LF images. Could work for a Date folder with many X.XX trials
# (but not a Date folder with fishX folders) or any double-nested folder heirarchy. 
# Where to put the stacks? New dir in the dir which contains the base folder?
def double_folder_stacker(pathFolder, warpFilePath, linearWarp, emailInterval, magnification, NA, pixLens, numSlices, umPerSlice):
    #*** make sure this is a doubly-nested heirarchy...
    
    # make save folder called stacks in the same folder as the folder to be processed...
    savePath = os.path.split(pathFolder)[0]+'/stacks/'+os.path.split(pathFolder)[1]     # not sure I like this
    if not os.path.isdir(savePath):
        print('creating focal stack directory: ' + savePath)
        os.makedirs(savePath)
    else:
        print('save path already exists: '+savePath)
    
    # get list of folders...
    folderList = sorted(os.listdir(pathFolder))
        
    print('Processing the folder of folders at '+pathFolder)
    for curFolder in folderList:
        # skip non-directory files...
        if not os.path.isdir(folderList):
            continue
        
        curFolderPath = pathFolder + '/' + curFolder
        curSavePath = savePath+'/'+curFolder
        
        #*** get warp... ***
        saveWarp(imagePath, pixLens, curSavePath, linearWarp)
        
        
        # call folder (seq) stacker...
        folder_stacker(curFolderPath
          , 0
          , warpFilePath
          , linearWarp
          , emailInterval
          , magnification
          , NA
          , pixLens
          , numSlices
          , umPerSlice
          , curSavePath)

#-----------------------------------------------------------------------------
#*** special case for my experiment folder heirarchy.
def date_stacker(pathDate, warpFilePath, linearWarp, emailInterval, magnification, NA, pixLens, numSlices, umPerSlice, savePath):
    print('\n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ')
    print('          - - - - - - - - - - - - - - - - - - - - - - - - - - - -          ')
    print('                 - - - - - - - - - - - - - - - - - - - - -               ')
    print('                        - - - - - - - - - - - - - -                      ')
    print('                               - - - - - - -                           ')
    print('                                    - -                               ')
    print('                                     -                             \n')
    
    # find date folder...
    dateStart = pathDate.find('200')            # YEAR 2010 BUG!!!!!  Only finds dates in the 2000's
    if dateStart == -1:
        dateStart = pathDate.find('201')        # find 2010 dates too
    date = pathDate[dateStart:dateStart+8]
    
    # get initial time...
    dateStartTime = time.time()
        
    # loop on trial folders...
    trialList = sorted(os.listdir(pathDate))
    print(trialList)
    for curTrial in trialList:
        curTrialPath = pathDate + '/' + curTrial
        
        # skip non-directory files...
        if not os.path.isdir(curTrialPath):
            continue
        
        print('\nProcessing the folder ' + curTrial)
        
        # get save folder...
        trialSavePath = savePath+'/'+curTrial
        print('saving to '+trialSavePath)
        
        # choose image on which to calc warp file...
        if warpFilePath=='':
            print('single_stacker will calculate warpfile')
        else:
            print('using supplied warp file: '+warpFilePath)
        
        # process folder...
        folder_stacker(curTrialPath, 0, warpFilePath, linearWarp, emailInterval, magnification, NA, pixLens, numSlices, umPerSlice, trialSavePath)
    
    # get elapsed time...
    dateEndTime = time.time()
    dateDiff = dateEndTime - dateStartTime
    #session = smtplib.SMTP('smtp.stanford.edu')
    #smtpresult = session.sendmail('tanders@stanford.edu', 'todicus@gmail.com', '[python] Finished trial '+curTrial+' [Batch_FocalStacker on '+socket.gethostname()+']')
    print('<<<< date '+date +' took '+'%0.1f' % (dateDiff/3600) +' hours >>>>')
