import os
#import socket		# for getting system name
#from focal_stacker import findNonStackedRaws
from focal_stacker import saveWarp, applyWarp


# stimulation record processer
#f = '/Users/todicus/Documents/LF/Data/20100603 - ogb1/visual stimulation/'
#name = '1.01 response_tester.txt (23-12-23).txt'


rectTiff = '/Users/logangrosenick/Documents/Data/20100720_ogb1/Raw_60x/AVG_2.01.tif'
rectTiffName = os.path.split(rectTiff)[1]
rectTiffFolder = os.path.split(rectTiff)[0]
warpFileName = rectTiffName[:-len(os.path.splitext(rectTiffName)[1])]+'.warp'
print(rectTiffFolder+'  |  '+warpFileName)
saveWarp(rectTiff, 9, rectTiffFolder+'/'+warpFileName)

