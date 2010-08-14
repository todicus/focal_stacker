import os
pathDate = '/Users/todicus/Documents/LF/Data/20100629/Raw_60x'
savePath = '/Users/todicus/Documents/LF/Data/20100629/stacks'

rawFolders = os.listdir(pathDate)
stackFolders = os.listdir(savePath)

# find last folder in raw which is also in stack...
lastCommonFolder = rawFolders.index(stackFolders[-1])

# check if last folder needs more processing...
lastRawFiles = os.listdir(pathDate + lastCommonFolder)
lastStackFiles = os.listdir(savePath + lastCommonFolder)

