import lfrecon
import os.path
import sys

settings = lfrecon.LFReconGlobals()

settings.verbose = {'i': 'INIT',
                 'c': 'CG',
                 'n': 'NEWTON',
                 'b': 'BARRIER'}

magnification = float(sys.argv[7])
NA = float(sys.argv[8])
mediumRI = 1.3333333    # for water dipping
microscopeRI = 1.0
umPerLenslet = 125.
lensletFocalLen = float(sys.argv[11]) # for f/30 array; MAKE THIS AN INPUT
numCores = int(sys.argv[1])
numCols = int(sys.argv[2])
numRows = int(sys.argv[3])
numSlices = int(sys.argv[4])
numU = int(sys.argv[6])
numV = int(sys.argv[6])
centerZ = 0.0
spacingZ = float(sys.argv[5])
filePath = sys.argv[9]
pathTemp = sys.argv[10]

settings.threads = numCores
settings.logfile = open(pathTemp+'produce_focal.log', 'wb')
settings.focal_alg = ''

recon = lfrecon.LFRecon(numCols, numRows, numSlices, numU, numV,
                        magnification, NA, mediumRI,
                        umPerLenslet, lensletFocalLen, microscopeRI,
                        centerZ, spacingZ)


print 'Getting light field from rectified.tmp'

lf = recon.load_lf(open(pathTemp+'rectified.tmp','rb').read())
#~ lf = recon.load_lf(open(filePath,'rb').read())		# for batch processing tmp files, filePath could be input variable

print('Computing focal stack')

focalStack = recon.blank_volume()   # make blank focal stack
recon.compute_focal_stack(lf, focalStack)
open(pathTemp+'focal_stack.tmp','wb').write(recon.dump_volume(focalStack))
