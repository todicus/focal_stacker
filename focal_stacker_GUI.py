# Ideas:
# - Have an option to examine the list of files in the input & output folders, and process
#   just the files missing in the output (rather than or in addition to a startTP 
#   input variable).
# - Create objective parameter sets which are saved as a text file, menu allows selecting.
#   Selection could also be accomplished through parsing the folder name, so 1.01_60x
#   would auto-select the 60x objective parameters.  Eventually each tiff should probably
#   have a header tag indicating the parameters underwhich it was recorded.
# - Estimate time to finishing, or at least display how many files are left to process. If
#   you keep track of average time per stack thus far, you can estimate time remaining.

import sys
import os
import time
from PyQt4 import QtGui, QtCore
from focal_stacker import single_stacker, folder_stacker, saveWarp, date_stacker

# paths:
pathHome = os.path.expanduser("~")
inputDefault = pathHome+'/Documents/LF/Data/'

# objective definitions (these will be selected by a menu, and should probably go in an external file):
obj_10x = {"mag":10, "NA":0.3, "umPerSlice":25}
obj_20x = {"mag":20, "NA":0.5, "umPerSlice":13}
obj_40x = {"mag":40, "NA":0.8, "umPerSlice":5}
obj_60x = {"mag":60, "NA":1.0, "umPerSlice":3}
obj_63x = {"mag":63, "NA":0.9, "umPerSlice":3}

# set defaults...
curObj = obj_40x
defaultMag = curObj["mag"]
defaultNA = curObj["NA"]
defaultUm = curObj["umPerSlice"]
defaultPx = 17           # automatically set this using image processing (find number of circles in a row)
defaultNumSlices = 25


# subclass the QMainWindow...
class MainWindow(QtGui.QMainWindow):
    def __init__(self):
        # check if mac or linux, and call main window accordingly (then this script should work on either platform)
        if os.uname()[0]=='Darwin':
            QtGui.QMainWindow.__init__(self, parent = None)
        else:
            QtGui.QMainWindow.__init__(self)
        self.setGeometry(760, 440, 200, 200)
        self.setWindowTitle('Focal Stacker')
        
        self.create_widgets()
        print('--------------------------------------------------------------')
        
    #--------------------------------------------------------------------------
    # Initialize the GUI objects and layout...
    def create_widgets(self):
        # create menu...
        self.menubar = self.menuBar()
        menu = self.menubar.addMenu('&Stuff')
        entry = menu.addAction("item")
        menu.addAction(entry)
        
        # create layouts...
#        input_group     = QtGui.QGroupBox('Input')
        input_layout    = QtGui.QVBoxLayout()
        input_field     = QtGui.QHBoxLayout()
        input_spins     = QtGui.QHBoxLayout()
        warp_field      = QtGui.QHBoxLayout()
        
        output_layout   = QtGui.QVBoxLayout()
        output_field    = QtGui.QHBoxLayout()
        output_spins    = QtGui.QHBoxLayout()
        
        buttons_layout  = QtGui.QHBoxLayout()
        
        # combine layouts...
        input_layout.addLayout(input_field)
        input_layout.addLayout(input_spins)
        input_layout.addLayout(warp_field)
#        input_group.setLayout(input_layout)
        output_layout.addLayout(output_field)
        output_layout.addLayout(output_spins)
        
        div_layout      = QtGui.QVBoxLayout()
        div_layout.addLayout(input_layout)
        div_layout.addLayout(output_layout)
        div_layout.addLayout(buttons_layout)
        div_layout.addStretch(.1)   # float doesn't seem to matter
        
        central_widget = QtGui.QWidget()
        central_widget.setLayout(div_layout)
        self.setCentralWidget(central_widget)
        
        # create input file field...
        self.inputFile_label = QtGui.QLabel("Input file:")
        self.input_field = QtGui.QLineEdit(self)
        self.input_field.setText(inputDefault)          # set default path
        self.input_but = QtGui.QPushButton(self)
        self.input_but.setText('choose')
        self.connect(self.input_but
            , QtCore.SIGNAL("clicked()")
            , self.input_clicked)
        self.startTP_label = QtGui.QLabel("start TP")
        self.startTP_box = QtGui.QSpinBox(self)
        self.startTP_box.setValue(0)
        self.startTP_box.setRange(0,9999)
        
        # create input spin boxes...
        self.mag_label = QtGui.QLabel("magnification")
        self.mag_box = QtGui.QSpinBox(self)
        self.mag_box.setRange(4,100)
        self.mag_box.setValue(defaultMag)
        
        self.na_label = QtGui.QLabel("numerical aperture")
        self.na_box = QtGui.QDoubleSpinBox(self)
        self.na_box.setRange(0.05,1.5)
        self.na_box.setValue(defaultNA)
        self.na_box.setSingleStep(0.05)
        
        self.pix_label = QtGui.QLabel("pix/lenslet")
        self.pix_box = QtGui.QSpinBox(self)
        self.pix_box.setRange(1,100)
        self.pix_box.setValue(defaultPx)
        self.pix_box.setToolTip('depends on binning; on Retiga, should be 5 for 4x4, 9 for 2x2, and 17 for 1x1')
        
        self.warp_check = QtGui.QCheckBox('linear warp', self)
        #self.warp_check.toggle()    # default to checked
        self.warp_check.setToolTip('if checked, use a linear model in rectification')
        
        # create warp field...
        self.warpFile_label = QtGui.QLabel("Warp file:")
        self.warp_field = QtGui.QLineEdit(self)
        self.warpFile_but = QtGui.QPushButton(self)
        self.warpFile_but.setText('choose')
        self.connect(self.warpFile_but
            , QtCore.SIGNAL("clicked()")
            , self.warpFile_clicked)
        self.warpFile_but.setToolTip('select a .warp file to use in rectifying the LF image(s)')
        self.warpCreate_but = QtGui.QPushButton(self)
        self.warpCreate_but.setText('create new')
        self.connect(self.warpCreate_but
            , QtCore.SIGNAL("clicked()")
            , self.warpCreate_clicked)
        self.warpCreate_but.setToolTip('calculates a new .warp file from a choosen LF image')
            
        
        # create output save field...
        self.outputFile_label = QtGui.QLabel("Output file:")
        self.output_field = QtGui.QLineEdit(self)
        self.output_field.setText(self.input_field.displayText())
        #self.output_field.setText(pathHome)          # set default path to pathHome
        self.output_but = QtGui.QPushButton(self)
        self.output_but.setText('choose')
        self.connect(self.output_but
            , QtCore.SIGNAL("clicked()")
            , self.output_clicked)
        self.copy_but = QtGui.QPushButton(self)
        self.copy_but.setText('copy input')
        self.connect(self.copy_but
            , QtCore.SIGNAL("clicked()")
            , self.copy_clicked)
        
        # create output spin boxes...
        self.numSlices_label = QtGui.QLabel("number of slices")
        self.numSlices_box = QtGui.QSpinBox(self)
        self.numSlices_box.setRange(1,499)
        self.numSlices_box.setValue(defaultNumSlices)
        
        self.umSlice_label = QtGui.QLabel("um per slice")
        self.umSlice_box = QtGui.QDoubleSpinBox(self)
        self.umSlice_box.setRange(0.1,50.)
        self.umSlice_box.setValue(defaultUm)
        self.umSlice_box.setSingleStep(0.25)
        
        self.cores_label = QtGui.QLabel("num cores")
        self.cores_box = QtGui.QSpinBox(self)
        self.cores_box.setRange(1,64)
        self.cores_box.setValue(2)
        self.cores_box.setToolTip('not yet implemented')
        
        # create buttons...
        self.single_but = QtGui.QPushButton(self)
        self.single_but.setText('Single')
        self.connect(self.single_but
            , QtCore.SIGNAL("clicked()")
            , self.single_clicked)
        self.single_but.setToolTip('process a single LF image')
            
        self.seq_but = QtGui.QPushButton(self)
        self.seq_but.setText('Sequence')
        self.connect(self.seq_but
            , QtCore.SIGNAL("clicked()")
            , self.seq_clicked)
        self.seq_but.setToolTip('process a folder full of LF images using a single .warp file; currently must all have same pixels/lenslet value')
        
        self.dub_folder_but = QtGui.QPushButton(self)
        self.dub_folder_but.setText('Folders')
        self.connect(self.dub_folder_but
            , QtCore.SIGNAL("clicked()")
            , self.dub_folder_clicked)
        self.dub_folder_but.setToolTip('process a folder of folders using .warp file calculated from second LF image in each folder; currently must all have same pixels/lenslet value')
        
        self.date_but = QtGui.QPushButton(self)
        self.date_but.setText('Date')
        self.connect(self.date_but
            , QtCore.SIGNAL("clicked()")
            , self.date_clicked)
        self.date_but.setToolTip('process an entire experiment Date folder which contains folders of LF images -- not yet implemented')
        

        # add input widgets to layout...
        input_field.addWidget(self.inputFile_label)
        input_field.addWidget(self.input_field)
        input_field.addWidget(self.input_but)
        input_field.addWidget(self.startTP_label)
        input_field.addWidget(self.startTP_box)
        
        input_spins.addWidget(self.mag_label)
        input_spins.addWidget(self.mag_box)
        input_spins.addWidget(self.na_label)
        input_spins.addWidget(self.na_box)
        input_spins.addWidget(self.pix_label)
        input_spins.addWidget(self.pix_box)
        input_spins.addWidget(self.warp_check)
        
        warp_field.addWidget(self.warpFile_label)
        warp_field.addWidget(self.warp_field)
        warp_field.addWidget(self.warpFile_but)
        warp_field.addWidget(self.warpCreate_but)

        # add output widgets to layout...        
        output_field.addWidget(self.outputFile_label)
        output_field.addWidget(self.output_field)
        output_field.addWidget(self.output_but)
        output_field.addWidget(self.copy_but)
        
        output_spins.addWidget(self.numSlices_label)
        output_spins.addWidget(self.numSlices_box)
        output_spins.addWidget(self.umSlice_label)
        output_spins.addWidget(self.umSlice_box)
        output_spins.addWidget(self.cores_label)
        output_spins.addWidget(self.cores_box)

        # add buttons to layout...
        buttons_layout.addWidget(self.single_but)
        buttons_layout.addWidget(self.seq_but)
        buttons_layout.addWidget(self.dub_folder_but)
        buttons_layout.addWidget(self.date_but)
        
    #--------------------------------------------------------------------------
    # Button actions:
    #--------------------------------------------------------------------------
    
    #--------------------------------------------------------------------------
    # The 'choose' button for input folder.
    def input_clicked(self):
        input_file = QtGui.QFileDialog.getExistingDirectory(self
            , 'Select folder'
            , str(self.input_field.displayText()).strip())   # default to current display 
            												 #path, and watch out for leading/trailing whitespace
        self.input_field.setText(input_file)
        
        # try to find objective from filename...
        filename = os.path.split(str(input_file))[1]
        self.setObjective(filename)
        
        # if output field is empty, fill with this folder (minus one branch?)
        if len(self.output_field.displayText())==0:
            self.output_field.setText(input_file)
        # set output to this folder
        #self.output_field.setText(input_file)
        
    #--------------------------------------------------------------------------
    # Finds objective power from file name, and sets parameters accordingly
    def setObjective(self, filename):
        if filename.find('10x') > 0:
            curObj = obj_10x
        elif filename.find('20x') > 0:
            curObj = obj_20x
        elif filename.find('40x') > 0:
            curObj = obj_40x
        elif filename.find('60x') > 0:
            curObj = obj_60x
        elif filename.find('63x') > 0:
            curObj = obj_63x
        else:
            print('no objective found, defaulting to 40x')
            curObj = obj_40x
        
        print('setting parameters to those of the '+str(curObj["mag"])+'x objective, aka:')
        print(curObj)
        self.mag_box.setValue(curObj["mag"])
        self.na_box.setValue(curObj["NA"])
        self.umSlice_box.setValue(curObj["umPerSlice"])
    
    #--------------------------------------------------------------------------
    # The 'choose' button for output folder.
    def output_clicked(self):
        output_file = QtGui.QFileDialog.getExistingDirectory(self
            , 'Choose location for saved tif stack'
            , self.output_field.text())
        self.output_field.setText(output_file)
        # create folder if doesn't exist? right now this is done in focal_staker.py
        
    #--------------------------------------------------------------------------
    # The 'choose' button for output folder.
    def copy_clicked(self):
        self.output_field.setText(self.input_field.displayText())

    #--------------------------------------------------------------------------
    # The 'choose' button for warp file.
    def warpFile_clicked(self):
        # set default folder...
        defDir = '~/Documents/LF/'
        if not len(self.input_field.displayText())==0:
            defDir = self.input_field.displayText()
        
        # ask user for .warp file...
        warp_file = QtGui.QFileDialog.getOpenFileName(self, 'Choose warp file', defDir, '*.warp')
        self.warp_field.setText(warp_file)
        
    #--------------------------------------------------------------------------
    # The 'create' button for warp file.
    def warpCreate_clicked(self):
        print('creating a new warp file')
        #*** maybe if no tif file in input field, popup dialog
        # ask for image on which to calc warp file...
        warpImage = str(QtGui.QFileDialog.getOpenFileName(self
            , 'Choose LF image to calc warp'
            , self.input_field.displayText()    # default to input field path
            , '*.tif'))
        print('calculating warp based on: '+warpImage)
        
        # name warp file based on the LF image used to calculate...
        warpFileName = os.path.split(warpImage[:-len(os.path.splitext(warpImage)[1])]+'.warp')[1]
        warpFilePath = str(self.input_field.displayText())+'/'+warpFileName
        print('saving to '+warpFilePath)
        
        # get the warp file...
        if self.warp_check.isChecked(): # using linear model
            linearWarp = 1
        else:
            linearWarp = 0
        
        save_return_flag = saveWarp(warpImage
                                    , self.pix_box.value()
                                    , warpFilePath
                                    , int(self.warp_check.isChecked()))

        # put path into the warp field...
        self.warp_field.setText(warpFilePath)
    
    #--------------------------------------------------------------------------
    # Process a single raw LF: asks for a tiff file path w/ file browser dialog, saves to the
    #  location in the output field. ***Should have a radio selector to choose same folder / one folder up
    #  / ./stacks?  
    def single_clicked(self):
        # get single raw LF tif file path...
        imagePath = str(QtGui.QFileDialog.getOpenFileName(self
            , 'Select a single raw LF image...'
            , self.input_field.displayText()))
        self.input_field.setText(os.path.split(imagePath)[0])   # put folder containing file in input field
        print('Processing: '+imagePath)
        
        # get warp file: calculate warp based on warp field in main window...        
        if len(str(self.warp_field.displayText()))==0:
            print('no warp file selected, please choose an LF image to use in calculating warp...')
            self.warpCreate_clicked()
        else:
            print('using warp file from GUI: '+warpFilePath)
        
        warpFilePath = str(self.warp_field.displayText())


        # get path to save folder...        
        savePath = str(self.output_field.displayText())
        if len(savePath) == 0:
            print('no save path chosen in GUI, please choose now...')
            self.output_clicked()   # trigger output button, which brings up input dialog box
        else:
            # make sure save folder exists...
            if not os.path.isdir(savePath):
                print('creating focal stack directory: ' + savePath)
                os.makedirs(savePath)
        
        # process the single stack...
        single_stacker(self.mag_box.value()
            , float(self.na_box.value())
            , imagePath      # path to raw LF file
            , int(self.pix_box.value())
            , warpFilePath
            , int(self.warp_check.isChecked())
            , int(self.numSlices_box.value())
            , float(self.umSlice_box.value())
            , str(self.output_field.displayText()))    # path to save stack file
        
    #--------------------------------------------------------------------------
    # Process a folder of raw LF images:
    def seq_clicked(self):
        # get path to folder to be processed...
        folderPath = str(self.input_field.displayText())
        if len(folderPath) == 0:
            print('no folder path chosen in GUI, please choose now...')
            self.input_clicked()
        
        # get warp file... (DEPENDS ON WARP CHECKBOX)
        if self.warp_check.isChecked(): # calculate warp based on the selected folder path...
            warpFilePath = ''    
        else: # calculate warp based on warp field in main window...
            warpFilePath = str(self.warp_field.displayText())
            warp_flag = 1
            print('using warp file from GUI: '+warpFilePath)        
        
        # get path to save folder...
        savePath = str(self.output_field.displayText())
        if len(savePath) == 0:
            print('no save path chosen in GUI, please choose now...')
            self.output_clicked()   # trigger output button, which brings up input dialog box
        
        folder_stacker(folderPath
            , int(self.startTP_box.value())     # start time point (NOT YET GUI'd)
            , warpFilePath
            , int(self.warp_check.isChecked())
            , 10000 # email interval (NOT YET GUI'd, needs menu item probably)
            , self.mag_box.value()
            , float(self.na_box.value())
            , int(self.pix_box.value())
            , int(self.numSlices_box.value())
            , float(self.umSlice_box.value())
            , str(self.output_field.displayText()))    # path to save folder of processed stacks 
    
    #--------------------------------------------------------------------------
    def dub_folder_clicked(self):
        print('processing nested folders...')
#        double_folder_stacker(pathFolder, warpFilePath, emailInterval, magnification, NA, pixLens, numSlices, umPerSlice):
    
    #--------------------------------------------------------------------------
    # should this coordinate a bunch of calls to folderstack, or just call date_Stacker in focal_stacker.py? If the latter, 
    # I'd have to write a separate controller for non-GUI processing, but I could check for obj/bin parameters  here.
    # I guess there could be a check box for auto-setting parameters from the file name, rather than running all folders
    # on the same parameters. Single/folder/date could be a drop-down or toggle list, which then alters the other options.
    # Date would have auto-check.
    def date_clicked(self):
        print('processing date folder (Toddy style!)...')
        
        # decide how to find warp files...
        if self.warp_check.isChecked(): # auto-calc warp from each run folder
            warpFilePath = ''
        else:   # use GUI-selected warp file for all folders
            warpFilePath = str(self.warp_field.displayText())
        
        date_stacker(str(self.input_field.displayText())
            , warpFilePath
            , int(self.warp_check.isChecked())
            , 10000 # email interval (NOT YET GUI'd, needs menu item probably)
            , self.mag_box.value()
            , float(self.na_box.value())
            , int(self.pix_box.value())
            , int(self.numSlices_box.value())
            , float(self.umSlice_box.value())
            , str(self.output_field.displayText()))    # path to save folder of processed stacks 

#------------------------------------------------------------------------------
if __name__ == "__main__":
    # Someone is launching this directly
    
    app = QtGui.QApplication(sys.argv)  # create the QApplication
    
    # show main window...
    main = MainWindow()
    main.show()
    main.raise_()
    sys.exit(app.exec_())

# Building Instructions:
# add "main.raise_()"   
# python pyinstaller/Makespec.py --out=StackMonsterTmp StackMonster/focal_stacker_GUI.py
# add stuff to bottom of spec file
# python pyinstaller/Build.py StackMonsterTmp/focal_stacker_GUI.spec
# cp -rv StackMonsterTmp/dist/focal_stacker_GUI/* StackMonster.app/Contents/MacOS
# cp -rv /Library/Frameworks/QtGui.framework/Versions/4/Resources/qt_menu.nib StackMonster.app/Contents/Resources