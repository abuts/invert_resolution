from mantid import config
from ReduceMARI_2016_1 import *
import time
import os
import datetime
import sys
try:
    #Note: due to the mantid-python implementation, one needs to run this 
    #script in Mantid  script window  TWICE!!!  to deploy the the changes made to MARIReduction_2016_1.py file.
    reload(sys.modules['ReduceMARI_2016_1'])
except:
    print "*** WARNING can not reload MARIReduction_2016_1 file"
    pass

# Run number and Ei
runno=21461 #'live'
sum_runs=False
ei=50

# White vanadium run number
wbvan=21334
# Default save directory
#config['defaultsave.directory'] = '/instrument/MARI/RBNumber/RB1610190' #data_dir 
data_dir = r'd:\Data\Python_proj\MantidMain\ReduceMARI'
    #config.appendDataSearchDir(map_mask_dir)
config.appendDataSearchDir(data_dir)

# Absolute normalisation parameters
#monovan=21803
#sam_mass=41.104
#sam_rmm=398.9439
monovan=0
sam_mass=0
sam_rmm=0

# Set to true to remove the constant ToF background from the data.
remove_bkg = True

# If necessary, add any sequence of reduction paramerters defined in MARIParameters.xml file 
# to the end ot the illiad string using the form: property=value 
# (e.g.:  iliad_mari(runno,ei,wbvan,monovan,sam_mass,sam_rmm,sum_runs,check_background=False)
iliad_mari(runno,ei,wbvan,monovan,sam_mass,sam_rmm,sum_runs,check_background=remove_bkg,save_file_name='tt')

# Empty can measurement from other user. Ei=50meV 350Hz
#iliad_mari(21435,ei,wbvan,monovan,sam_mass,sam_rmm,sum_runs,check_background=remove_bkg)

#Rebin2D(InputWorkspace='MAR0Reduced#49.81_SQW', OutputWorkspace='MAR0Reduced#49.81_SQW_line', Axis1Binning='4,5,9', Axis2Binning='-30,0.1,50', UseFractionalArea=True, Transpose=True)
