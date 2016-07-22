"""
sample Direct inelastic reduction for MERLIN performed in absolute units
"""
from qtiGenie import *
#from PySlice2 import *

# Calculates 

inst='MER'
iliad_setup(inst)

#mapfile = 'rings_123.map'

# where to save resutls (usually specified in Mantid, data search directories)

config['defaultsave.directory']=r"c:\Users\wkc26243\Documents\Visual Studio 2012\Projects\reduceMERLIN" 
save_dir = config.getString('defaultsave.directory')
if len(save_dir) ==0 :
    config['defaultsave.directory']=os.getcwd()
    save_dir = config.getString('defaultsave.directory')
    
#print "Data will be saved into: ",save_dir
# map mask and cal file, again the values from Mantid, data search directories can be modified here
#config.appendDataSearchDir('/home/merlin/mprogs/InstrumentFiles/merlin') 
config.appendDataSearchDir(r'c:\Users\wkc26243\Documents\work\Libisis\InstrumentFiles\merlin');
# data (raw or nxs) run files -- values from data search directories can be modified here
#config.appendDataSearchDir('/isisdatar55/NDXMERLIN/Instrument/data/cycle_13_3') 
#config.appendDataSearchDir('/archive/NDXMERLIN/Instrument/data/cycle_13_3') 
config.appendDataSearchDir(r'd:\Data\Mantid_Testing\14_04_28') 

#load vanadium file    
whitebeamfile="14338"
LoadRaw(Filename=whitebeamfile,OutputWorkspace="wb_wksp",LoadLogFiles="0")
MonoVanWB="wb_wksp"

# Mandatory positional parameters 
ei=15
rebin_params=[-20, 0.1, 14]
MonoVanRun=[14746]
monows=LoadRaw(Filename=str(MonoVanRun[0]),OutputWorkspace="monovan_wksp",LoadLogFiles="0",LoadMonitors='Include')
#RenameWorkspace(InputWorkspace="monovan_wksp_Monitors",OutputWorkspace="monovan_wksp_monitors")


#params['monovan_mapfile']='rings_125.map' 
#   Other positional  parameters
mapfile='rings_125' # ring map file is used for powder.  if absend idf file value is used instead
# key-coded parameters
params={}
params['vanadium-mass']=32.62   #   7.85
params['monovan_integr_range']=[-15,13.5]
params['norm_method']='current'
params['det_cal_file']='det_corr_125.dat'  #det_cal_file must be specified if the reduction sends out put to a workpsace
params['sample_mass']=3.56   #CePdSi3
params['sample_rmm']=330.750  #CePdSi3
params['monovan_mapfile']='rings_125.map' # good 
#params['hardmaskPlus']='mask6.msk'
#params['hardmaskOnly']='mask6.msk'
#following two lines need to be uncomented tohad  maks work ok with "False"
#params['use_hard_mask_only']=False
#params['hard_mask_file']='mask6v2.msk'
#params['hardmaskOnly']='' # no mask
params['hardmaskOnly']='Bjorn_mask.msk' 
#CePdSi3 15meV runs
runs=[14695, 14696, 14697,14698,14799,14702]
############## normal reduction####################

#save .nxspe file
for runfile in runs:
    save_file=inst+str(runfile)+'_abs.spe'
    LoadRaw(Filename=str(runfile),OutputWorkspace="run_wksp",LoadLogFiles="0")

  #  w1=iliad("wb_wksp","run_wksp",ei,rebin_params,mapfile,**params)
    
    w1=iliad("wb_wksp","run_wksp",ei,rebin_params,mapfile,monows,**params)

    SaveNXSPE('w1',save_file)
print "All done"
