##import DirectEnergyConversion as direct
##mono_sample = direct.DirectEnergyConversion('MARI')
##mono_sample.prop_man.log_to_mantid = True
##mono_sample.prop_man.normalise_method = 'none'
##mono_sample.prop_man.background = True
##mono_sample.prop_man.background_range = [18000, 19500]
##mono_sample.prop_man.fix_ei = False
##mono_sample.prop_man.energy_bins = '-11,0.1,11'
##mono_sample.prop_man.save_format = 'nxs'

##mono_sample.prop_man.motor_name = None
##mono_sample.prop_man.motor_offset = None
##mono_sample.prop_man.sum_runs = True
##mono_sample.prop_man.save_file_name = r'MAR11001.spe'
##mono_sample.convert_to_energy(r'd:\Data\MantidSystemTests\Data\MAR11060.raw', [r'd:\Data\MantidSystemTests\Data\MAR11001.raw'], 12, None, None, None, None)

#import DirectEnergyConversion as direct
#mono_sample = direct.DirectEnergyConversion('MARI')
#mono_sample.prop_man.log_to_mantid = True
#mono_sample.prop_man.normalise_method = 'none'
#mono_sample.prop_man.background = True
#mono_sample.prop_man.background_range = [18000, 19500]
#mono_sample.prop_man.fix_ei = False
#mono_sample.prop_man.energy_bins = '0,0.05,10'
#mono_sample.prop_man.map_file = r'd:\Data\MantidSystemTests\Data\mari_res.map'
#mono_sample.prop_man.monovan_mapfile = r'd:\Data\MantidSystemTests\Data\mari_res.map'
#mono_sample.prop_man.monovan_integr_range=[float(-0.8),float(0.8)]
#mono_sample.prop_man.sample_mass = 1
#mono_sample.prop_man.sample_rmm = 1
#mono_sample.prop_man.van_mass = 32.58
#mono_sample.prop_man.save_format = 'nxs','nxspe'

#mono_sample.prop_man.motor_name = r'mo'
#mono_sample.prop_man.motor_offset = 12
#mono_sample.prop_man.sum_runs = True
#mono_sample.prop_man.save_file_name = r'MAR11001.spe'
#mono_sample.convert_to_energy(r'd:\Data\MantidSystemTests\Data\MAR11060.raw', [r'MAR11001.raw'], 12, None, None, [r'd:\Data\MantidSystemTests\Data\MAR11015.raw'], None)



from qtiGenie import *
from PySlice2 import *
from autoEi import *
inst='mar'
iliad_setup(inst)
ext='.raw'
#mapfile='mari_res2012'

mapfile='mari_res2013'


#det_cal_file must be specified if the reduction sends out put to a workspace
#cal_file='MAR17578.raw'
cal_file='MAR18622.raw'
#load vanadium file

run='17690'
WB='18622'
ei,rebin_params,BkgdRange=autoEi(str(run),BkgdGen=True,monspecin=2,BkgdLevel=.4)
print "ei = ",ei
print "rebin_params = ",rebin_params
#print "BkgdRange = ",BkgdRange,bkgd_range = BkgdRange

w1=iliad(WB,run,ei,rebin_params,mapfile,det_cal_file=cal_file,norm_method='current',bkgd_range = [14500,18000])
