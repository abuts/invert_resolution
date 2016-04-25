""" Sample MARI reduction scrip used in testing ReductionWrapper """ 
from ReductionWrapper import *

from ISIS_MERLINReduction import ReduceMERLIN
from ISIS_MariReduction import ReduceMARIFromFile,ReduceMARIFromWorkspace,ReduceMARIMon2Norm,MARIReductionSum
from ISIS_LETReduction import ReduceLET_OneRep,ReduceLET_MultiRep2014


if __name__=="__main__":
     maps_dir = 'd:/Data/MantidSystemTests/Data'
     data_dir ='d:/Data/Mantid_Testing/14_11_27'
     ref_data_dir = 'd:/Data/MantidSystemTests/SystemTests/AnalysisTests/ReferenceResults' 
     config.setDataSearchDirs('{0};{1};{2}'.format(data_dir,maps_dir,ref_data_dir))
     #config.appendDataSearchDir('d:/Data/Mantid_GIT/Test/AutoTestData')
     config['defaultsave.directory'] = data_dir # folder to save resulting spe/nxspe files. Defaults are in

     # execute stuff from Mantid
     #rd = MARIReductionSum();
     rd = ReduceLET_MultiRep2014()
     rd.def_advanced_properties()
     rd.def_main_properties()


     using_web_data = False;
     if not using_web_data:
        run_dir=os.path.dirname(os.path.realpath(__file__))
        file = os.path.join(run_dir,'reduce_vars.py');
        rd.export_changed_values(file);

     rd.main(); 


