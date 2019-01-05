'''
This is provided for testing purposes. It allows you to run the single system test
within Mantid environent
It can be useful for debugging because the errors do not alway 'get out' of
the sub-process used for running the tests in the regular way
'''

import sys
import os
os.environ["PATH"] = r"c:\\mprogs\\MantidNightlyInstall\\bin;c:\\mprogs\\MantidNightlyInstall\\lib\\paraview-5.4\\site-packages;c:\\mprogs\\MantidNightlyInstall\\lib\\paraview-5.4\\site-packages\\vtk"+os.environ["PATH"]
sys.path.append('c:\\mprogs\\MantidNightlyInstall\\bin')
sys.path.append('c:\\Users\\abuts\\Documents\\Developing_soft\\Mantid\\Mantid_dev\\Testing\\SystemTests\\lib\\systemtests')
sys.path.append('c:\\Users\\abuts\\Documents\\Developing_soft\\Mantid\\Mantid_dev\\Testing\\SystemTests\\tests\\analysis')
sys.path.append('c:\\Users\\abuts\\Documents\\developing_soft\\Mantid\Mantid_dev\\scripts\\Inelastic')
print sys.path
from mantid.simpleapi import *
from mantid import config
import inspect

#this_dir = sys.path[0]
#print "Script resides in : ",this_dir
#python directories
#os.chdir(r'd:\Data\MantidSystemTests\SystemTests')
#stressmodule_dir = r'd:\Data\MantidSystemTests\StressTestFramework'
#tests_dir = r'd:/Data/MantidSystemTests/SystemTests/AnalysisTests'


#sys.path.insert(0,tests_dir)
#sys.path.append(stressmodule_dir)


data_dir='c:/Users/abuts/Documents/developing_soft/Mantid/Mantid_TestData'
#data_dir2=''C:/Backup/Backup_folder1/work/code/Mantid/git/mantid/Test/systemtests/Data/LOQ'

#maps_dir = 'd:/Data/MantidSystemTests/Data'
#data_dir ='d:/Data/isis/Let/June2013'
#data_dir ='d:/Data/Mantid_Testing'
#data_dir = 'd:/ttData'
#ref_data_dir = 'd:/Data/MantidSystemTests/SystemTests/AnalysisTests/ReferenceResults' 
#config.setDataSearchDirs('{0};{1};{2}'.format(data_dir,maps_dir,ref_data_dir))
config.appendDataSearchDir(data_dir)
config['defaultsave.directory'] = data_dir # folder to save resulting spe/nxspe files. Defaults are in


MARITests = ['MARIReductionFromFile','MARIReductionAutoEi','MARIReductionFromFileCache',
             'MARIReductionFromWorkspace','MARIReductionMon2Norm','MARIReductionMonSeparate',
             'MARIReductionSum','MARIReductionWaitAndSum']
MARIComponentTests = ['ISIS_ReductionWebLike','ISIS_ReductionWrapperValidate']

# data search directories
#reference_dir = 'C:/Backup/Backup_folder1/work/code/Mantid/git/mantid/Test/systemTests/SystemTests/AnalysisTests/ReferenceResults'
#mtdpy_header_dir = 'C:/Backup/Backup_folder1/work/code/Mantid/mantid/Code/Mantid/debug'
#sys.path.append(r'C:\Backup\Backup_folder1\work\code\Mantid\builds\all\bin\Debug')

# Find these first
#modlToRun = ['POLDIDataAnalysisTest']
modlToRun = ['ISISDirectInelastic']
#modlToRun = ['DirectInelasticDiagnostic']
#modlToRun = ['ISISDirectReductionComponents']
#modlToRun = ['PowderDiffProfileCalibrateTest'];
#testToRun = ['VulcanSeqRefineProfileFromScratch'];
#modlToRun = ['DirectInelasticDiagnostic'];
#testToRun = ['DirectInelasticDiagnostic2'];
#testToRun = ['DirectInelasticDiagnosticSNS','DirectInelasticDiagnostic']
#modlToRun=['ISISIndirectAbsCorTest','ISISIndirectBayesTest','ISISIndirectInelastic','ISISIndirectLoadAsciiTest'];

#modlToRun = ['Diffraction_Workflow_Test'];
#testToRun = ['Diffraction_Workflow_Test']
#testToRun = ['MARIReductionAutoEi']
#testToRun = ['POLDIDataAnalysisTestSiIndividual']
#testToRun = ['ISISLoadFilesLET']
#testToRun = ['DirectInelasticDiagnostic2']
#testToRun = ['MARIReductionMonSeparate']
#testToRun = ['MARIReductionFromFile','MARIReductionSum']
#testToRun = ['MARIReductionFromWorkspace','MARIReductionFromFile','MARIReductionSum','MAPSDgreduceReduction','LETReduction','LETReductionEvent2015Multirep','MERLINReduction']
#testToRun = ['LETReductionEvent2015Multirep']
testToRun = ['MARIReductionFromFileCache','MARIReductionWaitAndSum']
#testToRun = ['DirectInelasticDiagnostic']
#,
#testToRun = ['ISIS_ReductionWebLike','ISIS_ReductionWrapperValidate']




for mod_name in modlToRun:

    module = __import__(mod_name )
    reload(module)
    clear_test_list = False

    if len(testToRun) == 0:
        clear_test_list = True
        for name,class_inst in inspect.getmembers(module):
            if inspect.isclass(class_inst):
                if name.endswith('Test'):
                    testToRun.append(name)
                #if hasattr(class_inst, 'runTest') and callable(getattr(class_inst, 'runTest')) and not hasattr(class_inst, "__abstractmethods__"):
                #    testToRun.append(name);
    #reload(sys.modules['isis_reduction_steps'])

    for className in testToRun:
        try:
            testClass = getattr(module, className)()
            testClass.runTest()
            outcome = testClass.doValidation()
            print 'Test result: ' + str(outcome)
        except:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print 'Test: {0} thrown exception: {1} '.format(className,exc_type)
            print 'Reason: {0}'.format(exc_value)
        #os.chdir(this_dir)
        #raise

    if clear_test_list:
        testToRun=[]

#os.chdir(this_dir)
