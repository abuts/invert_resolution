from mantid.simpleapi import *
from mantid import config


def print_tableWS(pTWS,nRows):
    ''' Method to print part of the table workspace '''
    tab_names=pTWS.keys();
    
    for name in tab_names:
        if len(name)>8:
           name= name[0:8];
        print "| {0:8} ".format(name),
    print "|\n",

    for i in xrange(0,nRows):
        for name in tab_names:
              col = pTWS.column(name);
              data2pr=col[i]
              if type(data2pr) is float:
                   print "| {0:8.3f} ".format(data2pr),
              else:
                  print "| {0:8} ".format(data2pr),   
        print "|\n",

#########################################################            
data_dir="d:/Data/MantidSystemTests/Data/";
maps_dir = "d:/Data/MantidSystemTests/Data"
ref_data_dir = 'd:/Data/MantidSystemTests/SystemTests/AnalysisTests/ReferenceResults' 
config.setDataSearchDirs('{0};{1};{2}'.format(data_dir,maps_dir,ref_data_dir))
#config.appendDataSearchDir('d:/Data/Mantid_GIT/Test/AutoTestData')
#########################################################                


# load test workspace
Load(Filename=r'TOPAZ_3132_event.nxs',OutputWorkspace='TOPAZ_3132_event',LoadMonitors='1')
   
# build peak workspace necessary for IntegrateEllipsoids algorithm to work
ConvertToMD(InputWorkspace='TOPAZ_3132_event',QDimensions='Q3D',dEAnalysisMode='Elastic',Q3DFrames='Q_sample',LorentzCorrection='1',OutputWorkspace='TOPAZ_3132_md',\
MinValues='-25,-25,-25',MaxValues='25,25,25',SplitInto='2',SplitThreshold='50',MaxRecursionDepth='13',MinRecursionDepth='7')
FindPeaksMD(InputWorkspace='TOPAZ_3132_md',PeakDistanceThreshold='0.37680000000000002',MaxPeaks='50',DensityThresholdFactor='100',OutputWorkspace='TOPAZ_3132_peaks')   
FindUBUsingFFT(PeaksWorkspace='TOPAZ_3132_peaks',MinD='3',MaxD='15',Tolerance='0.12')
IndexPeaks(PeaksWorkspace='TOPAZ_3132_peaks',Tolerance='0.12')

# integrate Ellipsoids   
result=IntegrateEllipsoids(InputWorkspace='TOPAZ_3132_event',PeaksWorkspace='TOPAZ_3132_peaks',RegionRadius='0.25',PeakSize='0.2',BackgroundInnerSize='0.2',BackgroundOuterSize='0.25',OutputWorkspace='TOPAZ_3132_peaks')

# print 10 rows of resulting table workspace
print_tableWS(result,10)
