from mantid.simpleapi import *
from mantid import config


def print_WSDifference(pTWS1,pTWS2,nRows):
       ''' Method to print difference between two table workspaces before and after applying CentroidPeaks '''
       # columns to compare
       tab_names=['RunNumber','DetID','Energy','DSpacing','QLab','QSample']      
       common = tab_names[0:2];
       long = tab_names[-2:];


       for name in tab_names:
           if name in common :
                print "| {0:>10} ".format(name),
           else:
               if name in long:
                    print "|FindPeaksMD found (old):{0:>7} |IntegrEllipsoids (new): {0:>7} ".format(name),
               else:
                   ntp = name;
                   if len(ntp )>6:
                       ntp = ntp[0:6]
                   print "| old {0:>6} | new {0:>6} ".format(ntp),
       print "|\n",
   
       for i in xrange(0,nRows):
           for name in tab_names:
                 col1 = pTWS1.column(name);
                 data1_toPr=col1[i]
                 col2 = pTWS2.column(name);
                 data2_toPr=col2[i]
                 if name in common :
                      print "| {0:>10} ".format(data1_toPr),
                 else:
                      if name in long:
                         print "| {0:>30} | {1:>30} ".format(data1_toPr,data2_toPr),
                      else:
                         print "| {0:>10.2f} | {1:>10.2f} ".format(data1_toPr,data2_toPr),

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
peaks=FindPeaksMD(InputWorkspace='TOPAZ_3132_md',PeakDistanceThreshold='0.37680',MaxPeaks='50',DensityThresholdFactor='100',OutputWorkspace='TOPAZ_3132_peaks')   


peaks2=CentroidPeaksMD(InputWorkspace='TOPAZ_3132_md', PeaksWorkspace='TOPAZ_3132_peaks', OutputWorkspace='TOPAZ_3132_peaks2')

print_WSDifference(peaks,peaks2,10)

    
    

