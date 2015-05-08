from mantid.simpleapi import *
from mantid import config
import numpy as np
import sys
import os


maps_dir = 'C:/Users/wkc26243/Documents/work/InstrumentFiles/let/'
#data_dir ='d:/Data/isis/Let/June2013'
data_dir ='d:/Data/Mantid_Testing/14_10_17'
#data_dir =maps_dir
#data_dir = 'd:/ttData'
ref_data_dir = 'd:/Data/MantidSystemTests/SystemTests/AnalysisTests/ReferenceResults' 
config.setDataSearchDirs('{0};{1};{2}'.format(data_dir,maps_dir,ref_data_dir))
#config.appendDataSearchDir('d:/Data/Mantid_GIT/Test/AutoTestData')
config['defaultsave.directory'] = data_dir # folder to save resulting spe/nxspe files. Defaults are in



filename = 'LET0007438'
groupedFilename = filename+'rings';
Ei= 3.4
e_min = -3.4
e_max = 3.2
dE = 0.1
#
Ei= 25
e_min = -20
e_max = 20
dE = 0.1
bgRange = [15000,18000]

if not("Tgrid" in mtd):

    if not(groupedFilename in mtd):
        Load(Filename=filename+'.nxs', OutputWorkspace=filename, LoadMonitors=True)
        GroupDetectors(InputWorkspace=filename, OutputWorkspace=groupedFilename , MapFile='LET_one2one_123.map', Behaviour='Average')

    wsParent = mtd[groupedFilename];
    
    nHist = wsParent.getNumberHistograms();
    print "Parent workspace contains {0:10} histograms".format(nHist)
    # Get the energy binning correspondent to the binning produced by rebin function (not to re-implement the same function)
    ws1s = ExtractSingleSpectrum(wsParent,0);
    ws1s = ConvertUnits(ws1s,'DeltaE','Direct',Ei);
    ws1s = Rebin(ws1s,Params=[e_min,dE,e_max]);
    e_bins = ws1s.dataX(0);
    nBins =e_bins.size;

    x=[e_bins[i] for i in xrange(0,nBins)]
    y=[0 for xx in xrange(0,len(x)-1)]*nHist
    x = x*nHist
    DeleteWorkspace(ws1s);
    
    eGrid = CreateWorkspace(DataX=x,DataY=y,UnitX='DeltaE',NSpec=nHist,VerticalAxisUnit='SpectraNumber',ParentWorkspace=wsParent)
    
    Tgrid=ConvertUnits(eGrid,'TOF',Emode='Direct',EFixed=Ei)
    
else:
    Tgrid = mtd['Tgrid'];
    eGrid = mtd['eGrid'];
    nHist = Tgrid.getNumberHistograms();	
    nBins = Tgrid.dataX(0).size;

if not('Bg' in mtd):
    Bg=Rebin(InputWorkspace=groupedFilename,  Params=[bgRange[0],bgRange[1]-bgRange[0],bgRange[1]],PreserveEvents=False)
    #Bg=CalculateFlatBackground(InputWorkspace=groupedFilename, StartX=bgRange[0], EndX=bgRange[1], Mode='Mean', OutputMode='Return Background', SkipMonitors=True)
else:
    Bg = mtd['Bg']
    
# Assign constant background to the Time grid workspace, minding different time bin width
for nspec in xrange(0,nHist):	
    bg            = Bg.dataY(nspec)
    if bg[0]>0:
    bgT           = Bg.dataX(nspec)		    
    TimeScale     = Tgrid.dataX(nspec);
    # Jacobian for the unit conversion
    Jac           = (TimeScale[1:nBins]-TimeScale[0:nBins-1])*(bg[0]/(bgT[1]-bgT[0]));  
    error         = np.sqrt(Jac);
    eGrid.setY(nspec, Jac)
    eGrid.setE(nspec, error)
    #print " bg at spectra {0} equal to : {1}".format(nspec,bg[0])

        
background = eGrid;
resultEt   = ConvertUnits(groupedFilename,'DeltaE',Emode='Direct',EFixed=Ei)
result     = Rebin(InputWorkspace=resultEt, Params=[e_min,dE,e_max],PreserveEvents=False)
fr         = result-background;
#
sourceSum  = SumSpectra(result,0,nHist);
bckgrdSum  = SumSpectra(background ,0,nHist);
finSum     = SumSpectra(fr ,0,nHist);

    



