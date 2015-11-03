from mantid.simpleapi import *
from mantid import config
import sys
import os
import inspect


ws=Load(Filename=r'D:\users\abuts\SVN\Mantid\Mantid_testing\Sample_dataISIS\MAR11001.raw', OutputWorkspace='MAR11001', LoadMonitors='Separate')
RenameWorkspace('MAR11001','REM11001',RenameMonitors=True)
RenameWorkspace(InputWorkspace='REM11001',OutputWorkspace='REM02',RenameMonitors=True)
RenameWorkspace('REM02','REM03',True)
rem4=RenameWorkspace('REM03')
rem5=RenameWorkspace('rem4',True)
ws = RenameWorkspace(rem5,OutputWorkspace='rem6',RenameMonitors=True)
print "Renamed workspace ",ws.name(),' and monitor workspace: ',ws.getMonitorWorkspace().name()