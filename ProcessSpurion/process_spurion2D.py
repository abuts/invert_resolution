import math
import os
#os.environ["PATH"] = r"c:\Mantid\Code\builds\br_master\bin\Release;"+os.environ["PATH"]
from mantid.simpleapi import *
from mantid import config,api

def filter_spurion(ws,Emin,Emax,filterAccuracy,filterLimit,visualize_result=True,save_filter_log=False):
    """ Filter workspace from spurion appearing in specified energy range
        Inputs:
        ws -- the name or pointer to workspace to filter
        Emin -- minimal (elastic) energy range to fitler within
        Emax -- maximal (elastic) energy range to fitler within
        filterLimit -- the intensity limit, after exceeding which the frame would be considered a spurion
        number of time steps to divide experiment into to identify spurion

        Returns:
        workspace filtered according to filtering rules, e.g. the workspace without spurion
    """
    default_save_dir = config['defaultsave.directory']
    # get workspace to process
    if isinstance(ws,str):
        ws_name = ws
        ws2proc = mtd[ws]
    else: # Assume we got the workspace pointer
        ws_name = ws.name()
        ws2proc = ws
    #end
    # Get time of experiment
    t_start = 0
    t_end = int(ws2proc.getRun().getLogData('duration').value)
    # define energy limits we are interested in
    Units = 'Energy'
    Xmin = Emin
    Xmax = Emax
    stepX = abs(Emax - Emin) / 199.
    # Define number of time steps one needs to do to integrate data
    # due to some Mantid oddities number of time-steps we are going to get is defined as the
    # square of this value
    core = int(math.sqrt(filterAccuracy)) + 1
    if save_filter_log:
        file_to_save = os.path.join(default_save_dir,ws_name + '_spurion_time.txt')
        fs = open(file_to_save,'w')
    #
    if visualize_result:
        # Create workspace to store spectra data
        sws = CreateSampleWorkspace(WorkspaceType='Histogram',NumBanks=1,\
            BankPixelWidth=core + 1,XUnit=Units ,XMin=Xmin,XMax=Xmax,BinWidth=stepX ,NumEvents=0)

    # define time-filtering steps
    n_steps = core * core
    step = (t_end - t_start) / (n_steps)
    # define experiment's time scale we are going to loop through
    x = range(t_start,t_end,step)
    filter_ws = SumSpectra(ws2proc)

    # variables to identify filtering limits
    inSpurion = False
    signal_time_range = []
    signal_time_start = t_start
    signal_time_end = t_start
    frame_time_end = 0
    nSpurions = 0
    current_spurion = 'none'
    # loop through experiment time to see when spurion occurs
    for ind,bin in enumerate(x):
        # get time chunk
        frame_time_end = bin + step
        chunk_ws = FilterByTime(filter_ws,StartTime=float(bin),StopTime=float(frame_time_end))
        try:
            chunk_ws = NormaliseByCurrent(chunk_ws)
        except RuntimeError as rer:
            filter_result = "frame: {0:5}; time: {1:10}sec; ---- proton charge is zero. Frame ignored".format(ind,bin)
            print filter_result  
            if save_filter_log:
                fs.write(filter_result + '\n')
            continue
            
        # get current, the accelerator was producing during this time
        # for simplicity, convert our signal received during the specified time interval into energy
        chunk_wsE = ConvertUnits(chunk_ws,Units,EMode="Elastic") 
        # calculate integrated neutron counts, counted during this time frame
        ws_value = Rebin(chunk_wsE,Params='{0},{1},{2}'.format(Xmin,(Xmax - Xmin) * 1.0001,Xmax),PreserveEvents=False)
        #------------------------------------------------------------------------------
        if visualize_result:    
            # calculate spectral density to see where spurion occurs
            t1 = Rebin(chunk_wsE,Params='{0},{1},{2}'.format(Xmin,stepX,Xmax),PreserveEvents=False) 
            # store this spectral density in 2D workspace for future analysis & visual control
            t_y = t1.readY(0)
            y_val = sws.dataY(ind)
            y_len = len(y_val)
            for y_in in xrange(0,y_len):       
                y_val[y_in] = t_y[y_in]
            #y +=t_y.tolist()
        #------------------------------------------------------------------------------
        total_signal = ws_value.readY(0)
        if  total_signal > filterLimit:
            is_signal = '--- spurion'
        else:
            is_signal = '**** signal'
        filter_result = "frame: {0:5}; time: {1:10}sec; total signal {2:>10.2f}; {3}".format(ind,bin,total_signal[0],is_signal)
        print filter_result
        if save_filter_log:
             fs.write(filter_result + '\n')
        #------------------------------------------------------------------------------
        if total_signal > filterLimit:
            if visualize_result:
                 if current_spurion in mtd:
                     DeleteWorkspace(current_spurion)
                 current_spurion = 'spurion_{0}'.format(nSpurions)
                 RenameWorkspace(InputWorkspace=t1,OutputWorkspace=current_spurion)
                 plot_spectrum(current_spurion ,0)
            nSpurions+=1
            if not inSpurion:
                signal_time_end = bin
            supurion_time_end = bin + step
            inSpurion = True
            continue
        if inSpurion:
          signal_time_range.append((signal_time_start,signal_time_end))
          signal_time_start = supurion_time_end
          inSpurion = False
    # END CHUNK LOOP.  Close last frame
    if inSpurion:
        signal_time_range.append((signal_time_start,signal_time_end))
    else:
        signal_time_range.append((signal_time_start,frame_time_end))
    #
    # FINALIZE CALUALATIONS
    if visualize_result:
        # plot final 2D signal to see spurion and signal area
        plotSlice(sws)
        pass
    RunNumber = ws2proc.getRunNumber()
    if nSpurions > 0:
       for time_range in signal_time_range:
           filtered_chunk = FilterByTime(ws2proc,StartTime=float(time_range[0]),StopTime=float(time_range[1]))
           if 'filtered_ws' in mtd:
                filtered_ws +=filtered_chunk
           else:  # Prepare output workspace
                RenameWorkspace(InputWorkspace='filtered_chunk',OutputWorkspace='filtered_ws')
                RenameWorkspace(InputWorkspace=ws_name + '_monitors',OutputWorkspace='filtered_ws_monitors')
                filtered_ws = mtd['filtered_ws']
    else:
        RenameWorkspace(InputWorkspace=ws_name,OutputWorkspace='filtered_ws')
        RenameWorkspace(InputWorkspace=ws_name + '_monitors',OutputWorkspace='filtered_ws_monitors')
    #
    if save_filter_log:
        fs.close()
        file_to_save = os.path.join(default_save_dir,ws_name + '_spurion_filter_signal_time.txt')
        fs = open(file_to_save,'w')
        fs.write('RUN{0}|'.format(RunNumber))
        for time_range in signal_time_range:
            fs.write('{0};{1}|'.format(time_range[0],time_range[1]))
        fs.write('\n')
        fs.close()
        if visualize_result:
            SaveNexus(sws,Filename=ws_name + 'spurion_filter.nxs')
        
    return mtd['filtered_ws']

if __name__ == '__main__':
        fname = r'd:\users\abuts\SVN\ISIS\Python\ProcessSpurion\LET00014493.nxs'
        ws_name = 'LET00014493'
        LoadEventNexus(Filename=fname,OutputWorkspace=ws_name,SingleBankPixelsOnly='0',LoadMonitors='1',MonitorsAsEvents='0')
        filtrered_ws = filter_spurion(ws_name,0.1,8.,40,25000.,True,True)
