import math
import os
import datetime, time
import numpy
os.environ["PATH"] = r"c:\mprogs\MantidInstall\bin;"+os.environ["PATH"]
from mantid.simpleapi import *
from mantid import config,api

def calc_count_rate(ws,Emin,Emax,NBins,visualize_result=True,use_proton_chare_log=True,trust_good_frames=True,log_name = "block_count_rate"):
    """ Filter workspace from spurion appearing in specified energy range
        Inputs:
        ws -- the name or pointer to workspace to filter
        Emin -- minimal (elastic) energy range to fitler within
        Emax -- maximal (elastic) energy range to fitler within

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
    rf = ws2proc.getRun()
    # remove count run if it is not already deleted
    if rf.hasProperty(log_name): 
        DeleteLog(ws2proc,log_name)
    
    gfl  = rf.getLogData('good_frames')
    n_filt_steps = gfl.size()
    n_good_frames = gfl.value
    
    run_start = ws2proc.getRun().getLogData('start_time').value
    time_format = '%Y-%m-%dT%H:%M:%S'
    run_sttime = datetime.datetime.strptime(run_start,time_format)
    #t_start = time.mktime(run_sttime.timetuple())
    # define energy limits we are interested in
    Units = 'Energy'
    Xmin = Emin
    Xmax = Emax
    stepX = abs(Emax - Emin) /(NBins-1)
    # define experiment's time scale we are going to loop through as time scale good frames log is recorded on
    frame_time = gfl.times;
    t_start = frame_time[0].total_nanoseconds()/1.e+9
    to_sec = lambda x : x.total_nanoseconds()/1.e+9-t_start
    frame_time  = map(to_sec,frame_time)
    t_step = (gfl.times[1].total_nanoseconds()-gfl.times[0].total_nanoseconds())/1.e+9
    #
    if visualize_result:
        # Create workspace to store spectra data
        #
        # Y-axis will be the time of the experiment
        fax = lambda x: str(x+0.5*t_step)
        AxLab = map(fax,frame_time)
        # Create X-data
        ix = range(0,NBins-1)
        fx = lambda x: Xmin+stepX*x
        dataX = map(fx,ix)
        dataX = dataX*n_filt_steps
        # And fake Y-data
        dataY = map(lambda x: 0,ix)
        dataY = dataY*n_filt_steps
        #
        sws = CreateWorkspace(OutputWorkspace = ws_name+'_count_in_time',DataX=dataX, DataY=dataY, NSpec=n_filt_steps,UnitX="Energy",\
                              VerticalAxisValues=AxLab,VerticalAxisUnit = 'Time')

    # define time-filtering steps
    filter_ws = SumSpectra(ws2proc)
    filter_wsE = ConvertUnits(filter_ws,Units,EMode="Elastic") 


    # variables to identify filtering limits
    total_signal = 0;
    # loop through experiment time to see if/when spurion occurs
    for ind,t_bin in enumerate(frame_time[:-1]):
        t_bin_sec = t_bin
        frame_time_end = frame_time[ind+1]
        bin_time = 0.5*(t_bin_sec+frame_time_end) # time in seconds associated with the timestamp
        # convert to the time log format
        ts = run_sttime+datetime.timedelta(microseconds=int(bin_time*1e+6))
        time_string  = ts.strftime(time_format)
        n_chunks = 0

        # Get frame data
        data_present = True
        if n_good_frames[ind] == n_good_frames[ind+1] and trust_good_frames: # data have not been recorded
            frame_signal = 0
            data_present = False
        else:
            # get time chunk        
            chunk_ws = FilterByTime(filter_wsE ,StartTime=float(t_bin),StopTime=float(frame_time_end))
            if use_proton_chare_log:
                try:
                    chunk_ws = NormaliseByCurrent(chunk_ws)
                except RuntimeError as rer:
                    filter_result = "frame: {0:5}; time: {1:10}sec; ---- proton charge is zero. Frame ignored".format(ind,bin_)
                    print filter_result  
                    frame_signal = 0
                    data_present = False
                    #AddTimeSeriesLog(ws, Name=log_name, Time=time_string, Value=0.)
                    #continue            
            # get current, the accelerator was producing during this time
            if data_present:
                # calculate integrated neutron counts, counted during this time frame
                ws_value = Rebin(chunk_ws,Params='{0},{1},{2}'.format(Xmin,(Xmax - Xmin) * 1.0001,Xmax),PreserveEvents=False)
                frame_signal =  ws_value.readY(0)[0]
        #endif
        total_signal =  total_signal +frame_signal
        if data_present:
            n_chunks +=1

        AddTimeSeriesLog(ws, Name=log_name, Time=time_string, Value=frame_signal)
        #------------------------------------------------------------------------------
        if visualize_result and data_present:    
            # calculate spectral density to see where spurion occurs
            s1 = Rebin(chunk_ws,Params='{0},{1},{2}'.format(Xmin,stepX,Xmax),PreserveEvents=False) 
            # store this spectral density in 2D workspace for future analysis & visual control
            t_y = s1.readY(0)
            y_val = sws.dataY(ind)
            numpy.copyto(y_val,t_y)
        #------------------------------------------------------------------------------
    # END CHUNK LOOP.  Close last frame
    #
    # FINALIZE CALUALATIONS
    if visualize_result:
        # plot final 2D signal to see spurion and signal area
        plotSlice(sws)
        if 's1' in mtd:
            DeleteWorkspace(s1)
    DeleteWorkspace(chunk_ws)
    DeleteWorkspace(filter_ws)
    DeleteWorkspace(filter_wsE)
    if 'ws_value' in mtd:    
        DeleteWorkspace(ws_value)
    if n_chunks == 0:
        n_chunks = 1
        total_signal = 0
    # Return average signal
    return total_signal/n_chunks

if __name__ == '__main__':
        #fname = r'd:\Data\Python_proj\ProcessSpurion\MER30314.nxs'
        ws_name = 'LET00024415'
        fname = r'd:\Data\Python_proj\ProcessSpurion\LET00024415.nxs'
        #ws_name = 'LET00026647'

        LoadEventNexus(Filename=fname,OutputWorkspace=ws_name,SingleBankPixelsOnly='0',LoadMonitors='1',MonitorsAsEvents='0')
        #ws_name = 'MER30314'

        av_signal = calc_count_rate(ws_name,1.,12.,100,True,True,True)
        print "avrg_signal : ",av_signal 
