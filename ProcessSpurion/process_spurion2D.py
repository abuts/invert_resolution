import math
import os
import datetime
import numpy
os.environ["PATH"] = r"c:\mprogs\MantidInstall\bin;"+os.environ["PATH"]
from mantid.simpleapi import *
from mantid import config,api

def filter_spurion(ws,Emin,Emax,n_filt_steps,visualize_result=True,use_proton_chare_log=True,log_name = "block_count_rate"):
    """ Filter workspace from spurion appearing in specified energy range
        Inputs:
        ws -- the name or pointer to workspace to filter
        Emin -- minimal (elastic) energy range to fitler within
        Emax -- maximal (elastic) energy range to fitler within
        n_filt_steps -- number of time steps to divide run into to identify spurion

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
    t_start = 0.
    t_end = ws2proc.getRun().getLogData('duration').value
    run_start = ws2proc.getRun().getLogData('start_time').value
    time_format = '%Y-%m-%dT%H:%M:%S'
    run_sttime = datetime.datetime.strptime(run_start,time_format)
    # define energy limits we are interested in
    Units = 'Energy'
    Xmin = Emin
    Xmax = Emax
    NBins = 200;
    stepX = abs(Emax - Emin) /(NBins-1)
    # define experiment's time scale we are going to loop through
    t_step = (t_end-t_start)/(n_filt_steps)
    it = range(0,n_filt_steps)
    ft = lambda x : t_start+x*t_step
    spec = map(ft,it)
    #
    if visualize_result:
        # Create workspace to store spectra data
        #
        # Y-axis will be the time of the experiment
        fax = lambda x: str(t_start+x*t_step)
        AxLab = map(fax,it)
        # Create X-data
        ix = range(0,NBins)
        fx = lambda x: Xmin+stepX*x
        dataX = map(fx,ix)
        dataX = dataX*n_filt_steps
        # And fake Y-data
        it = range(1,NBins)
        fy = lambda x: x
        dataY = map(fy,it)
        dataY = dataY*n_filt_steps
        #
        sws = CreateWorkspace(OutputWorkspace = ws_name+'_count_in_time',DataX=dataX, DataY=dataY, NSpec=n_filt_steps,UnitX="Energy",\
                              VerticalAxisValues=AxLab,VerticalAxisUnit = 'Time')

    # define time-filtering steps
    filter_ws = SumSpectra(ws2proc)
    filter_wsE = ConvertUnits(filter_ws,Units,EMode="Elastic") 


    # variables to identify filtering limits
    total_signal = 0;
    # loop through experiment time to see when spurion occurs
    for ind,t_bin in enumerate(spec):
        # get time chunk
        frame_time_end = t_step*(ind+1)
        bin_time = 0.5*(t_bin+frame_time_end) # time in seconds associated with the timestamp
        chunk_ws = FilterByTime(filter_wsE ,StartTime=float(t_bin),StopTime=float(frame_time_end))
        ts = run_sttime+datetime.timedelta(microseconds=int(bin_time*1000000))
        time_string  = ts.strftime(time_format)
        if use_proton_chare_log:
            try:
                chunk_ws = NormaliseByCurrent(chunk_ws)
            except RuntimeError as rer:
                filter_result = "frame: {0:5}; time: {1:10}sec; ---- proton charge is zero. Frame ignored".format(ind,bin)
                print filter_result  
                AddTimeSeriesLog(ws, Name=log_name, Time=time_string, Value=0.)
                continue
            
        # get current, the accelerator was producing during this time
        # for simplicity, convert our signal received during the specified time interval into energy
        # calculate integrated neutron counts, counted during this time frame
        ws_value = Rebin(chunk_ws,Params='{0},{1},{2}'.format(Xmin,(Xmax - Xmin) * 1.0001,Xmax),PreserveEvents=False)
        frame_signal =  ws_value.readY(0)
        total_signal =  total_signal +frame_signal
        AddTimeSeriesLog(ws, Name=log_name, Time=time_string, Value=frame_signal[0])
        #------------------------------------------------------------------------------
        if visualize_result:    
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
        DeleteWorkspace(s1)
        pass
    DeleteWorkspace(chunk_ws)
    DeleteWorkspace(filter_ws)
    DeleteWorkspace(filter_wsE)
    DeleteWorkspace(ws_value)

    return total_signal/(len(spec))

if __name__ == '__main__':
        #fname = r'd:\Users\wkc26243\Data\Python_proj\ProcessSpurion\MER30303.nxs'
        #ws_name = 'MER30302'
        #LoadEventNexus(Filename=fname,OutputWorkspace=ws_name,SingleBankPixelsOnly='0',LoadMonitors='1',MonitorsAsEvents='0')
        ws_name = 'MER30314'
        av_signal = filter_spurion(ws_name,0.1,100.,200,True,False)
        print "avrg_signal : ",av_signal 
