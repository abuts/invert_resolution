""" Sample LET reduction script """ 
import os,sys
#os.environ["PATH"] = r"c:/Mantid/Code/builds/br_master/bin/Release;"+os.environ["PATH"]
from mantid import *
from Direct.ReductionWrapper import *
# this exact name is necessary for web services to work until we implement factory
class LETReduction(ReductionWrapper):
#------------------------------------------------------------------------------------#++++++
   @MainProperties
   def def_main_properties(self):
       """ Define main properties used in reduction. These are the property 
           a user usually wants to change
       """ 
       prop = {}
#---------------------------------------------------------------------------#
#------------- To change for users ---------------------------------------#
#---------------------------------------------------------------------------#
      # multiple energies provided in the data file
       prop['incident_energy'] = [12,5,2.7]
       # if energy is specified as a list (even with single value e.g. ei=[81])
       # The numbers are treated as a fraction of ei [from ,step, to ]. If energy is 
       # a number, energy binning assumed to be absolute (e_min, e_step,e_max)
       #
       prop['energy_bins'] = [-0.25,0.0025,0.85] #binning of the energy for the spe file. 
       #prop['energy_bins'] = [-2.5,0.01,4.9] #binning of the energy for the spe file. 
       #
       # the range of files to reduce.
       #
       # the range of files to reduce. This range ignored when deployed from autoreduction,
       # unless you going to sum these files. 
       # The range of numbers or run number is used when you run reduction from PC.
       prop['sample_run'] = [24414, 24415, 24416]#range(24409,24412) # 'LET18547.n001'
       prop['wb_run'] = 24363   # monovan run number 
       #
       prop['sum_runs'] = True # set to true to sum everything provided to sample_run
       #   
       # Absolute units reduction properties. Set prop['monovan_run']=None to do arbitrary units
       prop['monovan_run'] = None # vanadium run in the same configuration as your sample
       #prop['sample_mass'] = 100 #  # mass of your sample
       #prop['sample_rmm'] = 200 #  molecular weight of scatterers in your sample
       return prop
#------------------------------------------------------------------------------------#
   @AdvancedProperties
   def def_advanced_properties(self):
      """  Set up advanced properties, describing reduction.
           These are the properties, usually provided by an instrument 
           scientist
            
           separation between simple and advanced properties depends
           on scientist, experiment and user.   All are necessary for reduction 
           to work properly
      """
      prop = {}
      prop['map_file'] = 'LET_rings_153.map'
      prop['det_cal_file'] = 'det_corrected_cycle153.dat'
      prop['bleed'] = False
      prop['norm_method']='current'
      prop['detector_van_range']=[4.8,5.2]
      prop['background_range'] = [92000,98000] # TOF range for the calculating flat background
      prop['hardmaskOnly']='hard_153.msk' # Use diag (hardmaskPlus option) to enhance hard masks
      prop['check_background']=True
      prop['save_format'] = 'nxspe'
      # if two input files with the same name and  different extension found, what to prefer. 
      prop['data_file_ext']='.nxs' # for LET it may be choice between event and histo mode if 
      # raw file is written in histo, and nxs -- in event mode
      # Absolute units: map file to calculate monovan integrals
      prop['monovan_mapfile'] = 'LET_rings_153.map'
#
      prop['load_monitors_with_workspace']=False
      # change this to correct value and verify that motor_log_names refers correct and existing 
      # log name for crystal rotation to write correct psi value into nxspe files
      prop['motor_offset']=None
      # vanadium  mass valid from cycle 2013/5
      prop['vanadium-mass']=8.47
      #
      prop['monovan_lo_frac']=-0.4
      prop['monovan_hi_frac']= 0.4
      #RAE added
      prop['bleed_maxrate']=0.005
      return prop
      #
#------------------------------------------------------------------------------------#
   @iliad
   def reduce(self,input_file=None,output_directory=None):
      """ Method executes reduction over single file

          Overload only if custom reduction is needed or 
          special features are requested
      """
      print 'Input file: ',input_file
      if input_file:
          self.reducer.prop_man.sample_run = input_file
      ws = PropertyManager.sample_run.get_workspace()
      calc_count_rate(ws,Emin,Emax,NBins)
      FilterByLogValue(ws,wsf,"block_count_rate",1.2e+5,1.8e+5);
      self.reducer.prop_man.sample_run = wsf
      PropertyManager.sample_run.synchronize_ws()
      
      results = ReductionWrapper.reduce(self,input_file,output_directory)
      #
      return results
   #
   def set_custom_output_filename(self):
      """ define custom name of output files if standard one is not satisfactory 
          In addition to that, example of accessing reduction properties 
          Changing them if necessary
      """ 
      def custom_name(prop_man):
          """ sample function which builds filename from 
              incident energy and run number and adds some auxiliary information 
              to it.
          """ 
          # Note -- properties have the same names  as the list of advanced and 
          # main properties
          ei = PropertyManager.incident_energy.get_current()
          # sample run is more then just list of runs, so we use 
          # the formalization below to access its methods
          run_num = PropertyManager.sample_run.run_number()
          name = "LET{0}_{1:<3.2f}meV_Rings".format(run_num ,ei)
          return name
       
      # Uncomment this to use custom filename function
      # Note: the properties are stored in prop_man class accessed as
      # below. 
      return lambda : custom_name(self.reducer.prop_man)
      # use this method to use standard file name generating function
      #return None
    #
   #
   def validation_file_place(self):
      """Redefine this to the place, where validation file, used in conjunction with
         'validate_run' property, located. Here it defines the place to this script folder.
          but if this function is disabled, by default it looks for/places it 
          in a default save directory"""
      return os.path.split(os.path.realpath(__file__))[0]

    
   def __init__(self,web_var=None):
       """ sets properties defaults for the instrument with Name"""
       ReductionWrapper.__init__(self,'LET',web_var)
       
def calc_count_rate(ws,Emin,Emax,NBins,visualize_result=True,use_proton_chare_log=True,trust_good_frames=True,log_name = "block_count_rate"):
    """ Filter workspace from spurion appearing in specified energy range
        Inputs:
        ws -- the name or pointer to workspace to filter
        Emin -- minimal (elastic) energy range to filter within
        Emax -- maximal (elastic) energy range to filter within

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
       
#
if __name__ == "__main__":
#------------------------------------------------------------------------------------#
# SECTION USED TO RUN REDUCTION FROM MANTID SCRIPT WINDOW #
#------------------------------------------------------------------------------------#
##### Here one sets up folders where to find input data and where to save results ####
    # It can be done here or from Mantid GUI:
    #      File->Manage user directory ->Browse to directory
    # Folder where map and mask files are located:
    map_mask_dir = '/usr/local/mprogs/InstrumentFiles/merlin'
    config.appendDataSearchDir(map_mask_dir)
    # folder where input data can be found
    #data_dir = 'd:/Data/Mantid_Testing/15_01_27/LET/data'
    # auxiliary folder with results
    #rez_dir = 'd:/Data/Mantid_Testing/15_01_27/LET'
    # Set input search path to values, specified above
    #config.setDataSearchDirs('{0};{1};{2}'.format(data_dir,map_mask_dir,rez_dir))
    # use appendDataSearch directory to add more locations to existing Mantid 
    # data search path
    #config.appendDataSearchDir('d:/Data/Mantid_GIT/Test/AutoTestData')
    # folder to save resulting spe/nxspe files.
    #config['defaultsave.directory'] = rez_dir

###### Initialize reduction class above and set up reduction properties.        ######
######  Note no web_var in constructor.(will be irrelevant if factory is implemented)
    rd = LETReduction()
    # set up advanced and main properties74,28275,28276,28277,28278
    rd.def_advanced_properties()
    rd.def_main_properties()

#### uncomment rows below to generate web variables and save then to transfer to   ###
    ## web services.
    run_dir = os.path.dirname(os.path.realpath(__file__))
    file = os.path.join(run_dir,'reduce_vars.py')
    rd.save_web_variables(file)

#### Set up time interval (sec) for reducer to check for input data file.         ####
    #  If this file is not present and this value is 0,reduction fails 
    #  if this value >0 the reduction waits until file appears on the data 
    #  search path checking after time specified below.
    rd.wait_for_file = 1200  # waiting time interval in seconds

### Define a run number to validate reduction against future changes    #############
    # After reduction works well and all settings are done and verified, 
    # take a run number with good reduced results and build validation
    # for this result. 
    # Then place the validation run together with this reduction script.
    # Next time, the script will run reduction and compare the reduction results against
    # the results obtained earlier.
    #rd.validate_run_number = 21968  # Enabling this property disables normal reduction
    # and forces reduction to reduce run specified here and compares results against
    # validation file, processed earlier or calculate this file if run for the first time.
    #This would ensure that reduction script have not changed,
    #allow to identify the reason for changes if it was changed 
    # and would allow to recover the script,used to produce initial reduction
    #if changes are unacceptable.

####get reduction parameters from properties above, override what you want locally ###
   # and run reduction. Overriding would have form:
   # rd.reducer.prop_man.property_name (from the dictionary above) = new value e.g. 
   # rd.reducer.prop_man.energy_bins = [-40,2,40]
   # or 
   ## rd.reducer.prop_man.sum_runs = False
   # 
   # 
    rd.run_reduction()

