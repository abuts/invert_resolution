import os
import numpy as np
import math
import time
import datetime as dtt

def process_mantid(filename):

    f=open(filename)
    k_of_i={'LoadNXSPE':[],'ConvertToMD':[],'SaveMD':[]}
    keys = k_of_i.keys();
    for line in f:
        if 'successful' in line:
            for key in keys:
                if key in line:
                    cont = line.split();
                    k_of_i[key].append(float(cont[3]));


    f.close();
    return k_of_i;


def process_horace(filename):

    f=open(filename)
    k_of_i={'convert from spe':[],'save sqw data to file':[]}
    keys = k_of_i.keys();
    process_next_line = False;
    for line in f:
        if not(process_next_line):
            for key in keys:
                 if key in line:
                     rez_key = key;
                     process_next_line = True;
        else: # process next line and it is next line
            process_next_line = False;
            cont = line.split();
            k_of_i[rez_key].append(float(cont[3])); # Add time to log


    f.close();

    return k_of_i;

def convert_hor_usage_block(raw_block,t0,size):
    t_start = t0.mktime(t0);
    final_block=[];
    t0=0;

    for tp in raw_block:
        pers = tp[0];
        dt   = tp[1]-t0;
        t0   = tp[1];
        speed= size*pers/dt;
        final_block.append((t0+t_start,speed));
    return final_block;

def convertTimeToSeconds(data_table,t0,col_to_convert):
    """ 
    """
    all_col_names = data_table.keys();

    t0 = time.mktime(t0);
    for name in all_col_names:
        if name in col_to_convert:
            colmn=data_table[name];
            nrows = len(colmn)
            for i in xrange(0,nrows):
                t1=time.mktime(colmn[i])
                colmn[i]=(t1-t0);

    return data_table;

def process_HorWrite_withTime(filename):
    ''' Process Horace write sqw/combine tmp files statistics log written with time marks''' 
    f=open(filename)

    perf_data={'Reading header(s)':[],'accumulating binning':[],'Writing to output file':[]}
    keys = perf_data.keys();

    process_next_lines = False;
    in_conv_time = False;
    BlockKey='';
    for line in f:
        if not(process_next_lines):
            for key in keys:
                 if key in line:
                     process_next_lines=True;
                     in_conv_time      =True;
                     BlockKey = key;
                     break;
            continue
        else:
            cont = line.split();
            if in_conv_time:
                in_conv_time = False;
                t0 = time.strptime(cont[3]+' '+cont[4],'%Y/%m/%d %H:%M:%S')
                t = time.mktime(t0);
                perf_data[BlockKey].append((t,0.));
                continue
            else:
                if cont[0] == 'Completed':
                    t =time.mktime(t0);
                    t = t+float(cont[5]);
                    cmplString = cont[1];
                    completed=float(cmplString[0:-1]);
                    perf_data[BlockKey].append((t,completed));
                else:
                    process_next_lines=False;
                    BlockKey = '';
                    continue


    
    f.close();
    return perf_data;


def process_horace_withTime(filename):
    ''' Process Horace convert spe to tmp statistics logs written with time marks''' 

    f=open(filename)
    proc_keys=['convert from spe','save sqw data to file','System time is', 'Processing spe file','read spe and detector'];
    k_of_i={'Conv_Time':[],'Conv_tLength':[],'Saved_Time':[],'Save_tLength':[],'Completed':[],'read_time':[]}


    process_next_line = False;
    in_conv_time = True;
    in_processing_spe=False;
    for line in f:
        if not(process_next_line):
            for key in proc_keys:
                 if key in line:
                     if key is proc_keys[0]: # Save data to file block
                         in_conv_time = True;
                         rez_key = 'Conv_tLength';
                         process_next_line = True;
                     elif key is proc_keys[1]:  # convert from spe block
                         in_conv_time = False;
                         rez_key = 'Save_tLength';
                         process_next_line = True;
                     elif key is 'System time is':
                         if in_processing_spe:
                            cont = line.split();
                            theTime = time.strptime(cont[3]+' '+cont[4],'%Y/%m/%d %H:%M:%S')
                            if in_conv_time:
                                k_of_i['Conv_Time'].append(theTime);
                            else:
                                k_of_i['Saved_Time'].append(theTime);
                                in_processing_spe = False;
                     elif key is 'read spe and detector':
                         rez_key = 'read_time';
                         process_next_line = True;
                     elif key is 'Processing spe file':
                         in_processing_spe = True;
                         cont = line.split();
                         den  = cont[5].split(':');
                         rez = 100*float(cont[3])/float(den[0])
                         k_of_i['Completed'].append(rez);


                     # end key in line
                     break;
                        
        else: # process next line and it is the next line
            cont = line.split();
            if process_next_line:
                process_next_line = False;
                k_of_i[rez_key].append(float(cont[3]));
                continue;

    f.close();

    if len(k_of_i['Conv_Time'])==0:
        rt = 0;
        read_time = k_of_i['read_time']
        conv_time = k_of_i['Conv_tLength']
        save_time = k_of_i['Save_tLength']
        t0=0;
        for i in xrange(0,len(read_time)):
            theTime=t0+read_time[i]+conv_time[i];
            k_of_i['Conv_Time'].append(theTime);
            theTime = theTime+save_time[i];
            k_of_i['Saved_Time'].append(theTime);
            t0=theTime;

    return k_of_i;

def calc_statistics(k_of_i):
    ncols = len(k_of_i.keys());
    nrows = len(k_of_i.itervalues().next());


    data = np.zeros((nrows,ncols));
    ncol = 0;
    #total_time=0;
    for key in k_of_i.iterkeys():
        rez = np.array(k_of_i[key]);
        nel = rez.shape;
        nruns = nel[0];
        if nruns == 0:
            continue

        rez.reshape(nruns,1);
        nr = min(nruns,nrows);
        data[0:nr,ncol]=rez[0:nr];
        ncol+=1;

        sum = np.sum(rez);
        #total_time = total_time+sum
        Avrg = float(sum/nruns);
        min_val=rez.min();
        max_val=rez.max();
        sdev = np.subtract(rez,Avrg);
        sigma= math.sqrt(np.dot(sdev,sdev)/nruns);


        print " Value {0}\t: N elements: {1:<17}\t Average: {2:<5},\t Sigma = {3:<5},\t min Val={4:<5},\t max Val={5:<} ".format(key,nruns,Avrg,sigma,min_val,max_val)
        #print " Time to run {0} instances of {1}".format(nruns,key)
    #print " Total conversion time: {0} sec or {1} min".format(total_time,total_time/60);
    return data;

def write_matl_table(filename,rezDict):
    f=open(filename,'w')

    ncols = len(rezDict.keys());
    col_names=['Conv_Time','Conv_tLength','Completed','Saved_Time','Save_tLength']

    for name in col_names:
        f.write(" {0:10}  ".format(name));
    f.write("\n");

    nrows = len(rezDict.itervalues().next());

    for i in xrange(0,nrows):
        for name in col_names:
            colmn=rezDict[name];
            result = colmn[i];
            f.write(' {0:10}  '.format(result))
        f.write("\n");

    f.close();


def process_iostat_log(filename):
    f=open(filename)
    result={'Time':[],'CPU_user':[],'CPU_system':[],'IOwait':[],'Idle':[],'sda:tps':[]}
    inBlock=False;
    process_cpu = False;
    nempty=0;
    for line in f:
        components = line.split();
        if len(components) == 0:
            continue

        if not inBlock:                    
            try:
                if len(components) == 3:
                    theTime = time.strptime(components[0]+' '+components[1]+' '+components[2],'%m/%d/%Y %I:%M:%S %p')
                elif len(components) == 2:
                    theTime = time.strptime(components[0]+' '+components[1],'%m/%d/%Y %I:%M:%S')
                else:
                    raise ValueError;

                result['Time'].append(theTime);
                inBlock = True;
                continue
            except ValueError:
                continue

        if 'avg-cpu' in components[0]:
            process_cpu = True;
            continue;
        if process_cpu:
            process_cpu = False;
            result['CPU_user'].append(components[0]);
            result['CPU_system'].append(components[2]);
            result['IOwait'].append(components[3]);
            result['Idle'].append(components[5]);
            continue;

        if 'sdi' == components[0]:
            result['sda:tps'].append(components[1]);
            inBlock = False;
            continue;



    f.close();
    return result;

def write_iostat_stats(filename,stats,t0):

    ncols = len(stats.keys());

    col_names=['Time','CPU_user','CPU_system','IOwait','Idle','sda:tps']


    f=open(filename,'w' )
    for name in col_names:
        f.write(" {0:10}  ".format(name));
    f.write("\n");

    nrows = len(stats.itervalues().next());
    for i in xrange(0,nrows):
        for name in col_names:
            colmn=stats[name];
            result = colmn[i];
            if name is col_names[0] :
               t1=time.mktime(result);
               result=t1-t0;
            f.write(' {0:10}'.format(result))

        f.write("\n");



    f.close();
def  write_matlwr_perf_table(filename,write_log,t0,size):
    f=open(filename,'w' )

    col_names=write_log.keys();

    for name in col_names:
        f.write("# {0} ------->\n".format(name))
        data = write_log[name];
        nrows = len(data);
        t_prev=0;
        for i in xrange(0,nrows):
            t1,compl = data[i];
            t = t1-t0;
            if i==0:
                speed=0;
            else:
                speed = (compl/100)*(size/(1024*1024))/(t-t_prev);
            t_prev=t;

            f.write('{0:10}, {1:10}, {2:10}\n'.format(t,compl,speed));


    f.close();
    pass


if __name__ == '__main__':
    #print "===========    MANTID: ============>"
    #manrez=process_mantid('mantid_log.txt')
    #man_stat=calc_statistics(manrez);
    #manrez=process_mantid('mantid_log1.txt')
    #man_stat=calc_statistics(manrez);
    #print "===========    HORACE: ============>"
    #horrez = process_horace('Horace_log.txt')
    #hor_stat=calc_statistics(horrez);
    #horrez = process_horace('Horace_log2.txt')
    #hor_stat=calc_statistics(horrez);
    #horrez = process_horace('Horace_runTimeMex.txt')
    #hor_stat=calc_statistics(horrez);

    #print "===================================>"

    #matlab1   = process_horace_withTime('Horace_perfJob1b_2014_06_07.log');
    #write1_log= process_HorWrite_withTime('Horace_perfJob1b_2014_06_07.log');
    #t0=(matlab1['Conv_Time'])[0];

    #t0=time.mktime(t0);
    #size1 = 156921310384#/231 # 
    #write_matl_table('Horace_RUN1_STATS_2014_06_07.txt',matlab1,t0)
    #write_matlwr_perf_table('Horace_SAVE1_STATS_2014_06_07.txt',write1_log,t0,size1);

    #size2 = 156242437492#/230 # 
    #matlab2 = process_horace_withTime('Horace_perfJob2b_2014_06_07.log');
    #write2_log= process_HorWrite_withTime('Horace_perfJob2b_2014_06_07.log');
    #write_matl_table('Horace_RUN2_STATS_2014_06_07.txt',matlab2,t0)
    #write_matlwr_perf_table('Horace_SAVE2_STATS_2014_06_07.txt',write2_log,t0,size2);


    #syslog=process_iostat_log('melehan_iostatPerf_2Jobsb2014_06_07.log')
    #write_iostat_stats('syslogStats2014_06_07.txt',syslog,t0)

    file1='HoraceWin_perfJob1_2014_06_07.log'
    file2='HoraceWin_perfJob2_2014_06_07.log'
    file1='Chadwick_Horace_perfJob1_2014_05_22.log'
    file1='hor_perfJob1b_150G_2014_06_16.log'

    file1='Horace_perfJob2_2014_06_12_br1.log'
    file1='Horace_perf300GbJob_Melehan_2014_06_11.log'
    file1='Horace_perfJob1_2014_06_08nomex.log'
    file2='SCARF_perf150GJob2.log'
    file1='horRHEL7_perfJob1b_150G_2014_06_16.log'
    file2='horRHEL7_perfJob2b_150G_2014_06_16.log'
    file1='horRHEL7_perfJob1_dirty_300G_2014_06_16.log'
    file2='horRHEL7_perfJob1_clean300G_2014_06_16.log'
    file1='RHEL7_hor_perfJob1_150G_2014_07_03a.log'
    file2='RHEL7_hor_perfJob2_150G_2014_07_03a.log '
    file1='NewMERLINSmallHeap_150GbJob128-Jul-2014_15-44-20.log'
    file2='NewMERLINSmallHeap_150GbJob228-Jul-2014_15-45-06.log'
    file0='NewMERLIN_150GbJob1_25-Jul-2014_13-37-00.log' 

    #fileM="Mantid_150Gb_convToMDMergeMD_Job1Performance2014_06_11.log"
    #mantid = process_mantid(fileM);
    #rez = calc_statistics(mantid);
    horrez = process_horace(file0)
    hor_stat=calc_statistics(horrez);
    horrez = process_horace(file1)
    hor_stat=calc_statistics(horrez);


    #write1_log= process_HorWrite_withTime(file1);

    matlab1   = process_horace_withTime(file1);
    t_start=(matlab1['Conv_Time'])[0];
    matlab1   = convertTimeToSeconds(matlab1,t_start,['Conv_Time','Saved_Time']);
    rez = calc_statistics(matlab1);
    write1_log= process_HorWrite_withTime(file1);

    t0=time.mktime(t_start);
    size1 = 156921310384#/231 # 
    write_matl_table('NewMERLINSmallHeap_job1run_150Gb1.txt',matlab1)
    write_matlwr_perf_table('NewMERLINSmallHeap_job1save_150Gb1.txt',write1_log,t0,size1);


    size2 = 156242437492#/230 # 
    matlab2 = process_horace_withTime(file2);
    #
    t_start=(matlab2['Conv_Time'])[0];
    t0=time.mktime(t_start);
    #
    matlab2   = convertTimeToSeconds(matlab2,t_start,['Conv_Time','Saved_Time']);
    write2_log= process_HorWrite_withTime(file2);
    write_matl_table('NewMERLINSmallHeap_job2run_150Gb1.txt',matlab2)
    write_matlwr_perf_table('NewMERLINSmallHeap_job2save_150Gb1.txt',write2_log,t0,size2);


   # syslog=process_iostat_log('rhel7_iostat_perf_2014_06_16.log')
   # write_iostat_stats('rhel7_syslogStats2014_06_16.txt',syslog,t0)



    #f=open('time_torun.txt','w')
    #
    #f.close();

