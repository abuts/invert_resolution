from collections import OrderedDict
from  math import ceil,floor,sqrt
import numpy as np

def process_file(filename):
    """ Process the file containing semicolon (;) separated list of detectors information
        imported from MANTID
       
        MANTID has to export detectors list as ASCII file with semicolon-separated option
        and header retained and procedure convert this file  into map with keys as column names,
        and the map values are the lists retrieved from column values.
    """

    f=open(filename)
    time=[]
    signal = []
    error = []
    for line in f:
        if line[0] == '#':
            continue
        cont = line.split(' ')
        time.append(float(cont[0]))
        signal.append(float(cont[1]))
        error.append(float(cont[2]))
    f.close()
    return (time,signal,error)
#
def read_det_table(filename):
    """ Read the detector table to find spectra"""
    f=open(filename)
    spec_num=[]
    q = []
    for line in f:
        if line[0] == '#':
            continue
        cont = line.split(',')
        spec_num.append(int(cont[0]))
        q.append(float(cont[4]))
    f.close()
    return (spec_num,q)

def read_s_of_tof(filename):
    """ Read the exported Mantid workspace, containing signal vs TOF"""
    f=open(filename)
    tof =[]
    sp_num = [];
    signal = {};
    num = 0;
    for line in f:
        if line[0] == '#':
            continue
        lin = line.rstrip(',\n')
        cont = lin.split(',')
        y = np.array(cont[1:])
        if num == 0:
            tof = y.astype(np.float)
        else:
            sp_num.append(int(cont[0]))
            signal[num-1] = y.astype(np.float)
        # end

        num +=1
    f.close()
    return (tof,sp_num,signal)



def write_mod_table(ouf_filename,filenums):

    tfp = open(ouf_filename,'w');
    #for fil,wl in zip(files,wavelength):
    for ind in fn:
        fil = 'mod/mod{0}.m'.format(ind)
        en = 0.5*(ind+ind+1);
        wl = 9.058/sqrt(en); # lambda in A 
        time,signal,error = process_file(fil);
        for t,s,e in zip(time,signal,error):
            tfp.write('{0},{1},{2},{3}\n'.format(wl,t,s,e))
    tfp.close();

def write_signal(ouf_filename,tof_line,signal):

    tfp = open(ouf_filename,'w');
    for t,s in zip(tof_line,signal):
        tfp.write('{0},{1}\n'.format(t,s))
    tfp.close();




if __name__=="__main__":
    #files =['lam=.4.m','lam=.5.m','lam=.7.m','lam=1.m','lam=1.5.m','lam=2.m','lam=3.m','lam=4.m']
    #wavelength =[0.4,0.5,0.7,1.,1.5,2,3,4]
    #fn = range(9,161);
    #write_mod_table('mod_e_table',fn);
    spn,q = read_det_table('tof_pulse/SR_MER22346_150_mevSTRIP-Detectors-1.txt')
    tof,spn1,sig_dic = read_s_of_tof('tof_pulse/SR_MER#150mEv#022346_tf.csv')
    qi_min = 1.5
    qi_max = 1.7
    #sp_in = [];
    tof_line = 0.5*(tof[:-1]+tof[1:])
    signal = np.zeros(tof_line.size)
    for sp,Q in zip(spn,q):
        if Q>=qi_min and Q<= qi_max:
            signal += sig_dic[sp]
            #sp_in.append(sp)
    write_signal('sig150mEv_Q1_5.csv',tof_line,signal)
    

