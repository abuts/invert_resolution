from collections import OrderedDict
from  math import ceil,floor

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



if __name__=="__main__":
    files =['lam=.4.m','lam=.5.m','lam=.7.m','lam=1.m','lam=1.5.m','lam=2.m','lam=3.m','lam=4.m']
    wavelength =[0.4,0.5,0.7,1.,1.5,2,3,4]
    tfp = open('table','w');
    for fil,wl in zip(files,wavelength):
        time,signal,error = process_file(fil);
        for t,s,e in zip(time,signal,error):
            tfp.write('{0},{1},{2},{3}\n'.format(wl,t,s,e))
    tfp.close();
