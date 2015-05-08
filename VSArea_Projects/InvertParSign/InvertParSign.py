from mantid.kernel import funcreturns
from mantid.simpleapi import *
from mantid import config
from DirectEnergyConversion import *
import sys
import time as time
import numpy
import nxs
import inspect as insp        


def invert_par_sign():
    f=open('4to1_095.par')
    ft=open('4to1_095_inv.par','w');

    for line in f:
        cont = line.split()
        if len(cont) > 1:
            rez = [];
            for num in cont:
                rez.append(float(num));
            ft.write(' {0: 9.5f} {1: 9.5f} {2: 9.5f} {3: 9.5f} {4: 9.5f} {5: 9d}\n'.format(rez[0],rez[1],-rez[2],rez[3],rez[4],int(rez[5])))
        else:
            ft.write(line);


    f.close();
    ft.close();

def rough_log_file(filename,out_file):
        f=open(filename)
        ft=open(out_file,'w')
        for line in f:
            cont = line.split(',');
            for num in cont:
                mn = num.translate(None,'[]');
                try:
                    fn = float(mn)
                    ft.write('{0:6.2f},'.format(fn))
                except:
                    ft.write(line.rstrip())
                    break

            ft.write('\n');
                            
        f.close()
        ft.close()

def export_nxspe2phx(filename):

    file = nxs.open(filename,'r');
    entries = file.getentries()
    if len(entries) == 0:
        raise IOError(" Can not find nexus entries in the file: {0}".format(filename));
    nxspe_data = entries.keys();
    if len(nxspe_data)> 1:
        raise IOError(" Currently do not support more then one nexus enrty in the file : {0}".format(filename));

        
    file.opengroup(nxspe_data[0])
    file.opengroup('data')
    phx_par = {"polar":[],"polar_width":[],"azimuthal":[],"azimuthal_width":[]};

    for key in phx_par:
        file.opendata(key)
        phx_par[key]=file.getdata();
        file.closedata();

    file.closegroup();
    file.closegroup();
    file.close();

    out_file = os.path.splitext(filename)[0]+'.phx'
    if os.path.exists(out_file) :
       os.unlink(out_file);

    ft = open(out_file,'w');
    numdet = len(phx_par["polar"]);
    ft.write("{0}\n".format(numdet))
    for i in xrange(0,numdet):
        ft.write("1\t0\t{0: 5.2f}\t{1: 5.2f}\t{2: 5.2f}\t{3: 5.2f}\t0\n".format(phx_par['polar'][i],phx_par['azimuthal'][i],phx_par['polar_width'][i],phx_par['azimuthal_width'][i]));
    ft.close()


    pass

def convert_det2phx(det_filename):
    f=open(det_filename)
    out_file = os.path.splitext(det_filename)[0]+'.phx'
    if os.path.exists(out_file) :
       raise LookupError(" File {0} already exist".format(out_file))   


    ft=open(out_file,'w');
    polar
    for line in f:
        cont = line.split()
        if len(cont) > 1:
            rez = [];
            for num in cont:
                rez.append(float(num));
            ft.write(' {0: 9.5f} {1: 9.5f} {2: 9.5f} {3: 9.5f} {4: 9.5f} {5: 9d}\n'.format(rez[0],rez[1],-rez[2],rez[3],rez[4],int(rez[5])))
        else:
            ft.write(line);


    f.close();
    ft.close();
def build_map_from_mantid_dets(det_text_file):

 
    angles={};
    f=open(det_text_file)
    for line in f:
        cont = line.split('\t')
        # column 4 is polar angle
        angl_key = "{0: 5.2f}".format(float(cont[4]));
        # column 1 is a spectra ID. 
        if angl_key in angles:
            angles[angl_key].append(int(cont[1]))
        else:
            angles[angl_key] = [int(cont[1])]

    f.close();

    selected_angles = angles.keys();
    selected_angles.sort(key=float)

    out_file = os.path.splitext(det_text_file)[0]+'.map'
    f = open(out_file,'w');    

    n_angles=len(angles);        
    f.write("    {0}\n".format(n_angles));
    
    for i in xrange(0,n_angles):
        det_angle =selected_angles[i]; 
        spectra = angles[det_angle];   
        
        f.write("{0} \n".format(i+1))
        f.write("{0} \n".format(len(spectra)));

        for spectra_id in spectra :
            f.write("{0}  ".format(spectra_id))
        f.write("\n\n");

    f.close();

def write_1to2map(targ_text_file,n_spectra):

 

    out_file = os.path.splitext(targ_text_file)[0]+'.map'

    f = open(out_file,'w');    

    f.write("    {0}\n".format(n_spectra));
    
    for i in xrange(0,n_spectra):       
        f.write("{0} \n".format(i+1))
        f.write("{0} \n".format(1));

        f.write("{0}  ".format(i+1))
        f.write("\n");

    f.close();

if __name__ == '__main__':

    #invert_par_sign();

    #rough_log_file(r'd:\Data\Mantid_Testing\13_11_22\TableHorace.txt','TableHorRough.txt')
    #rough_log_file(r'd:\Data\Mantid_Testing\13_11_22\TableMantid.txt','TableMantidRough.txt')

    #build_map_from_mantid_dets('MAR18622-Detectors-1.DAT')
    #export_nxspe2phx('mar18939.nxspe')
    write_1to2map('mar285_one2oneRings',285)
    pass