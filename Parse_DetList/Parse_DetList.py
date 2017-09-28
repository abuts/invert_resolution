from collections import OrderedDict
from  math import ceil,floor

def process_detectors(filename):
    """ Process the file containing semicolon (;) separated list of detectors information
        imported from MANTID
       
        MANTID has to export detectors list as ASCII file with semicolon-separated option
        and header retained and procedure convert this file  into map with keys as column names,
        and the map values are the lists retrieved from column values.
    """

    f=open(filename)
    det_col={}
    guide_lines = []
    ic = 0
    for line in f:
        cont = line.split(';')
        if ic==0:
            # parse header
            for key in cont:
                kr = key.strip()
                det_col[kr]=[]
                guide_lines.append(kr)
            ic+=1
        else:
          if line.find('yes')!=-1:
             continue
          for ind,keyl in enumerate(guide_lines):
                det_col[keyl].append(cont[ind])
    f.close()
    return det_col


#
def find_params(det_col):
    """ Procedure to find min-max values for detectors"""
    min_R = min(det_col['R'])
    max_R = max(det_col['R'])
    min_ind = det_col['R'].index(min_R)
    max_ind = det_col['R'].index(max_R)
    min_specID = det_col['Spectrum No'][min_ind]
    max_specID = det_col['Spectrum No'][max_ind]
    min_detID = det_col['Detector ID(s)'][min_ind]
    max_detID = det_col['Detector ID(s)'][max_ind]

    return (min_R,min_specID,min_detID,max_R,max_specID,max_detID)

def build_rings_spectra_map(det_col,step=0.5):
    """ Procedure to group spectra in lists according to Theta column values
        with accuracy 0.5 in Theta.
    """
    det_theta=det_col['Theta']
    theta_min = 360;
    for i,th in enumerate(det_theta):
        TH = float(th)
        if TH !=0 and TH<theta_min:
            theta_min = TH
        det_theta[i] = TH

    spec_num = det_col['Spectrum No']
    distance = det_col['R']

    theta_max = max(det_theta)
    n_steps = (ceil(theta_max) - floor(theta_min))/step
    #indexes = range(0,n_steps)
    tml = 0
    ring_list = {};
    for det in zip(det_theta,spec_num,distance):
        if(float(det[2]) < 5 and  float(det[1]) > 0): #Skip monitors && missing spectra
            theta = float(det[0])
            if theta > tml:
                tml = theta
            #theta0 = round(theta,0)
            #theta1 = round(theta+0.5,0)
            #key = '{0:.1f}'.format(0.5*(theta0+theta1))
            key = floor((theta-floor(theta_min))/step)
            spec_id = det[1].replace(',','')
            if key in ring_list:
                ring_list[key].append(spec_id)
            else:
                ring_list[key] = [spec_id]
    return ring_list

def save_rings_map(ring_dic,filename):
    """ Save ISIS ring map calculated by build_rings_spectra_map procedure
        in ISIS map format. First string of the script contains 
        the assumptions about map accuracy (accuracy is less or equal then 0.1degree in Theta)
    """

    ring_map = OrderedDict(sorted(ring_dic.iteritems(), key = lambda key_value : int(10*float(key_value[0]))))

    f = open(filename,'w')
    f.write('{0:>8}\n'.format(len(ring_map)))
    ic = 1
    for key,val in ring_map.iteritems():
        f.write('{0}\n'.format(ic))
        f.write('{0}\n'.format(len(val)))
        for kw in val:
            f.write('  {0}'.format(kw))
        f.write('\n\n')
        ic+=1
    f.close()



if __name__=="__main__":
    print '--------------   MAP  ----------------------------------'
    i=10
    det_col =process_detectors('HET05240-Detectors-1.txt')
    print "The detector file contains the following columns: ", det_col.keys()
    ring_list = build_rings_spectra_map(det_col,0.2)
    print 'Identified ',len(ring_list),' rings'
    save_rings_map(ring_list,'HET_Rings_NoPSD_2017.map')
    #det_col =process_detectors('MAP21385-Detectors-1.txt')

    #def mapper(x):
    #    if float(x)<0.1:
    #        return '6.05'
    #    else:
    #        return x
    #det_col['R'] = map(lambda x : mapper(x),det_col['R'])

    #minR,min_specID,min_detID,maxR,max_specID,max_detID = find_params(det_col)
    #print 'Min_det: R={0}, Spec_ID={1}, Det_ID={2}'.format(minR,min_specID,min_detID)
    #print 'Max_det: R={0}, Spec_ID={1}, Det_ID={2}'.format(maxR,max_specID,max_detID)


    #print '--------------   LET  ----------------------------------'
    #det_col =process_detectors('LET00005545-Detectors-1.txt')
    #minR,min_specID,min_detID,maxR,max_specID,max_detID = find_params(det_col)
    #print 'Min_det: R={0}, Spec_ID={1}, Det_ID={2}'.format(minR,min_specID,min_detID)
    #print 'Max_det: R={0}, Spec_ID={1}, Det_ID={2}'.format(maxR,max_specID,max_detID)

    #print '--------------   MER  ----------------------------------'
    #det_col =process_detectors('MER18962-Detectors-1.txt')
    #minR,min_specID,min_detID,maxR,max_specID,max_detID = find_params(det_col)
    #print 'Min_det: R={0}, Spec_ID={1}, Det_ID={2}'.format(minR,min_specID,min_detID)
    #print 'Max_det: R={0}, Spec_ID={1}, Det_ID={2}'.format(maxR,max_specID,max_detID)


    #det_col =process_detectors('MAR11001-Detectors-1.txt')
    #minR,min_specID,min_detID,maxR,max_specID,max_detID = find_params(det_col)
    #print '--------------   MAR  ----------------------------------'
    #print 'Min_det: R={0}, Spec_ID={1}, Det_ID={2}'.format(minR,min_specID,min_detID)
    #print 'Max_det: R={0}, Spec_ID={1}, Det_ID={2}'.format(maxR,max_specID,max_detID)

