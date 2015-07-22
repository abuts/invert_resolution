from collections import OrderedDict

def process_detectors(filename):

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
    min_R = min(det_col['R'])
    max_R = max(det_col['R'])
    min_ind = det_col['R'].index(min_R)
    max_ind = det_col['R'].index(max_R)
    min_specID = det_col['Spectrum No'][min_ind]
    max_specID = det_col['Spectrum No'][max_ind]
    min_detID = det_col['Detector ID(s)'][min_ind]
    max_detID = det_col['Detector ID(s)'][max_ind]

    return (min_R,min_specID,min_detID,max_R,max_specID,max_detID)

def build_rings_spectra_map(det_col):
    det_theta=det_col['Theta']
    spec_num = det_col['Spectrum No']
    distance = det_col['R']

    ring_list = {};
    for det in zip(det_theta,spec_num,distance):
        if(float(det[2]) > 5.5):
            theta = float(det[0])
            theta0 = round(theta,0)
            theta1 = round(theta+0.5,0)
            key = '{0:.1f}'.format(0.5*(theta0+theta1))
            spec_id = det[1].replace(',','')
            if key in ring_list:
                ring_list[key].append(spec_id)
            else:
                ring_list[key] = [spec_id]
    return ring_list

def save_rings_map(ring_dic,filename):    
    """ Save ISIS ring map """
    def numeric_compare(x, y):
        return int(10*(float(x) - float(y)))

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
    det_col =process_detectors('MAP_SpectraList.txt')
    print det_col.keys()
    ring_list = build_rings_spectra_map(det_col)
    print 'Identified ',len(ring_list),' rings'
    save_rings_map(ring_list,'MAP_rings2015.map')
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

