import os
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
                k_name = key.strip()
                det_col[k_name ]=[]
                guide_lines.append(k_name)
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

if __name__=="__main__":
    print '--------------   MAP  ----------------------------------'
    det_col =process_detectors('MAP21385-Detectors-1.txt')
    def mapper(x):
        if float(x)<0.1:
            return '6.05'
        else:
            return x
    det_col['R'] = map(lambda x : mapper(x),det_col['R'])

    minR,min_specID,min_detID,maxR,max_specID,max_detID = find_params(det_col)
    print 'Min_det: R={0}, Spec_ID={1}, Det_ID={2}'.format(minR,min_specID,min_detID)
    print 'Max_det: R={0}, Spec_ID={1}, Det_ID={2}'.format(maxR,max_specID,max_detID)

    print '--------------   MER  ----------------------------------'
    det_col =process_detectors('MER18962-Detectors-1.txt')
    minR,min_specID,min_detID,maxR,max_specID,max_detID = find_params(det_col)
    print 'Min_det: R={0}, Spec_ID={1}, Det_ID={2}'.format(minR,min_specID,min_detID)
    print 'Max_det: R={0}, Spec_ID={1}, Det_ID={2}'.format(maxR,max_specID,max_detID)


    print '--------------   LET  ----------------------------------'
    det_col =process_detectors('LET00018482-Detectors-1.txt')
    minR,min_specID,min_detID,maxR,max_specID,max_detID = find_params(det_col)
    print 'Min_det: R={0}, Spec_ID={1}, Det_ID={2}'.format(minR,min_specID,min_detID)
    print 'Max_det: R={0}, Spec_ID={1}, Det_ID={2}'.format(maxR,max_specID,max_detID)



    det_col =process_detectors('MAR11001-Detectors-1.txt')
    minR,min_specID,min_detID,maxR,max_specID,max_detID = find_params(det_col)
    print '--------------   MAR  ----------------------------------'
    print 'Min_det: R={0}, Spec_ID={1}, Det_ID={2}'.format(minR,min_specID,min_detID)
    print 'Max_det: R={0}, Spec_ID={1}, Det_ID={2}'.format(maxR,max_specID,max_detID)

