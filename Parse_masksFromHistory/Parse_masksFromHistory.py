def process_detectors(filename):

    f=open(filename)
    for line in f:
        pos = line.find('DetectorList=')
        if pos != -1:
            cont = line[pos+len('DetectorList='):]
            cont = filter(lambda x : x in '0123456789,',cont)
            data = cont.split(',')
            break
    f.close()
    for ind,val in enumerate(data):
        data[ind] = int(val)
    return data

if __name__ == '__main__':
    
    det1=process_detectors('reduce_right.txt')
    det2=process_detectors('reduce_wrong.txt')
    det1.sort()
    det2.sort()
    rez = zip(det1,det2)

    f=open('comparison.txt','w')
    f.write('masks comparison:\n')
    for ind,val in enumerate(rez):
        if val[0] != val[1]:
            f.write('{0} : {1}:\n'.format(ind,val))
    f.write('Completed\n')
    f.close()