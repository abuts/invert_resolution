#/usr/bin/python
from multiprocessing import Pool

def process_file(filename):
    """ """
    fh = fopen(filename,'rb')
    if fh<0:
        raise ValueError("Can not open file: "+filename)


if __name__ == '__main__':
    argi 