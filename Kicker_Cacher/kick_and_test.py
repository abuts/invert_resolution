#!/usr/bin/python
#
import multiprocessing
import threading
import time
import sys
import argparse
import numpy
# Emulate multithreaded Horace workflow to check CEPH validity

class arr_holder(threading.Thread):
    """ Class to generate test data"""

    def __init__(self,n_chunk=0,chunk_size=1024,n_chunks=None):
        threading.Thread.__init__(self)

        self.rc = 0
        self.data_lock = threading.Condition()
        self.write_lock= threading.Condition()
        self.generating = True
        self.ready_to_write = False


        self.reset_data(n_chunk,chunk_size,n_chunks)

    def reset_data(self,n_chunk,chunk_size,n_chunks):
        self._chunk_size = chunk_size
        self._arr_holder = []

        if n_chunks is None:
            n_chunks = 1
        self._n_chunks = n_chunks
        #
        self.gen_test_data(666)


    def gen_test_data(self,n_chunk):
        numpy.random.seed(n_chunk)
        tmp = self._arr_holder
        #self._arr_holder = numpy.random.normal(-1,1,self._chunk_size);
        self._arr_holder = numpy.ones(self._chunk_size)*n_chunk;

        return tmp;

    def swap(self,other_arr=None):

        tmp = self._arr_holder
        if other_arr:
            self._arr_holder = other_arr
        else:
            self._arr_holder = []
        #self.data_lock.notify()
        return tmp

    def run(self):
       for nch in xrange(0,self._n_chunks):
           while self.generating:
                self.data_lock.acquire()
                self.gen_test_data(nch)
                #print('generated ChN: {0}'.format(nch))
                self.generating = False
                self.data_lock.release()
           self.write_lock.acquire()
           self.ready_to_write = True
           self.write_lock.release()

           while not self.generating :
               try:
                   self.data_lock.wait()
               except RuntimeError:
                   pass
       self.write_lock.acquire()
       self.ready_to_write = True
       self.write_lock.release()

               # try:
               #     self.write_lock.wait()
               # except RuntimeError:
                #    pass



#
input_data = arr_holder()


class progress_rep:
    """ report progress of the chunk-reading job"""

    def __init__(self,all_size,n_threads=None):
        """ Constructor, initiating progress reporting
            Usage:
            prog = progress_rep(total_size, n_threads)
            where:
            total_size -- the data size to be read
            n_threads  -- number of threads will be used
                          to read the file.
        """
        self._all_size = all_size
        self._tick_size = all_size/100
        if n_threads is None:
            self._n_threads = 1
        else:
            self._n_threads = n_threads

        self._tick_cnt = 0
        self._tick_barrier = self._tick_size
        self._prev_size = 0
        self._start_time = time.time()
        self._prev_time = self._start_time
        self._run_time  = self._start_time
    #
    def check_progress(self,cur_size):
        """ check current progress and report if achieved the limit for reporting.

            Adjust the reporting limit to the next position if current limit was achieved.
            Usage:
            prog.check_progress(size)
            where:
            size -- the amount of data read by current thread.
        """

        self._tick_cnt += 1
        if cur_size >=self._tick_barrier:
            self._report_progress(cur_size)
            self._tick_barrier += self._tick_size
    #
    def _report_progress(self,cur_size):
        """ Internal method used to check time and print progress message"""

        pers = 100*float(cur_size)/float(self._all_size);
        self._prev_time  = self._run_time
        self._run_time   = time.time()
        tot_size = float(cur_size*self._n_threads)/(1024*1024)
        block_size = float(cur_size - self._prev_size)*self._n_threads/(1024*1024)
        tot_time = self._run_time-self._start_time
        block_time = self._run_time-self._prev_time
        Av_speed = tot_size / tot_time
        Loc_speed = block_size / block_time
        sys.stdout.write(\
           "Read: {0:.1f}MB Completed: {1:3.1f}% Av Speed: {2:3.2f}MB/s: Loc speed: {3:3.2f}MB/s\r"\
           .format(tot_size,pers,Av_speed,Loc_speed))
        sys.stdout.flush()
        self._prev_size = cur_size




def write_test_file(filename,fielsize,chunk_size):
    """ Write test data file"""

    fd = open(filename,'wb')
    if fd<0:
        raise RuntimeError("Can not open test file %s".format(filename))
    n_chunks = fielsize/chunk_size

    input_data.reset_data(0,chunk_size,n_chunks)
    input_data.generating = True
    input_data.start()
    #input_data.run()
    while input_data.generating:
       pass

    time1= time.time()
    size1 = 0
    tot_size = float(n_chunks*chunk_size*8)/(1024*1024)

    for n_ch in xrange(0,n_chunks):
        input_data.data_lock.acquire()
        test_data = input_data.swap()
        input_data.generating = True
        input_data.data_lock.release()
        try:
          input_data.data_lock.notify()
        except RuntimeError:
          pass

        while not input_data.ready_to_write:
            pass
        fd.write(test_data)
        input_data.write_lock.acquire()
        input_data.ready_to_write = False
        input_data.write_lock.release()

        #print 'nch=',n_ch,' td=',test_data[0:5]
        cur_size = float(n_ch*chunk_size*8)/(1024*1024)
        pers     = 100*cur_size /tot_size
        cur_time = time.time()
        ds = cur_size - size1
        dt = cur_time - time1
        Wr_speed = ds / dt
        sys.stdout.write(\
           "Wrote: {0:.1f}MB Completed: {1:3.1f}% Write Speed: {2:3.2f}MB/s\r"\
           .format(cur_size,pers,Wr_speed))
        sys.stdout.flush()
        size1 = cur_size
        time1 = cur_time


    input_data.generating = True
    fd.close()
    input_data.join()
    sys.stdout.write("\n")
    sys.stdout.flush()



def chunk_wrapper(pout,filename,chunk_num,chunk_size,n_chunks):
    """ multiprocessing piped wrapper around read_and_test_chunk"""
    ok,n_faling_chunk = read_and_test_chunk(filename,chunk_num,chunk_size,n_chunks)
    pout.send([ok,n_faling_chunk])
    pout.close()

def read_and_test_chunk(filename,chunk_num,chunk_size,n_chunks):
    """ read and verify chunk of binary data.

        Used for read data in a single thread.
        Inputs:
        filename -- name of the binary file to read
        start    -- initial position to read data from
        size     -- the number of bytes to read from the file
        buf_size -- the size of the buffer to use while reading the data
        progress -- True if progress messages should be printed and False if not
        n_workers -- number of threads to read file. Used to estimate the progress
                    of multithreaded job.
    """
    success=False
    checker = arr_holder(0,chunk_size)

    fh = open(filename,'rb')
    if fh<0:
        raise ValueError("Can not open file: "+filename)
    start = chunk_num*chunk_size*8L
    fh.seek(start,0)
    with open(filename, "rb") as f:
        for nch in xrange(chunk_num,chunk_num+n_chunks+1):
            data = numpy.frombuffer(fh.read(chunk_size*8));
            checker.gen_test_data(nch)
            sample = checker.swap()
            ok = numpy.array_equal(sample,data)
            if not ok:
                return (False,nch)
    return (True,0)


def check_file(filename,file_size,chunk_size,n_threads):
    """ Read special file using multiple threads and check file integrity 

        Input:
        dictionary, generated by ArgumentParser, containing the
        name of the file to read and some auxiliary information
        about reading job parameters.
    """

    # Estimate the file size
    fh = open(filename,'rb')
    if fh<0:
        raise ValueError("Can not open file: "+filename)
    fh.close()

    # Evaluate the parameters of the file reading jobs.
    n_chunks = file_size/chunk_size
    n_chunks_per_thread = int(n_chunks/n_threads)
    nch_per_thread = [n_chunks_per_thread]*n_threads
    n_chunks_tot_distr = numpy.sum(nch_per_thread)

    for i in xrange(0,n_threads):
        if n_chunks_tot_distr < n_chunks :
            nch_per_thread[i] = nch_per_thread[i]+1
            n_chunks_tot_distr= n_chunks_tot_distr+1
        else:
            break
    start_chunk_num = numpy.append([0],numpy.cumsum(nch_per_thread))

    #---------------------------------------------------------------------------------------------
    # Start parallel jobs:
    job_list  = []
    result_p  = []
    #read_chunk(filename,chunk_beg[0],chunk_size[0],buf_size,True,n_threads)
    for nthr in xrange(n_threads):
        if nthr == 0:
            log=True
        else:
            log=False
        parent_conn, child_conn = multiprocessing.Pipe()
        p = multiprocessing.Process(target=chunk_wrapper, args=(child_conn,filename,start_chunk_num[nthr],chunk_size,nch_per_thread[nthr]-1,))
        p.start()
        result_p.append(parent_conn)
        job_list.append(p)
    # Wait for jobs to finish.
    ok = True
    for p,proc_out in zip(job_list,result_p):
        out = proc_out.recv()
        if not out[0]:
            ok = False
            print('Error reading chunk N {0}'.format(out[1]))
        p.join()
    return ok
#------------------------------------------------

if __name__ == '__main__':
    """test io operations over large files on CEPH"""

    parser = argparse.ArgumentParser(add_help=True, version='0.1',description='test Horace-like IO operations on CEPHs')
    parser.add_argument('-n',action='store', dest='nthreads',type=int,default=16,help='number of threads to process file. Default is 16 threads')
    parser.add_argument('-b',action='store', dest='buffer',type=int,default=4096,help='Buffer size to read each chunk of data. Default is 4096 bytes.')

    args = vars(parser.parse_args())
    #process_file(args)
    #chunk_size = 5*10000000
    chunk_size =  5*1000000
    filesize = 100*chunk_size

    write_test_file("test_file.tmp",filesize,chunk_size)
#    ok=read_and_test_chunk("test_file.tmp",75,chunk_size,25)
    ok = check_file("test_file.tmp",filesize,chunk_size,4)
