import numpy, sharedmem, time
from load_data import *

def sortsnapshot(flagno,n2div,dirname1,preci):
    #precb = (preci + 1)*4
    snappid_argsorted = numpy.memmap(dirname1 + '/snappid_argsortednew',dtype=numpy.int,mode='r')
    n1 = len(snappid_argsorted)
    del snappid_argsorted
    pid_arg = loadsnapshot(4,dirname1,n2div,1,preci)
    if (flagno == 0):
       f = open(dirname1 + '/snappidsort','w')
       snap_tst = loadsnapshot(0,dirname1,n2div,1,preci)
    if (flagno == 1):
       f = open(dirname1 + '/snappossort','w')
       snap_tst = loadsnapshot(1,dirname1,n2div,0,preci)
    if (flagno == 2):
       f = open(dirname1 + '/snapptypesort','w')
       snap_tst = loadsnapshot(2,dirname1,n2div,1,preci)
    if (flagno == 3):
       f = open(dirname1 + '/snapmasssort','w')
       snap_tst = loadsnapshot(3,dirname1,n2div,1,preci)
    chunksize = numpy.int(n1/n2div) 
    slices = [slice(i, i+chunksize) for i in range(0, n1-chunksize,
chunksize)]
    for i in range(0, n2div + 1):
             snap_sorttst = snap_tst[pid_arg[slices[i]]]
             snap_sorttst.tofile(f)
    f.close()
    print 'done'
    del f
    return n1
