#Loads the processed data

import numpy, sharedmem, time

def loadsnapshot(flagno,dirname,n2div,SMbool,preci):

    precb = (preci + 1)*4
    snappid_notsorted = numpy.memmap(dirname+'/snappid_notsorted',dtype=numpy.uint64,mode='r')
    n1 = len(snappid_notsorted)
    del snappid_notsorted
    if (flagno == 0):
       #snap_notsorted = numpy.zeros(n1,dtype = numpy.uint64)
       shp1 = n1
       dtp1 = numpy.uint64
       dtp2 = numpy.uint64
       str1 = 'snappid_notsorted'
    if (flagno == 1):
       #snap_notsorted = numpy.zeros((n1,3),dtype = numpy.float64)
       shp1 = (n1,3)
       dtp1 = 'f'+str(precb)
       dtp2 = ('f'+str(precb),3)
       str1 = 'snappos_notsorted'
    if (flagno == 2):
       #snap_notsorted = numpy.zeros(n1,dtype = numpy.int32)
       shp1 = n1
       dtp1 = numpy.int32
       dtp2 = numpy.int32
       str1 = 'snapptype_notsorted'
    if (flagno == 3):
       #snap_notsorted = numpy.zeros(n1,dtype = numpy.float64)
       shp1 = n1
       dtp1 = 'f'+str(precb)
       dtp2 = 'f'+str(precb)
       str1 = 'snapmass_notsorted'
    if (flagno == 4):
       #snap_notsorted = numpy.zeros(n1,dtype = numpy.int)
       shp1 = n1
       dtp1 = numpy.int
       dtp2 = numpy.int
       str1 = 'snappid_argsortednew'
    if (SMbool == 0):
       snap_notsorted = numpy.zeros(shp1,dtype = dtp1)
    if (SMbool == 1):
       snap_notsorted = sharedmem.empty(shp1,dtp1) 
    chunksize = numpy.int(n1/n2div)
    slices = [slice(i, i+chunksize) for i in range(0, n1-chunksize,chunksize)]
    slices.append(slice(chunksize*n2div,n1))
    f = open(dirname+'/'+str1)
    for i in range(0,n2div):
          snap_notsorted[slices[i]] = numpy.fromfile(f,dtype = dtp2,count = chunksize)
          print i
    snap_notsorted[slices[n2div]] = numpy.fromfile(f,dtype = dtp2,count = n1 - (n2div*chunksize))
    del slices, chunksize, n1, shp1, dtp1, dtp2, str1
    return snap_notsorted

def loadsnapshot_sorted(flagno,dirname,n2div,SMbool,preci):

    precb = (preci + 1)*4
    snappidsort = numpy.memmap(dirname+ '/snappidsort',dtype=numpy.uint64,mode='r')
    n1 = len(snappidsort)
    del snappidsort
    if (flagno == 0):
       #snap_notsorted = numpy.zeros(n1,dtype = numpy.uint64)
       shp1 = n1
       dtp1 = numpy.uint64
       dtp2 = numpy.uint64
       str1 = 'snappidsort'
    if (flagno == 1):
       #snap_notsorted = numpy.zeros((n1,3),dtype = numpy.float64)
       shp1 = (n1,3)
       dtp1 = 'f'+str(precb)
       dtp2 = ('f'+str(precb),3)
       str1 = 'snappossort'
    if (flagno == 2):
       #snap_notsorted = numpy.zeros(n1,dtype = numpy.int32)
       shp1 = n1
       dtp1 = numpy.int32
       dtp2 = numpy.int32
       str1 = 'snapptypesort'
    if (flagno == 3):
       #snap_notsorted = numpy.zeros(n1,dtype = numpy.float64)
       shp1 = n1
       dtp1 = 'f'+str(precb)
       dtp2 = 'f'+str(precb)
       str1 = 'snapmasssort'
    if (flagno == 4):
       shp1 = n1
       dtp1 = numpy.int32
       dtp2 = numpy.int32
       str1 = 'snap_DMmask'
    if (flagno == 6):
        snap_DMpos = numpy.memmap(dirname+'/snapDMpos',dtype = ('f'+str(precb),3),mode='r')
        n1 = len(snap_DMpos)
        shp1 = (n1,3)
        dtp1 = 'f'+ str(precb)
        dtp2 = ('f' + str(precb),3)
        str1 = 'snapDMpos'
    print n1
    if (SMbool == 0):
       snap_sort = numpy.zeros(shp1,dtype = dtp1)
    if (SMbool == 1):
       snap_sort = sharedmem.empty(shp1,dtp1) 
    chunksize = numpy.int(n1/n2div)
    slices = [slice(i, i+chunksize) for i in range(0, n1-chunksize,chunksize)]
    slices.append(slice(chunksize*n2div,n1))
    f = open(dirname+'/'+str1)
    for i in range(0,n2div):
          snap_sort[slices[i]] = numpy.fromfile(f,dtype = dtp2,count = chunksize)
          print i
    snap_sort[slices[n2div]] = numpy.fromfile(f,dtype = dtp2,count = n1 - (n2div*chunksize))
    del slices, chunksize, n1, shp1, dtp1, dtp2, str1
    return snap_sort

#############

def loadgrpsubgrp_data(subgrpbool,flagno,dirname,n2div,SMbool,preci):

    precb = (preci + 1)*4
    if (subgrpbool == 0):
       gstr2 = 'grpids_'
    if (subgrpbool == 1):
       gstr2 = 'subgrpids_'
    ids_correct = numpy.memmap(dirname+'/'+gstr2+'idscorrect',dtype=numpy.uint64)
    n2 = len(ids_correct)
    del ids_correct
    if (flagno == 0):
       gshp1 = n2
       gdtp1 = numpy.uint64
       gdtp2 = numpy.uint64
       gstr1 = 'idscorrect'
    if (flagno == 1):
       gshp1 = (n2,3)
       gdtp1 = 'f'+str(precb)
       gdtp2 = ('f'+str(precb),3)
       gstr1 = 'poscorrect'
    if (flagno == 2):
       gshp1 = n2
       gdtp1 = numpy.int32
       gdtp2 = numpy.int32
       gstr1 = 'ptypecorrect'
    if (flagno == 3):
       gshp1 = n2
       gdtp1 = 'f'+str(precb)
       gdtp2 = 'f'+str(precb)
       gstr1 = 'masscorrect'
    if (flagno == 4):
       gshp1 = n2
       gdtp1 = numpy.int32
       gdtp2 = numpy.int32
       gstr1 = 'DMmask'
    if (flagno == 5):
       gshp1 = n2
       gdtp1 = numpy.int32
       gdtp2 = numpy.int32
       gstr1 = 'starmask'
    if (SMbool == 0):
       gsnap_corrected = numpy.zeros(gshp1,dtype = gdtp1)
    if (SMbool == 1):
       gsnap_corrected = sharedmem.empty(gshp1,gdtp1)
    chunksize = numpy.int(n2/n2div)
    slices = [slice(i, i+chunksize) for i in range(0, n2-chunksize,chunksize)]
    slices.append(slice(chunksize*n2div,n2)) 
    f = open(dirname+'/'+gstr2+gstr1)
    for i in range(0,n2div):
          gsnap_corrected[slices[i]] = numpy.fromfile(f,dtype = gdtp2,count = chunksize)
          print i
    gsnap_corrected[slices[n2div]] = numpy.fromfile(f,dtype = gdtp2,count = n2 - (n2div*chunksize))
    del slices, chunksize, n2, gshp1, gdtp1, gdtp2, gstr1, gstr2
    return gsnap_corrected

#############

def loadgrpsubgrp_idspiv(subgrpbool,dirname,n2div):

    if (subgrpbool == 0):
       str1 = 'grpids_piv'
       str2 = 'grpids_'
    if (subgrpbool == 1):
       str1 = 'subgrpids_piv'
       str2 = 'subgrpids_' 
    ids_piv = numpy.memmap(dirname+'/'+str1,dtype=numpy.int)
    n2 = len(ids_piv)
    del ids_piv
    ids_piv  = sharedmem.empty(n2,numpy.int) 
    f = open(dirname+'/'+str1)
    chunksize = numpy.int(n2/n2div)
    slices = [slice(i, i+chunksize) for i in range(0, n2-chunksize,chunksize)]
    slices.append(slice(chunksize*n2div,n2))   
    for i in range(0,n2div):
          ids_piv[slices[i]] = numpy.fromfile(f,dtype = numpy.int,count = chunksize)   
          print i
    ids_piv[slices[n2div]] = numpy.fromfile(f,dtype = numpy.int,count = n2 - (n2div*chunksize))  
    return ids_piv

#############
def loadfinalshared(subgrpbool,flagno,ids_piv,n2div,SMbool,preci):

    n2 = len(ids_piv)
    precb = (preci + 1)*4
    if (subgrpbool == 0):
       str1 = 'grpids_piv'
       str2 = 'grpids_'
    if (subgrpbool == 1):
       str1 = 'subgrpids_piv'
       str2 = 'subgrpids_'
    if (flagno == 0):
       gshp1 = n2
       gdtp1 = numpy.uint64
       gdtp2 = numpy.uint64
       gstr1 = 'pidcorrect'
    if (flagno == 1):
       gshp1 = (n2,3)
       gdtp1 = 'f'+str(precb)
       gdtp2 = ('f'+str(precb),3)
       gstr1 = 'poscorrect'
    if (flagno == 2):
       gshp1 = n2
       gdtp1 = numpy.int32
       gdtp2 = numpy.int32
       gstr1 = 'ptypecorrect'
    if (flagno == 3):
       gshp1 = n2
       gdtp1 = 'f'+str(precb)
       gdtp2 = 'f'+str(precb)
       gstr1 = 'masscorrect'
    grpids_psome  = sharedmem.empty(gshp1,gdtp1)
    return grpids_psome, str2, gstr1
