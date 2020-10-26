import numpy as np

def groupcoord(fname,preci) :
    
    precb = (preci + 1)*4
    header = np.dtype([
                       ('Npart',(np.uint32,6)),
                       ('Massarr',(np.double,6)),
                       ('Time',np.double),
                       ('Redshift',np.double),
                       ('FlagSfr',np.uint32),
                       ('FlagFeedback',np.uint32),
                       ('Nall',(np.uint32,6)),
                       ('FlagCooling',np.uint32),
                       ('NumFiles',np.uint32),
                       ('BoxSize',np.double),
                       ('Omega0',np.double),
                       ('OmegaLambda',np.double),
                       ('HubbleParam',np.double),
                       ('Na11HW',(np.uint32,6)),
                       ('fill',(np.int8,72))])
    assert header.itemsize == 256
    f = open(fname)
    int1 = np.fromfile(f, dtype=np.uint32, count=1)[0]
    h = np.fromfile(f, dtype=header, count = 1)[0]
    int2 = np.fromfile(f, dtype=np.uint32, count=1)[0]
    assert int1 == int2
    totNpart = h['Npart'].sum()
    int3 = np.fromfile(f, dtype=np.uint32, count=1)[0]
    pos = np.fromfile(f, dtype = ('f'+str(precb),3), count = totNpart)
    int4 = np.fromfile(f, dtype=np.uint32, count=1)[0]
    assert int3 == int4    
    f.seek(totNpart*3*precb + 8,1)
    int5 = np.fromfile(f, dtype=np.uint32, count=1)[0]
    pid = np.fromfile(f, dtype = np.uint64, count = totNpart)
    int6 = np.fromfile(f, dtype=np.uint32, count=1)[0]
    assert int5 == int6
    Nm = 0
    for i in range(0,6):
          if (h['Massarr'][i] == 0):
             Nm = Nm + (h['Npart'][i])
    if (Nm > 0):
       int7 = np.fromfile(f,dtype=np.uint32,count=1)[0]
       vpmass = np.fromfile(f, dtype='f'+str(precb), count = Nm)
       int8 = np.fromfile(f,dtype=np.uint32,count=1)[0]
       assert int7 == int8
    ptype = np.zeros(totNpart,dtype = np.uint32)
    pmass = np.zeros(totNpart,dtype = 'f'+str(precb))
    count = 0
    vpcount = 0
    for i in range(len(h['Npart'])):
          ptype[count : count + h['Npart'][i]] = i
          if (h['Massarr'][i] == 0):
             pmass[count : count + h['Npart'][i]] = (vpmass[vpcount : vpcount + h['Npart'][i]]).astype('f'+str(precb))
             vpcount = vpcount + h['Npart'][i]
          if (h['Massarr'][i] > 0):
             pmass[count : count + h['Npart'][i]] = h['Massarr'][i]
          count = count + h['Npart'][i]
    print count, totNpart, vpcount, Nm  
    return pid, pos, ptype, h, totNpart, pmass

def spacecount_snapshot(dirname,dirno,n):
    
    header = np.dtype([
                       ('Npart',(np.int32,6)),
                       ('Massarr',(np.double,6)),
                       ('Time',np.double),
                       ('Redshift',np.double),
                       ('FlagSfr',np.int32),
                       ('FlagFeedback',np.int32),
                       ('Nall',(np.int32,6)),
                       ('FlagCooling',np.int32),
                       ('NumFiles',np.int32),
                       ('BoxSize',np.double),
                       ('Omega0',np.double),
                       ('OmegaLambda',np.double),
                       ('HubbleParam',np.double),
                       ('Na11HW',(np.int32,6)),
                       ('fill',(np.int8,72))])
    assert header.itemsize == 256
    Nidstot = 0
    for i in range(0,n):
          f = open(str(dirname)+'.%d'%i)
          int1 = np.fromfile(f, dtype=np.int32, count=1)[0]
          h = np.fromfile(f, dtype=header, count = 1)[0]
          int2 = np.fromfile(f, dtype=np.int32, count=1)[0]
          assert int1 == int2
          Nidstot = Nidstot + h['Npart'].sum()
          print [i, h['Npart'].sum()]
    print 'particle count completed : ' + str(Nidstot)
    snappid_notsorted = np.memmap('/home/vat/physics/snapdir_'+str(dirno)+'data/snappid_notsorted',shape = Nidstot,dtype=np.uint64,mode='w+')
    snappos_notsorted = np.memmap('/home/vat/physics/snapdir_'+str(dirno)+'data/snappos_notsorted',shape = (Nidstot,3),dtype=np.float64,mode='w+')
    snapptype_notsorted = np.memmap('/home/vat/physics/snapdir_'+str(dirno)+'data/snapptype_notsorted',shape = Nidstot,dtype=np.int32,mode='w+')
    snapmass_notsorted = np.memmap('/home/vat/physics/snapdir_'+str(dirno)+'data/snapmass_notsorted',shape = Nidstot,dtype=np.float64,mode='w+')
    idpos = 0
    for i in range(0,n):
             print [i, 'bfr']
             pid, pos, ptype, h, totNpart, pmass = groupcoord(str(dirname)+'.%d'%i)
             print [i, 'afr']
             snappid_notsorted[idpos : idpos + totNpart] = pid
             snappos_notsorted[idpos : idpos + totNpart] = pos
             snapptype_notsorted[idpos : idpos + totNpart] = ptype
             snapmass_notsorted[idpos : idpos + totNpart] = pmass
             idpos = idpos + totNpart
             print [i, idpos, Nidstot]
             del pid, pos, ptype, pmass, h, totNpart
    del snappid_notsorted, snappos_notsorted, snapptype_notsorted, snapmass_notsorted
    snappid_notsorted = np.memmap('/home/vat/physics/snapdir_'+str(dirno)+'data/snappid_notsorted',shape = Nidstot,dtype=np.uint64,mode='r')
    snappos_notsorted = np.memmap('/home/vat/physics/snapdir_'+str(dirno)+'data/snappos_notsorted',shape = (Nidstot,3),dtype=np.float64,mode='r')
    snapptype_notsorted = np.memmap('/home/vat/physics/snapdir_'+str(dirno)+'data/snapptype_notsorted',shape = Nidstot,dtype=np.int32,mode='r')
    snapmass_notsorted = np.memmap('/home/vat/physics/snapdir_'+str(dirno)+'data/snapmass_notsorted',shape = Nidstot,dtype=np.float64,mode='r') 
    return snappid_notsorted, snappos_notsorted, snapptype_notsorted, snapmass_notsorted, Nidstot, idpos

def groupcoord_header(fname) :

    header = np.dtype([
                       ('Npart',(np.int32,6)),
                       ('Massarr',(np.double,6)),
                       ('Time',np.double),
                       ('Redshift',np.double),
                       ('FlagSfr',np.int32),
                       ('FlagFeedback',np.int32),
                       ('Nall',(np.int32,6)),
                       ('FlagCooling',np.int32),
                       ('NumFiles',np.int32),
                       ('BoxSize',np.double),
                       ('Omega0',np.double),
                       ('OmegaLambda',np.double),
                       ('HubbleParam',np.double),
                       ('Na11HW',(np.int32,6)),
                       ('fill',(np.int8,72))])
    assert header.itemsize == 256
    f = open(fname)
    int1 = np.fromfile(f, dtype=np.int32, count=1)[0]
    h = np.fromfile(f, dtype=header, count = 1)[0]
    int2 = np.fromfile(f, dtype=np.int32, count=1)[0]
    assert int1 == int2
    totNpart = h['Npart'].sum()
    return h, totNpart

def groupcoord_dmpidsum(dirname,dirno,ns):
    dmcc = 0
    for i in range(0,ns):
          h, totNpart1 = groupcoord_header(dirname+'/snapdir_' + dirno + '/snapshot_'+dirno+'.'+str(i))
          dmcc = dmcc + h['Npart'][1]
          print [i,dmcc]
    return dmcc

def groupcoord_ids(fname,Ncnt) :

    f = open(fname)
    f.seek(256 + 8 + (2*Ncnt*3*8) + (2*8),1)
    int5 = np.fromfile(f, dtype=np.int32, count=1)[0]
    pid = np.fromfile(f, dtype = np.uint64, count = Ncnt)
    int6 = np.fromfile(f, dtype=np.int32, count=1)[0]
    assert int5 == int6
    totNpart = len(pid)
    return pid, totNpart 

def spacecount_snapshotonlyids(dirname,dirno,n):

    Nidstot = 0
    Nids_ind = np.zeros(n,dtype=np.uint64)
    for i in range(0,n):
          h, totNpart = groupcoord_header(str(dirname)+'.%d'%i)
          Nidstot = Nidstot + h['Npart'].sum()
          Nids_ind[i] = totNpart
          print [i, h['Npart'].sum(), totNpart]
    print 'particle count completed : ' + str(Nidstot)
    snappid_notsorted = np.memmap('/home/vat/physics/snapdir_'+str(dirno)+'data/snappid_notsorted_ids',shape = Nidstot,dtype=np.uint64,mode='w+')
    idpos = 0
    del totNpart
    for i in range(0,n):
             print [i, 'bfr']
             pid, totNpart = groupcoord_ids(str(dirname)+'.%d'%i,Nids_ind[i])
             print [i, 'afr']
             snappid_notsorted[idpos : idpos + totNpart] = pid
             idpos = idpos + totNpart
             print [i, idpos, Nidstot]
    del snappid_notsorted
    snappid_notsorted = np.memmap('/home/vat/physics/snapdir_'+str(dirno)+'data/snappid_notsorted_ids',shape = Nidstot,dtype=np.uint64,mode='r')
    return snappid_notsorted, Nidstot, idpos
    
def snapconcids(dir,n):
   
    pidConc = []
    posConc = []
    ptypeConc = []
    pmassConc = []
    Ntotal = 0
    for i in range(n):
          result = groupcoord(str(dir)+'.%d'%i)
          pidConc.append(result[0])
          Ntotal = Ntotal + result[4]
    pidConc = array(np.concatenate(pidConc))
    numConc = array([i for i in range(Ntotal)])
    par = np.empty(shape=Ntotal,
                   dtype = [
                     ('pid', np.uint64),
                     ('num',np.uint64)])
    par['pid'] = pidConc
    par['num'] = numConc
    par.sort(order='pid')
    numsort = par['num'].copy()
    pidsort = par['pid'].copy()
    del par, pidConc
    print ' Sorting Completed '
    possort = array(zeros((Ntotal,3),dtype=np.float64))
    Ntotal = 0
    for i in range(n):
          result = groupcoord(str(dir)+'.%d'%i)
          posConc.append(result[1])
          Ntotal = Ntotal + result[4]
    posConc = array(np.concatenate(posConc))
    for i in range(Ntotal):
          possort[i] = (posConc[numsort[i]]).copy()
    del posConc, result
    ptypesort = array(zeros(Ntotal,dtype=np.uint64))
    pmasssort = array(zeros(Ntotal,dtype=np.float64))
    Ntotal = 0
    for i in range(n):
          result = groupcoord(str(dir)+'.%d'%i)
          ptypeConc.append(result[2])
          pmassConc.append(result[5])
          Ntotal = Ntotal + result[4]
    ptypeConc = array(np.concatenate(ptypeConc))
    pmassConc = array(np.concatenate(pmassConc))
    for i in range(Ntotal):
             ptypesort[i] = (ptypeConc[numsort[i]])
             pmasssort[i] = (pmassConc[numsort[i]])
    del ptypeConc, pmassConc, result
    
    return pidsort, possort, ptypesort, pmasssort, Ntotal

def snapconcpos(dir,n) :
     
    pidConc = []
    posConc = []
    Ntotal = 0
    for i in range(n):
          result = groupcoord(str(dir)+'.%d'%i)
          pidConc.append(result[0])
          posConc.append(result[1])
          Ntotal = Ntotal + result[4]
    print 'Ntotal = ' + str(Ntotal)
    pidConc = np.concatenate(pidConc)
    posConc = np.concatenate(posConc)
    
    par = np.empty(shape=Ntotal,
               dtype=[
                 ('pid', np.uint64),
                 ('pos', (np.float32, 3))])
    par['pid'] = pidConc
    par['pos'] = posConc
    par.sort(order='pid')
    pidConc = []
    posConc = []
    return par, Ntotal

def snapconcptype(dir,n) :
     
    pidConc = []
    ptypeConc = []
    Ntotal = 0
    for i in range(n):
          result = groupcoord(str(dir)+'.%d'%i)
          pidConc.append(result[0])
          ptypeConc.append(result[2])
          Ntotal = Ntotal + result[4]
    pidConc = np.concatenate(pidConc)
    ptypeConc = np.concatenate(ptypeConc) 

    par = np.empty(shape=Ntotal,
               dtype=[
                 ('pid', np.uint64),
                 ('ptype',np.int32)])
    par['pid'] = pidConc
    par['ptype'] = ptypeConc
    par.sort(order='pid')
    pidConc = []
    ptypeConc = []
    return par, Ntotal

def snapconcpmass(dir,n) :
     
    pidConc = []
    pmassConc = []
    Ntotal = 0
    for i in range(n):
          result = groupcoord(str(dir)+'.%d'%i)
          pidConc.append(result[0])
          pmassConc.append(result[5])
          Ntotal = Ntotal + result[4]
    pidConc = np.concatenate(pidConc)
    pmassConc = np.concatenate(pmassConc) 

    par = np.empty(shape=Ntotal,
               dtype=[
                 ('pid', np.uint64),
                 ('pmass',double)])
    par['pid'] = pidConc
    par['pmass'] = pmassConc
    par.sort(order='pid')
    pidConc = []
    pmassConc = []
    return par, Ntotal
