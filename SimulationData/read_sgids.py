import numpy as np
import sharedmem
from read_sgtab import * 

def subgrpids(fname):

    header = np.dtype([
           ('Ngroups',np.int32),
           ('TotNgroups',np.int32),           
           ('Nids',np.int32),
           ('TotNids',np.int64),
           ('NTask',np.int32),
           ('Send_offset',np.int32)])
    f = open(fname)
    h = np.fromfile(f, dtype=header, count=1)[0]
    ids = np.fromfile(f, dtype=np.uint64,count=h['Nids'])
    return ids,h 

def subgrppidsconc(dir,n) :
     
    subidsconc = []
    Ngrptotal = 0;
    Nidstotal = 0;
    for i in range(n):
          result = subgrpids(str(dir)+'.%d'%i)
          subidsconc.append(result[0])
          Ngrptotal = Ngrptotal + result[1]['Ngroups']
          Nidstotal = Nidstotal + result[1]['Nids']
          print i
    subidsconc = np.concatenate(subidsconc)
    return subidsconc,Ngrptotal,Nidstotal   

def subgrpstructs(dirname,dirno,n,grno):

    subgrpConc = subtabsconc(dirname+'/groups_'+dirno+'/subhalo_tab_'+dirno,n)

    Nsubstr = subgrpConc[3][grno]

    if (grno == 0):
       Nsubbfr = 0        
    else:
       Nsubbfr = sum(subgrpConc[3][0:grno])

    if (Nsubstr > 0):
       lensubstr = subgrpConc[2][Nsubbfr:(Nsubbfr+Nsubstr)]
       suboffset = subgrpConc[4][Nsubbfr:(Nsubbfr+Nsubstr)]
       idMbndoffset = subgrpConc[8][Nsubbfr:(Nsubbfr+Nsubstr)].copy()
    else:
       lensubstr = []  
       suboffset = []

    subIDlist = subgrppidsconc(dirname+'/groups_'+dirno+'/subhalo_ids_'+dirno,n)
    subPOSlist = []
    substrID = [] 
    substrPOS = []
    idMbndPos = []
    for i in range(len(lensubstr)):
          substri = subIDlist[0][suboffset[i]:(suboffset[i]+lensubstr[i])]
          substrposi = []
          substri = substri.copy()
          substrposi = substrposi.copy()
          idMbndPos.append(substrposi[0].copy())        
          substrID.append(substri.copy())
          substrPOS.append(substrposi.copy())
    substrIDconc = np.concatenate(substrID)
    substrPOSconc = np.concatenate(substrPOS)
    return substrID,substrPOS, substrIDconc, substrPOSconc, idMbndoffset, idMbndPos

def subgrpstructsnew(dirname,dirno,n,grno):

    subgrpConc = subtabsconc(dirname+'/groups_'+dirno+'/subhalo_tab_'+dirno,n)

    Nsubstr = subgrpConc[3][grno]

    if (grno == 0):
       Nsubbfr = 0        
    else:
       Nsubbfr = sum(subgrpConc[3][0:grno])

    if (Nsubstr > 0):
       lensubstr = subgrpConc[2][Nsubbfr:(Nsubbfr+Nsubstr)]
       suboffset = subgrpConc[4][Nsubbfr:(Nsubbfr+Nsubstr)]
       idMbndoffset = subgrpConc[8][Nsubbfr:(Nsubbfr+Nsubstr)].copy()
    else:
       lensubstr = []  
       suboffset = []
    subIDlist = subgrppidsconc(dirname+'/groups_'+dirno+'/subhalo_ids_'+dirno,n)
    substrID = [] 
    print 'done'
    for i in range(len(lensubstr)):
          substri = subIDlist[0][suboffset[i]:(suboffset[i]+lensubstr[i])]
          substri = substri.copy()        
          substrID.append(substri.copy())
          print [i, len(lensubstr)]
    return substrID

def subgrpids_allgroups(dirname,dirno,n):

    subgrpConc = subtabsconc(dirname+'/groups_'+dirno+'/subhalo_tab_'+dirno,n)

    Nsubstr = (subgrpConc[3]).sum()

    Nsubbfr = 0        

    if (Nsubstr > 0):
       lensubstr = subgrpConc[2][Nsubbfr:(Nsubbfr+Nsubstr)]
       suboffset = subgrpConc[4][Nsubbfr:(Nsubbfr+Nsubstr)]
    else:
       lensubstr = []  
       suboffset = []
    subIDlist = subgrppidsconc(dirname+'/groups_'+dirno+'/subhalo_ids_'+dirno,n)
    substrID = [] 
    print 'done'   
    for i in range(len(lensubstr)):
          substri = subIDlist[0][suboffset[i]:(suboffset[i]+lensubstr[i])]
          substri = substri.copy()        
          substrID.append(substri.copy())
          print [i, len(lensubstr)]
    subIDconc = np.array(subIDlist[0])
    return substrID, subIDconc, subgrpConc, subIDlist

def subgrpids_allgroups_lowz(dirname_cc,dirname_nc,dirno,ng,nsg):

    subgrpConc = subtabsconc_lowz(dirname_cc,dirname_nc,dirno,ng,nsg)
    nsub_cc = subgrpConc[10]
    nsub_nc = subgrpConc[11]
    inigrp = subgrpConc[12]

    Nsubstr = (subgrpConc[3]).sum()

    Nsubbfr = 0        

    if (Nsubstr > 0):
       lensubstr = subgrpConc[2][Nsubbfr:(Nsubbfr+Nsubstr)]
       suboffset = subgrpConc[4][Nsubbfr:(Nsubbfr+Nsubstr)]
    else:
       lensubstr = []  
       suboffset = []
    subIDlist = subgrppidsconc(dirname_nc+'/groups_'+dirno+'/subhalo_ids_'+dirno,nsg)
    substrID = [] 
    print 'done'
    ncoff = 0
    for hnocc in range(0,inigrp):
              subidsconchno,totnidshno,totnsghno = subgrppidsconc_cc(dirname_cc,dirno,hnocc)
              assert totnidshno == len(subidsconchno) 
              subIDlist[0][suboffset[ncoff]:suboffset[ncoff]+totnidshno] = subidsconchno[:]
              assert suboffset[ncoff + totnsghno] == suboffset[ncoff] + totnidshno
              ncoff = ncoff + totnsghno
              print [hnocc, totnidshno, totnsghno]
              del subidsconchno, totnidshno, totnsghno
    assert ncoff == nsub_cc
    for i in range(0,len(lensubstr)):
          substri = subIDlist[0][suboffset[i]:(suboffset[i]+lensubstr[i])]
          substri = substri.copy()        
          substrID.append(substri.copy())
          print [i, len(lensubstr)]
    subIDconc = np.array(subIDlist[0])
    return substrID, subIDconc, subgrpConc, subIDlist

def subgrppidsconc_cc(dirname_cc,dirno_cc,hno):
 
    a1 = hno/10
    a2 = hno%10
    subgrph0 = subgrptab(dirname_cc + '/isolate_' + dirno_cc + '/0'+'%d'%a1 + '%d'%a2 + '/groups_' + dirno_cc + '/subhalo_tab_' + dirno_cc + '.0')
    h = subgrph0[0]
    totnids = h['TotNids']
    totnsg = h['TotNsubgroups']
    chknids = 0
    i1 = 0
    while (chknids != totnids):
          subgrphi = subgrptab(dirname_cc + '/isolate_' + dirno_cc + '/0'+'%d'%a1 + '%d'%a2 + '/groups_' + dirno_cc + '/subhalo_tab_' + dirno_cc+'.' + str(i1))
          chknids = chknids + subgrphi[0]['Nids']
          print [i1, hno]
          i1 = i1 + 1
    subidsconchno,Ngrptotal,Nidstotal = subgrppidsconc(dirname_cc + '/isolate_' + dirno_cc + '/0'+'%d'%a1 + '%d'%a2 + '/groups_' + dirno_cc + '/subhalo_ids_'+dirno_cc,i1)
    assert totnids == Nidstotal
    assert Ngrptotal == 1
    return subidsconchno,totnids,totnsg         
   
def subgrpidsorted_allgroups(dirname,dirno,n):

    subgrpConc = subtabsconc(dirname+'/groups_'+dirno+'/subhalo_tab_'+dirno,n)

    Nsubstr = (subgrpConc[3]).sum()

    Nsubbfr = 0        

    if (Nsubstr > 0):
       lensubstr = subgrpConc[2][Nsubbfr:(Nsubbfr+Nsubstr)]
       suboffset = subgrpConc[4][Nsubbfr:(Nsubbfr+Nsubstr)]
    else:
       lensubstr = []  
       suboffset = []
    subIDlist = subgrppidsconc(dirname+'/groups_'+dirno+'/subhalo_ids_'+dirno,n)
    substrID = [] 
    substr_pivID = []
    print 'done'   
    snappidsort = np.memmap('/home/snapdir_' + dirno + 'data/snappidsort',dtype=np.uint64,mode='r')
    for i in range(10000,len(lensubstr)):
          substri = np.array(subIDlist[0][suboffset[i]:(suboffset[i]+lensubstr[i])])    
          substrID.append(substri.copy())
          piv = np.searchsorted(snappidsort,substri,side='left')
          substr_pivID.append(piv.copy())
    subIDconc = np.array(subIDlist[0])
    return substrID, subIDconc, substr_pivID
