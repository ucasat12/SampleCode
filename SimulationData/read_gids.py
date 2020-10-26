import numpy as np
import sharedmem
from read_gtab import *

def grpids(fname):
    
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

def pidsconc(dir,n) :
     
    idsconc = []
    Ngrptotal = 0;
    Nidstotal = 0;
    for i in range(0,n):
          result = grpids(str(dir)+'.%d'%i)
          idsconc.append(result[0])
          Ngrptotal = Ngrptotal + result[1]['Ngroups']
          Nidstotal = Nidstotal + result[1]['Nids']
          print i
    idsconc = np.concatenate(idsconc)
    #idsconc.sort()
    return idsconc,Ngrptotal,Nidstotal   

def grpidsofgroup(dirname,n,grno):

    idsconc = []
    grptab0 = tabsconc(dirname,n)
    grplen = grptab0[2]
    grpoffset = grptab0[3]
    grplentype = grptab0[4]
    print len(grplen)
    Nlen = grplen[grno]
    Noffset = grpoffset[grno]
    Noffsetp1 = Noffset + Nlen
    Nlentype = grplentype[grno]
    i = 0
    chk = 0
    idsNids = 0
    while (idsNids < Noffset):
          idsNids = idsNids + grpids(dirname+'.'+str(int(i)))[1]['Nids']
          if (idsNids >= (Noffsetp1)):
             chk == 1
          i = i + 1
    if (chk == 1):
       i = i - 1
       idslist = grpids(dirname+'.'+str(int(i)))[0]
       kcount = Noffset - grpids(dirname+'.'+str(int(i)))[1]['Send_offset'] 
       for j in range(kcount,kcount + Nlen):
           idsconc.append(idslist[j])
    if (chk == 0):
       kcount = Noffset - idsNids
       idslist = []
       while (idsNids < Noffsetp1):
             result = grpids(dirname+'.'+str(int(i)))[0]
             idslist.append(result)
             idsNids = idsNids + grpids(dirnane+'.'+str(int(i)))[1]['Nids']
       idslist = np.concatenate(idslist)
       for j in range(kcount,kcount + Nlen):
           idsconc.append(idslist[j])
    return idsconc    

def grpidsofgroup1(dirname,dirno,n,grno):

    grptab0 = tabsconc('/home/groups_'+dirno+'/group_tab_'+dirno,n)
    grplen = grptab0[2]
    grpoffset = grptab0[3]
    grplentype = grptab0[4]

    Nlen = grplen[grno]
    Noffset = grpoffset[grno]
    Noffsetp1 = Noffset + Nlen
    Nlentype = grplentype[grno]
    i = 0
    chk = 0
    idsNids = 0
    while (idsNids < Noffset):
          idsNids = idsNids + grpids('/home/groups_'+dirno+'/group_ids_'+dirno+'.'+str(int(i)))[1]['Nids']
          if (idsNids >= (Noffsetp1)):
             chk == 1
          i = i + 1
    if (chk == 1):
       if (i > 0):
          i = i - 1
       idslist = np.array(grpids('/home/groups_'+dirno+'/group_ids_'+dirno+'.'+str(int(i)))[0])
       kcount = Noffset - grpids('/home/groups_'+dirno+'/group_ids_'+dirno+'.'+str(int(i)))[1]['Send_offset'] 
       idsconc = idslist[kcount:kcount + Nlen]
    if (chk == 0):
       if (i > 0):
          i = i - 1
       kcount = Noffset - grpids('/home/groups_'+dirno+'/group_ids_'+dirno+'.'+str(int(i)))[1]['Send_offset']
       idslist = []
       idslist.append(grpids('/home/groups_'+dirno+'/group_ids_'+dirno+'.'+str(int(i)))[0])
       while (idsNids < Noffsetp1):
             i = i + 1
             result = grpids('/home/groups_'+dirno+'/group_ids_'+dirno+'.'+str(int(i)))[0]
             idslist.append(result)
             idsNids = idsNids + grpids('/home/groups_'+dirno+'/group_ids_'+dirno+'.'+str(int(i)))[1]['Nids']
       idslist = np.array(np.concatenate(idslist))
       idsconc = idslist[kcount:kcount + Nlen]
    return idsconc  

def grpids_allgroups(dirname,dirno,n):

    grpconc = tabsconc(dirname + '/groups_'+dirno+'/group_tab_'+dirno,n)
    grplen = grpconc[2]
    grpoffset = grpconc[3]
    Nsubbfr = 0
    grpIDlist = pidsconc(dirname + '/groups_'+dirno+'/group_ids_'+dirno,n)
    grpID = []
    print 'done'
    for i in range(len(grplen)):
          grpstr = grpIDlist[0][grpoffset[i]:(grpoffset[i]+grplen[i])]
          grpstr = grpstr.copy()        
          grpID.append(grpstr.copy())
          print [i, len(grplen)]
    grpID = np.array(grpID)
    grpIDconc = np.array(grpIDlist[0])
    return grpID, grpIDconc, grpconc, grpIDlist

def grpids_allgroups_lowz(dirname_nc,dirno,ng):

    grpconc = tabsconc_lowz(dirname_nc + '/groups_'+dirno+'/group_tab_'+dirno,ng)
    grplen = grpconc[2]
    grpoffset = grpconc[3]
    Nsubbfr = 0
    grpIDlist = pidsconc(dirname_nc + '/groups_'+dirno+'/group_ids_'+dirno,ng)
    grpID = []
    print 'done'
    for i in range(len(grplen)):
          grpstr = grpIDlist[0][grpoffset[i]:(grpoffset[i]+grplen[i])]
          grpstr = grpstr.copy()        
          grpID.append(grpstr.copy())
          print [i, len(grplen)]
    grpID = np.array(grpID)
    grpIDconc = np.array(grpIDlist[0])
    return grpID, grpIDconc, grpconc, grpIDlist

