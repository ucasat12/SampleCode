import numpy
from read_gtab import *
from read_gids import *
from load_data import *

def getidslen_grp(dirname,dirno,dirname1,ng,n2div,SMbool,preci,ptype):
   
    grpids_piv = loadgrpsubgrp_idspiv(0,dirname1,n2div)
    grpConc = tabsconc(dirname+'/groups_'+dirno+'/group_tab_'+dirno,ng)
    cc = 0
    lenConc = grpConc[2]
    offsetConc = grpConc[3]
    Ngroups = len(grpConc[0])
    if (ptype == 1):
       mask_grp = loadgrpsubgrp_data(0,4,dirname1,n2div,SMbool,preci)
       strmask = 'grp_DMpiv'
    if (ptype == 4):
       mask_grp = loadgrpsubgrp_data(0,5,dirname1,n2div,SMbool,preci)
       strmask = 'grp_starpiv'
    tstma = numpy.ma.array(grpids_piv, mask = mask_grp)
    grplen = []
    offseti = 0
    f = open(dirname1 + '/' + strmask,'w' )
    for i in range(0,Ngroups): 
          li = lenConc[i]
          tst4ma = tstma[offseti : offseti + li]
          pivma = tst4ma[~tst4ma.mask]
          pivdata = pivma.data
          grplen.append([i, len(pivdata)])
          pivdata.tofile(f)
          offseti = offseti + li
          print i, Ngroups
    f.close()
    del f
    grplen = numpy.array(grplen)
    if (ptype == 1):
       grplen.tofile(dirname1+'/gsl_DMgrps')
    if (ptype == 4):
       grplen.tofile(dirname1+'/gsl_stargrps')
    return grplen

def getidslen_grp_lowz(dirname,dirno,dirname1,ng,n2div,SMbool,preci,ptype):
       
    grpids_piv = loadgrpsubgrp_idspiv(0,dirname1,n2div)
    grpConc = tabsconc_lowz(dirname+'/groups_'+dirno+'/group_tab_'+dirno,ng)
    cc = 0
    lenConc = grpConc[2]
    offsetConc = grpConc[3]
    Ngroups = len(grpConc[0])
    if (ptype == 1):
       mask_grp = loadgrpsubgrp_data(0,4,dirname1,n2div,SMbool,preci)
       strmask = 'grp_DMpiv'
    if (ptype == 4):
       mask_grp = loadgrpsubgrp_data(0,5,dirname1,n2div,SMbool,preci)
       strmask = 'grp_starpiv'
    tstma = numpy.ma.array(grpids_piv, mask = mask_grp)
    grplen = []
    offseti = 0
    f = open(dirname1 + '/' + strmask,'w' )
    for i in range(0,Ngroups): 
          li = lenConc[i]
          tst4ma = tstma[offseti : offseti + li]
          pivma = tst4ma[~tst4ma.mask]
          pivdata = pivma.data
          grplen.append([i, len(pivdata)])
          pivdata.tofile(f)
          offseti = offseti + li
          print i, Ngroups
    f.close()
    del f
    grplen = numpy.array(grplen)
    if (ptype == 1):
       grplen.tofile(dirname1+'/gsl_DMgrps')
    if (ptype == 4):
       grplen.tofile(dirname1+'/gsl_stargrps')
    return grplen

