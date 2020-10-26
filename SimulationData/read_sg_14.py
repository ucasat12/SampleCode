import numpy
from read_gtab import *
from read_gids import *
from load_data import *

def getidslen_subgrp(dirname,dirno,dirname1,nsg,n2div,SMbool,preci,ptype):
    
    subgrpids_piv = loadgrpsubgrp_idspiv(1,dirname1,n2div)
    subgrpConc = subtabsconc(dirname+'/groups_'+dirno+'/subhalo_tab_'+dirno,nsg)
    cc = 0
    grnrconc = subgrpConc[5]
    lensubstrConc = subgrpConc[2]
    suboffsetConc = subgrpConc[4]
    Nsubstr = len(grnrconc)
    pgrno = -1
    subgrpno = -1
    if (ptype == 1):
       mask_subgrp = loadgrpsubgrp_data(1,4,dirname1,n2div,SMbool,preci)
       strmask = 'subgrp_DMpiv'
    if (ptype == 4):
       mask_subgrp = loadgrpsubgrp_data(1,5,dirname1,n2div,SMbool,preci)
       strmask = 'subgrp_starpiv'
    tstma = numpy.ma.array(subgrpids_piv, mask = mask_subgrp)
    grpsubgrplen = []
    offseti = 0
    f = open(dirname1 + '/' + strmask,'w' )
    for i in range(0,Nsubstr): 
          li = lensubstrConc[i]
          grno = grnrconc[i]
          subgrpno = subgrpno + 1
          if (grno > pgrno):
             subgrpno = 0 
          pgrno = grno
          tst4ma = tstma[offseti : offseti + li]
          pivma = tst4ma[~tst4ma.mask]
          pivdata = pivma.data
          grpsubgrplen.append([grno, subgrpno, len(pivdata)])
          pivdata.tofile(f)
          offseti = offseti + li
          print i, Nsubstr
    f.close()
    del f
    grpsubgrplen = numpy.array(grpsubgrplen)
    # data type of grpsubgrplen : (numpy.int,3)
    if (ptype == 1):
       grpsubgrplen.tofile(dirname1+'/gsl_DMsubgrps')
    if (ptype == 4):
       grpsubgrplen.tofile(dirname1+'/gsl_starsubgrps')
    return grpsubgrplen

def getidslen_subgrp(dirname_cc,dirname_nc,dirno,dirname1,ng,nsg,n2div,SMbool,preci,ptype):
    
    subgrpids_piv = loadgrpsubgrp_idspiv(1,dirname1,n2div)
    subgrpConc = subtabsconc_lowz(dirname_cc,dirname_nc,dirno,ng,nsg)
    cc = 0
    grnrconc = subgrpConc[5]
    lensubstrConc = subgrpConc[2]
    suboffsetConc = subgrpConc[4]
    Nsubstr = len(grnrconc)
    pgrno = -1
    subgrpno = -1
    if (ptype == 1):
       mask_subgrp = loadgrpsubgrp_data(1,4,dirname1,n2div,SMbool,preci)
       strmask = 'subgrp_DMpiv'
    if (ptype == 4):
       mask_subgrp = loadgrpsubgrp_data(1,5,dirname1,n2div,SMbool,preci)
       strmask = 'subgrp_starpiv'
    tstma = numpy.ma.array(subgrpids_piv, mask = mask_subgrp)
    grpsubgrplen = []
    offseti = 0
    f = open(dirname1 + '/' + strmask,'w' )
    for i in range(0,Nsubstr): 
          li = lensubstrConc[i]
          grno = grnrconc[i]
          subgrpno = subgrpno + 1
          if (grno > pgrno):
             subgrpno = 0 
          pgrno = grno
          tst4ma = tstma[offseti : offseti + li]
          pivma = tst4ma[~tst4ma.mask]
          pivdata = pivma.data
          grpsubgrplen.append([grno, subgrpno, len(pivdata)])
          pivdata.tofile(f)
          offseti = offseti + li
          print i, Nsubstr
    f.close()
    del f
    grpsubgrplen = numpy.array(grpsubgrplen)
    if (ptype == 1):
       grpsubgrplen.tofile(dirname1+'/gsl_DMsubgrps')
    if (ptype == 4):
       grpsubgrplen.tofile(dirname1+'/gsl_starsubgrps')
    return grpsubgrplen
