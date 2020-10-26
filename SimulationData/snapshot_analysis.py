import numpy, sharedmem
from load_data import *
from snapshot_readsorted import *
from read_gtab import *
from read_gids import *
from read_sgtab import *
from read_sgids import *
from read_g_14 import *
from read_sg_14 import *

preci = 1
dirname = '/home/snapdir_' + dirno + '/groups_and_subgroups'
dirno = '100'
ns = 256
ng = 256
nsg = 256
dirname1 = '/home/snapdir_' + dirno + '/data'
n2div = 5000
SMbool = 1

snappid_notsorted = loadsnapshot(0,dirname1,n2div,SMbool,preci)
snappid_argsorted = numpy.argsort(snappid_notsorted)

n1 = len(snappid_argsorted)
chunksize = numpy.int(n1/n2div)
slices = [slice(i, i+chunksize) for i in range(0, n1-chunksize,chunksize)]
slices.append(slice(chunksize*n2div,n1))
f = open(dirname1 + '/snappid_argsortednew','w')
for i in range(0,len(slices)):
      (snappid_argsorted[slices[i]]).tofile(f)
f.close()

del snappid_notsorted, snappid_argsorted

tst1new = sortsnapshot(0,n2div,dirname1,preci)
tst1new = sortsnapshot(1,n2div,dirname1,preci)
tst1new = sortsnapshot(2,n2div,dirname1,preci)
tst1new = sortsnapshot(3,n2div,dirname1,preci)

grpID, grpIDconc, grpconc, grpIDlist = grpids_allgroups(dirname,dirno,ng)
f = open(dirname1 + '/grpids_idscorrect','w')
for i in range(0,len(grpID)):
      (grpID[i]).tofile(f)
f.close()

del grpID, grpIDconc, grpconc, grpIDlist

substrID, subIDconc, subgrpConc, subIDlist = subgrpids_allgroups(dirname,dirno,nsg)
f = open(dirname1 + '/subgrpids_idscorrect','w')
for i in range(0,len(substrID)):
      (substrID[i]).tofile(f)
f.close()
del substrID, subIDconc, subgrpConc, subIDlist

grpids_ids = loadgrpsubgrp_data(0,0,dirname1,n2div/2,SMbool,preci)
subgrpids_ids = loadgrpsubgrp_data(1,0,dirname1,n2div,SMbool,preci)
snappidsort = loadsnapshot_sorted(0,dirname1,n2div,SMbool,preci)
grpids_piv = numpy.searchsorted(snappidsort,grpids_ids,side='left')
subgrpids_piv = numpy.searchsorted(snappidsort,subgrpids_ids,side='left')
grpids_piv.tofile(dirname1+'/grpids_piv')
subgrpids_piv.tofile(dirname1+'/subgrpids_piv')

del grpids_ids, subgrpids_ids, snappidsort, grpids_piv, subgrpids_piv

grpids_piv = loadgrpsubgrp_idspiv(0,dirname1,n2div/2)
subgrpids_piv = loadgrpsubgrp_idspiv(1,dirname1,n2div)
snapptypesort = loadsnapshot_sorted(2,dirname1,n2div,SMbool,preci)

n2 = len(grpids_piv)
chunksize2 = numpy.int(n2*2/n2div)
slices2 = [slice(i, i+chunksize2) for i in range(0, n2-chunksize2,chunksize2)]
slices2.append(slice(chunksize2*n2div/2,n2))
f = open(dirname1 + '/grpids_ptypecorrect','w')
for i in range(0,n2div/2 + 1):
      grpptypes = snapptypesort[grpids_piv[slices2[i]]]
      (grpptypes).tofile(f)
f.close()

n2 = len(subgrpids_piv)
chunksize2 = numpy.int(n2/n2div)
slices2 = [slice(i, i+chunksize2) for i in range(0, n2-chunksize2,chunksize2)]
slices2.append(slice(chunksize2*n2div,n2))
f = open(dirname1 + '/subgrpids_ptypecorrect','w')
for i in range(0,n2div + 1):
      subgrpptypes = snapptypesort[subgrpids_piv[slices2[i]]]
      (subgrpptypes).tofile(f)
f.close()

grpids_ptype = loadgrpsubgrp_data(0,2,dirname1,n2div/2,SMbool,preci)
try1 = numpy.ma.masked_not_equal(grpids_ptype,1)
try4 = numpy.ma.masked_not_equal(grpids_ptype,4)
chk1 = try1.mask.astype(int32)
chk4 = try4.mask.astype(int32)
chk1.tofile(dirname1 + '/grpids_DMmask')
chk4.tofile(dirname1 + '/grpids_starmask')
del grpids_ptype, try1, try4, chk1, chk4 

subgrpids_ptype = loadgrpsubgrp_data(1,2,dirname1,n2div,SMbool,preci)
try1 = numpy.ma.masked_not_equal(subgrpids_ptype,1)
try4 = numpy.ma.masked_not_equal(subgrpids_ptype,4)
chk1 = try1.mask.astype(int32)
chk4 = try4.mask.astype(int32)
chk1.tofile(dirname1 + '/subgrpids_DMmask')
chk4.tofile(dirname1 + '/subgrpids_starmask')
del subgrpids_ptype, try1, try4, chk1, chk4

grplen1 = getidslen_grp(dirname,dirno,dirname1,ng,n2div/2,SMbool,preci,1)
grplen4 = getidslen_grp(dirname,dirno,dirname1,ng,n2div/2,SMbool,preci,4)
subgrplen1 = getidslen_subgrp(dirname,dirno,dirname1,ng,n2div,SMbool,preci,1)
subgrplen4 = getidslen_subgrp(dirname,dirno,dirname1,ng,n2div,SMbool,preci,4)
