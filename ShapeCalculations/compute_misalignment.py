import numpy, matplotlib
from read_sgtab import *

nmin = 10
dirno = '100'
mmin = 10.0
mmax = 15.0
subgrpbool = 1
ndim = 3
gdim = str(ndim) + 'd'
boxlen = 1000.0
permax = boxlen/2.0

if (subgrpbool == 1):
   
   rdstr1 = ''
   center_DM = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/centerdata_' + str(ndim) + 'dDM',dtype = (numpy.float,3))
   grplen_DM = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/grplen_' + str(ndim) + 'dDM',dtype = (numpy.int,3))
   pno_DM = grplen_DM[:,2]
   eigvec_DM = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/eigvecdata_' + str(ndim) + 'dDM',dtype = (numpy.float,(ndim,ndim)))
   mass_DM = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/massdata_' + str(ndim) + 'dDM',dtype = numpy.float)
   mass_DM = mass_DM*pow(10,10)
   q_DM = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/qval_' + str(ndim) + 'dDM',dtype = numpy.float)
   s_DM = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/sval_' + str(ndim) + 'dDM',dtype = numpy.float)
   idCent_DM = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/idcc_' + str(ndim) + 'dDM',dtype = (numpy.int))

   center_star = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/centerdata_' + str(ndim) + 'dstar',dtype = (numpy.float,3))
   grplen_star = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/grplen_' + str(ndim) + 'dstar',dtype = (numpy.int,3))
   pno_star = grplen_star[:,2]
   eigvec_star = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/eigvecdata_' + str(ndim) + 'dstar',dtype = (numpy.float,(ndim,ndim)))
   mass_star = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/massdata_' + str(ndim) + 'dstar',dtype = numpy.float)
   mass_star = mass_star*pow(10,10)
   q_star = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/qval_' + str(ndim) + 'dstar',dtype = numpy.float)
   s_star = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/sval_' + str(ndim) + 'dstar',dtype = numpy.float)
   idCent_star = numpy.fromfile('/home/subgroupshapes_' + rdstr1 + '/idcc_' + str(ndim) + 'dstar',dtype = (numpy.int))

   piv = numpy.searchsorted(idCent_DM,idCent_star,side='left')
   chk1 = ((piv < len(idCent_DM)))
   piv3d = piv[chk1]

   center_DM2 = center_DM[piv3d]
   grplen_DM2 = grplen_DM[piv3d]
   pno_DM2 = pno_DM[piv3d]
   eigvec_DM2 = eigvec_DM[piv3d]
   mass_DM2 = mass_DM[piv3d]
   idCent_DM2 = idCent_DM[piv3d]
   q_DM2 = q_DM[piv3d]
   s_DM2 = s_DM[piv3d]

   center_star2 = center_star[chk1]
   grplen_star2 = grplen_star[chk1]
   pno_star2 = pno_star[chk1]
   eigvec_star2 = eigvec_star[chk1]
   mass_star2 = mass_star[chk1]
   idCent_star2 = idCent_star[chk1]
   q_star2 = q_star[piv3d]
   s_star2 = s_star[piv3d]

   chk2 = ((grplen_DM2[:,2] > nmin) & (grplen_star2[:,2] > nmin) & (idCent_DM2 == idCent_star2) & (mass_DM2 >= pow(10,mmin)) & (mass_DM2 < pow(10,mmax)))

   center_DM3 = center_DM2[chk2]
   grplen_DM3 = grplen_DM2[chk2]
   pno_DM3 = pno_DM2[chk2]
   eigvec_DM3 = eigvec_DM2[chk2]
   mass_DM3 = mass_DM2[chk2]
   idCent_DM3 = idCent_DM2[chk2]
   q_DM3 = q_DM2[chk2]
   s_DM3 = s_DM2[chk2]

   center_star3 = center_star2[chk2]
   grplen_star3 = grplen_star2[chk2]
   pno_star3 = pno_star2[chk2]
   eigvec_star3 = eigvec_star2[chk2]
   mass_star3 = mass_star2[chk2]
   idCent_star3 = idCent_star2[chk2]
   q_star3 = q_star2[chk2]
   s_star3 = s_star2[chk2]

DMeigvec = eigvec_DM3[:,0]
stareigvec = eigvec_star3[:,0]
misalign_angle = numpy.array(numpy.arccos(numpy.abs(numpy.sum(DMeigvec[:]*stareigvec[:],axis = 1)))*(180/numpy.pi))
hmisalign,binmisalign = numpy.histogram(misalign_angle,bins=nbins,range=(0.0,90.0),density=True)   

fig1 = matplotlib.pyplot.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(0.5*(binmisalign[1:]+binmisalign[:-1]),double(hmisalign),lw=2.0,color = 'r',label= 'Angles')
ax1.set_xlabel('orientation (in degrees)')
ax1.set_ylabel('Normalized histogram')
ax1.legend()
fig1.show()

