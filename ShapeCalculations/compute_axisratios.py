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

misalign_angle = []

q_valDM = numpy.array(q_DM3) 
q_valstar = numpy.array(q_star3)
h1,b1 = numpy.histogram(q_valDM,bins=nbins,range=(0.0,1.0),density=True)
h3,b3 = numpy.histogram(q_valstar,bins=nbins,range=(0.0,1.0),density=True)
if (ndim == 3):
   s_valDM = numpy.array(s_DM3)
   s_valstar = numpy.array(s_star3)
   h2,b2 = numpy.histogram(s_valDM,bins=nbins,range=(0.0,1.0),density=True)
   h4,b4 = numpy.histogram(s_valstar,bins=nbins,range=(0.0,1.0),density=True)
hmisalign,binmisalign = numpy.histogram(misalign_angle,bins=nbins,range=(0.0,90.0),density=True)   

ma_starpos = numpy.array(ma_starpos)
ma_DMpos = numpy.array(ma_DMpos)

fig1 = matplotlib.pyplot.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(0.5*(b1[1:]+b1[:-1]),double(h1),lw=2.0,color = 'b',label= 'q (1)')
ax1.legend(loc = 0)
ax1.set_xlabel('q',fontsize=20)
ax1.set_ylabel('Normalized Probability Density')
fig1.show()

fig2 = matplotlib.pyplot.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(0.5*(b2[1:]+b2[:-1]),double(h2),lw=2.0,color = 'g',label= 's (1)')
ax2.legend(loc = 0)
ax2.set_xlabel('s',fontsize=20)
ax2.set_ylabel('Normalized Probability Density')
fig2.show()

fig3 = matplotlib.pyplot.figure()
ax3 = fig3.add_subplot(111)
ax3.plot(0.5*(b3[1:]+b3[:-1]),double(h3),lw=2.0,color = 'r',label= 'q (4)')
ax3.legend(loc = 0)
ax3.set_xlabel('q',fontsize=20)
ax3.set_ylabel('Normalized Probability Density')
fig3.show()

fig4 = matplotlib.pyplot.figure()
ax4 = fig4.add_subplot(111)
ax4.plot(0.5*(b4[1:]+b4[:-1]),double(h4),lw=2.0,color = 'c',label= 's (4)')
ax4.legend(loc = 0)
ax4.set_xlabel('q',fontsize=20)
ax4.set_ylabel('Normalized Probability Density')
fig4.show()

