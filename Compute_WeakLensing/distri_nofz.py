import numpy as np, scipy
from scipy.interpolate import interp1d

# n(z) parameters

def get_pars_nz():

    alpha = 1.27
    beta = 1.02
    z0 = 0.5
    zmin = 0.0
    zmax = 3.5
    ngal_dens = 26.0

    return alpha, beta, z0, zmin, zmax, ngal_dens

# photo-z parameters

def get_pars_photoz():

    alpha, beta, z0, zmin, zmax, ngal_dens = get_pars_nz()
    zmin_pz = zmin
    zmax_pz = zmax
    nbins_pz = 10
    Az_sigz = 0.05

    return zmin_pz, zmax_pz, nbins_pz, Az_sigz

# survey parameters
    
def get_pars_survey():

    survey_area = 18000.0
    shape_noise = 0.05

    return survey_area, shape_noise

def nz_speczi():

    nbins = 1000 ## Number of red-shift data points 
    alpha, beta, z0, zmin, zmax, ngal_dens = get_pars_nz()
    zrange = np.linspace(zmin,zmax,nbins+1)
    dz = np.diff(zrange)
    zmid = 0.5*(zrange[1:]+zrange[:-1])
    pzfunc = pow(zrange,alpha)*np.exp(-1*(pow(zrange/z0,beta)))
    psum = ((0.5*(pzfunc[1:]+pzfunc[:-1]))*dz).sum()
    pzfunc = pzfunc/psum

    return zrange, pzfunc

def func_pzi(zph,zsp,az):

    sigz_pzi = az*(1+zsp)
    fval_pzi = (1.0/(np.sqrt(2.0*np.pi)*sigz_pzi))*(np.exp(-1.0*(zph-zsp)**2.0/(2.0*(sigz_pzi**2.0))))

    return fval_pzi

def pz_tomographic_bins():

    zvals_pz, fpz_pz = nz_speczi()
    zmin_pz, zmax_pz, nbins_pz, Az_sigz = get_pars_photoz()
    func_zpz = scipy.interpolate.interp1d(zvals_pz, fpz_pz)
    fpz_cdf = np.zeros(len(zvals_pz))
    for i in range(0,len(zvals_pz)):
          fpz_cdf[i] = scipy.integrate.quad(func_zpz,0,zvals_pz[i])[0]
    zvals_pzbins = zvals_pz[np.searchsorted(fpz_cdf,np.linspace(0,1.0,nbins_pz+1))] ## Splitting bins to ensure equal number density of galaxies in each bin
    fpz_pzij = np.zeros((nbins_pz,len(zvals_pz)),dtype=np.float32)
    for i in range(0,nbins_pz):
          zmin_pzi = zvals_pzbins[i]
          zmax_pzi = zvals_pzbins[i+1]
          zvals_pzi = zvals_pz
          fpz_speczi = fpz_pz
          for j in range(0,len(fpz_speczi)):
                zij_pz = zvals_pzi[j]
                nij_pz = fpz_speczi[j]
                fpz_pzij[i,j] = nij_pz*scipy.integrate.quad(func_pzi,zmin_pzi,zmax_pzi,args=(zij_pz,Az_sigz))[0]/scipy.integrate.quad(func_pzi,zmin_pz,zmax_pz,args=(zij_pz,Az_sigz))[0]
                print i, j
    return fpz_pzij

def plot_nz():

    alpha, beta, z0, zmin, zmax, ngal_dens = get_pars_nz()
    zmin_pz, zmax_pz, nbins_pz, Az_sigz = get_pars_photoz()
    zvals_pz, fpz_pz = nz_speczi()
    fpz_pzij = pz_tomographic_bins()
    matplotlib.pyplot.figure(1)
    matplotlib.pyplot.plot(zvals_pz,ngal_dens*fpz_pz) 
    for i in range(0,nbins_pz):
          matplotlib.pyplot.plot(zvals_pz,ngal_dens*fpz_pzij[i,:],ls='--')
    matplotlib.pyplot.xlabel(r'$z$',fontsize=20)
    matplotlib.pyplot.ylabel(r'$\frac{dN}{d \Omega dz} [\rm{arcmin} ^{-2}]$',fontsize=20)
    matplotlib.pyplot.xlim([0,3.5])
    matplotlib.pyplot.show()

def get_cosmo_params():
    h0 = 0.6774
    H0 = 100*h0
    om0 = 0.3089
    ol0 = 0.6911
    c0 = 3.0*1e5  
    return h0, om0, ol0, c0

def Hzinv(z):
    h0, om0, ol0, c0 = get_cosmo_params()
    H0 = h0*100.0
    az = 1.0/(1+z)
    Hz = H0*(pow(om0*((1.0/az)**3.0) + ol0,0.5))
    Hzinv = 1.0*c0/Hz
    return Hzinv
   
def func_chi_z(z):
    chi_z = scipy.integrate.quad(Hzinv,0,z)[0]
    return chi_z

def get_chi_z():

    zvals_pz, fpz_pz = nz_speczi()
    chiz = np.zeros(len(zvals_pz))
    for i in range(0,len(zvals_pz)):
          chiz[i] = func_chi_z(zvals_pz[i])
    chi_func_z = interp1d(zvals_pz,chiz)
    chi_inf_z = scipy.integrate.quad(Hzinv,0,np.inf)[0]

    return zvals_pz, chi_func_z, chi_inf_z

def get_a_chi():
    
    zvals_pz, chi_func_z, chi_inf_z = get_chi_z()
    a_func_chi = scipy.interpolate.interp1d(chi_func_z(zvals_pz),1.0/(1.0 + zvals_pz))

    return zvals_pz, a_func_chi, chi_func_z, chi_inf_z

def gxi(chid,fpz_pz_chid,chimin):

    fval = fpz_pz_chid*(chid-chimin)/chid

    return fval

def gfunc_bin():

    h0, om0, ol0, c0 = get_cosmo_params()
    zvals_pz, a_func_chi, chi_func_z, chi_inf_z = get_a_chi()

    pre_fac = (3.0/2.0)*om0*((h0*100.0/c0)**2.0)
    
    fpz_pzij = pz_tomographic_bins()

    for i in range(len(fpz_pzij[:,0])):
          n_i = fpz_pzij[i,:]
          n_sum = 0.5*((n_i[1:]+n_i[:-1])*np.diff(zvals_pz[:])).sum()
          n_i = n_i/n_sum
          fpz_pzij[i,:] = n_i

    qis_chiz = np.zeros((len(fpz_pzij[:,0]),len(fpz_pzij[0,:])),dtype=np.float32)

    #nfunc_bin_i = nfunc_bin(bini,zvals,chiz_func,fpz,zvals_wlbins)
    #intp_zbins = zvals_pz #np.linspace(zvals[0],zvals[-1],101)
    intp_chiz = chi_func_z(zvals_pz)
    intp_a = 1.0/(1.0 + zvals_pz)
    for bini in range(0,len(fpz_pzij[:,0])):
             for i in range(0,len(zvals_pz)):
                   fs1 = gxi(intp_chiz[i:],fpz_pzij[bini,i:],intp_chiz[i])
                   fs2 = intp_chiz[i:]
                   qis_chiz[bini,i] = 0.5*((fs1[1:]+fs1[:-1])*np.diff(fs2)).sum()
                   qis_chiz[bini,i] = (intp_chiz[i]/intp_a[i])*qis_chiz[bini,i]
                   print bini, i
    qis_chiz = qis_chiz*pre_fac
    return zvals_pz, intp_chiz, fpz_pzij, qis_chiz

zvals_pz, intp_chiz, fpz_ij, qis_chiz = gfunc_bin()

