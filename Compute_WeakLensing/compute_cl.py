import numpy as np
import scipy


def cl_func(chiz,func_ni,func_gj,func_ij_chi_pk,func_az):

    fval = func_ni(chiz)*func_gj(chiz)*func_ij_chi_pk(chiz)/(func_az(chiz)**2.0)
    return fval

nbins_wl = 10
inp_chiz = np.array([chiz_func(0.06),chiz_func(1.0),chiz_func(2.0)])
norm_pm0 = 9.0*((100.0*0.702)**4)*(0.205**2)/(4.0*((3.0e5)**4))

for bin_i in range(0,nbins_wl):
          for bin_j in range(bin_i,nbins_wl):
                    zmin_pzi = zvals_pzbins[bin_i]
                    zmax_pzi = zvals_pzbins[bin_i+1]
                    zmin_pzj = zvals_pzbins[bin_j]
                    zmax_pzj = zvals_pzbins[bin_j+1]     
                    zvals_pz_ij = zvals_pz
                    fpz_specz_ij = fpz_pz 
                    fpz_pzi = np.zeros(len(fpz_specz_ij),dtype=np.float32) 
                    fpz_pzj = np.zeros(len(fpz_specz_ij),dtype=np.float32) 
                    for b_pz_ij in range(0,len(fpz_speczi)):
                                zij_pz = zvals_pzi[b_pz_ij]
                                nij_pz = fpz_speczi[b_pz_ij]
                                fpz_pzi[b_pz_ij] = nij_pz*scipy.integrate.quad(func_pzi,zmin_pzi,zmax_pzi,args=(zij_pz,Az_sigz))[0]/scipy.integrate.quad(func_pzi,zmin_pz,zmax_pz,args=(zij_pz,Az_sigz))[0]
                                fpz_pzj[b_pz_ij] = nij_pz*scipy.integrate.quad(func_pzi,zmin_pzj,zmax_pzj,args=(zij_pz,Az_sigz))[0]/scipy.integrate.quad(func_pzi,zmin_pz,zmax_pz,args=(zij_pz,Az_sigz))[0]
                    nfunc_bin_ij = gfunc_bin(bin_i,zvals_pz,chiz_func,fpz_pzi,zvals_pzbins)
                    gfunc_bin_ij = gfunc_bin(bin_j,zvals_pz,chiz_func,fpz_pzj,zvals_pzbins)
                    cls_ij = np.zeros(len(lvals),dtype=np.float32)
                    for  b_ij in range(0,len(lvals)):
                              ij_l = np.log10(lvals[b_ij])
                              log_l_pk_vals = np.array([func_i0l_logPk_intpz2(ij_l),func_i0l_logPk_intpz1(ij_l),func_i0l_logPk_intpz0(ij_l)]) 
                              func_chi_pk = interp1d(inp_chiz,10.0**(log_l_pk_vals))
                              cls_ij[b_ij] = scipy.integrate.quad(cl_func,chiz_func(0.06),chiz_func(2.0),args=(nfunc_bin_ij,gfunc_bin_ij,func_chi_pk,func_az))[0]
