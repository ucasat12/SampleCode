import numpy as np
from scipy.interpolate import interp1d

dirname = '/physics2/vat/cosmosis_install/cosmosis/ia_hmodel/ia_halo_model/ia_mock_inifiles/mock_barpk_euclid/save_path/shear_cl'
ntbin = 10

sige = 0.26
ngal = 20.0
ngal_ti = ngal*1.0/ntbin
omgs = 15000.0

def readdata(dirname,ntbin,omgs,sige,ngal_ti):

    ell = np.genfromtxt(dirname + '/ell.txt')
    #ell_diff = np.diff(ell) 
    lmin = 100
    lmax = 5000
    n_lbins = 12
    ell_bins = np.logspace(np.log10(lmin),np.log10(lmax),n_lbins+1)
    ell_diff = np.diff(ell_bins)
    ell_cen = 10.0**(0.5*(np.log10(ell_bins[1:]) + np.log10(ell_bins[:-1])))
    cov_mat = np.zeros((n_lbins*ntbin*(ntbin+1)/2,n_lbins*ntbin*(ntbin+1)/2),dtype=np.float32)
    data_vec_ij = np.zeros(n_lbins*ntbin*(ntbin+1)/2,dtype=np.float32)
    data_vec_kl = np.zeros(n_lbins*ntbin*(ntbin+1)/2,dtype=np.float32)
    ij_ini = 0
    ij_fin = n_lbins
    kl_ini = 0
    kl_fin = n_lbins
    for i in range(0,ntbin):
        for j in range(0,i+1):
            kl_ini = 0
            kl_fin = n_lbins
            data_cl_ij = np.zeros(n_lbins,dtype=np.float32)
            Cij = np.genfromtxt(dirname + '/bin_' + str(i+1) + '_' + str(j+1) + '.txt')
            cfunc_ij = interp1d(ell,Cij)
            for k in range(0,ntbin):
                for l in range(0,k+1):
                    cov_ijkl = np.zeros((len(ell_cen),len(ell_cen)))
                    data_cl_kl = np.zeros(n_lbins,dtype=np.float32)
                    Ckl = np.genfromtxt(dirname + '/bin_' + str(k+1) + '_' + str(l+1) + '.txt')
                    cfunc_kl = interp1d(ell,Ckl)                    
                    deltaik = 0
                    deltajl = 0
                    deltail = 0
                    deltajk = 0
                    if (i == k): 
                       deltaik = 1
                    if (j == l):
                       deltajl = 1
                    if (i == l):
                       deltail = 1
                    if (j == k):
                       deltajk = 1
                    if (i >= k):
                       Cik = np.genfromtxt(dirname + '/bin_' + str(i+1) + '_' + str(k+1) + '.txt')
                    if (i < k):
                       Cik = np.genfromtxt(dirname + '/bin_' + str(k+1) + '_' + str(i+1) + '.txt')
                    if (j >= l):
                       Cjl = np.genfromtxt(dirname + '/bin_' + str(j+1) + '_' + str(l+1) + '.txt')
                    if (j < l):
                       Cjl = np.genfromtxt(dirname + '/bin_' + str(l+1) + '_' + str(j+1) + '.txt')
                    if (i >= l):
                       Cil = np.genfromtxt(dirname + '/bin_' + str(i+1) + '_' + str(l+1) + '.txt')
                    if (i < l):
                       Cil = np.genfromtxt(dirname + '/bin_' + str(l+1) + '_' + str(i+1) + '.txt')
                    if (j >= k):
                       Cjk = np.genfromtxt(dirname + '/bin_' + str(j+1) + '_' + str(k+1) + '.txt')
                    if (j < k):
                       Cjk = np.genfromtxt(dirname + '/bin_' + str(k+1) + '_' + str(j+1) + '.txt')
                    cfunc_ik = interp1d(ell,Cik)
                    cfunc_jl = interp1d(ell,Cjl)
                    cfunc_il = interp1d(ell,Cil)
                    cfunc_jk = interp1d(ell,Cjk)
                    for l12 in range(0,len(ell_cen)):
                            cov_ijkl[l12,l12] = (2.0*np.pi/(omgs*ell_cen[l12]*ell_diff[l12]))*((cfunc_ik(ell_cen[l12]) + deltaik*sige*sige/(2.0*ngal_ti))*(cfunc_jl(ell_cen[l12]) + deltajl*sige*sige/(2.0*ngal_ti)) + (cfunc_il(ell_cen[l12]) + deltail*sige*sige/(2.0*ngal_ti))*(cfunc_jk(ell_cen[l12]) + deltajk*sige*sige/(2.0*ngal_ti)))
                            data_cl_ij[l12] = cfunc_ij(ell_cen[l12])
                            data_cl_kl[l12] = cfunc_kl(ell_cen[l12])
                    cov_mat[ij_ini:ij_fin,kl_ini:kl_fin] = cov_ijkl[:,:]
                    print(data_vec_kl.shape, data_cl_kl.shape)
                    data_vec_kl[kl_ini:kl_fin] = data_cl_kl
                    kl_ini = kl_ini + n_lbins
                    kl_fin = kl_fin + n_lbins 
                    print(i,j,k,l)
            data_vec_ij[ij_ini:ij_fin] = data_cl_ij
            ij_ini = ij_ini + n_lbins
            ij_fin = ij_fin + n_lbins
    return cov_mat, data_vec_ij, data_vec_kl
fcov, fij, fkl = readdata(dirname,ntbin,omgs,sige,ngal_ti)
fcov_inv = np.linalg.inv(fcov)

np.savetxt('./data_vec_cls_barpk_euclid.txt',fij)

