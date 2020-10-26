import numpy as np
import matplotlib

def grptab(fnam1):

    headstruct = np.dtype([
              ('Ngroups',np.int32),
              ('TotNgroups',np.int32),           
              ('Nids',np.int32),
              ('TotNids',np.int64),
              ('NTask',np.int32)])
    f = open(fnam1)
    h = np.fromfile(f,dtype=headstruct,count=1)[0]
    grplen = np.fromfile(f,dtype=np.int32,count=h['Ngroups'])
    grpoffset = np.fromfile(f,dtype=np.uint32,count=h['Ngroups'])
    mass = np.fromfile(f,dtype=np.float32,count=h['Ngroups'])
    pos = np.fromfile(f,dtype=(np.float32,3),count=h['Ngroups'])
    f.seek(h['Ngroups']*3*4,1)
    grplentype = np.fromfile(f,dtype=(np.int32,6),count=h['Ngroups'])
    grpmasstype = np.fromfile(f,dtype=(np.float32,6),count=h['Ngroups'])
    return mass,pos,h,grplen,grpoffset,grplentype,grpmasstype

def tabsconc(dir,n) :
     
    msList = [];
    posList = [];
    grplenList = [];
    grpoffsetList = [];
    grplentypeList = [];
    grpmasstypeList = []
    for i in range(0,n):
          result = grptab(str(dir)+'.%d'%i)
          msList.append(result[0])
          posList.append(result[1])
          grplenList.append(result[3])
          grpoffsetList.append(result[4])
          grplentypeList.append(result[5])
          grpmasstypeList.append(result[6])
          print i
    massconc = np.concatenate(msList)
    posconc = np.concatenate(posList)
    grplenconc = np.concatenate(grplenList)
    grpoffsetconc = np.concatenate(grpoffsetList)
    grplentypeconc = np.concatenate(grplentypeList)
    grpmasstypeconc = np.concatenate(grpmasstypeList)
    return massconc,posconc,grplenconc,grpoffsetconc,grplentypeconc, grpmasstypeconc

def tabsconc_lowz(dir_nc,ng):

    massconc,posconc,grplenconc,grpoffsetconc,grplentypeconc, grpmasstypeconc = tabsconc(dir_nc,ng)
    grpoffset = grpoffsetconc
    tst_grpoff = (grpoffset[:-1] > grpoffset[1:]).astype(np.int32)
    grp_corr = np.zeros(len(tst_grpoff) + 1,dtype = np.int32)
    cc1 = 0
    grpoff1 = np.zeros(len(grpoffset),dtype=np.uint64)
    for i in range(0,len(tst_grpoff)):
          cc1 = cc1 + tst_grpoff[i]
          grp_corr[i+1] = cc1 
    grpoff1 = np.uint64(grpoffset + 4294967296*grp_corr) 
    return massconc,posconc,grplenconc,grpoff1,grplentypeconc, grpmasstypeconc

def grpmassfunc_obs(dir,n,ini,rang) :

    tabsMassPos = tabsconc(dir,n)
    tabsMassconc = tabsMassPos[0]/(7.02033421204296e-11)
    tabsPosconc = tabsMassPos[1]
    massbin = np.zeros(pow(10,rang),np.int32)
    densfunc = np.zeros(pow(10,rang),np.double)
    massp = np.zeros(pow(10,rang),np.double)
    Bsize = 50000
    for i in range(len(tabsMassconc)):
          testm = tabsMassconc[i]
          cno = np.int(testm/(pow(10,ini)))
          if (((cno-1) >= 0) & ((cno-1) < pow(10,rang))):
             massbin[cno-1] = massbin[cno-1] + 1
    densfunc = np.double(massbin)/(pow(Bsize,3))
    for i in range(pow(10,rang)):
          massp[i] = pow(10,ini) + pow(10,ini)*i
    return massbin,densfunc,tabsMassconc,tabsPosconc,massp
    
def grpmassfunc(dirname,dirno,boxlen,n):

    final = tabsconc('/home/groups_' + dirno + '/group_tab_' + dirno,ng)
    sm = final[0]/(1e-10)
    h,bins = histogram(sm,bins=logspace(8,16,n))
    ncount = np.zeros(n - 1)
    t1 = len(sm)
    t2 = sm.min()
    t3 = sm.max()
    V = (boxlen)**3
    Mcenter = zeros(n - 1)
    hdndm = zeros(n - 1)
    for i in range(n - 1):
        Mcenter[i] = 0.5*(bins[i] + bins[i+1])
        hdndm[i] = h[i]/(bins[i+1]-bins[i])
        print i
    hdndm = hdndm/V
    matplotlib.pyplot.plot(Mcenter,hdndm, label = ' ' + strRS)
    matplotlib.pyplot.legend()
    matplotlib.pyplot.xscale('log')
    matplotlib.pyplot.yscale('log')
    matplotlib.pyplot.xlabel('M (in Solar Mass)')
    matplotlib.pyplot.ylabel('dn/dM')
    matplotlib.pyplot.title('Halo Mass function from Group Tab')
    matplotlib.pyplot.show()
    return Mcenter,hdndm,bins,ncount,sm,V

def grpstellarmassfunc(dirname,dirno,boxlen,n):

    final = tabsconc('/home/groups_' + dirno + '/group_tab_' + dirno,ng)
    sm = np.array(final[5])/(1e-10)
    h,bins = histogram(sm[:,4],bins=logspace(7,13,n))
    ncount = np.zeros(n - 1)
    t1 = len(sm)
    t2 = sm.min()
    t3 = sm.max()

    V = (boxlen)**3
    Mcenter = zeros(n - 1)
    hdndm = zeros(n - 1)
    for i in range(n - 1):
        Mcenter[i] = 0.5*(bins[i] + bins[i+1])
        hdndm[i] = h[i]/(bins[i+1]-bins[i])
        print i
    hdndm = hdndm/V
    matplotlib.pyplot.plot(Mcenter,hdndm, label = ' ' + strRS)
    matplotlib.pyplot.legend()
    matplotlib.pyplot.xscale('log')
    matplotlib.pyplot.yscale('log')
    matplotlib.pyplot.xlabel('M (in Solar Mass)')
    matplotlib.pyplot.ylabel('dn/dM')
    matplotlib.pyplot.title('Halo Mass function from Group Tab')
    matplotlib.pyplot.show()
    return Mcenter,hdndm,bins,ncount,sm,V

def finalplot():
    
    chk1 = grptab('/home/groups_' + dirno + '/group_tab_' + dirno + '.31')
    grphead32 = chk1[2]
    Ntot = grphead32['TotNids']
    grphist = tenplot()
    grpH = grphist[3]
    bins = grphist[2]
    center = grphist[1]
    nH = np.zeros(len(grpH))
    for i in range(len(center)):
         nH[i] = grpH[i]*(bins[i+1] - bins[i])*Ntot    
    mf = np.zeros(len(grpH))
    for i in range(len(center)):
         mf[i] = (nH[i]/center[i])    
    Vol = (50) ** 3
    dlogM = diff(log10(bins))[0]
    dndM = mf/(Vol*dlogM)
    return center, dndM

