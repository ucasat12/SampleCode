import numpy as np
import matplotlib

def subgrptab(fnam):
    
    headstruct = np.dtype([
                           ('Ngroups',np.int32),
                           ('TotNgroups',np.int32),           
                           ('Nids',np.int32),
                           ('TotNids',np.int64),
                           ('NTask',np.int32),
                           ('Nsubgroups',np.int32),
                           ('TotNsubgroups',np.int32)])
    f = open(fnam)
    h = np.fromfile(f,dtype=headstruct,count=1)[0]
    f.seek(h['Ngroups']*4,1)
    grpidoffset = np.fromfile(f,dtype=np.int32,count=h['Ngroups'])
    f.seek(h['Ngroups']*4,1)
    grpMinPot = np.fromfile(f,dtype=(np.float32,3),count=h['Ngroups'])
    f.seek(h['Ngroups']*4,1)
    f.seek(h['Ngroups']*4,1)
    f.seek(h['Ngroups']*4,1)
    f.seek(h['Ngroups']*4,1)
    f.seek(h['Ngroups']*4,1)
    f.seek(h['Ngroups']*4,1)
    f.seek(h['Ngroups']*4,1)
    f.seek(h['Ngroups']*4,1)
    f.seek(h['Ngroups']*4,1)
    contlen = np.fromfile(f,dtype=np.int32,count=h['Ngroups'])
    contmass = np.fromfile(f,dtype=np.float32,count=h['Ngroups'])
    Nsubstruct = np.fromfile(f,dtype=np.int32,count=h['Ngroups'])
    firsubstruct = np.fromfile(f,dtype=np.int32,count=h['Ngroups'])
    lensubstruct = np.fromfile(f,dtype=np.int32,count=h['Nsubgroups'])
    suboffset = np.fromfile(f,dtype=np.uint32,count=h['Nsubgroups'])  
    parent_sub = np.fromfile(f,dtype=np.int32,count=h['Nsubgroups'])
    subgrp_mass = np.fromfile(f,dtype=np.float32,count=h['Nsubgroups'])
    subgrp_pos = np.fromfile(f,dtype=(np.float32,3),count=h['Nsubgroups'])
    subgrp_vel = np.fromfile(f,dtype=(np.float32,3),count=h['Nsubgroups'])
    subgrp_CM = np.fromfile(f,dtype=(np.float32,3),count=h['Nsubgroups'])
    subgrp_spin = np.fromfile(f,dtype=(np.float32,3),count=h['Nsubgroups'])
    f.seek(h['Nsubgroups']*4,1)
    f.seek(h['Nsubgroups']*4,1)
    f.seek(h['Nsubgroups']*4,1)
    f.seek(h['Nsubgroups']*4,1)
    idMbnd = np.fromfile(f,dtype=np.uint64,count=h['Nsubgroups'])
    grnr = np.fromfile(f,dtype=np.int32,count=h['Nsubgroups'])  
    masstab = np.fromfile(f,dtype=(np.float32,6),count = h['Nsubgroups'])  
    return h,grpidoffset,subgrp_mass,subgrp_pos,Nsubstruct,firsubstruct,contlen,contmass,suboffset,parent_sub,grnr,lensubstruct,subgrp_CM,grpMinPot,idMbnd,masstab,subgrp_vel,subgrp_spin
    
def subtabsconc(dir,n) :
   
    massList = []
    posList = []
    NsubstrList = []
    suboffsetList = []
    lensubstrList = []
    grnrList = []
    CMlist = []
    grpMinPotList = []
    idMbndList = []
    masstab = []
    velList = []
    spinList = []
    for i in range(0,n):
          result = subgrptab(str(dir)+'.%d'%i)
          massList.append(result[2])
          posList.append(result[3])
          NsubstrList.append(result[4])
          suboffsetList.append(result[8])
          lensubstrList.append(result[11])
          grnrList.append(result[10])
          CMlist.append(result[12])
          grpMinPotList.append(result[13])
          idMbndList.append(result[14])
          masstab.append(result[15])
          velList.append(result[16])
          spinList.append(result[17])
          print i
    massConc = np.concatenate(massList)
    posConc = np.concatenate(posList)
    NsubstrConc = np.concatenate(NsubstrList)
    suboffsetConc = np.concatenate(suboffsetList)
    lensubstrConc = np.concatenate(lensubstrList)
    grnrConc = np.concatenate(grnrList)
    subgrpCMconc = np.concatenate(CMlist)
    grpMinPotConc = np.concatenate(grpMinPotList)
    idMbndListConc = np.concatenate(idMbndList)
    masstabConc = np.concatenate(masstab)
    velConc = np.concatenate(velList)
    spinConc = np.concatenate(spinList)
    return massConc,posConc,lensubstrConc,NsubstrConc,suboffsetConc,grnrConc,subgrpCMconc,grpMinPotConc,idMbndListConc,masstabConc,velConc,spinConc

def grpradii(fnam):

    headstruct = np.dtype([
                           ('Ngroups',np.int32),
                           ('TotNgroups',np.int32),           
                           ('Nids',np.int32),
                           ('TotNids',np.int64),
                           ('NTask',np.int32),
                           ('Nsubgroups',np.int32),
                           ('TotNsubgroups',np.int32)])
    f = open(fnam)
    h = np.fromfile(f,dtype=headstruct,count=1)[0]
    f.seek(h['Ngroups']*4,1)
    grpidoffset = np.fromfile(f,dtype=np.int32,count=h['Ngroups'])
    grpMass = np.fromfile(f,dtype=np.float32,count=h['Ngroups'])
    grpMinPot = np.fromfile(f,dtype=(np.float32,3),count=h['Ngroups'])
    Mmean200 = np.fromfile(f,dtype=np.float32,count=h['Ngroups'])
    Rmean200 = np.fromfile(f,dtype=np.float32,count=h['Ngroups'])
    Mcrit200 = np.fromfile(f,dtype=np.float32,count=h['Ngroups'])
    Rcrit200 = np.fromfile(f,dtype=np.float32,count=h['Ngroups'])
    Mtop200 = np.fromfile(f,dtype=np.float32,count=h['Ngroups'])
    Rtop200 = np.fromfile(f,dtype=np.float32,count=h['Ngroups'])
    f.seek(h['Ngroups']*4,1)
    f.seek(h['Ngroups']*4,1)
    f.seek(h['Ngroups']*4,1)
    contlen = np.fromfile(f,dtype=np.int32,count=h['Ngroups'])
    contmass = np.fromfile(f,dtype=np.float32,count=h['Ngroups'])
    Nsubstruct = np.fromfile(f,dtype=np.int32,count=h['Ngroups'])
    firsubstruct = np.fromfile(f,dtype=np.int32,count=h['Ngroups'])
    lensubstruct = np.fromfile(f,dtype=np.int32,count=h['Nsubgroups'])
    suboffset = np.fromfile(f,dtype=np.int32,count=h['Nsubgroups'])  
    parent_sub = np.fromfile(f,dtype=np.int32,count=h['Nsubgroups'])
    subgrp_mass = np.fromfile(f,dtype=np.float32,count=h['Nsubgroups'])
    subgrp_pos = np.fromfile(f,dtype=(np.float32,3),count=h['Nsubgroups'])
    f.seek(h['Nsubgroups']*4*3,1)
    subgrp_CM = np.fromfile(f,dtype=(np.float32,3),count=h['Nsubgroups'])
    f.seek(h['Nsubgroups']*4*3,1)
    f.seek(h['Nsubgroups']*4,1)
    f.seek(h['Nsubgroups']*4,1)
    f.seek(h['Nsubgroups']*4,1)
    f.seek(h['Nsubgroups']*4,1)
    idMbnd = np.fromfile(f,dtype=np.uint64,count=h['Nsubgroups'])
    grnr = np.fromfile(f,dtype=np.int32,count=h['Nsubgroups'])   
    return Rmean200, Rcrit200, Rtop200, Mmean200, Mcrit200, Mtop200, grpMass 

def concgrpradii(dir,n):
    
    Rmean200list = []
    Rcrit200list = []
    Rtop200list = []
    Mmean200list = []
    Mcrit200list = []
    Mtop200list = []
    grpMasslist = [] 
    for i in range(n):
          result = grpradii(str(dir)+'.%d'%i)
          Rmean200list.append(result[0])
          Rcrit200list.append(result[1])
          Rtop200list.append(result[2])
          Mmean200list.append(result[3])
          Mcrit200list.append(result[4])
          Mtop200list.append(result[5])
          grpMasslist.append(result[6])
    Rmean200list = np.concatenate(Rmean200list)
    Rcrit200list = np.concatenate(Rcrit200list)
    Rtop200list = np.concatenate(Rtop200list)
    Mmean200list = np.concatenate(Mmean200list)
    Mcrit200list = np.concatenate(Mcrit200list)
    Mtop200list = np.concatenate(Mtop200list)
    grpMasslist = np.concatenate(grpMasslist)
    Rvir = zeros(len(Rmean200list))
    for i in range(len(Rmean200list)):
          mfr = (grpMasslist[i]*1000)/3.0
          Rvir[i] = 7.84*(pow(mfr,1.0/3.0))
    return Rmean200list,Rcrit200list,Rtop200list,Rvir,Mmean200list,Mcrit200list,Mtop200list,grpMasslist 

def subgrpmassfunc(dirname,dirno,boxlen,n):
    
    final = subtabsconc('/home/groups_' + dirno + '/subhalo_tab_' + dirno,nsg)
    sm = final[0]/(1e-10)
    h,bins = histogram(sm,bins=logspace(7,15,n))
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
    pyplot.plot(Mcenter,hdndm, label = ' ' + strRS)
    pyplot.legend()
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.xlabel('M (in Solar Mass)')
    pyplot.ylabel('dn/dM')
    pyplot.title('Halo Mass function from Subgroup Tab')
    pyplot.show()
    return Mcenter,hdndm,bins,ncount,sm,V   

def subgrpstellarmassfunc(dirname,dirno,boxlen,n):
    
    final = subtabsconc('/home/groups_' + dirno + '/subhalo_tab_' + dirno,nsg)

    sm = np.array(final[9])/(1e-10)
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
    pyplot.plot(Mcenter,hdndm, label = ' ' + strRS)
    pyplot.legend()
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.xlabel('M (in Solar Mass)')
    pyplot.ylabel('dn/dM')
    pyplot.title('Halo Mass function from Subgroup Tab')
    pyplot.show()
    return Mcenter,hdndm,bins,ncount,sm,V 

