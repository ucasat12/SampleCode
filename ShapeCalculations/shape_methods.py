import numpy

def shapecalc_uw3d(sgrno,r200lj,cen,dmpos1,boxlen,permax,nmin,errtol,ntol):

    rdiff = dmpos1 - cen
    Rmask = 1 - (numpy.abs(rdiff/permax)).astype(numpy.int32)
    Roffset = numpy.ma.masked_array(boxlen - 2.0*(dmpos1),Rmask)
    DMpos = numpy.array(dmpos1 + Roffset.filled(0)).copy()
                                  
    rdiff = DMpos - cen
    rdiff2 = pow(rdiff,2.0)
    rmeas = pow(rdiff2[:,0] + rdiff2[:,1] + rdiff2[:,2],1.0/2.0)
    piv1 = numpy.argsort(rmeas)
    rdiff = rdiff[piv1]
    rmeas = rmeas[piv1]
    chk_rsq = (rmeas > 0)
    rdiff = rdiff[chk_rsq]
    rmeas = rmeas[chk_rsq]
    rin = r200lj*(rmeas.max())    
    chkrin = (rmeas <= rin)
    redDMpos = rdiff[chkrin]
                    
    qin = 1
    sin = 1
    ain = rin/pow(qin*sin,1.0/3.0)
    bin = rin*qin/pow(qin*sin,1.0/3.0)
    cin = rin*sin/pow(qin*sin,1.0/3.0)
    nin = len(redDMpos)
    if (nin > 0):
       Rwght_sq = pow(redDMpos[:,0]/ain,2.0) + pow(redDMpos[:,1]/bin,2.0) + pow(redDMpos[:,2]/cin,2.0)
    errmax = 1.0
    nc = 0              
    pareig = numpy.empty(shape=3,
                         dtype = [
                                  ('eigno',numpy.float32),
                                  ('eigvec',(numpy.float32,3))])
    qfin = 0
    sfin = 0

    eigvaldata_fin = numpy.zeros(3)
    eigvecdata_fin = numpy.zeros((3,3))     
    eigvecdata_rotfin = numpy.zeros((3,3))
    nIeigvec0t = numpy.array([[1,0,0],[0,1,0],[0,0,1]])
    Itens = numpy.zeros((3,3))
    if (len(redDMpos) >= nmin):

              Itens = numpy.dot(redDMpos.T,redDMpos)            
              Ieigs = numpy.linalg.eig(Itens)
              Ieigval = numpy.array(Ieigs[0])
              Ieigvec = numpy.array(Ieigs[1]).T

              pareig['eigno'] = Ieigval
              pareig['eigvec'] = Ieigvec 
              pareig.sort(order='eigno')
              nIeigval = numpy.array((pareig['eigno'][::-1]).copy())
              nIeigvec = numpy.array((pareig['eigvec'][::-1]).copy())

              qp = qin
              sp = sin
              qin = pow(nIeigval[1]/nIeigval[0],1.0/2.0)
              sin = pow(nIeigval[2]/nIeigval[0],1.0/2.0)
              ain = rin/pow(qin*sin,1.0/3.0)
              bin = rin*qin/pow(qin*sin,1.0/3.0)
              cin = rin*sin/pow(qin*sin,1.0/3.0)
              nIeigvec0t = numpy.dot(nIeigvec,nIeigvec0t)
              rdiff = (numpy.dot(nIeigvec,rdiff.T)).T
              rmeas2 = pow(rdiff[:,0]/ain,2.0) + pow(rdiff[:,1]/bin,2.0) + pow(rdiff[:,2]/cin,2.0)
              chkrin = (rmeas2 <= 1.0)
              redDMpos = rdiff[chkrin]
              Rwght_sq = rmeas2[chkrin]
              
              errmax = numpy.max([(qin-qp)**2/qp**2,(sin-sp)**2/sp**2,(1 - numpy.abs(nIeigvec[0,0]))**2,(1 - numpy.abs(nIeigvec[1,1]))**2,(1 - numpy.abs(nIeigvec[2,2]))**2])
              nin = len(redDMpos)
              nc = nc + 1               
    if (len(redDMpos) >= nmin):
       qfin = pow(nIeigval[1]/nIeigval[0],1.0/2.0)
       sfin = pow(nIeigval[2]/nIeigval[0],1.0/2.0)
       eigvaldata_fin = nIeigval
       eigvecdata_fin = nIeigvec0t
       eigvecdata_rotfin = nIeigvec
    return qfin, sfin, eigvaldata_fin, eigvecdata_fin, errmax, nin, nc, eigvecdata_rotfin, Itens

def shapecalc_rd3d(sgrno,r200lj,cen,dmpos1,boxlen,permax,nmin,errtol,ntol):

    rdiff = dmpos1 - cen
    Rmask = 1 - (numpy.abs(rdiff/permax)).astype(numpy.int32)
    Roffset = numpy.ma.masked_array(boxlen - 2.0*(dmpos1),Rmask)
    DMpos = numpy.array(dmpos1 + Roffset.filled(0)).copy()
                                  
    rdiff = DMpos - cen
    rdiff2 = pow(rdiff,2.0)
    rmeas = pow(rdiff2[:,0] + rdiff2[:,1] + rdiff2[:,2],1.0/2.0)
    piv1 = numpy.argsort(rmeas)
    rdiff = rdiff[piv1]
    rmeas = rmeas[piv1]
    chk_rsq = (rmeas > 0)
    rdiff = rdiff[chk_rsq]
    rmeas = rmeas[chk_rsq]

    rin = r200lj*(rmeas.max())
    chkrin = (rmeas <= rin)
    redDMpos = rdiff[chkrin]
                    
    qin = 1
    sin = 1
    ain = rin/pow(qin*sin,1.0/3.0)
    bin = rin*qin/pow(qin*sin,1.0/3.0)
    cin = rin*sin/pow(qin*sin,1.0/3.0)
    nin = len(redDMpos)
    if (nin > 0):
       Rwght_sq = pow(redDMpos[:,0]/ain,2.0) + pow(redDMpos[:,1]/bin,2.0) + pow(redDMpos[:,2]/cin,2.0)
    errmax = 1.0
    nc = 0              
    pareig = numpy.empty(shape=3,
                         dtype = [
                                  ('eigno',numpy.float32),
                                  ('eigvec',(numpy.float32,3))])
    qfin = 0
    sfin = 0

    eigvaldata_fin = numpy.zeros(3)
    eigvecdata_fin = numpy.zeros((3,3))     
    eigvecdata_rotfin = numpy.zeros((3,3))
    nIeigvec0t = numpy.array([[1,0,0],[0,1,0],[0,0,1]])    
    Itens = numpy.zeros((3,3))
    if (len(redDMpos) >= nmin):
             
              redDMpos1 = numpy.zeros((len(redDMpos),3))
              redDMpos1[:,0] = redDMpos[:,0]/Rwght_sq[:]
              redDMpos1[:,1] = redDMpos[:,1]/Rwght_sq[:]
              redDMpos1[:,2] = redDMpos[:,2]/Rwght_sq[:]

              Itens = numpy.dot(redDMpos.T,redDMpos1)            
              Ieigs = numpy.linalg.eig(Itens)
              Ieigval = numpy.array(Ieigs[0])
              Ieigvec = numpy.array(Ieigs[1]).T

              pareig['eigno'] = Ieigval
              pareig['eigvec'] = Ieigvec 
              pareig.sort(order='eigno')
              nIeigval = numpy.array((pareig['eigno'][::-1]).copy())
              nIeigvec = numpy.array((pareig['eigvec'][::-1]).copy())

              qp = qin
              sp = sin
              qin = pow(nIeigval[1]/nIeigval[0],1.0/2.0)
              sin = pow(nIeigval[2]/nIeigval[0],1.0/2.0)
              ain = rin/pow(qin*sin,1.0/3.0)
              bin = rin*qin/pow(qin*sin,1.0/3.0)
              cin = rin*sin/pow(qin*sin,1.0/3.0)
              nIeigvec0t = numpy.dot(nIeigvec,nIeigvec0t)
              rdiff = (numpy.dot(nIeigvec,rdiff.T)).T
              rmeas2 = pow(rdiff[:,0]/ain,2.0) + pow(rdiff[:,1]/bin,2.0) + pow(rdiff[:,2]/cin,2.0)
              chkrin = (rmeas2 <= 1.0)
              redDMpos = rdiff[chkrin]
              Rwght_sq = rmeas2[chkrin]
              
              errmax = numpy.max([(qin-qp)**2/qp**2,(sin-sp)**2/sp**2,(1 - numpy.abs(nIeigvec[0,0]))**2,(1 - numpy.abs(nIeigvec[1,1]))**2,(1 - numpy.abs(nIeigvec[2,2]))**2])
              nin = len(redDMpos)
              nc = nc + 1             
    if (len(redDMpos) >= nmin):
       qfin = pow(nIeigval[1]/nIeigval[0],1.0/2.0)
       sfin = pow(nIeigval[2]/nIeigval[0],1.0/2.0)
       eigvaldata_fin = nIeigval
       eigvecdata_fin = nIeigvec0t
       eigvecdata_rotfin = nIeigvec
    return qfin, sfin, eigvaldata_fin, eigvecdata_fin, errmax, nin, nc, eigvecdata_rotfin, Itens

def shapecalc_uw3diter(sgrno,r200lj,cen,dmpos1,boxlen,permax,nmin,errtol,ntol):

    rdiff = dmpos1 - cen
    Rmask = 1 - (numpy.abs(rdiff/permax)).astype(numpy.int32)
    Roffset = numpy.ma.masked_array(boxlen - 2.0*(dmpos1),Rmask)
    DMpos = numpy.array(dmpos1 + Roffset.filled(0)).copy()
                                  
    rdiff = DMpos - cen
    rdiff2 = pow(rdiff,2.0)
    rmeas = pow(rdiff2[:,0] + rdiff2[:,1] + rdiff2[:,2],1.0/2.0)
    piv1 = numpy.argsort(rmeas)
    rdiff = rdiff[piv1]
    rmeas = rmeas[piv1]
    chk_rsq = (rmeas > 0)
    rdiff = rdiff[chk_rsq]
    rmeas = rmeas[chk_rsq]

    rin = r200lj    
    chkrin = (rmeas <= rin)
    redDMpos = rdiff[chkrin]
                    
    rin = r200lj
    qin = 1
    sin = 1
    ain = rin/pow(qin*sin,1.0/3.0)
    bin = rin*qin/pow(qin*sin,1.0/3.0)
    cin = rin*sin/pow(qin*sin,1.0/3.0)
    nin = len(redDMpos)
    if (nin > 0):
       Rwght_sq = pow(redDMpos[:,0]/ain,2.0) + pow(redDMpos[:,1]/bin,2.0) + pow(redDMpos[:,2]/cin,2.0)
    errmax = 1.0
    nc = 0              
    pareig = numpy.empty(shape=3,
                         dtype = [
                                  ('eigno',numpy.float32),
                                  ('eigvec',(numpy.float32,3))])
    qfin = 0
    sfin = 0

    eigvaldata_fin = numpy.zeros(3)
    eigvecdata_fin = numpy.zeros((3,3))     
    eigvecdata_rotfin = numpy.zeros((3,3))
    nIeigvec0t = numpy.array([[1,0,0],[0,1,0],[0,0,1]])
    Itens = numpy.zeros((3,3))
    while ((nin >= nmin) & (errmax >= errtol) & (nc <= ntol)):
    
          if (len(redDMpos) >= nmin):

              Itens = numpy.dot(redDMpos.T,redDMpos)            
              Ieigs = numpy.linalg.eig(Itens)
              Ieigval = numpy.array(Ieigs[0])
              Ieigvec = numpy.array(Ieigs[1]).T

              pareig['eigno'] = Ieigval
              pareig['eigvec'] = Ieigvec 
              pareig.sort(order='eigno')
              nIeigval = numpy.array((pareig['eigno'][::-1]).copy())
              nIeigvec = numpy.array((pareig['eigvec'][::-1]).copy())

              qp = qin
              sp = sin
              qin = pow(nIeigval[1]/nIeigval[0],1.0/2.0)
              sin = pow(nIeigval[2]/nIeigval[0],1.0/2.0)
              ain = rin/pow(qin*sin,1.0/3.0)
              bin = rin*qin/pow(qin*sin,1.0/3.0)
              cin = rin*sin/pow(qin*sin,1.0/3.0)
              nIeigvec0t = numpy.dot(nIeigvec,nIeigvec0t)
              rdiff = (numpy.dot(nIeigvec,rdiff.T)).T
              rmeas2 = pow(rdiff[:,0]/ain,2.0) + pow(rdiff[:,1]/bin,2.0) + pow(rdiff[:,2]/cin,2.0)
              chkrin = (rmeas2 <= 1.0)
              redDMpos = rdiff[chkrin]
              Rwght_sq = rmeas2[chkrin]
              
              errmax = numpy.max([(qin-qp)**2/qp**2,(sin-sp)**2/sp**2,(1 - numpy.abs(nIeigvec[0,0]))**2,(1 - numpy.abs(nIeigvec[1,1]))**2,(1 - numpy.abs(nIeigvec[2,2]))**2])
              nin = len(redDMpos)
              nc = nc + 1             
    if ((nin >= nmin) & (nc <= ntol)):
       qfin = pow(nIeigval[1]/nIeigval[0],1.0/2.0)
       sfin = pow(nIeigval[2]/nIeigval[0],1.0/2.0)
       eigvaldata_fin = nIeigval
       eigvecdata_fin = nIeigvec0t
       eigvecdata_rotfin = nIeigvec
    return qfin, sfin, eigvaldata_fin, eigvecdata_fin, errmax, nin, nc, eigvecdata_rotfin, Itens

def shapecalc_rd3diter(sgrno,r200lj,cen,dmpos1,boxlen,permax,nmin,errtol,ntol):

    rdiff = dmpos1 - cen
    Rmask = 1 - (numpy.abs(rdiff/permax)).astype(numpy.int32)
    Roffset = numpy.ma.masked_array(boxlen - 2.0*(dmpos1),Rmask)
    DMpos = numpy.array(dmpos1 + Roffset.filled(0)).copy()
                                  
    rdiff = DMpos - cen
    rdiff2 = pow(rdiff,2.0)
    rmeas = pow(rdiff2[:,0] + rdiff2[:,1] + rdiff2[:,2],1.0/2.0)
    piv1 = numpy.argsort(rmeas)
    rdiff = rdiff[piv1]
    rmeas = rmeas[piv1]
    chk_rsq = (rmeas > 0)
    rdiff = rdiff[chk_rsq]
    rmeas = rmeas[chk_rsq]

    rin = r200lj    
    chkrin = (rmeas <= rin)
    redDMpos = rdiff[chkrin]
                    
    rin = r200lj
    qin = 1
    sin = 1
    ain = rin/pow(qin*sin,1.0/3.0)
    bin = rin*qin/pow(qin*sin,1.0/3.0)
    cin = rin*sin/pow(qin*sin,1.0/3.0)
    nin = len(redDMpos)
    if (nin > 0):
       Rwght_sq = pow(redDMpos[:,0]/ain,2.0) + pow(redDMpos[:,1]/bin,2.0) + pow(redDMpos[:,2]/cin,2.0)
    errmax = 1.0
    nc = 0              
    pareig = numpy.empty(shape=3,
                         dtype = [
                                  ('eigno',numpy.float32),
                                  ('eigvec',(numpy.float32,3))])
    qfin = 0
    sfin = 0

    eigvaldata_fin = numpy.zeros(3)
    eigvecdata_fin = numpy.zeros((3,3))     
    eigvecdata_rotfin = numpy.zeros((3,3))
    nIeigvec0t = numpy.array([[1,0,0],[0,1,0],[0,0,1]])
    Itens = numpy.zeros((3,3))
    while ((nin >= nmin) & (errmax >= errtol) & (nc <= ntol)):
    
          if (len(redDMpos) >= nmin):
             
              redDMpos1 = numpy.zeros((len(redDMpos),3))
              redDMpos1[:,0] = redDMpos[:,0]/Rwght_sq[:]
              redDMpos1[:,1] = redDMpos[:,1]/Rwght_sq[:]
              redDMpos1[:,2] = redDMpos[:,2]/Rwght_sq[:]
              Itens = numpy.dot(redDMpos.T,redDMpos1)            
              Ieigs = numpy.linalg.eig(Itens)
              Ieigval = numpy.array(Ieigs[0])
              Ieigvec = numpy.array(Ieigs[1]).T

              pareig['eigno'] = Ieigval
              pareig['eigvec'] = Ieigvec 
              pareig.sort(order='eigno')
              nIeigval = numpy.array((pareig['eigno'][::-1]).copy())
              nIeigvec = numpy.array((pareig['eigvec'][::-1]).copy())

              qp = qin
              sp = sin
              qin = pow(nIeigval[1]/nIeigval[0],1.0/2.0)
              sin = pow(nIeigval[2]/nIeigval[0],1.0/2.0)
              ain = rin/pow(qin*sin,1.0/3.0)
              bin = rin*qin/pow(qin*sin,1.0/3.0)
              cin = rin*sin/pow(qin*sin,1.0/3.0)
              nIeigvec0t = numpy.dot(nIeigvec,nIeigvec0t)
              rdiff = (numpy.dot(nIeigvec,rdiff.T)).T
              rmeas2 = pow(rdiff[:,0]/ain,2.0) + pow(rdiff[:,1]/bin,2.0) + pow(rdiff[:,2]/cin,2.0)
              chkrin = (rmeas2 <= 1.0)
              redDMpos = rdiff[chkrin]
              Rwght_sq = rmeas2[chkrin]
              
              errmax = numpy.max([(qin-qp)**2/qp**2,(sin-sp)**2/sp**2,(1 - numpy.abs(nIeigvec[0,0]))**2,(1 - numpy.abs(nIeigvec[1,1]))**2,(1 - numpy.abs(nIeigvec[2,2]))**2])
              nin = len(redDMpos)
              nc = nc + 1  
             
    if ((nin >= nmin) & (nc <= ntol)):
       qfin = pow(nIeigval[1]/nIeigval[0],1.0/2.0)
       sfin = pow(nIeigval[2]/nIeigval[0],1.0/2.0)
       eigvaldata_fin = nIeigval
       eigvecdata_fin = nIeigvec0t
       eigvecdata_rotfin = nIeigvec
    return qfin, sfin, eigvaldata_fin, eigvecdata_fin, errmax, nin, nc, eigvecdata_rotfin, Itens
