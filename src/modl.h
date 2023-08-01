      parameter (NKNTS=222)
      parameter (MXDSC=30)
      real lcon,ncon
      common/modl1/n,nic,noc,moho,nsl,ifanis,r(NKNTS)
     1            ,rho(NKNTS),qrho(3,NKNTS),g(NKNTS),ell(NKNTS),eta(NKNTS)
      common/modl2/acon(NKNTS),qacon(3,NKNTS),ccon(NKNTS),qccon(3,NKNTS)
     1            ,lcon(NKNTS),qlcon(3,NKNTS),ncon(NKNTS),qncon(3,NKNTS)
     2            ,fcon(NKNTS),qfcon(3,NKNTS)
      common/modl3/qshear(NKNTS),qkappa(NKNTS)
      common/modl4/ndisc,ndsc(MXDSC)
      common/modl5/pres(NKNTS)
      common/modl6/crn,cvn,crb
