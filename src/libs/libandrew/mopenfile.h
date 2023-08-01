      parameter (MXOPEN=16)
      parameter (MXFILS=2048)
      character*80 mopnames
      common/mopnprms/mopnfls,mopnames(MXFILS),moplus(MXFILS)
     1     ,mopnact,mopknt,mopifls(40),mopknts(40)
     1   ,mopiap(MXFILS),moplrec(MXFILS),mopinew(MXFILS)

