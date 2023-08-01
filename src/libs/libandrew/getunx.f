c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      character*80 function getunx(ident,nni,nbyts)
      character*(*) ident
      character*80 getgnl
      include 'getgnl.h'
      include 'getunx.h'
      getunx=getgnl(ident,nni,nbyts
     1  ,swnam,cline,deflt
     1  ,nswt,lswnam,ldeflt,ifreq,nmin,nmax
     1  ,icnt,iopn,iptr,itable,lcom)

      return
      end
