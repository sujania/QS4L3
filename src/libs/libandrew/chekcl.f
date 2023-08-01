c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine chekcl(string)
      character*(*) string

      include 'getgnl.h'
      include 'getunx.h'
      external myiargc,mygetarg
      call chekgl(string,myiargc,mygetarg,ierr
     1  ,swnam,cline,deflt
     1  ,nswt,lswnam,ldeflt,ifreq,nmin,nmax,icnt,iopn,iptr,itable,lcom)
      if(ierr.ne.0) call exit(2)
      return
      end

      integer function myiargc()
      myiargc=iargc()
      return
      end

      subroutine mygetarg(iarg,arg)
      character*(*) arg
      call getarg(iarg,arg)
      return
      end

