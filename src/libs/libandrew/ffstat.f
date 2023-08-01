c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine ffstat(lufl,lent,iftyp,isize)
      include 'openfile.h'
      !Get size (bytes),mode, and record length of file

c     all cfstat(jchn(lufl),isize,iftyp,ierrno)
c        if(ierrno.ne.0) call check('cfstat on ffstat')
      lent=jrecl(lufl)
      return
      end
