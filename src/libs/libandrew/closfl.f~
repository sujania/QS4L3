c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine closfl(lufl,istat)
      include 'openfile.h'
      common/optdr/ifopt(nlu)
      if(jchn(lufl).le.0) return

      if(ifopt(lufl).ne.0) then
        call cmtio(jchn(lufl),5,1,ires,ierrno)
        if(ierrno.ne.0) call cperror('cmtio:rewind in closfl')
        ifopt(lufl)=0
        call oasopen()
        call oassend('offline',mess,lmess,-1)
        call oassend('dismount',mess,lmess,-1)
        call oasclose()
      endif
      call cclose(jchn(lufl),ires,ierrno)
      if(ierrno.ne.0) then
        call cperror('cclose in closfl')
        print*,'logical unit ',lufl,' channel ',jchn(lufl)
      endif

      jchn(lufl)=0
      return
      end
