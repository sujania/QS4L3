
c---------------------------------------------------------------
      subroutine mreopen(ifl,lu,istat)

      include 'mopenfile.h'

      if(mopnact.ge.MXOPEN) then
        lmin=0
        minknt=z'7fffffff'
        do i=1,40
          if(mopifls(i).ne.0) then
            if(mopknts(i).lt.minknt) then
              lmin=i
              minknt=mopknts(i)
            endif
          endif
        enddo
        if(lmin.eq.0) stop 'yyyy'
        call closfl(lmin,istat)
        mopknts(lmin)=0
        moplus(mopifls(lmin))=-1
        mopifls(lmin)=0
        mopnact=mopnact-1
      endif
      call opnflc(moplus(ifl),mopnames(ifl),mopiap(ifl)
     1     ,0,0,istat,moplrec(ifl),mopinew(ifl))
      if(istat.ne.0) then
c       write(6,*) 'mreopen: open error:'//mopnames(ifl)(1:istlen(mopnames(ifl))),mopiap(ifl)
        istat=9
        moplus(ifl)=-1
        lu=moplus(ifl)
        return
      endif
      mopinew(ifl)=0
      mopifls(moplus(ifl))=ifl
      mopnact=mopnact+1
      lu=moplus(ifl)
c     write(6,*) 'mreopen: ifl=',ifl,'  lu=',lu
      return
      end
