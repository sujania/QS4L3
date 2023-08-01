
c---------------------------------------------------------------

      subroutine mbffi(ifl,ifbin,ibuf,nbytes,istat,nread,irec)
      dimension ibuf(*)

      include 'mopenfile.h'


      if(irec.le.0) stop 'mbffi: only direct access supported'
      lu=moplus(ifl)
      if(lu.gt.0) then
      else if(lu.eq.0) then
        stop 'mbffi: lu=0'
      else
        call mreopen(ifl,lu,istat)
      endif
      if(lu.le.0) stop 'mbffo: lu le 0'
      call bffi(lu,ifbin,ibuf,nbytes,istat,nread,irec)
      mopknt=mopknt+1
      mopknts(lu)=mopknt
      return
      end
