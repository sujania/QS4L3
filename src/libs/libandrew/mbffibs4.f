c--------------------------------------------------------
      subroutine mbffibs4(lu,ifbin,ibuf,nbyts,j,m,irec)
      dimension ibuf(*)
      call mbffi(lu,ifbin,ibuf,nbyts,j,m,irec)
      if(j.ne.2) return
      call byswap4(ibuf,(m+3)/4)
      return
      end

