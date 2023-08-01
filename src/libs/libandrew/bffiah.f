c--------------------------------------------------------

      subroutine bffiah(lu,ifbin,iah,nbyts,jah,m,irec)
      dimension iah(0:*)
      include 'ahheader.h'
      call bffi(lu,ifbin,iah,nbyts,jah,m,irec)
      if(jah.ne.2) return
      call byswap4(iah,LAHHDR)
      call byswap4(iah(AAHACD),2)
      call byswap4(iah(AAHACH),2)
      call byswap4(iah(AAHATP),2)
      call byswap4(iah(AAHAEC),MAHCOM)
      call byswap4(iah(AAHADC),MAHCOM)
      call byswap4(iah(AAHALG),MAHLOG)
      if(m.gt.LAHHDR*4) call byswap4(iah(LAHHDR),(m+3)/4-LAHHDR)
      return
      end

