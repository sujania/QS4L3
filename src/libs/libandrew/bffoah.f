c--------------------------------------------------------

      subroutine bffoah(lu,ifbin,iah,nbyts,jah,irec)
      dimension iah(0:*)
      include 'ahheader.h'
      call byswap4(iah,LAHHDR)
      call byswap4(iah(AAHACD),2)
      call byswap4(iah(AAHACH),2)
      call byswap4(iah(AAHATP),2)
      call byswap4(iah(AAHAEC),MAHCOM)
      call byswap4(iah(AAHADC),MAHCOM)
      call byswap4(iah(AAHALG),MAHLOG)
      if(nbyts.gt.LAHHDR*4)  call byswap4(iah(LAHHDR),(nbyts+3)/4-LAHHDR)
      call bffo(lu,ifbin,iah,nbyts,jah,irec)
      call byswap4(iah,LAHHDR)
      call byswap4(iah(AAHACD),2)
      call byswap4(iah(AAHACH),2)
      call byswap4(iah(AAHATP),2)
      call byswap4(iah(AAHAEC),MAHCOM)
      call byswap4(iah(AAHADC),MAHCOM)
      call byswap4(iah(AAHALG),MAHLOG)
      if(nbyts.gt.LAHHDR*4)  call byswap4(iah(LAHHDR),(nbyts+3)/4-LAHHDR)
      return
      end
