      subroutine fetch(nord,iq,lord,kbuf,ifexst)
c     nord = over tone number of the mode
c     iq = mode type, 't' for toroidal, and 's' for spheroidal
c     lord = angular degree of the mode
c     kbuf = 1 or 2, tells the routine which common block
c            to store the eigenfunction in
c     ifexist = 1 if the input mode is contained in the file
c
c     Contents of the eigenfunction common block:
c     
c     n     = overtone number
c     jcom  = mode type, 1 = radial, 2 = toroidal, 3 = spheroidal
c     l     = angular order 
c     om    = eigenfrequency in (normalized) rad/s
c     cpg   = group velocity of the mode (units?)
c     avert = vertical acceleration at the sea floor
c     ahor  = horizontal acceleration at the sea floor
c     phis  = gravitational acceleration at the sea floor
c     u     = vertical eigenfunction
c     up    = radial derivative of vertical eigenfunction
c     v     = tangential eigenfunction (V for spheroidals, W 
c             for spheroidals)
c     vp    = radial derivative of vertical eigenfunction
c     ph    = gravitational eigenfunction
c     php   = radial derivative of gravitational eigenfunction
c
c
      character*1 iq,it,is
      integer*2 nrec,maxn,maxl
      common/eig1/n1,jcom1,l1,om1,q1,cgp1,avert1,ahor1,phis1
     1        ,u1(222),up1(222),v1(222),vp1(222),ph1(222),php1(222)
      common/eig2/n2,jcom2,l2,om2,q2,cgp2,avert2,ahor2,phis2
     1        ,u2(222),up2(222),v2(222),vp2(222),ph2(222),php2(222)
      common/get/nrec(500,2),maxn(500,2),maxl(2),llu
      data it,is/'t','s'/
c     check that the eigenfunctions are being stored in 
c     one of the possible common blocks
      if(kbuf.ne.1.and.kbuf.ne.2) stop 'error 7 in fetch'
      ifexst=0
c     is the mode toroidal (itors = 1) or spheroidal (itors = 2)
      itors=0
      if(iq.eq.it) itors=1
      if(iq.eq.is) itors=2
      if(itors.ne.1.and.itors.ne.2) pause 'error 5 in fetch'
c     get the maximum value of l for the appropriate mode type
      itemp=maxl(itors)
c     if lord is greater than maxl, then return with 
c     ifexst = 0
      if(lord.gt.itemp) goto 99
c     work out how many overtone numbers there are 
c     for the given value of l using the array
c     maxn.
      ind=lord+1
      itemp=maxn(ind,itors)
c     if nord is greater than maxn for the given l return
c     with ifexst = 0
      if(nord.gt.itemp) goto 99
c     get the record number of the mode
      nr=nrec(ind,itors)+nord
c     call sysio(ipblk,92,llu,n,532,nr)
c     call ioerr(ipblk,ier)
c     read the eigenfunction record, and out the 
c     data into the appropriate common block

      if(kbuf.eq.1) then
        call bffi(llu,1,n1,1341*4,jstat,nread,nr)
        call byswap4(n1,1341)
        jcom=jcom1
        l=l1
        n=n1
      else
        call bffi(llu,1,n2,1341*4,jstat,nread,nr)
        call byswap4(n2,1341)
        jcom=jcom2
        l=l2
        n=n2
      endif
      if(n.ne.nord) goto 99
      if(iq.eq.it.and.jcom.ne.2) goto 99
      if(iq.eq.is.and.l.eq.0.and.jcom.ne.1) goto 99
      if(iq.eq.is.and.l.ne.0.and.jcom.ne.3) goto 99
      if(lord.ne.l) goto 99
      ifexst=1
c     write(6,1234)n1,jcom1,l1,om1,q1,cgp1,avert1,ahor1,phis1
c    1            ,n2,jcom2,l2,om2,q2,cgp2,avert2,ahor2,phis2
 1234 format(3i3,6e10.3)
c     write(6,1233)
 1233 format('?')
c     read(5,1232)ijk
 1232 format(i2)
c     if(ijk.eq.0)goto 99

c     write(6,1235)(u1(i),up1(i),v1(i),vp1(i),ph1(i),php1(i),
c    1              u2(i),up2(i),v2(i),vp2(i),ph2(i),php2(i),i=1,222)
c1235 format(6e12.5)

   99 return
      end
