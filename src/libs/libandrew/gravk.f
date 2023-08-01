      subroutine gravk(mins,maxs)
      parameter(MLL=20)
      parameter(MM=2*MLL+1)

      common/grav/grv,qgrv

      real lcon,ncon
      common/modl1/n,nic,noc,moho,nsl,ifanis,r(222)
     1            ,rho(222),qrho(3,222),g(222),ell(222),eta(222)
      common/modl2/acon(222),qacon(3,222),ccon(222),qccon(3,222)
     1            ,lcon(222),qlcon(3,222),ncon(222),qncon(3,222)
     2            ,fcon(222),qfcon(3,222)
      common/nond/rn,wn,vn,gn,rhobar
      common/eig1/nord1,jcom1,lord1,wcom1,qbar1,cgp1,avert1,ahor1,phis1
     1          ,eif1(222,6)
      common/eig2/nord2,jcom2,lord2,wcom2,qbar2,cgp2,avert2,ahor2,phis2
     1          ,eif2(222,6)

      dimension u1(222),up1(222),v1(222),vp1(222),p1(222),pp1(222),
     1          u2(222),up2(222),v2(222),vp2(222),p2(222),pp2(222)
      dimension q1(3),qp1(3),q2(3),qp2(3)


      equivalence (eif1(1,1),u1(1)),(eif1(1,2),up1(1)),(eif1(1,3),v1(1)),
     1            (eif1(1,4),vp1(1)),(eif1(1,5),p1(1)),(eif1(1,6),pp1(1)),
     2            (eif2(1,1),u2(1)),(eif2(1,2),up2(1)),(eif2(1,3),v2(1)),
     3            (eif2(1,4),vp2(1)),(eif2(1,5),p2(1)),(eif2(1,6),pp2(1))
      equivalence (q1(1),uri),(q1(2),vri),(q1(3),pri),
     1            (qp1(1),upri),(qp1(2),vpri),(qp1(3),ppri),
     2            (q2(1),urii),(q2(2),vrii),(q2(3),prii),
     3            (qp2(1),uprii),(qp2(2),vprii),(qp2(3),pprii)

      double precision bd(5)
      real q(3,222),work(3,222),b(5)
      real grv(222,0:MLL),qgrv(3,222,0:MLL),ga(222,0:MLL),gb(222,0:MLL)
      real gaint(222,0:MLL),gbint(222,0:MLL)

      data n670/180/,n220/202/,nmidc/216/
      data capom/7.292115e-05/
      data coe/1.05688 76266 5/
c             coe=(sqpai/sqrt(5.))*4/3.
      data nint/2/

ca      write(6,*) 'n is ',n

      do j=0,MLL
         do i=1,222
            grv(i,j)=0.d0
            ga(i,j)=0.d0
            gb(i,j)=0.d0
            gaint(i,j)=0.d0
            gbint(i,j)=0.d0
         enddo
      enddo

      do j=0,MLL
         do i=1,222
            do k=1,3
               qgrv(k,i,j)=0.0
            enddo
         enddo
      enddo

      if(jcom1.eq.0.or.jcom2.eq.0) return

      if(jcom1.eq.1)jcom1=3
      if(jcom2.eq.1)jcom2=3
      jcom=3
      if(jcom1.eq.2.and.jcom2.eq.2) jcom=2

      if(jcom.eq.2) return

      fli=float(lord1)
      flii=float(lord2)
      fl3i=fli*(fli+1.)
      fl3ii=flii*(flii+1.)

      capon2=(capom/wn)**2
      delv=1000./vn
      deldis=1000./rn
      wcom=sqrt(wcom1*wcom2)
      omn2=wcom1*wcom2/wn**2
      w0n2=(w0/wn)**2

      j11=1
      j12=3
      if(jcom1.eq.2) then
        j11=2
        j12=2
      endif

      j21=1
      j22=3
      if(jcom2.eq.2) then
        j21=2
        j22=2
      endif

      nstart=1
      if(jcom1.eq.2.or.jcom2.eq.2) nstart=noc

      !APV changed upper limit from n to n-1
      do 1000 iq=nstart,n-1

      call corfac(iq,wcom,jcom,xac,xf,xln)

      iq1=iq+1
      r1=r(iq)
      r2=r(iq1)
ca      hn=r2-r1
ca      if(hn.lt.1.e-4) write(6,*) 'discontinuity at iq is ',iq

      rr=r1
      rr2=rr*rr

      do 1101 i=j11,j12
      i1=2*i-1
      i2=i1+1
      q1(i)=eif1(iq,i1)
 1101 qp1(i)=eif1(iq,i2)*rr

      do 1102 i=j21,j22
      i1=2*i-1
      i2=i1+1
      q2(i)=eif2(iq,i1)
 1102 qp2(i)=eif2(iq,i2)*rr

      rrho=rho(iq)

      if(jcom1.eq.2) goto 1105
      fi=2.*uri-fl3i*vri
 1105 if(jcom2.eq.2) goto 1110
      fii=2.*urii-fl3ii*vrii

 1110 if(jcom.eq.2)goto 1106
      uivii=uri*vrii
      uiivi=urii*vri
      fivii=fi*vrii
      fiivi=fii*vri
      upivii=upri*vrii
      upiivi=uprii*vri
      vpiuii=vpri*urii
      vpiiui=vprii*uri

 1106 if(jcom1.ne.jcom2) goto 1120

      rkp1=0.
      rkpp1=0.

      if(jcom1.eq.2) goto 1150

      ur2=uri*urii
      fiuii=fi*urii
      fiiui=fii*uri

      rkp0a=rrho*ur2
      rkp1a=.5*rrho*(vpiiui+uivii-upivii-2.*fivii)
      rkp1b=.5*rrho*(vpiuii+uiivi-upiivi-2.*fiivi)
      rkpp0=rrho*(fiuii+fiiui)
      rkpp1a=.5*rrho*uivii
      rkpp1b=.5*rrho*uiivi

      goto 1150

 1120 if(jcom1.eq.2)then
        rkp1=-.5*rrho*(vpiuii+uiivi-upiivi-2.*fiivi)
        rkpp1=-.5*rrho*uiivi
      else
        rkp1=.5*rrho*(vpiiui+uivii-upivii-2.*fivii)
        rkpp1=.5*rrho*uivii
      endif

 1150 continue

 
c
c   compute kernels G(a)s and G(b)s
c

      minsn=abs(fli-flii)
      maxsn=fli+flii

      if(minsn.gt.maxs) return
      if(maxsn.gt.maxs) maxsn=maxs

      if(minsn.lt.mins) minsn=mins

      do 1170 i=mins,maxs

      fsi=float(i)
      fs2i=1.+fsi
      fs3i=fsi*(1.+fsi)

      if(jcom1.eq.jcom2) then
        add1=.5*(rkpp1a*fs2i-rkp1a)*(fl3ii+fs3i-fl3i)
        add2=.5*(rkpp1b*fs2i-rkp1b)*(fl3i+fs3i-fl3ii)
        add3=-1.*(rkpp0*fs2i+rkp0a*fs3i)
        ga(iq,i)=add1+add2+add3
        add1=.5*(rkpp1a*fsi+rkp1a)*(fl3ii+fs3i-fl3i)
        add2=.5*(rkpp1b*fsi+rkp1b)*(fl3i+fs3i-fl3ii)
        add3=-1.*(rkpp0*fsi-rkp0a*fs3i)
        gb(iq,i)=add1+add2+add3
      endif

      if(jcom1.eq.2.or.jcom2.eq.2) then
        ga(iq,i)=(rkpp1*fs2i-rkp1)
        gb(iq,i)=(rkpp1*fsi+rkp1)
      endif
         
 1170 continue

 1000 continue 

c
c Now integrate kernels
c

      do is=mins,maxs

      fsi=float(is)
      call rspln(1,n,r,gb(1,is),q,work)

      gbint(1,is)=gb(1,is)/(fsi+1.0)
      gbint(1,is)=0.0
      gaint(1,is)=0.
      gaint(n,is)=0.

      do iq=nstart,n-1

      r1=r(iq)
      r2=r(iq+1)
      dr=(r2-r1)/2.
      if(dr.lt.1.e-4) then
        gbint(iq+1,is)=gbint(iq,is)
        goto 245
      endif
      rr=r1+dr

      afac=gb(iq,is)-q(1,iq)*r1+q(2,iq)*r1**2-q(3,iq)*r1**3
      bfac=q(1,iq)-2.*q(2,iq)*r1+3.*q(3,iq)*r1**2
      cfac=q(2,iq)-3.*q(3,iq)*r1
      dfac=q(3,iq)


      fac=gb(iq,is)+dr*(q(1,iq)+dr*(q(2,iq)+dr*q(3,iq)))
      gdum=afac+(bfac*rr)+(cfac*(rr**2))+(dfac*(rr**3))

      gbint(iq+1,is)=grintb(afac,bfac,cfac,dfac,r2,fsi)
     .             -((r1/r2)**(is+1))*grintb(afac,bfac,cfac,dfac,r1,fsi)
     .             +((r1/r2)**(is+1))*gbint(iq,is)
ca      write(6,*) iq+1,gbint(iq+1,is)
245   continue

      enddo

      fsi=float(is)
      call rspln(1,n,r,ga(1,is),q,work)

      gaint(n,is)=0.

      do iq=n-1,nstart+1,-1

      r1=r(iq)
      r2=r(iq+1)
      dr=(r2-r1)/2.
      rr=r1+dr

      afac=ga(iq,is)-q(1,iq)*r1+q(2,iq)*r1**2-q(3,iq)*r1**3
      bfac=q(1,iq)-2.*q(2,iq)*r1+3.*q(3,iq)*r1**2
      cfac=q(2,iq)-3.*q(3,iq)*r1
      dfac=q(3,iq)

      fac=ga(iq,is)+dr*(q(1,iq)+dr*(q(2,iq)+dr*q(3,iq)))
      gdum=afac+(bfac*rr)+(cfac*(rr**2))+(dfac*(rr**3))

      gaint(iq,is)=((r1/r2)**(is))*grinta(afac,bfac,cfac,dfac,r2,fsi)
     .             -grinta(afac,bfac,cfac,dfac,r1,fsi)
     .             +((r1/r2)**(is))*gaint(iq+1,is)
      if(abs(gaint(iq,is)).le.1.e-35) gaint(iq,is)=0.
ca      write(6,*) iq,gaint(iq,is)

      enddo

      gaint(nstart,is)=gaint(nstart+1,is)

      do iq=nstart,n
         grv(iq,is)=(4./(2.0*fsi+1.))*(gaint(iq,is)-gbint(iq,is))
         if(abs(grv(iq,is)).le.1.e-30) grv(iq,is)=0.
      enddo

      call rspln(1,n,r,grv(1,is),qgrv(1,1,is),work)

      enddo
      
      return
      end

c--------------------------------------------------------------------
