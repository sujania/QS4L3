c     subroutine for reading in and processing the PREM
c     reference model
      subroutine modl(lumod,lupr,rn,wn,vn,gn,rhobar)
      include 'modl.h'
      common/work/wk(3,NKNTS),coef(5),y(4),y1(4),y2(4),y3(4),ydot(4)
c
      data bigg,capom/6.6723e-11,7.292115e-05/
      data pi,tau,rhon/3.1415926,1000.0,5514.3/
c
      read(lumod,7001) ifdeck,ifanis
 7001 format(49x,i1,8x,i1)
      backspace lumod
      if(ifdeck.eq.0) goto 200
      read(lumod,1000) n,nic,noc,moho
 1000 format(/,4i5)
c
      if(ifanis.eq.0) read(lumod,1001)(r(i),rho(i)
     1     ,ccon(i),lcon(i),qshear(i),qkappa(i),i=1,n)
 1001 format(6f10.2)
      if(ifanis.ne.0) read(lumod,1002)(r(i),rho(i)
     1     ,ccon(i),lcon(i),qshear(i),qkappa(i),acon(i)
     2     ,ncon(i),fcon(i),i=1,n)
 1002 format(f8.0,8f9.2)

      goto 300
c
  200 read(lumod,1005) nreg,nic,noc,rx,moho
 1005 format(3i5,e15.8,i5)
      print *, "nreg", nreg, "nic", nic, "noc", noc, "rx", rx
      knt=0
      knt1=0
      if(ifanis.eq.0) npar=5
      if(ifanis.ne.0) npar=8
c
      do 10 nn=1,nreg
      read(lumod,1101) nlay,ifs,r1,r2
 1101 format(1x,i4,i5,2e15.8)
      if(nn.ne.1) r(knt)=r1*tau
      dr=(r2-r1)/float(nlay-1)
      do 2 i=1,nlay
      knt=1+knt
    2 r(knt)=r1+dr*float(i-1)
      do 1 j=1,npar
      read(lumod,3002) (coef(i),i=1,5)
 3002 format(5e16.5)
      do 3 i=1,nlay
      ind=knt1+i
      rt=r(ind)/rx
      if(j.eq.npar) r(ind)=r(ind)*tau
      val=coef(5)
      do 4 ii=1,4
    4 val=val*rt+coef(5-ii)
      if(j.eq.1) rho(ind)=val*tau
      if(j.eq.2) ccon(ind)=val*tau
      if(j.eq.3) lcon(ind)=val*tau
      if(j.eq.4) qshear(ind)=val
      if(j.eq.5) qkappa(ind)=val
      if(j.eq.6) acon(ind)=val*tau
      if(j.eq.7) ncon(ind)=val*tau
      if(j.eq.8) fcon(ind)=val
    3 continue
    1 continue
   10 knt1=knt1+nlay
      n=knt
c
  300 rn=r(n)
      if(ifanis.ne.0) goto 400
      do 25 i=1,n
      acon(i)=ccon(i)
      ncon(i)=lcon(i)
   25 fcon(i)=1.
c
  400 nsl=n
  420 if(lcon(nsl).ne.0.) goto 410
      nsl=nsl-1
      goto 420
c
  410 do 26 i=1,n
      r(i)=r(i)/rn
   26 rho(i)=rho(i)/rhon

      call rspln(1,n,r,rho,qrho,wk)
c
c  y(3) is .75*g/r
c  y(4) is int_0^r rho*g
c
      y(1)=1.
      y(2)=0.
      y(3)=rho(1)
      y(4)=0.
      rr=1.e-4
      do 100 i=2,n
c         print *, "i,n", i,n
      rstep=r(i)-r(i-1)
      if(abs(rstep).gt.1.e-4) goto 110
      ell(i)=ell(i-1)
      eta(i)=eta(i-1)
      g(i)=g(i-1)
      pres(i)=pres(i-1)
      goto 100
c
  110 do 101 ii=1,4
  101 y1(ii)=y(ii)
      rr1=rr
      if(rr+rstep.gt.r(i)) rstep=r(i)-rr
      iback=1
      goto 1100
c
 1201 do 201 ii=1,4
  201 y2(ii)=y(ii)
      rr2=rr
  901 rstep=rstep*.5
      do 321 ii=1,4
  321 y(ii)=y1(ii)
      rr=rr1
      iback=2
      goto 1100
c
 1202 do 401 ii=1,4
  401 y3(ii)=y(ii)
      rr3=rr
      iback=3
      goto 1100
c
 1203 do 501 ii=1,4
      if(abs(y(ii)-y2(ii)).gt..5e-4) goto 601
  501 continue
      goto 701
c
  601 do 801 ii=1,4
  801 y2(ii)=y3(ii)
      rr2=rr3
      goto 901
  701 if(abs(rr-r(i)).lt.1.e-5) goto 1300
      rstep=4.*rstep
      goto 110
c
 1300 ell(i)=y(1)
      eta(i)=rr*y(2)/y(1)
      g(i)=4.*rr*y(3)/3.
      pres(i)=y(4)
      goto 100
c
 1100 k=krunge(4,y,ydot,rr,rstep)
      if(k.ne.1) goto (1201,1202,1203),iback
c
      ydot(1)=y(2)
      t=rr-r(i-1)
      im1=i-1
      rhot=rho(im1)+t*(qrho(1,im1)+t*(qrho(2,im1)+t*qrho(3,im1)))
      ydot(3)=3.*(rhot-y(3))/rr
c      print *, y(3)
      ydot(2)=-2.*(ydot(3)*y(1)+3.*rhot*y(2))/(y(3)*rr)
      ydot(4)=rhot*rr*1.33333333*y(3)
      goto 1100
c
  100 continue
c
      goto 1400
c
c
 1400 factr=.75*g(n)/r(n)
      factr2=factr**2
      rhobar=rhon*factr
      wn=(pi*bigg*rhobar)**.5
      vn=rn*wn
      gn=rn*wn**2
      fac=2.5*(capom/wn)**2/(4.*ell(n)*(eta(n)+2.)/3.)
      g(1)=0.
      ell(1)=1.
      pres(1)=0.
      do 1401 i=1,n
      rho(i)=rho(i)/factr
      g(i)=g(i)/factr
      pres(i)=(pres(n)-pres(i))/factr2
      ell(i)=ell(i)*fac
      do 1401 j=1,3
 1401 qrho(j,i)=qrho(j,i)/factr
c
c
c
      do 7 i=1,n
      ccon(i)=rho(i)*(ccon(i)/vn)**2
      lcon(i)=rho(i)*(lcon(i)/vn)**2
      acon(i)=rho(i)*(acon(i)/vn)**2
      ncon(i)=rho(i)*(ncon(i)/vn)**2
    7 fcon(i)=fcon(i)*(acon(i)-2.*lcon(i))
c
c
      call rspln(1,n,r,ccon,qccon,wk)
      call rspln(1,n,r,lcon,qlcon,wk)
      call rspln(1,n,r,acon,qacon,wk)
      call rspln(1,n,r,ncon,qncon,wk)
      call rspln(1,n,r,fcon,qfcon,wk)
c
      ndisc=0.
      do 56 i=2,n
      if(abs(r(i)-r(i-1)).gt.1.e-4) goto 56
      ndisc=1+ndisc
      ndsc(ndisc)=i-1
   56 continue
      ndisc=1+ndisc
      ndsc(ndisc)=n
c
      if(lupr.eq.0) goto 99
      write(lupr,909) rhobar,1./ell(n),nic,noc,moho,nsl,n
  909 format(///' mean density =',f10.3,'     surface ellipticity ='
     1   ,' 1/',f7.3/' nic =',i4,'   noc =',i4,'   moho =',i4
     2   ,'   nsl =',i4,'   n =',i4///)
      i=0
      rlas=-10000.
      lines=9
  903 write(lupr,902)
  902 format(1x,'level',
     1 3x,'radius',7x,'rho',8x,'vpv',8x,'vph',7x,'vsv',
     2 7x,'vsh',6x,'eta',7x,'g',7x,'1/qm',5x,'1/qk'
     3  ,'      1/ell    etaell'/)
  906 i=1+i
      if(i.gt.n) goto 99
      rr=r(i)*rn
      rrho=rho(i)*rhobar
      vpv=sqrt(ccon(i)/rho(i))*vn
      vph=sqrt(acon(i)/rho(i))*vn
      vsv=sqrt(lcon(i)/rho(i))*vn
      vsh=sqrt(ncon(i)/rho(i))*vn
      eeta=fcon(i)/(acon(i)-2.*lcon(i))
      gg=g(i)*gn
      if(abs(rlas-rr).lt.500.) write(lupr,912)
  912 format(1x,132('-'))
      if(rr.eq.rlas) lines=1+lines
      rlas=rr
      write(lupr,9021) i,rr,rrho,vpv,vph,vsv,vsh,eeta,gg
     1   ,qshear(i),qkappa(i),ell(i),eta(i),pres(i)*rhobar*vn*vn
 9021 format(2x,i3,f11.1,f10.3,2f11.3,2f10.3,4f9.5,f11.7,f9.5,1pe10.3)
      lines=1+lines
      if(lines.lt.57) goto 906
      lines=1
      write(lupr,911)
  911 format(1h1)
      goto 903
   99 continue
      crn=rn
      cvn=vn
      crb=rhobar

      end
