      program compuderiv

      parameter (MMODES=500)

      character*80 getunx,mantlemodel
      character*2 ftype
      character*1 iq(MMODES)
      integer lord(MMODES),nord(MMODES),ity(MMODES),ntot(MMODES)

      real time(2),tcom

      common/modes/nord,iq,lord,ity,ntot,nsum,mtot

      call chekcl('|-lu2:o:1:[/eejit/home/talavera/dst2model/github/dta/foanis05.222]'
     1          //'|-lu7:o:1:[mdcpl.out]'
     1          //'|-lu3:o:1:[/eejit/home/talavera/dst2model/github/dta/m1084x2.htm] model on unit 3 (rdmdl)'
     1          //'|-model:o:1:[/eejit/home/talavera/dst2model/github/dta/mzero_4_l3.in] Model'
     1          //'|-pc:o:1:[default.pc] startup plotting commands'
     1          //'|')

      mantlemodel=getunx('-model',1,ll)
      write(6,'(a)') mantlemodel

cs    Number of layers
      kstart=1
      kend=3

ca      write(6,*) 'Input ldiff,mdiff'
ca      read(5,*) ldiff,mdiff

      nsum=0
      open(10,file='modes.in')
      rewind(10)
      read(10,*) mtot

      do i=1,mtot
         read(10,222) nord(i),iq(i),lord(i)
         if(iq(i).eq.'s') ity(i)=538976339
         if(iq(i).eq.'t') ity(i)=538976340
         if(iq(i).eq.'S') then
                ity(i)=538976339
                iq(i)='s'
         else if(iq(i).eq.'T') then
                ity(i)=538976340
                iq(i)='t'
         endif
         ntot(i)=(lord(i)*2)+1
         nsum=nsum+ntot(i)
         write(6,222) nord(i),iq(i),lord(i)
      enddo
  222 format(i3,1x,a1,1x,i3)
      close(10)


      tcom=dtime(time)

      do kdiff=kstart,kend ! looping through layers

        do ldiff=0,4,2 ! looping through degrees
        do mdiff=0,2*ldiff ! looping through each coeff in a degree

        write(6,*) 'parameters are ', ldiff,mdiff,kdiff

        call mdcplmdelta(ldiff,mdiff,kdiff)

        enddo
        enddo

        tcom=dtime(time)
        write(6,*) 'computation time is ',tcom

      enddo

      end

c --------------------------------------------------------------------------
      subroutine mdcplmdelta(ldiff,mdiff,kdiff)

c ML  is the maximum angular order of the modes
c MLL is the maximum degree of the heterogeneity
c MMODES is the maximum number of modes

c In this version of the program, w0=0.5*(om1+om2)
c and wsum is not substracted from the diagonal elements.

      parameter (ML=75)
      parameter (MDIM=2500)
      parameter (MMODES=500)
      complex*16 aa(MDIM,MDIM),ev(MDIM),ef(MDIM,MDIM),wk(2*MDIM),
     1           w(MDIM),z(-ML:ML,-ML:ML),A(MDIM,MDIM)
      complex    hh(MDIM*MDIM)
      complex    w0,w1,w2,dsw1,dsw2,wsum
      integer lord(MMODES),nord(MMODES),ity(MMODES),ntot(MMODES)
      real*8 wr(MDIM),wq(MDIM)

      character*1 iq(MMODES)
      character*80 file
      character*80 getunx

c common blocks for the eigenfunctions U, U', V, V',ph1,ph1'
c jcom=1 --> radial mode
c jcom=2 --> toroidal mode
c jcom=3 --> spheroidal mode
c n = radial order
c l  = angular order
c om = angular frequency in rad/s
c cgp = group velocity in km/s
c avert, ahor = vertical and horizontal accelerations at the ocean floor
c               -- used to calculate seismograms.


      common/eig1/n1,jcom1,l1,om1,q1,cgp1,avert1,ahor1,phis1
     1        ,u1(222),up1(222),v1(222),vp1(222),ph1(222),php1(222)
      common/eig2/n2,jcom2,l2,om2,q2,cgp2,avert2,ahor2,phis2
     1        ,u2(222),up2(222),v2(222),vp2(222),ph2(222),php2(222)
      common/matr/w0,z
      common/modes/nord,iq,lord,ity,ntot,nsum,mtot

c conversion factors returned by subroutine model
c rn =     length conversion to m (=6371000.00)
c vn =     velocity conversion to m/s
c wn =     frequency conversion to s**(-1)
c gn =     acceleration conversion to m/s/s
c rhobar = density conversion to kg/m**3

      common/nond/rn,wn,vn,gn,rhobar

c default models are:
c foanis05.222 = PREM
c m1084x2.htm  = M84C
c lmnthet.dta  = L02.56


      file=getunx('-lu2',1,ll)
      open(2,file=file,status='old')
      file=getunx('-lu7',1,ll)
      open(7,file=file)
      file=getunx('-lu3',1,ll)
      open(3,file=file)

c open the PREM mode catalogue

c1      write(6,*) 'openfl'
      call openfl(1,'/eejit/home/talavera/dst2model/github/dta/PREM222.BIN',1,0,0,istat,5364)

c1      write(6,*) 'seteig'
      call seteig(1)

c read in the PREM model and set it up in common blocks modl1, modl2
c dimensional conversion factors rn, wn, vn, gn, rhobar
c are returned. These are conversion factors from a dimensionless
c system of units such that
c   earth's mean density = 1, pi*G = 1, earth radius=1.
c The model is stored in dimensionless units.

c1      write(6,*) 'modl'
      call modl(2,0,rn,wn,vn,gn,rhobar)

c     I guess fus, fup, flp, fls
c     are factors associated with
c     upper mantle S and P, lower mantle S and P
c     They could be automatically set to 1.0
c
c  -> these values are not needed anymore, replaced by
c     alfau,alfal,betau,betal.
c     fus is still needed in some subroutines, I think it just
c     has to be set to 1.0
c
      fus=1.0
ca    9 write(6,"('Input alfau,alfal,betau,betal')")
ca      read(5,*) alfau,alfal,betau,betal
ca      write(6,*) 'alfau,alfal,betau,betal',alfau,alfal,betau,betal
      call intpltnnew(getunx('-model',1,ll),xper,ldiff,mdiff,kdiff)


cL    8 continue
ca      write(6,"('Input mins & maxs  eg:0 8')")
ca      read(5,*)mins, maxs
ca      write(6,*) 'Minimum and maximum heterogeity are ',mins,maxs
ca      write(6,"('do you want include rotation? yes=1 no=0')")
ca      read(5,"(i1)")isw1
ca      write(6,"('do you want include ellip.? yes=1 no=0')")
ca      read(5,"(i1)")isw2
      isw1=0
      isw2=0

c get the eigenfunctions from the mode catalogue
c and compute the coupling matrices

      do 200 im=1,mtot
c
c Self-coupling
c
         call fetch(nord(im),iq(im),lord(im),1,ifexst)
ca         write(6,*) nord(im),iq(im),lord(im),1,ifexst
         call fetch(nord(im),iq(im),lord(im),2,ifexst)
ca         write(6,*) nord(im),iq(im),lord(im),2,ifexst
         if(ifexst.ne.1) then
          write(6,"('mode does not exist in catalogue')")
          write(6,*) 'In self-coupling', nord(im),iq(im),lord(im),ifexst
          stop
         endif

         w1=cmplx(om1,.5*om1*q1)
c         w0=w1
         w0=om1

c1         write(6,*) 'Self-coupling ',im,' mode with omega = ',om1

         call memdcpl(isw1,isw2,fus,ldiff,ldiff,mdiff,kdiff)

         if(im.eq.mtot) go to 12
c
c Cross-coupling
c
         do 400 jm=im+1,mtot
            call fetch(nord(jm),iq(jm),lord(jm),2,ifexst)
c            write(6,*) nord(jm),iq(jm),lord(jm),2,ifexst
            if(ifexst.ne.1) then
             write(6,"('mode does not exist in catalogue')")
             write(6,*) 'In cross coupling', nord(jm),iq(jm),lord(jm),ifexst
             stop
            endif

            w2=cmplx(om2,.5*om2*q2)
c            w0=0.5*(w1+w2)
            w0=0.5*(om1+om2)
c            w0=sqrt(om1*om2)

c         write(6,*) 'Cross-coupling ',im,' mode with ',jm,' mode'
c         write(6,*) 'omega = ',w2

            call memdcpl(isw1,isw2,fus,ldiff,ldiff,mdiff,kdiff)

400      continue
200    continue

12     continue
       close(1)
       close(2)
       close(3)
       close(7)
      return

      end
c
c-----------------------------------------------
cc   to calculate the matrix elements for mode coupling, (k'm'1z1km)
cc         for -lord2 < m' <lord2 and -lord1 < m <lord1


      subroutine  memdcpl(isw1,isw2,fus,minll,maxll,mdiff,kdiff)

      parameter(MLL=20)
      parameter(ML=75)

      common/eig1/nord1,jcom1,lord1,wcom1,qbar1,cgp1,avert1,ahor1,phis1
     1          ,eif1(222,6)
      common/eig2/nord2,jcom2,lord2,wcom2,qbar2,cgp2,avert2,ahor2,phis2
     1          ,eif2(222,6)


      common/int/rot,asp ! asp is htg in intmdcpl
      common/matr/w0,z

      complex    iunt,w0
      complex*16 zr(-ML:ML,-ML:ML),z(-ML:ML,-ML:ML),tt(0:MLL,-MLL:MLL),t

      double precision rot,asp(0:MLL,2*MLL+1),sj(-ML:ML,-ML:ML),bfac,wig0

      data iunt/(0.,1.)/
      data capom/7.292115e-05/       ! siderial angular velocity
      data pai/3.14159 62535 9/

      w0r=cabs(w0)

c do the radial integrations

      call intmdcpl(w0r,isw2,fus,minll,maxll,kdiff) ! this gives asp and rot

      if(jcom1.eq.1)jcom1=3
      if(jcom2.eq.1)jcom2=3

ca      write(6,*) 'In intmdcpl'
ca      write(6,*) 'mdiff is ',mdiff

cc  coupling due to rotation
cc     z will be 0 unless:mord1=mord2 and
cc                         jcom1.eq.jcom2 and lord1=lord2 , or
cc                         jcom1.ne.jcom2 and lord1=lord2(+,-)1

cc  coupling due to asphericity
      open(666,file='smax.dat')
      read(666,*) smax
      close(666)
      open(25,file='A.dat',access='append')

cs    converting complex to real SH
      do 200 ll=0,4,2

      tt(ll,0)=asp(ll,1) ! asphericity
      if(lord1.eq.lord2) then
         if (ll.le.smax) then
ca      write(6,*) 'rot ',lord1,rot,wcom1,qbar1
ca      write(6,*) lord1,ll,0,1.e6*dble(tt(ll,0))/(2.0*wcom1*2.0*pai)
         if(mdiff.eq.0) write(25,*) 1.e6*dble(tt(ll,0))/(2.0*wcom1*2.0*pai)
         if(mdiff.ne.0) write(25,*) dble(0.000000000000000)
ca        if(mdiff.eq.0) write(6,*) 1.e6*dble(tt(ll,0))/(2.0*wcom1*2.0*pai)
         else
         if(mdiff.eq.0) write(25,*) dble(0.000000000000000)
         if(mdiff.ne.0) write(25,*) dble(0.000000000000000)
         endif
      endif
      do 202 mm=1,ll
      ma=2*mm
      tt(ll,mm)=0.5*(asp(ll,ma)-asp(ll,ma+1)*iunt)
      if(lord1.eq.lord2) then
       if (ll.le.smax) then
ca      write(6,*) lord1,ll,(2*mm-1),1.e6*dble(tt(ll,mm))/(2.0*wcom1*2.0*pai)
ca      write(6,*) lord1,ll,2*mm,1.e6*dimag(tt(ll,mm))/(2.0*wcom1*2.0*pai)
       mcheck=2*mm-1
       if(mdiff.eq.mcheck) write(25,*) 1.e6*dble(tt(ll,mm))/(2.0*wcom1*2.0*pai)
       if(mdiff.ne.mcheck) write(25,*) dble(0.000000000000000)
       mcheck=2*mm
       if(mdiff.eq.mcheck) write(25,*) 1.e6*dimag(tt(ll,mm))/(2.0*wcom1*2.0*pai)
       if(mdiff.ne.mcheck) write(25,*) dble(0.000000000000000)
       else
       mcheck=2*mm-1
       if(mdiff.eq.mcheck) write(25,*) dble(0.000000000000000)
       if(mdiff.ne.mcheck) write(25,*) dble(0.000000000000000)
       mcheck=2*mm
       if(mdiff.eq.mcheck) write(25,*) dble(0.000000000000000)
       if(mdiff.ne.mcheck) write(25,*) dble(0.000000000000000)
      endif
      endif
  202 continue

  200 continue

      close(25)

  400 return
      end
c-----------------------------------------------
      subroutine intmdcpl(w0,isw,fus,mins,maxs,kdiff)
      parameter(MLL=20)
      parameter(MM=2*MLL+1)

      common/int/rot,htg
      common/grav/grv,qgrv

      real lcon,ncon,llp,nnp,ll,nn,kappa,mu
      common/modl1/n,nic,noc,moho,nsl,ifanis,r(222)
     1            ,rho(222),qrho(3,222),g(222),ell(222),eta(222)
      common/modl2/acon(222),qacon(3,222),ccon(222),qccon(3,222)
     1            ,lcon(222),qlcon(3,222),ncon(222),qncon(3,222)
     2            ,fcon(222),qfcon(3,222)
      common/modl3/qshear(222),qkappa(222)
      common/nond/rn,wn,vn,gn,rhobar
      common/eig1/nord1,jcom1,lord1,wcom1,qbar1,cgp1,avert1,ahor1,phis1
     1          ,eif1(222,6)
      common/eig2/nord2,jcom2,lord2,wcom2,qbar2,cgp2,avert2,ahor2,phis2
     1          ,eif2(222,6)

c  common block for harmonic coefficients of the model
c  perturbations -- perturbations in A, C, L, N, F, Rho, Phi0, Phi0'
c  In this program Rho, Phi0, Phi0' are not perturbed -- upgrade in the future.

      common/vioem/da(0:MLL,MM,222,5),dc(0:MLL,MM,222,5),df(0:MLL,MM,222,
     1          5),dl(0:MLL,MM,222,5),dn(0:MLL,MM,222,5),dr(0:MLL,MM,222,
     2          5),dp(0:MLL,MM,222,5),dpp(0:MLL,MM,222,5),
     3             rhoin(222,5),grin(222,5),elin(222,5),etin(222,5)
      common/hetmdl/iw(5),crt(81)

      dimension u1(222),up1(222),v1(222),vp1(222),p1(222),pp1(222),
     1          u2(222),up2(222),v2(222),vp2(222),p2(222),pp2(222)
      dimension dta(0:MLL,MM),dtc(0:MLL,MM),dtf(0:MLL,MM),dtl(0:MLL,MM),
     1          dtn(0:MLL,MM),dtr(0:MLL,MM),dtp(0:MLL,MM),dtpp(0:MLL,MM),
     2          h(0:MLL,MM)
      dimension q1(3),qp1(3),q2(3),qp2(3)
      dimension xi(5),wt(5),a1(3),b1(3),a2(3),b2(3)

      equivalence (eif1(1,1),u1(1)),(eif1(1,2),up1(1)),(eif1(1,3),v1(1)),
     1            (eif1(1,4),vp1(1)),(eif1(1,5),p1(1)),(eif1(1,6),pp1(1)),
     2            (eif2(1,1),u2(1)),(eif2(1,2),up2(1)),(eif2(1,3),v2(1)),
     3            (eif2(1,4),vp2(1)),(eif2(1,5),p2(1)),(eif2(1,6),pp2(1))
      equivalence (q1(1),uri),(q1(2),vri),(q1(3),pri),
     1            (qp1(1),upri),(qp1(2),vpri),(qp1(3),ppri),
     2            (q2(1),urii),(q2(2),vrii),(q2(3),prii),
     3            (qp2(1),uprii),(qp2(2),vprii),(qp2(3),pprii)

      double precision temp(2,7),acum(2,7),add,rot,b(5),
     1                 asp(0:MLL,MM,7),ast(0:MLL,MM,7),htg(0:MLL,MM)
      real grv(222,0:MLL),qgrv(3,222,0:MLL),gint(0:MLL)

c Gauss-Legendre absicae and weights from Abromowitz & Stegun

      data xi/-0.90617 98459 38664,-0.53846 93101 05683
     1       , 0.00000 00000 00000, 0.53846 93101 05683
     2       , 0.90617 98459 38664/
      data wt/ .23692 68850 56189, 0.47862 86704 99366
     1       , .56888 88888 88889, 0.47862 86704 99366
     2       , .23692 68850 56189/


      data n670/180/,n220/202/,nmidc/216/
      data capom/7.292115e-05/
      data coe/1.05688 76266 5/
c             coe=(sqpai/sqrt(5.))*4/3.
      data nint/2/
      data drn,drnsl,drmidc,drmoho
     1   ,droc,druc,drlc
     2   ,dvpoc,dvpuc,dvplc
     3   ,dvsoc,dvsuc,dvslc/
     4     3.5028e-02, 2.1017e-01,-3.7886e-01,-7.8982e-01
     5   ,-3.5028e-04, 0.        , 0.
     6   , 0.        ,-2.1017e-02, 3.5028e-03
     7   , 0.        ,-5.2542e-03, 1.7514e-03/


      call gravk(mins,maxs)

      if(jcom1.eq.0.or.jcom2.eq.0) return

      if(jcom1.eq.1)jcom1=3
      if(jcom2.eq.1)jcom2=3
      jcom=3
      if(jcom1.eq.2.and.jcom2.eq.2)jcom=2

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

      do 101 i=1,nint
      do 101 j=1,7
  101 acum(i,j)=0.

      do 102 i=0,MLL
      do 102 j=1,2*i+1
      do 102 k=1,7
 102     asp(i,j,k)=0.

      do 103 i=0,MLL
      do 103 j=1,2*i+1
 103     htg(i,j)=0.


      nstart=1
      if(jcom1.eq.2.or.jcom2.eq.2) nstart=noc

      do 1000 iq=nstart,n ! loop over PREM layers

      do 1001 i=1,nint
      do 1001 j=1,7
 1001 temp(i,j)=0.

      do 1002 i=mins,maxs
      do 1002 j=1,2*i+1
      do 1002 k=1,7
 1002 ast(i,j,k)=0.

      if(iq.eq.n) goto 1200

      call corfac(iq,wcom,jcom,xac,xf,xln) ! calculating fac. for poly. extrap.

      iq1=iq+1
      r1=r(iq)
      r2=r(iq1)
      hn=r2-r1
      hnh=hn*.5
      if(hn.lt.1.e-4) goto 1200
      hr=1./hn
      hsq=hr*hr
      hcu=hr*hsq

c      do i=1,3
c         a1(i)=0.
c         b1(i)=0.
c         a2(i)=0.
c         b2(i)=0.
c      enddo

      do 1003 i=j11,j12
      i1=2*i-1
      i2=i1+1
      a1(i)=(eif1(iq,i2)+eif1(iq1,i2))*hsq
     1     +2.*(eif1(iq,i1)-eif1(iq1,i1))*hcu
 1003 b1(i)=-(2.*eif1(iq,i2)+eif1(iq1,i2))*hr
     1     -3.*(eif1(iq,i1)-eif1(iq1,i1))*hsq

      do 1004 i=j21,j22
      i1=2*i-1
      i2=i1+1
      a2(i)=(eif2(iq,i2)+eif2(iq1,i2))*hsq
     1     +2.*(eif2(iq,i1)-eif2(iq1,i1))*hcu
 1004 b2(i)=-(2.*eif2(iq,i2)+eif2(iq1,i2))*hr
     1     -3.*(eif2(iq,i1)-eif2(iq1,i1))*hsq


      do 1100 il=1,5 ! loop through sub layers

      t=.5*hn*(xi(il)+1.)
      rr=r1+t
      rr2=rr*rr

c      do i=1,3
c         q1(i)=0.
c         qp1(i)=0.
c         q2(i)=0.
c         qp2(i)=0.
c      enddo

      do 1101 i=j11,j12
      i1=2*i-1
      i2=i1+1
      q1(i)=eif1(iq,i1)+t*(eif1(iq,i2)+t*(b1(i)+t*a1(i)))
 1101 qp1(i)=(eif1(iq,i2)+t*(2.*b1(i)+t*3.*a1(i)))*rr

      do 1102 i=j21,j22
      i1=2*i-1
      i2=i1+1
      q2(i)=eif2(iq,i1)+t*(eif2(iq,i2)+t*(b2(i)+t*a2(i)))
 1102 qp2(i)=(eif2(iq,i2)+t*(2.*b2(i)+t*3.*a2(i)))*rr

ca      write(6,*) 'eig func are ',iq,eif1(iq,1),uri,vri,upri,vpri

      el=elin(iq,il)
      et=etin(iq,il)
      gr=grin(iq,il)

      if(jcom1.eq.2.or.jcom2.eq.2) goto 1103 ! S mode
      aa=xac*(acon(iq)+t*(qacon(1,iq)+t*(qacon(2,iq)+t*qacon(3,iq))))
      cc=xac*(ccon(iq)+t*(qccon(1,iq)+t*(qccon(2,iq)+t*qccon(3,iq))))
      ff=xf*(fcon(iq)+t*(qfcon(1,iq)+t*(qfcon(2,iq)+t*qfcon(3,iq))))
      aap=xac*(qacon(1,iq)+t*(2.*qacon(2,iq)+t*3.*qacon(3,iq)))
      ccp=xac*(qccon(1,iq)+t*(2.*qccon(2,iq)+t*3.*qccon(3,iq)))
      ffp=xf*(qfcon(1,iq)+t*(2.*qfcon(2,iq)+t*3.*qfcon(3,iq)))
 1103 ll=xln*(lcon(iq)+t*(qlcon(1,iq)+t*(qlcon(2,iq)+t*qlcon(3,iq)))) ! T mode
      nn=xln*(ncon(iq)+t*(qncon(1,iq)+t*(qncon(2,iq)+t*qncon(3,iq))))
      rrho=rho(iq)+t*(qrho(1,iq)+t*(qrho(2,iq)+t*qrho(3,iq)))
      llp=xln*(qlcon(1,iq)+t*(2.*qlcon(2,iq)+t*3.*qlcon(3,iq)))
      nnp=xln*(qncon(1,iq)+t*(2.*qncon(2,iq)+t*3.*qncon(3,iq)))
      rhop=qrho(1,iq)+t*(2.*qrho(2,iq)+t*3.*qrho(3,iq))
      etan=ff/(aa-2.*ll)
      kappa=(cc+4.*aa-4.*nn+4.*ff)/9.
      mu=(cc+aa+6.*ll+5.*nn-2.*ff)/15.
      qmu=qshear(iq)
      qkap=qkappa(iq)

cs    No need to multiply by qmu, the 3D qmu model read in intpltnnew is in dqmu
      do 1104 i=mins,maxs ! loop through degrees
      do 1104 j=1,2*i+1 ! loop though coeffs in degree
      dta(i,j)=mu*da(i,j,iq,il)
      dtc(i,j)=mu*dc(i,j,iq,il)
      dtl(i,j)=mu*dl(i,j,iq,il)
      dtn(i,j)=mu*dn(i,j,iq,il)
      dtf(i,j)=mu*df(i,j,iq,il)
      dtr(i,j)=0.!dr(i,j,iq,il)
      dtp(i,j)=0.!dp(i,j,iq,il)
ca      write(6,*) i,j,da(i,j,iq,il)

 1104 dtpp(i,j)=0.!dpp(i,j,iq,il)

 6000 if(isw.eq.0)goto 5001
      dta(2,1)=dta(2,1)+coe*rr*el*aap
      dtc(2,1)=dtc(2,1)+coe*rr*el*ccp
      dtf(2,1)=dtf(2,1)+coe*rr*el*ffp
      dtl(2,1)=dtl(2,1)+coe*rr*el*llp
      dtn(2,1)=dtn(2,1)+coe*rr*el*nnp
cs      dtr(2,1)=dtr(2,1)+coe*rr*el*rhop
cs      dtp(2,1)=dtp(2,1)+coe*0.5*capon2*rr2
cs      dtpp(2,1)=dtpp(2,1)+coe*capon2*rr


 5001 continue

      do i=mins,maxs ! gravity
      gint(i)=0.!grv(iq,i)+t*(qgrv(1,iq,i)+t*(qgrv(2,iq,i)+t*qgrv(3,iq,i)))
      enddo

      xxi=vpri-vri
      xxii=vprii-vrii
      vr2=vri*vrii

      if(jcom1.eq.2) goto 1105
      fi=2.*uri-fl3i*vri
      xxi=xxi+uri
 1105 if(jcom2.eq.2) goto 1110
      fii=2.*urii-fl3ii*vrii
      xxii=xxii+urii

 1110 xx2=xxi*xxii
      if(jcom.eq.2)goto 1106
      uivii=uri*vrii
      uiivi=urii*vri
      fivii=fi*vrii
      fiivi=fii*vri
      upivii=upri*vrii
      upiivi=uprii*vri
      vpiuii=vpri*urii
      vpiiui=vprii*uri

 1106 if(jcom1.ne.jcom2)goto 1120

      ccc=vr2
      rkl1=xx2
      rkn2=vr2
      rkr1=-w0n2*vr2*rr2
c ---- I think the following two values have to be set to zero here -----
      rkp1=0.
      rkpp1=0.

      if(jcom1.eq.2) goto 1150

      ur2=uri*urii
      fiuii=fi*urii
      fiiui=fii*uri

      ccc=vr2+uivii+uiivi
      rka0=fi*fii
      rkc0=upri*uprii
      rkf0=uprii*fi+upri*fii
      rkn0=-rka0
      rkr0=((8.*rrho-w0n2)*ur2*rr+ppri*urii+pprii*uri
     1     -.5*gr*(4.*ur2+fiuii+fiiui))*rr
      rkr1=rkr1+(pri*vrii+prii*vri+.5*gr*(uivii+uiivi))*rr
      rkp0a=rrho*ur2
      rkp1a=.5*rrho*(vpiiui+uivii-upivii-2.*fivii)
      rkp1b=.5*rrho*(vpiuii+uiivi-upiivi-2.*fiivi)
      rkpp0=-rrho*(fiuii+fiiui)*rr
      rkpp1a=.5*rrho*rr*uivii
      rkpp1b=.5*rrho*rr*uiivi

      goto 1150

 1120 rkl1=xx2
      rkn2=vr2
      rkr1=-w0n2*vr2*rr2
      if(jcom1.eq.2)then
        ccc=.5*((fl3i+fl3ii-2.)*vr2-(fl3i-fl3ii+2.)*uiivi)
        rkr1=rkr1+(prii+.5*gr*urii)*vri*rr
        rkp1=-.5*rrho*(vpiuii+uiivi-upiivi-2.*fiivi)
        rkpp1=-.5*rrho*uiivi*rr
      else
        ccc=-.5*((fl3ii+fl3i-2.)*vr2-(fl3ii-fl3i+2.)*uivii)
        rkl1=-rkl1
        rkn2=-rkn2
        rkr1=-rkr1-(pri+.5*gr*uri)*vrii*rr
        rkp1=.5*rrho*(vpiiui+uivii-upivii-2.*fivii)
        rkpp1=.5*rrho*uivii*rr
      endif

 1150 continue

c
c   normalization integral
c
      if(jcom1.ne.jcom2) goto  1160
      if(lord1.ne.lord2) goto  1160
      t1=vr2*fl3i
      if(jcom1.ne.2) t1=t1+ur2
      add=wt(il)*rr2*t1*rrho*hnh
      temp(1,1)=temp(1,1)+add
 1160 continue

c
c   rotational splitting
c
      add=wt(il)*rr2*ccc*rrho*hnh
      temp(1,2)=temp(1,2)+add

c
c   aspherical perturbation (including ellipticity)
c

      do 1170 i=mins,maxs
      do 1170 j=1,2*i+1
      add=wt(il)*(rkl1*dtl(i,j)+rkr1*dtr(i,j)+rkp1*dtp(i,j)+rkpp1*dtpp(
     1    i,j))*hnh
      ast(i,j,2)=ast(i,j,2)+add
      if(jcom1.ne.jcom2) then
ca        if(i.eq.2.and.j.eq.1) write(6,*) rkl1,rkr1,rkp1,rkpp1
      endif
      add=wt(il)*rkn2*dtn(i,j)*hnh
      ast(i,j,3)=ast(i,j,3)+add
      add=wt(il)*gint(i)*dtr(i,j)*rr2*hnh
      ast(i,j,7)=ast(i,j,7)+add

      if(jcom1.eq.2.or.jcom2.eq.2)goto 1170
      add=wt(il)*(rka0*dta(i,j)+rkc0*dtc(i,j)+rkf0*dtf(i,j)+rkn0*dtn(i,j
     1         )+rkr0*dtr(i,j)+rkpp0*dtpp(i,j))*hnh
      ast(i,j,1)=ast(i,j,1)+add
      add=wt(il)*rkp0a*dtp(i,j)*hnh
      ast(i,j,4)=ast(i,j,4)+add
      add=wt(il)*(rkp1a*dtp(i,j)+rkpp1a*dtpp(i,j))*hnh
      ast(i,j,5)=ast(i,j,5)+add
      add=wt(il)*(rkp1b*dtp(i,j)+rkpp1b*dtpp(i,j))*hnh
      ast(i,j,6)=ast(i,j,6)+add
 1170 continue

c
c    crustal correction
c
 5003 if(iq.lt.moho) goto 1184
      vs=sqrt(ll/rrho)
      vvp=sqrt(aa/rrho)
      if(iq.gt.nmidc) goto 1181
      ddrho=drlc
      ddvs2=2.*vs*dvslc*delv
      ddvp2=2.*vvp*dvplc*delv
      goto 1183
 1181 if(iq.gt.nsl) goto 1182
      ddrho=druc
      ddvs2=2.*vs*dvsuc*delv
      ddvp2=2.*vvp*dvpuc*delv
      goto 1183
 1182 ddrho=droc
      ddvs2=2.*vs*dvsoc*delv
      ddvp2=2.*vvp*dvpoc*delv
 1183 add=wt(il)*hnh*(rrho*ddvs2*rkl1+(ll*rkl1+rrho*rkr1)*ddrho/rrho)
      temp(2,2)=temp(2,2)+add
ca      if(jcom1.ne.jcom2) write(6,*) 'a',iq,temp(2,2),add
      add=wt(il)*hnh*(rrho*ddvs2+nn*ddrho/rrho)*rkn2
      temp(2,3)=temp(2,3)+add
      if(jcom1.eq.2.or.jcom2.eq.2)goto 1184
      etan=ff/(aa-2.*ll)
      add=wt(il)*hnh*(rrho*(ddvp2*(rka0+rkc0+etan*rkf0)+ddvs2*(rkn0-2.*
     1 etan*rkf0))+(aa*rka0+cc*rkc0+ff*rkf0+nn*rkn0+rrho*rkr0)*ddrho/rrho)
      temp(2,1)=temp(2,1)+add
 1184 continue

 1100 continue

      goto  1300


c
c  discontinuity contribution
c

 1200 rr=r(iq)
      rr2=rr*rr

      do 1299 idis=1,3,2

      iqt=iq+(idis-1)/2
      if(iqt.gt.n) goto 1299

      call corfac(iqt,wcom,jcom,xac,xf,xln)

c      do i=1,3
c         q1(i)=0.
c         qp1(i)=0.
c         q2(i)=0.
c         qp2(i)=0.
c      enddo

      do 1201 i=j11,j12
      i1=2*i-1
      i2=i1+1
      q1(i)=eif1(iqt,i1)
 1201 qp1(i)=eif1(iqt,i2)*rr

      do 1202 i=j21,j22
      i1=2*i-1
      i2=i1+1
      q2(i)=eif2(iqt,i1)
 1202 qp2(i)=eif2(iqt,i2)*rr

      el=ell(iqt)

      if(jcom.eq.2) goto 1203
      if(jcom1.ne.2.and.jcom2.ne.2)  aa=xac*acon(iqt)
      gr=g(iqt)
      cc=xac*ccon(iqt)
      ff=xf*fcon(iqt)
 1203 ll=xln*lcon(iqt)
      nn=xln*ncon(iqt)
      rrho=rho(iqt)

      do i=mins,maxs
      gint(i)=grv(iqt,i)
      enddo

      vr2=vri*vrii
ca      xxi=vpri-vri+uri
ca      xxii=vprii-vrii+urii
ca      xx2=xxi*xxii

      xxi=vpri-vri
      xxii=vprii-vrii

      if(jcom1.eq.2) goto 4405
      xxi=xxi+uri
      fi=2.*uri-fl3i*vri
 4405 if(jcom2.eq.2) goto 4440
      xxii=xxii+urii
      fii=2.*urii-fl3ii*vrii

4440  xx2=xxi*xxii

      if(jcom1.ne.jcom2)goto 1220

      rkl1=xx2-vpri*xxii-vprii*xxi
      rkn2=vr2
      rkr1=-w0n2*vr2*rr2

      if(jcom1.eq.2)goto 1230

      fi=2.*uri-fl3i*vri
      fii=2.*urii-fl3ii*vrii
      ur2=uri*urii

      rka0=fi*fii
      rkc0=-upri*uprii
      rkc11=upri*vrii
      rkc12=uprii*vri
      rkf11=fi*vrii
      rkf12=fii*vri
      rkn0=-rka0
      rkr0=((8.*rrho-w0n2)*ur2*rr+(ppri*urii+pprii*uri)
     1     -.5*gr*(4.*ur2+fi*urii+fii*uri))*rr
      rkr1=rkr1+((pri*vrii+prii*vri)+.5*gr*(uri*vrii+urii*vri))*rr

      goto 1230

 1220 rkr1=-w0n2*vr2*rr2
      if(jcom1.eq.2)then
        rkc1=-uprii*vri
        rkf1=-fii*vri
        rkl1=xx2-vprii*xxi-vpri*xxii
        rkn2=vr2
        rkr1=rkr1+(prii+.5*gr*urii)*vri*rr
      else
        rkc1=upri*vrii
        rkf1=fi*vrii
        rkl1=-xx2+vpri*xxii+xxi*vprii
ca        write(6,*) 'rkl1 in if-statment',iq,rkl1
        rkn2=-vr2
        rkr1=-rkr1-(pri+.5*gr*uri)*vrii*rr
      endif


c
c    aspherical perturbation(only ellipticity)
c


 1230 if(isw.eq.0)goto 1240
c 1230 continue
cl -------------------------
cl Discontinuity topography
cl -------------------------

      do i=mins,maxs
         do j=1,2*i+1
            h(i,j)=0.
         enddo
      enddo

      h(2,1)=-coe*rr*el

      do 1239 i=mins,maxs
      do 1239 j=1,2*i+1

      ttt=float(idis-2)*h(i,j)

      add=-ttt*(ll*rkl1+rrho*rkr1)
      ast(i,j,2)=ast(i,j,2)+add
ca      if(i.eq.2.and.j.eq.1) write(6,*) 'ast(2,1,2) = ',ast(i,j,2)
      add=-ttt*nn*rkn2
      ast(i,j,3)=ast(i,j,3)+add
      add=-ttt*gint(i)*rrho*rr2
      ast(i,j,7)=ast(i,j,7)+add

      if(jcom1.eq.2.or.jcom2.eq.2)goto 1231
      add=-ttt*(aa*rka0+cc*rkc0+nn*rkn0+rrho*rkr0)
      ast(i,j,1)=ast(i,j,1)+add
      add=-ttt*(cc*rkc11+ff*rkf11)
      ast(i,j,5)=ast(i,j,5)+add
      add=-ttt*(cc*rkc12+ff*rkf12)
      ast(i,j,6)=ast(i,j,6)+add

 1231 if(jcom1.eq.jcom2)goto 1239
      add=-ttt*(cc*rkc1+ff*rkf1)
      ast(i,j,2)=ast(i,j,2)+add
ca      write(6,*) 'a', ast(i,j,2), add

 1239 continue

c
c  crustal thickness
c
 1240 if(iq.ne.nmidc.and.iq.ne.moho.and.iq.ne.nsl.and.iq.ne.n) goto 1249
      if(iq.eq.n) delh=drn*deldis
      if(iq.eq.nsl) delh=drnsl*deldis
      if(iq.eq.nmidc) delh=drmidc*deldis
      if(iq.eq.moho) delh=drmoho*deldis
      ttt=float(idis-2)*delh

      add=-ttt*(ll*rkl1+rrho*rkr1)
      temp(2,2)=temp(2,2)+add
ca      if(jcom1.ne.jcom2) write(6,*) 'b',iq,temp(2,2),add

      add=-ttt*nn*rkn2
      temp(2,3)=temp(2,3)+add

      if(jcom1.eq.2.or.jcom2.eq.2)goto 1241
      add=-ttt*(aa*rka0+cc*rkc0+nn*rkn0+rrho*rkr0)
      temp(2,1)=temp(2,1)+add
      add=-ttt*(cc*rkc11+ff*rkf11)
      temp(2,5)=temp(2,5)+add
      add=-ttt*(cc*rkc12+ff*rkf12)
      temp(2,6)=temp(2,6)+add

 1241 if(jcom1.eq.jcom2)goto 1249
      add=-ttt*(cc*rkc1+ff*rkf1)
      temp(2,2)=temp(2,2)+add
ca      write(6,*) 'c',iq,temp(2,2),add
 1249 continue

 1299 continue

 1300 do 1301 i=1,nint
      do 1301 j=1,7
 1301 acum(i,j)=acum(i,j)+temp(i,j)

      do 1302 i=mins,maxs
      do 1302 j=1,2*i+1
      do 1302 k=1,7
 1302 asp(i,j,k)=asp(i,j,k)+ast(i,j,k)


 1000 continue

ca      write(6,*) 'Write out acum'
ca      do 1551 i=1,nint
ca      do 1551 j=1,6
ca 1551 write(6,"(2e12.4)") acum(i,j)

ca      write(6,*) '   '
ca      write(6,*) 'Write out asp'
ca      do 1552 i=mins,maxs
ca      do 1552 j=1,2*i+1
ca      do 1552 k=1,6
ca 1552 write(6,"(2e12.4)") asp(i,j,k)

      anorm=acum(1,1)*omn2
      rot=acum(1,2)*omn2

      do 2000 i=mins,maxs

      call bcoff1a(lord1,i,lord2,b)

      if(jcom1.ne.jcom2)goto 2020
      fl3=float(i*(i+1))
      do 2010 j=1,2*i+1 ! htg is the asphericity in common/int
      htg(i,j)=((asp(i,j,1)+asp(i,j,4)*fl3+.5*(asp(i,j,2)*(fl3i+fl3ii-
     1   fl3)+asp(i,j,5)*(fl3ii+fl3-fl3i)+asp(i,j,6)*(fl3i+fl3-fl3ii))
     2            )*b(1)+asp(i,j,3)*b(3)+asp(i,j,7)*b(1))*omn2*wn**2
 2010 continue
      goto 2000

 2020 do 2021 j=1,2*i+1
 2021 htg(i,j)=(asp(i,j,2)*b(4)+asp(i,j,3)*b(5)+asp(i,j,7)*b(4))*omn2*wn**2

 2000 continue
      return
      end

c-----------------------------------------------
      subroutine bcoff1a(l1,is,l2,b)
      implicit double precision (a-h,o-z)
      double precision b
      dimension b(5)
      dfloat(n)=n
      ls=l1+l2+is
      l1c=l1*(l1+1)
      l2c=l2*(l2+1)
      isc=is*(is+1)
      if((ls/2)*2.ne.ls) goto 10
ca      b(1)=wig0(l1,is,l2)
ca      write(6,*) 'bfac is ',l1,is,l2,b(1)
      b(1)=1.d0
      ld=l1c+l2c-isc
      b(2)=dfloat(ld)*b(1)*0.5d0
      b(3)=0.5d0*dfloat((ld-2)*ld-2*l1c*l2c)*b(1)
      b(4)=0.d0
      b(5)=0.d0
      return
10    b(1)=0.d0
      b(2)=0.d0
      b(3)=0.d0
      b(4)=wig0(l1+1,l2+1,is+1)
      if(is.lt.iabs(l1-l2).or.is.gt.(l1+l2))then
        t1=0.
        goto 11
      endif
      t1=dfloat((ls+2)*(ls+4))*dfloat((ls+1-2*l2)*(ls+1-2*l1))
     1     *dfloat(ls+1-2*is)/dfloat(ls+3)
11    b(4)=b(4)*0.5d0*dsqrt(t1)
      b(5)=dfloat(l1c+l2c-isc-2)*b(4)
      return
      end

c-----------------------------------------------
c This routine evaluates the harmonic coefficients
c of delta(A), delta(C) etc. using the input heterogeneous
c models. It will need to be changed to accomodate different
c model parameterisations etc.
c

      subroutine intpltnnew(hetmodel,xper,ldiff,mdiffo,kdiff)

      parameter(MLL=20)
      parameter(MM=2*MLL+1)
      parameter (MXPARM=30)
      parameter (MXLENY=(MLL+1)**2)
      parameter (MXMDLL=MXLENY*MXPARM)

      real lcon,ncon,ll,nn
      double precision rknots(3),vals(MLL,MM,21),rho_x
      common/modl1/n,nic,noc,moho,nsl,ifanis,r(222)
     1            ,rho(222),qrho(3,222),g(222),ell(222),eta(222)
      common/modl2/acon(222),qacon(3,222),ccon(222),qccon(3,222)
     1            ,lcon(222),qlcon(3,222),ncon(222),qncon(3,222)
     2            ,fcon(222),qfcon(3,222)
      common/nond/rn,wn,vn,gn,rhobar
      common/vioem/da(0:MLL,MM,222,5),dc(0:MLL,MM,222,5),df(0:MLL,MM,222,
     1          5),dl(0:MLL,MM,222,5),dn(0:MLL,MM,222,5),dr(0:MLL,MM,222,
     2          5),dp(0:MLL,MM,222,5),dpp(0:MLL,MM,222,5),
     3             rhoin(222,5),grin(222,5),elin(222,5),etin(222,5)
      common/lmmdl/iz(4),p(588),b1(588)   ! L02.56
      common/hetmdl/iw(5),s(588),b2(588)  ! M84C

      character*80 file

      dimension dvs2(0:MLL,MM,4),dvp(0:MLL,MM,5),xi(5)

      data xi/-0.90617 98459 38664,-0.53846 93101 05683
     1       , 0.00000 00000 00000, 0.53846 93101 05683
     2       , 0.90617 98459 38664/
      data thrd,fot/.33333 33333 33333,1.33333 33333 33333/
      data n670/180/,n220/202/,nmidc/216/

      character*(*) hetmodel
      character*2 ftype
      integer smax,nspl

csc Legendre polynomials
cs
cs      p1(x)=sqrt(1.5)*x
cs      p2(x)=sqrt(5./8.)*(3.*x*x-1.)
cs      p3(x)=sqrt(7./8.)*(5.*x**3-3.*x)
cs      p4(x)=sqrt(9./128.)*(35.*x**4-3.*x**2+3.)
cs      p0=sqrt(.5)

ca      write(6,*) 'Input l(0-20), m(0-40), k(4-24)'
c1      read(5,*) ldiff,99623293046617489,kdiff
      mdiff=mdiffo+1
ca      write(6,*) 'Input xper'
ca      read(5,*) xper
      xper=1.0  ! pert. to get the model parameters derivatives (linear inv.)

cs      rknots(7)=24.0
cs      rknots(6)=502.0
cs      rknots(5)=980.0
cs      rknots(4)=1458.0
cs      rknots(3)=1936.0
cs      rknots(2)=2414.0
cs      rknots(1)=2892.0

      rknots(3)=24.0
      rknots(2)=1458.0
      rknots(1)=2892.0
      rknots(:) = (6371.0 - rknots(:)) / 6371.0 ! normalize
cs      open(13,file='bsplines-test.dat',access='append')

      open(27,file=hetmodel)
      rewind(27)
      read(27,*) smax, nspl
cs      do k=1,nspl
cs        do i=0,smax
csc          do j=1,2*i+1
cs             read(27,*) (vals(i,j,k),j=1,2*i+1)
csc          enddo
cs        enddo
cs        if(smax.lt.MLL) then
cs          do i=smax+1,MLL
cs            do j=1,2*i+1
cs              vals(i,j,k)=0.0
cs            enddo
cs          enddo
cs        endif
cs      enddo
      close(27)

      fac=1000./vn

ca      write(6,*) 'model parameter is ',dvsr(ldiff,mdiff,kdiff)

      do 1000 iq=1,n  ! do loop over PREM levels iq=1 (centre) to iq=n (surface)
      iq1=iq+1
      r1=r(iq)
      r2=r(iq1)
      hn=r2-r1
      hnh=hn*.5
      if(hn.lt.1.e-4)goto 1000
      hr=1./hn
      hsq=hr*hr
      hcu=hr*hsq

      elld=0
      if(iq.ne.1) elld=ell(iq)*eta(iq)/r(iq)
      elld1=ell(iq1)*eta(iq1)/r(iq1)
      ae=(elld+elld1)*hsq+2.*(ell(iq)-ell(iq1))*hcu
      be=-(2.*elld+elld1)*hr-3.*(ell(iq)-ell(iq1))*hsq

      gd=fot*rho(iq)
      if(iq.ne.1) gd=4.*rho(iq)-2.*g(iq)/r(iq)
      gd1=4.*rho(iq1)-2.*g(iq1)/r(iq1)
      ag=(gd+gd1)*hsq+2.*(g(iq)-g(iq1))*hcu
      bg=-(2.*gd+gd1)*hr-3.*(g(iq)-g(iq1))*hsq


      do 1100 il=1,5

      t=.5*hn*(xi(il)+1.)
      rr=r1+t

      gr=g(iq)+t*(gd+t*(bg+t*ag))
      el=ell(iq)+t*(elld+t*(be+ae*t))
      et=(elld+t*(2.*be+t*3.*ae))*rr/el

      grin(iq,il)=gr
      elin(iq,il)=el
      etin(iq,il)=et

      aa=acon(iq)+t*(qacon(1,iq)+t*(qacon(2,iq)+t*qacon(3,iq)))
      cc=ccon(iq)+t*(qccon(1,iq)+t*(qccon(2,iq)+t*qccon(3,iq)))
      ff=fcon(iq)+t*(qfcon(1,iq)+t*(qfcon(2,iq)+t*qfcon(3,iq)))
      ll=lcon(iq)+t*(qlcon(1,iq)+t*(qlcon(2,iq)+t*qlcon(3,iq)))
      nn=ncon(iq)+t*(qncon(1,iq)+t*(qncon(2,iq)+t*qncon(3,iq)))
      rrho=rho(iq)+t*(qrho(1,iq)+t*(qrho(2,iq)+t*qrho(3,iq)))
      rhoin(iq,il)=rrho

      etan=ff/(aa-2.*ll)

      if(iq.le.nic) then ! innner core
        do i=0,MLL
          do j=1,2*i+1
            da(i,j,iq,il)=0.
            dc(i,j,iq,il)=0.
            dl(i,j,iq,il)=0.
            dn(i,j,iq,il)=0.
            df(i,j,iq,il)=0.
            dr(i,j,iq,il)=0.
          enddo
        enddo
      else if(iq.le.noc) then ! outer core
        do i=0,MLL
          do j=1,2*i+1
            da(i,j,iq,il)=0.
            dc(i,j,iq,il)=0.
            dl(i,j,iq,il)=0.
            dn(i,j,iq,il)=0.
            df(i,j,iq,il)=0.
            dr(i,j,iq,il)=0.
          enddo
        enddo
      else if(iq.le.moho) then ! mantle
c
c
c        f1=2.*cc*alfau*sqrt(rrho/ll)
c        f2=2.*sqrt(ll*rrho)
c        f3=betau*rrho*sqrt(rrho/ll)
c        rx=(2.*rr-r(noc)-r(moho))/(r(moho)-r(noc))  ! -1 at cmb 1 at moho
        do i=0,MLL
          do j=1,2*i+1
            dqmurel=0.
cs          Gives a value between 0-1 for the value of a particular spline
cs          at a particular depth
            call spl(kdiff-1,nspl,rknots(1:nspl),dble(rr),rho_x)
cs          normalization for 3 bsplines w/ knots=24,1458,2892
cs            rho_x=rho_x/7.6096669308370943
            dqmurel=xper*rho_x
c           figure out adjustments to A,C,L,N,F
c           assuming (because it's built into splh format that
c           dvp/vp = .5*dvs/vs
c           drho/rho = 0
c
c            dvprel=.5*dvsrel   ! This may change
c            drrel=.3*dvsrel          ! This may change

            da(i,j,iq,il)=0.
            dc(i,j,iq,il)=0.
            dn(i,j,iq,il)=0.
            dl(i,j,iq,il)=0.
            df(i,j,iq,il)=0.
            dr(i,j,iq,il)=0.

            if(i.eq.ldiff.and.j.eq.mdiff) then
ca              write(6,*) diffrel
ca              write(6,*) 'model parameter is ',dvsr(ldiff,mdiff,kdiff+4)
            da(i,j,iq,il)=(4./3.)*dqmurel ! + kappa
            dc(i,j,iq,il)=(4./3.)*dqmurel ! + kappa
            dn(i,j,iq,il)=dqmurel
            dl(i,j,iq,il)=dqmurel
            df(i,j,iq,il)=(-2./3.)*dqmurel ! + kappa
            dr(i,j,iq,il)=0.0 ! rrho*drrel
cs          i = degree, j = each coeff of degree, iq PREM layer, il sub layer
            endif

c           t=(p0*dvs2(i,j,1)+p1(rx)*dvs2(i,j,2)+p2(rx)*dvs2(i,j,3)+
c    1         p3(rx)*dvs2(i,j,4))*0.5*sqrt(rrho/ll)
c           dl(i,j,iq,il)=t*f2
c           dn(i,j,iq,il)=t*f2
c           da(i,j,iq,il)=t*f1
c           dc(i,j,iq,il)=t*f1
c           df(i,j,iq,il)=etan*(da(i,j,iq,il)-2.*dl(i,j,iq,il))
c           dr(i,j,iq,il)=t*f3
          enddo
        enddo
      else    ! crust
      endif
 1100 continue
cs      the bspline of each of the params in the sph file is printed
cs      write(13,*) kdiff, 6371.0-(rr*6371.0), rho_x*7.6096669308370943, rho_x
 1000 continue

      return
      end

c-----------------------------
c       BOXCAR SHAPE FUNCTIONS
c-----------------------------
        subroutine bcfun(ord, nknots, knot, xi, rho_x)

        integer, intent(in)           :: ord, nknots
        integer, parameter            :: NKNOTS_MAX = 1000
        double precision, parameter   :: TOL = 1e-4
        double precision, dimension(1:nknots), intent(in) :: knot
        double precision, intent(in)  :: xi
        double precision, intent(out) :: rho_x

        integer          :: ii, Nx
c        double precision, dimension(1:NKNOTS_MAX) :: hh

        Nx = nknots - 1
        !! Compute vector hh of spacings
        !call fill_hh(hh, knot, Nx)

        !! Consistency checks
        if ((xi - TOL) .lt. knot(nknots)) then
           write(6,*) 'err1'
           stop
        elseif ((xi + TOL) .gt. knot(1)) then
           write(6,*) 'err2'
           stop
        elseif (ord .gt. Nx) then
           write(6,*) 'err3'
           stop
        endif

        if (ord .eq. 1) then ! LHS
           if (xi .le. knot(ord) .and. xi .ge. knot(ord+1)) then
              ! x0<=x<=x1
              rho_x = 1.0
           else
              ! x>x2
              rho_x = 0.0
           endif
        elseif (ord .eq. Nx) then ! RHS - 1
           if (xi .le. knot(ord) .and. xi .ge. knot(ord+1)) then
              ! x0<=x<=x1
              rho_x = 1.0
           else
              ! x>x4
              rho_x = 0.0
           endif
        else ! Away from borders
           if (xi .le. knot(ord-1) .and. xi .ge. knot(ord)) then
              ! x0<=x<=x1
              rho_x = 0.0
           elseif (xi .le. knot(ord) .and. xi .ge. knot(ord+1)) then
              ! x1<=x<=x2
              rho_x = 1.0
           else
              ! x>x4
              rho_x = 0.0
           endif
        endif

        return
        end

c--------------------------------
c       B-SPLINES SHAPE FUNCTIONS
c--------------------------------
! Returns the ord'th spline's value defined at nknots points provided in
! rknots array at xi.

        subroutine spl(ord, nknots, knot, xi, rho_x)

        integer, intent(in)            :: ord, nknots
        integer, parameter :: NKNOTS_MAX = 1000
        double precision, parameter    :: TOL = 1.5
        double precision, dimension(0:nknots-1), intent(in) :: knot
        double precision, intent(in)   :: xi
        double precision, intent(out)  :: rho_x

        integer             ::  ii, Nx
        double precision    ::  coefa, coefb, coefc, coefd, denom
        double precision    :: denomsum, dd, denom1, denom2
        double precision, dimension(0:NKNOTS_MAX) :: hh

        Nx = nknots - 1
        !! Compute vector hh of spacings */
        !! hh = farray1(0, Nx - 1); */
        call fill_hh(hh, knot, Nx)

        !! Consistency checks */
        if ((xi - TOL) .gt. knot(Nx)) then
            print*,"xi=",xi," / knot(",Nx,")=",knot(Nx)
            stop
            !return 0.0;
        else if ((xi + TOL) .lt. knot(0)) then
            print*,"xi=",xi," / knot(0)=", knot(0)
            stop
            !return 0.0;
        else if (ord .gt. Nx) then
            print*,"Warning: spl index ",ord," exceeds knot count", Nx
            stop
            !return 0.0;
        endif

        if (ord .eq. 0) then !          ! LHS */
            denom = 3. * hh(ord) * hh(ord) + 3. * hh(ord) * hh(ord + 1)
     1              + hh(ord + 1) * hh(ord + 1)
            if (xi .ge. knot(ord) .and. xi .le. knot(ord + 1)) then !! x0.le.x.le.x1 */
                coefa = 4. / (hh(ord) * (hh(ord) + hh(ord + 1)) * denom)
                coefb = 0.0
                coefc = -12 / denom
                coefd = 4 * (2 * hh(ord) + hh(ord + 1)) / denom
                rho_x = coefa * (xi - knot(ord)) * (xi - knot(ord)) *
     1                  (xi - knot(ord))
                rho_x = rho_x + coefb * (xi - knot(ord)) *
     1                  (xi - knot(ord))
                rho_x = rho_x + coefc * (xi - knot(ord))
                rho_x = rho_x + coefd
            else if (xi .gt. knot(ord + 1) .and. xi .le. knot(ord + 2)) then !  ! x1.le.x.le.x2 */
                coefa = -4. / (hh(ord + 1) * (hh(ord) + hh(ord + 1)) *
     1                  denom)
                coefb = 12 / ((hh(ord) + hh(ord + 1)) * denom)
                coefc = -12. * hh(ord + 1) / ((hh(ord) + hh(ord + 1)) *
     1                  denom)
                coefd = 4. * hh(ord + 1) * hh(ord + 1) / ((hh(ord) +
     1                  hh(ord + 1)) * denom)

                rho_x = coefa * (xi - knot(ord + 1)) * (xi - knot(ord +
     1                  1)) * (xi - knot(ord + 1))
                rho_x = rho_x + coefb * (xi - knot(ord + 1)) *
     1                  (xi - knot(ord + 1))
                rho_x = rho_x + coefc * (xi - knot(ord + 1))
                rho_x = rho_x + coefd
            else    ! x.gt.x2 */
                rho_x = 0.0
            endif
        else if (ord .eq. 1) then !     ! LHS+1 */
            denom = (3. * hh(ord - 1) * hh(ord - 1) + 4. * hh(ord - 1) *
     1               hh(ord) + hh(ord) * hh(ord) + 2. * hh(ord - 1) *
     2               hh(ord + 1) + hh(ord) * hh(ord + 1))
            denomsum = hh(ord - 1) + hh(ord) + hh(ord + 1)
            dd = denomsum * denom
            if (xi .ge. knot(ord - 1) .and. xi .le. knot(ord)) then !! x0.le.x.le.x1 */
                coefa = -4. * (3. * hh(ord - 1) + 2. * hh(ord) + hh(ord
     1                 +1))/(hh(ord - 1) * (hh(ord - 1) + hh(ord)) * dd)
                coefb = 0.
                coefc = 12. / denom
                coefd = 0.

                rho_x = coefa * (xi - knot(ord - 1)) * (xi - knot(ord -
     1                  1)) * (xi - knot(ord - 1))
                rho_x = rho_x + coefb * (xi - knot(ord - 1)) * (xi -
     1                  knot(ord - 1))
                rho_x = rho_x + coefc * (xi - knot(ord - 1))
                rho_x = rho_x + coefd
            else if (xi .ge. knot(ord) .and. xi .le. knot(ord + 1)) then!  ! x1.le.x.le.x2 */
                coefa = 4. * (2. * hh(ord - 1) * hh(ord - 1) + 6. *
     1                  hh(ord - 1) * hh(ord) + 3. * hh(ord) * hh(ord) +
     2                  3. * hh(ord - 1) * hh(ord + 1) + 3. * hh(ord) *
     3                  hh(ord + 1) + hh(ord + 1) * hh(ord + 1)) /
     4                  (hh(ord) * (hh(ord - 1) + hh(ord)) * (hh(ord) +
     5                  hh(ord + 1)) * dd)
                coefb = -12. * (3. * hh(ord - 1) + 2. * hh(ord) + hh(ord
     1                  + 1)) / ((hh(ord - 1) + hh(ord)) * dd)
                coefc = 12. * (-2. * hh(ord - 1) * hh(ord - 1) + hh(ord)
     1                  * hh(ord) + hh(ord) * hh(ord + 1)) /
     2                  ((hh(ord - 1) + hh(ord)) * dd)
                coefd = 4. * hh(ord - 1) * (4. * hh(ord - 1) * hh(ord) +
     1                  3. * hh(ord) * hh(ord) + 2. * hh(ord - 1) *
     2                  hh(ord + 1) + 3. * hh(ord) * hh(ord + 1)) /
     3                  ((hh(ord - 1) + hh(ord)) * dd)

                rho_x = coefa * (xi - knot(ord)) * (xi - knot(ord)) *
     1                  (xi - knot(ord))
                rho_x = rho_x +  coefb * (xi - knot(ord)) * (xi -
     1                  knot(ord))
                rho_x = rho_x +  coefc * (xi - knot(ord))
                rho_x = rho_x +  coefd
            else if (xi .ge. knot(ord + 1) .and. xi .le. knot(ord + 2)) then !  ! x2.le.x.le.x3 */
                dd = dd *  (hh(ord) + hh(ord + 1))
                coefa = -4. * (2. * hh(ord - 1) + hh(ord)) / (hh(ord +
     1                  1) * dd)
                coefb = 12. * (2. * hh(ord - 1) + hh(ord)) / dd
                coefc = -12. * (2. * hh(ord - 1) + hh(ord)) * hh(ord +
     1                  1) / dd
                coefd = 4. * (2. * hh(ord - 1) + hh(ord)) * hh(ord + 1)
     1                  * hh(ord + 1) / dd

                rho_x = coefa * (xi - knot(ord + 1)) * (xi - knot(ord +
     1                  1)) * (xi - knot(ord + 1))
                rho_x = rho_x +  coefb * (xi - knot(ord + 1)) * (xi -
     1                  knot(ord + 1))
                rho_x = rho_x +  coefc * (xi - knot(ord + 1))
                rho_x = rho_x +  coefd
            else ! x.gt.x3 */
                rho_x = 0.0
            endif
        else if (ord .eq. Nx - 1) then !        ! RHS-1 */
            denom = hh(ord - 2) * hh(ord - 1) + hh(ord - 1) * hh(ord -
     1              1) + 2. * hh(ord - 2) * hh(ord) + 4. * hh(ord - 1) *
     2              hh(ord) + 3. * hh(ord) * hh(ord)
            denomsum = hh(ord - 2) + hh(ord - 1) + hh(ord)
            dd = denomsum * denom
            if (xi .ge. knot(ord - 2) .and. xi .le. knot(ord - 1)) then !       ! x0.le.x.le.x1 */
                coefa = 4. * (hh(ord - 1) + 2. * hh(ord)) / (hh(ord - 2)
     1                  * (hh(ord - 2) + hh(ord - 1)) * dd)
                coefb = 0.0
                coefc = 0.0
                coefd = 0.0

                rho_x = coefa * (xi - knot(ord - 2)) * (xi - knot(ord -
     1                  2)) * (xi - knot(ord - 2))
                rho_x = rho_x +  coefb * (xi - knot(ord - 2)) * (xi -
     1                  knot(ord - 2))
                rho_x = rho_x +  coefc * (xi - knot(ord - 2))
                rho_x = rho_x +  coefd
            else if (xi .ge. knot(ord - 1) .and. xi .le. knot(ord)) then!       ! x1.le.x.le.x2 */
                coefa = -4. * (hh(ord - 2) * hh(ord - 2) + 3. * hh(ord -
     1                  2) * hh(ord - 1) + 3. * hh(ord - 1) * hh(ord -
     2                  1) + 3. * hh(ord - 2) * hh(ord) + 6. * hh(ord -
     3                  1) * hh(ord) + 2. * hh(ord) * hh(ord)) /
     4                  (hh(ord - 1) * (hh(ord - 2) + hh(ord - 1)) *
     5                  (hh(ord - 1) + hh(ord)) * dd)
                coefb = 12. * (hh(ord - 1) + 2. * hh(ord)) / ((hh(ord -
     1                  2) + hh(ord - 1)) * dd)
                coefc = 12. * hh(ord - 2) * (hh(ord - 1) + 2. * hh(ord))
     1                  / ((hh(ord - 2) + hh(ord - 1)) * dd)
                coefd = 4. * hh(ord - 2) * hh(ord - 2) * (hh(ord - 1) +
     1                  2. * hh(ord))/((hh(ord - 2) + hh(ord - 1)) * dd)

                rho_x = coefa * (xi - knot(ord - 1)) * (xi - knot(ord -
     1                  1)) * (xi - knot(ord - 1))
                rho_x = rho_x +  coefb * (xi - knot(ord - 1)) * (xi -
     1                  knot(ord - 1))
                rho_x = rho_x +  coefc * (xi - knot(ord - 1))
                rho_x = rho_x +  coefd
            else if (xi .ge. knot(ord) .and. xi .le. knot(ord + 1)) then!       ! x2.le.x.le.x3 */
                dd = dd *  (hh(ord - 1) + hh(ord))
                coefa = 4. * (hh(ord - 2) + 2. * hh(ord - 1) + 3. *
     1                  hh(ord)) / (hh(ord) * dd)
                coefb = -12. * (hh(ord - 2) + 2. * hh(ord - 1) + 3. *
     1                  hh(ord)) / dd
                coefc = 12. * (-hh(ord - 2) * hh(ord - 1) - hh(ord - 1)
     1                  * hh(ord - 1) + 2. * hh(ord) * hh(ord)) / dd
                coefd = 4. * hh(ord) * (3. * hh(ord - 2) * hh(ord - 1) +
     1                  3. * hh(ord - 1) * hh(ord - 1) + 2. * hh(ord -
     2                  2) * hh(ord) + 4. * hh(ord - 1) * hh(ord)) / dd

                rho_x = coefa * (xi - knot(ord)) * (xi - knot(ord)) *
     1                  (xi - knot(ord))
                rho_x = rho_x +  coefb * (xi - knot(ord)) * (xi -
     1                  knot(ord))
                rho_x = rho_x +  coefc * (xi - knot(ord))
                rho_x = rho_x +  coefd
            else ! x.gt.x4 */
                rho_x = 0.0
            endif
        else if (ord .eq. Nx) then !    ! RHS */
            denom = (hh(ord - 2) + hh(ord - 1)) * (hh(ord - 2) * hh(ord
     1              - 2) + 3. * hh(ord - 2) * hh(ord - 1) + 3. * hh(ord
     2              - 1) * hh(ord - 1))
            if (xi .ge. knot(ord - 2) .and. xi .le. knot(ord - 1)) then !! x0.le.x.le.x1 */
                coefa = 4. / (hh(ord - 2) * denom)
                coefb = 0.0
                coefc = 0.0
                coefd = 0.0
                rho_x = coefa * (xi - knot(ord - 2)) * (xi -
     1                  knot(ord - 2)) * (xi - knot(ord - 2))
                rho_x = rho_x +  coefb * (xi - knot(ord - 2)) * (xi -
     1                  knot(ord - 2))
                rho_x = rho_x +  coefc * (xi - knot(ord - 2))
                rho_x = rho_x +  coefd
            else if (xi .ge. knot(ord - 1) .and. xi .le. knot(ord)) then!       ! x1.le.x.le.x2 */
                coefa = -4. / (hh(ord - 1) * denom)
                coefb = 12 / denom
                coefc = 12 * hh(ord - 2) / denom
                coefd = 4. * hh(ord - 2) * hh(ord - 2) / denom

                rho_x = coefa * (xi - knot(ord - 1)) * (xi -
     1                  knot(ord - 1)) * (xi - knot(ord - 1))
                rho_x = rho_x +  coefb * (xi - knot(ord - 1)) * (xi -
     1                  knot(ord - 1))
                rho_x = rho_x +  coefc * (xi - knot(ord - 1))
                rho_x = rho_x +  coefd
            else ! x.gt.x2 */
                rho_x = 0.0
            endif
        else !          ! Away from borders */
            denom1 = hh(ord - 2) + hh(ord - 1) + hh(ord) + hh(ord + 1)
            if (xi .ge. knot(ord - 2) .and. xi .le. knot(ord - 1)) then !       ! x0.le.x.le.x1 */
                coefa = 4. / (hh(ord - 2) * (hh(ord - 2) + hh(ord - 1))*
     1                  (hh(ord - 2) + hh(ord - 1) + hh(ord)) * denom1)
                coefb = 0.0
                coefc = 0.0
                coefd = 0.0
                rho_x = coefa * (xi - knot(ord - 2)) * (xi -
     1                  knot(ord - 2)) * (xi - knot(ord - 2))
                rho_x = rho_x +  coefb * (xi - knot(ord - 2)) * (xi -
     1                  knot(ord - 2))
                rho_x = rho_x +  coefc * (xi - knot(ord - 2))
                rho_x = rho_x +  coefd
            else if (xi .ge. knot(ord - 1) .and. xi .le. knot(ord)) then!       ! x1.le.x.le.x2 */
                denom2 = (hh(ord - 2) + hh(ord - 1)) * (hh(ord - 2) +
     1                   hh(ord - 1) + hh(ord))
                denom = denom1 * denom2

                coefa = -4. * (hh(ord - 2) * hh(ord - 2) + 3. *
     1                  hh(ord - 2) * hh(ord - 1) + 3. * hh(ord - 1) *
     2                  hh(ord - 1) + 2. * hh(ord - 2) * hh(ord) +  4. *
     3                  hh(ord - 1) * hh(ord) + hh(ord) * hh(ord) +
     4                  hh(ord - 2) * hh(ord + 1) + 2. * hh(ord - 1) *
     5                  hh(ord + 1) + hh(ord) * hh(ord + 1)) /
     6                  (hh(ord - 1) * (hh(ord - 1) + hh(ord)) *
     7                  (hh(ord - 1) + hh(ord) + hh(ord + 1)) * denom)
                coefb = 12. / denom
                coefc = 12. * hh(ord - 2) / denom
                coefd = 4. * hh(ord - 2) * hh(ord - 2) / denom

                rho_x = coefa * (xi - knot(ord - 1)) * (xi -
     1                  knot(ord - 1)) * (xi - knot(ord - 1))
                rho_x = rho_x +  coefb * (xi - knot(ord - 1)) * (xi -
     1                  knot(ord - 1))
                rho_x = rho_x +  coefc * (xi - knot(ord - 1))
                rho_x = rho_x +  coefd
            else if (xi .ge. knot(ord) .and. xi .le. knot(ord + 1)) then!       ! x2.le.x.le.x3 */
                denom2 = (hh(ord - 1) + hh(ord)) * (hh(ord - 2) +
     1                   hh(ord - 1) + hh(ord)) * (hh(ord - 1) + hh(ord)
     2                   + hh(ord + 1))
                denom = denom1 * denom2

                coefa = 4. * (hh(ord - 2) * hh(ord - 1) + hh(ord - 1) *
     1                  hh(ord - 1) + 2. * hh(ord - 2) * hh(ord) + 4. *
     2                  hh(ord - 1) * hh(ord) + 3. * hh(ord) * hh(ord) +
     3                  hh(ord - 2) * hh(ord + 1) + 2. * hh(ord - 1) *
     4                  hh(ord + 1) + 3. * hh(ord) * hh(ord + 1) +
     5                  hh(ord + 1) * hh(ord + 1)) / (hh(ord) *
     6                  (hh(ord) + hh(ord + 1)) * denom)
                coefb = -12. * (hh(ord - 2) + 2. * hh(ord - 1) + 2. *
     1                  hh(ord) + hh(ord + 1)) / denom
                coefc = 12. * (-hh(ord - 2) * hh(ord - 1) - hh(ord - 1)
     1                  * hh(ord - 1) + hh(ord) * hh(ord) + hh(ord) *
     2                  hh(ord + 1)) / denom
                coefd = 4. * (2. * hh(ord - 2) * hh(ord - 1) * hh(ord) +
     1                  2. * hh(ord - 1) * hh(ord - 1) * hh(ord) +
     2                  hh(ord - 2) * hh(ord) * hh(ord) + 2. *
     3                  hh(ord - 1) * hh(ord) * hh(ord) + hh(ord - 2) *
     4                  hh(ord - 1) * hh(ord + 1) + hh(ord - 1) *
     5                  hh(ord - 1) * hh(ord + 1) + hh(ord - 2) *
     6                  hh(ord) * hh(ord + 1) + 2. * hh(ord - 1) *
     7                  hh(ord) * hh(ord + 1)) / denom
                rho_x = coefa * (xi - knot(ord)) * (xi - knot(ord)) *
     1                  (xi - knot(ord))
                rho_x = rho_x +  coefb * (xi - knot(ord)) * (xi -
     1                  knot(ord))
                rho_x = rho_x +  coefc * (xi - knot(ord))
                rho_x = rho_x +  coefd
            else if (xi .ge. knot(ord + 1) .and. xi .le. knot(ord + 2)) then !  ! x3.le.x.le.x4 */
                denom2 = (hh(ord) + hh(ord + 1)) * (hh(ord - 1) +
     1                   hh(ord) + hh(ord + 1))
                denom = denom1 * denom2

                coefa = -4. / (hh(ord + 1) * denom)
                coefb = 12 / denom
                coefc = -12 * hh(ord + 1) / denom
                coefd = 4. * hh(ord + 1) * hh(ord + 1) / denom

                rho_x = coefa * (xi - knot(ord + 1)) * (xi -
     1                  knot(ord + 1)) * (xi - knot(ord + 1))
                rho_x = rho_x +  coefb * (xi - knot(ord + 1)) * (xi -
     1                  knot(ord + 1))
                rho_x = rho_x +  coefc * (xi - knot(ord + 1))
                rho_x = rho_x +  coefd
            else ! x.gt.x4 */
                rho_x = 0.0
            endif
        endif

        return
        end


c-----------------------------------
c       DISTANCES BETWEEN KNOT RADII
c-----------------------------------
!       Compute the distance between the b-spline knot radii

        subroutine fill_hh(hh, knot, Nx)

        integer, intent(in)                       :: Nx
        integer, parameter :: NKNOTS_MAX = 1000
        double precision, dimension(Nx), intent(in) :: knot
        double precision, dimension(1:NKNOTS_MAX), intent(out):: hh

        !LOCAL VARIABLES
        integer     :: ii

        do ii = 1, Nx
            hh(ii) = knot(ii + 1) - knot(ii)
        enddo

        return
        end


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

c--------------------------------------------------------------------------

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

      do 1000 iq=nstart,n

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

      function grinta(afac,bfac,cfac,dfac,rr,fsi)

      if(fsi.eq.0.) then
        add1=afac*log(rr)
      else
        add1=afac/(-1.*fsi)
      endif
      if(fsi.eq.1) then
        add2=(bfac*rr)*log(rr)
      else
        add2=(bfac*rr)/(1.-fsi)
      endif
      if(fsi.eq.2.) then
        add3=(cfac*rr**2)*log(rr)
      else
        add3=(cfac*rr**2)/(2.-fsi)
      endif
      if(fsi.eq.3.) then
        add4=(dfac*rr**3)*log(rr)
      else
        add4=(dfac*rr**3)/(3.-fsi)
      endif
      grinta=add1+add2+add3+add4

      return
      end


c--------------------------------------------------------------------

      function grintb(afac,bfac,cfac,dfac,rr,fsi)

      add1=afac/(fsi+1.)
      add2=(bfac*rr)/(fsi+2.0)
      add3=(cfac*rr**2)/(fsi+3.)
      add4=(dfac*rr**3)/(fsi+4.)
      grintb=add1+add2+add3+add4

      return
      end

c ----------------------------------------------------------------------

      subroutine hereada(lu,file,x,nstr,lmax,mask,lask,inorm)
      save
      character*(*) file
      dimension x(*),mask(*),lask(*)
      character*132 abuf
      character*80 form
      sqth=sqrt(.5)
      open(lu,file=file,status='old')
      read(lu,'(a132)') abuf
      llen=132
      do 30 i=1,132
      if(abuf(llen:llen).ne.' ') goto 40
   30 llen=llen-1
   40 read(abuf,'(10x,i5)') lmax
      num=llen-15-1-(lmax+1)-1-3-1

      if(llen.gt.15) then
        idef=0
        write(form,'(''(16x,'',i3,''i1,1x,i3,1x,80i1)'')') lmax+1
ca        write(6,'(a30,2i8)') form(1:30),lmax,num
        read(abuf,form) (lask(i+1),i=0,lmax),nstr,(mask(i),i=1,num)
        if(num.ne.nstr) pause 'error 1 in heread'
      else
        idef=1
        do 50 i=0,lmax
        lask(i+1)=1
   50   continue
      endif

      leny=0
      do 60 l=0,lmax
      if(lask(l+1).ne.0) leny=leny+2*l+1
   60 continue

      ip=0
  100 read(lu,'(a132)',end=99) abuf
      backspace lu
      ip=1+ip
      if(idef.ne.0) mask(ip)=1

      ind=1+(ip-1)*leny
      do 704 l=0,lmax
      if(lask(l+1).ne.0) then
        ind1=ind+2*l
        read(lu,'(11e12.4)')(x(i),i=ind,ind1)
        if(inorm.eq.1) x(ind)=x(ind)/sqth
        ind=ind1+1
      endif
  704 continue
      goto 100

   99 if(idef.ne.0) then
        nstr=ip
        nstruc=ip
      else
        nstruc=0
        do i=1,nstr
          if(mask(i).ne.0) nstruc=nstruc+1
        enddo
c       if(nstruc.ne.ip) pause 'error 3 in heread'
      endif

      ii=istlen(file)
ca      write(6,'(''model from'',1x,a)') file(1:ii)
ca      write(6,'(i4,'' parameters:  mask '',80i1)') nstr,(mask(i),i=1,nstr)
ca    write(6,'(4x,'' lmax ='',i4,'':  mask '',80i1)') lmax,(lask(i+1),i=0,lmax)
      close(lu)
      return
      end
