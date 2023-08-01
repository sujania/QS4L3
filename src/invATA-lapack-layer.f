      program invATA

      parameter(MMOD=5000)
      parameter(MMOD2=500500)
      parameter(MDATA=100000)

      real*8 mzero(MMOD),mnew(MMOD),mtot,mzerotot,madd(MMOD),dTd
      double precision ATd(MMOD),mold(MMOD),ATAmnew(MMOD),xmisfit,varred
      double precision ATA(MMOD,MMOD),det,ATAold(MMOD,MMOD),Atest(MMOD,MMOD)
      double precision Rmatrix(MMOD,MMOD),Cmatrix(MMOD,MMOD),saveatd(MMOD)
      double precision identitym(MMOD,MMOD)

ca parameters for ahouse

      real*8 Y(MMOD),A(MMOD),B(MMOD),P(MMOD),TA(MMOD),TB(MMOD)
      real*8 eta(MMOD),etafac
      real*8 W(MMOD),V(MMOD),EV(MMOD),C(MMOD2),eva(MMOD),evec(MMOD,MMOD)

      double precision ATAinv(MMOD,MMOD),effeig
      integer info, lwork
      integer, allocatable, dimension(:) :: ipiv
      double precision, allocatable, dimension(:) :: ain, work

ca      write(6,*) 'Input eta'
      read(5,*) etafac
      write(6,*) 'Damping parameter eta is ',etafac
      write(6,*) '10*d damping for how many parameters, give 1st and last loc?'
      read(5,*) ndamp1
      read(5,*) ndamp2


      open(70,file='ATd.dat')
      rewind(70)
      maxmod=0
      do i=1,MMOD
      read(70,*,END=111) ATd(i)
      saveatd(i)=ATd(i)
      maxmod=maxmod+1
      enddo
      write(6,*) 'ERROR: increase model vector size'
 111  continue
      close(70)

      open(80,file='dTd.dat')
      rewind(80)
      read(80,*) dTd
      close(80)

      do i=1,maxmod
         eta(i)=dble(etafac)
      enddo

      if(ndamp1.ne.0) then
        do i=ndamp1,ndamp2
           eta(i)=dble(etafac)*10 !0.d0
        enddo
      endif

      do i=1,maxmod
         write(6,*) i, eta(i)
      enddo

      open(88,file='inversion.out')
      rewind(88)
      open(89,file='inversion-simplemisfit.dat',access='append')

ca Read in starting model

      open(40,file='mzero.dat')
      open(50,file='mold.dat')
      rewind(40)
      rewind(50)
      do i=1,maxmod
         read(40,*) mzero(i)
         read(50,*) mold(i)
      enddo
      close(40)
      close(50)

ca Read in ATA

      open(60,file='ATAmatrix.dat')
      rewind(60)
      do j=1,maxmod
      do i=1,maxmod
      read(60,*) ATA(i,j)
      enddo
      enddo
      close(60)


cl Invert ATA+damp*I

      N=maxmod
      allocate(ain(N*N))
      icount=0
      do i=1,N
         do j=1,N
            icount=icount+1
            if(i.eq.j) then ! only diag elements aare damped
            ain(icount)=ATA(i,j)+eta(i) !etafac
            else
            ain(icount)=ATA(i,j)
            endif
         enddo
      enddo

      allocate(ipiv(N))

cl ---- STEP 1: Compute LU factorization
      call dgetrf(N,N,ain,N,ipiv,info)

      if (info > 0) then
       write(*,*) "(dgetrf) The algorithm failed to compute eigenvalues"
      else if (info < 0) then
       write(*,*) "(dgetrf) The algorithm failed due to illegal value"
      endif

cl ---- STEP 2: Invert matrix through LU factorization
cl    get optimum size of work-array first:
      lwork=-1
      allocate(work(1))
      call dgetri(N,ain,N,ipiv,work,lwork,info)
      lwork=work(1)
      deallocate(work)
      allocate(work(lwork))

      call dgetri(N,ain,N,ipiv,work,lwork,info)
      deallocate(work,ipiv)

      icount=0
      do i=1,N
         do j=1,N
            icount=icount+1
            ATAinv(i,j)=ain(icount)
         enddo
      enddo
      deallocate(ain)

      do i=1,maxmod
         ATd(i)=ATd(i)-eta(i)*(mold(i)-mzero(i))
      enddo

      open(60,file='ATAtot.dat')
      do j=1,maxmod
      do i=1,maxmod
      write(60,*) ATA(i,j)
      ATAold(i,j)=ATA(i,j)
      if(i.eq.j) ATAold(i,i)=ATAold(i,i)+eta(i)
      enddo
      enddo
      close(60)

      do i=1,maxmod
      mnew(i)=mold(i)
      madd(i)=0.d0
      do j=1,maxmod
         mnew(i)=mnew(i)+sngl(ATAinv(i,j))*ATd(j)
         madd(i)=madd(i)+sngl(ATAinv(i,j))*ATd(j)
      enddo
      enddo

      effeig = 0.d0
      open(55,file='Rmatrix.dat')
      do i=1,maxmod
      do j=1,maxmod
         Atest(i,j)=0.d0
         Rmatrix(i,j)=0.d0
      do k=1,maxmod
         Atest(i,j)=Atest(i,j)+ATAold(i,k)*ATAinv(k,j)
         Rmatrix(i,j)=Rmatrix(i,j) + ATAinv(i,k)*ATA(k,j)
      enddo
      write(55,*) Rmatrix(i,j)
      if(i.eq.j) effeig = effeig + Rmatrix(i,j)
      enddo
      enddo
      close(55)

      open(66,file='Cmatrix.dat')
      do j=1,maxmod
      do i=1,maxmod
      Cmatrix(i,j)=ATAinv(i,j)
      write(66,*) Cmatrix(i,j)
      enddo
      enddo
      close(66)

      open(44,file='mnew.dat')
      rewind(44)
      mtot=0.0
      mzerotot=0.0
      do i=1,maxmod
      write(44,*) mnew(i),Rmatrix(i,i),Cmatrix(i,i),madd(i)
      mtot=mtot+(mnew(i)**2)
      mzerotot=mzerotot+(mzero(i)**2)
      enddo
      close(44)

      eps=0.0
      do i=1,maxmod
      if(mold(i).ne.0) then
      eps=eps+((mnew(i)-mold(i))/mold(i))**2
      endif
      if(mold(i).eq.0) then
      eps=eps+1.0
      endif
      enddo
      eps=sqrt(eps/maxmod)

ca Calculate variance reduction

      do i=1,maxmod
         ATAmnew(i)=0.d0
         do j=1,maxmod
            ATAmnew(i)=ATAmnew(i)+ATA(i,j)*madd(j)
         enddo
      enddo
      xmisfit=0.d0
      do i=1,maxmod
         xmisfit=xmisfit+(madd(i)*ATAmnew(i))-(2.0*madd(i)*saveatd(i))
      enddo
      varred=(xmisfit+dTd)/dTd

ca      write(6,*) ' '
      write(6,*) 'Epsilon is ',eps
      write(6,*) 'Effective # eigenvalues is ',effeig
      write(6,*) 'Mzero tot is ',mzerotot

      write(88,*) etafac,eps,effeig,mtot,xmisfit
      close(88)
      write(89,*) etafac,eps,effeig,mtot,xmisfit,varred
      close(89)


      end
