      program buildATA

      parameter(MMOD=4000)
      parameter(MDATA=7000)

      real*8 A(MDATA,MMOD),dobs(MDATA),dsyn(MDATA),weight(MDATA)
      real*8 ATd(MMOD)
      real*8 ATA(MMOD,MMOD),dTd

      write(6,*) 'Input number of model parameters'
      read(5,*) maxmod
      write(6,*) 'Input scaling factor'
      read(5,*) xfac

      do j=1,MMOD
         ATd(j)=0.0
      do i=1,MMOD
         ATA(i,j)=0.d0
      enddo
      enddo

      do iblock=1,1

ca Read in data

      if(iblock.eq.1) then
      open(30,file='cstobs.dat')
      rewind(30)
ca      read(30,*) maxdata
      maxdata=0
      do i=1,MDATA
         read(30,*,END=333) dobs(i)
         maxdata=maxdata+1
      enddo
333   continue
      write(6,*) 'maxdata is ',maxdata
      close(30)
      endif

ca Read in data weighting

      open(30,file='cstweight.dat')
      rewind(30)
      do i=1,maxdata
         read(30,*) weight(i)
      enddo
      close(30)

ca Read in derivatives

      if(iblock.eq.1) then
      open(20,file='A.dat')
      rewind(20)
ca      read(20,*) maxdata
      do j=1,maxmod
      do i=1,maxdata
         read(20,*) A(i,j)
         A(i,j)=xfac*A(i,j)
      enddo
      enddo
      close(20)
      endif

ca Read in synthetic data

      if(iblock.eq.1) then
      open(30,file='cstsyn.dat')
      rewind(30)
      do i=1,maxdata
         read(30,*) dsyn(i)
      enddo
      close(30)
      endif

      do i=1,maxdata
         dsyn(i)=-1.0*dsyn(i)+dobs(i)
         dTd=dTd+(dsyn(i)*dsyn(i))
      enddo

ca Compute (AT*A)

      do j=1,maxmod
      do i=1,maxmod
      do k=1,maxdata
ca         ATA(i,j)=ATA(i,j)+dble(A(k,i)*A(k,j)/(sigmad**2))
         ATA(i,j)=ATA(i,j)+dble(A(k,i)*A(k,j)/(weight(k)))
      enddo
      enddo
      enddo


      do i=1,maxmod
      do j=1,maxdata
ca         ATd(i)=ATd(i)+(A(j,i)*dsyn(j)/(sigmad**2))
         ATd(i)=ATd(i)+(A(j,i)*dsyn(j)/(weight(j)))
      enddo
      enddo

      enddo

      open(60,file='ATA.dat')
      rewind(60)
      do j=1,maxmod
      do i=1,maxmod
      write(60,*) ATA(i,j)
      enddo
      enddo
      close(60)

      open(70,file='ATd.dat')
      rewind(70)
      do i=1,maxmod
      write(70,*) ATd(i)
      enddo
      close(70)

      open(80,file='dTd.dat')
      rewind(80)
      write(80,*) dTd
      close(80)

      
      end
