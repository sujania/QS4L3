      program addATAm

      parameter(MMOD=1000)
      parameter(MDATA=7000)

      double precision ATd(MMOD),ATdadd
      double precision ATA(MMOD,MMOD),ATAadd,dTd,dTdadd
      character*2048 ATAfile,ATdfile,dTdfile

ca      write(6,*) 'Input sigma data'
ca      read(5,*) sigmad
      write(6,*) 'Input number of model paraaddmeters'
      read(5,*) maxmod
ca      write(6,*) 'Input scaling factor'
ca      read(5,*) xfac

      write(6,*) 'Input number of ATA, ATd and dTd files'
      read(5,*) maxfiles

      dTd=0.d0
      do j=1,MMOD
         ATd(j)=0.d0
      do i=1,MMOD
         ATA(i,j)=0.d0
      enddo
      enddo

      open(60,file='ATA.dat')
      rewind(60)
      do j=1,maxmod
      do i=1,maxmod
      read(60,*) ATA(i,j)
      enddo
      enddo
      close(60)

      open(70,file='ATd.dat')
      rewind(70)
      do i=1,maxmod
      read(70,*) ATd(i)
      enddo
      close(70)

      open(80,file='dTd.dat')
      rewind(80)
      read(80,*) dTd
      close(80)

      do ifile=1,maxfiles

      read(5,'(a)') ATAfile
      read(5,'(a)') ATdfile
      read(5,'(a)') dTdfile

      open(60,file=ATAfile)
      rewind(60)
      open(70,file=ATdfile)
      rewind(70)
      open(80,file=dTdfile)
      rewind(80)

ca Add ATA, ATd and dTd

      read(80,*) dTdadd
      dTd=dTd+dTdadd
      do j=1,maxmod
      read(70,*) ATdadd
      ATd(j)=ATd(j)+ATdadd
      do i=1,maxmod
         read(60,*) ATAadd
         ATA(i,j)=ATA(i,j)+ATAadd
      enddo
      enddo

      close(60)
      close(70)
      close(80)

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
