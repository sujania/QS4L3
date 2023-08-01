      program avmisfit
      character*80 dum
      real dum1,dum2,dum3,dum4
      integer ndat,cnt
      real*8 mf,xmisfit,totmf
      logical :: exist

      open(12,file='misfit.dat')
      totmf=0.0
      cnt=0
      ndattot=0
      do i=1,1000
         read(12,'(a)',END=100) dum
         read(12,*) mf, ndat, dum3
         if (mf.ge.0.000) then
         totmf=totmf+mf
         ndattot=ndattot+ndat
         endif
         cnt=cnt+1
      enddo

100   continue

      xmisfit=totmf/ndattot
      write(6,*) 'total misfit is',totmf
      write(6,*) 'total amount of modes is',cnt
      inquire(file='inversion.out', exist=exist)
      if (exist) then
        write(6,*) 'Average misfit is',xmisfit
        open(17,file='inversion.out')
        rewind(17)
        read(17,*) dum1,dum2,dum3,dum4
        close(17)

        open(18,file='inversion.dat',access='append')

        write(18,*) dum1,dum2,dum3,dum4,xmisfit
        close(18)
      else
        write(6,*) 'Average misfit is',xmisfit
        open(18,file='inversion.dat',access='append')
        write(18,*) '-    ','-    ','-    ','-    ',xmisfit
        close(18)
      endif
      end
