      program m2sph

      parameter (MXLH=40)
      parameter (MXLENY=(MXLH+1)**2)
      parameter (MXDEP=21)
      parameter (MXSPL=MXDEP+3)
      parameter (MXMSZ=MXSPL*MXLENY)
      parameter (MAXP=1024)
      dimension wk1(MXLENY),wk2(MXLENY),wk3(MXLENY)
      dimension d0(MXLENY)
      dimension spl(MXDEP)
      dimension x(MXMSZ)
      integer smax,nspl
      character*80 mfl

c    read sph model
      read(5,111) mfl
111   format(a80)
      open(10,file=mfl,status='old')
      read(10,*) smax,nspl
      if(smax.gt.MXLH) stop'lmx.gt.MXLH'
      if(nspl.gt.MXSPL) stop'ndmx.gt.MXSPL'
      natd=(smax+1)**2

      open(20,file='mnew.dat',status='old')
      
      ind=1
      do n=1,nspl
      do i=0,smax
         ind1=ind+2*i
         do j=ind,ind1,1 
         if(mod(i,2)==0) read(20,*) x(j)
         if(mod(i,2)==1) x(j)=0.0000000
         enddo
         ind=ind1+1
      enddo
      enddo

      call wmod(10,x,smax,nspl)
      
      end

c ------------------------------------------------------------------------

      subroutine wmod(ifile,x,smax,nspl)


      character*80 infl
      dimension x(*)
      integer smax,nspl

      natd=(smax+1)**2
      ind=1
      do i=1,nspl
       do j=0,smax
        ind1=ind+2*j
        write(ifile,'(11e12.4)')(x(k),k=ind,ind1)
        ind=ind1+1
       enddo
      enddo
      end
