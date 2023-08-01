      program sph2m

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
      character*80 mfl
      integer smax,nspl

c    read sph model
      read(5,111) mfl
111   format(a80)
      call rmod(mfl,x,smax,nspl)
ca      write(6,*) lmx,MXLH,ndmx,MXSPL
      if(smax.gt.MXLH) stop'lmx.gt.MXLH'
      if(nspl.gt.MXSPL) stop'ndmx.gt.MXSPL'
      natd=(smax+1)**2

      ind=1
      do n=1,nspl
      do j=0,smax
         ind1=ind+2*j
         do i=ind,ind1,1 
         if(mod(j,2)==0) write(6,'(e12.4)') x(i)
         enddo
         ind=ind1+1
      enddo
      enddo
      end

c ------------------------------------------------------------------------

      subroutine rmod(infl,x,smax,nspl)


      character*80 infl
      dimension x(*)
      integer smax,nspl

      open(10,file=infl,status='old')
      read(10,*) smax,nspl
     
      natd=(smax+1)**2
      ind=1
      do i=1,nspl
       do j=0,smax
        ind1=ind+2*j
        read(10,'(11e12.4)',end=100)(x(k),k=ind,ind1)
        ind=ind1+1
       enddo
      enddo

      goto 200

 100  stop 'incompatible sph header'

 200  continue
      
      end
