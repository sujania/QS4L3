      program misfit_permode

      double precision obs(5000),new(5000),syn(5000),wt(1000)
      double precision avwt,nf,datasum,diffsum,diff,wttot,mftot,mf
      integer smax,n,ndata

      open(10,file='cstobs.dat')
      open(30,file='cstnew.dat')
      open(31,file='cstobs_S20+CRUST.dat')
      open(20,file='smax.dat')
      open(44,file='cstweight.dat')
      open(40,file='../misfit.dat', access='append')

      read(20,*) smax
      ndata=0
      do i=1,200
        read(31,*,END=222) dum
        ndata=ndata+1
      enddo
222   continue

      n=15 !(smax+1)*(smax+2)/2 ! n=15 for smax= 4
      do i=1,n
         read(10,*) obs(i)
         read(30,*) new(i)
         read(44,*) wt(i)
      enddo

      diffsum=0.0
      datasum=0.0
      ndata=0
      do i=1,n
         if(wt(i).ne.1.000000) then
         diff=sqrt((new(i)-obs(i))*(new(i)-obs(i)))
         dat=sqrt(obs(i)*obs(i))
         diffsum=diffsum+diff
         datasum=datasum+dat
         ndata=ndata+1
         endif
      enddo
      if(ndata>n) ndata=n
      mftot=diffsum/datasum
      
       write(40,*) mftot*ndata, ndata, mftot

      close(40)
      close(44)
      close(31)
      close(30)
      close(20)
      close(10)
      end
      
