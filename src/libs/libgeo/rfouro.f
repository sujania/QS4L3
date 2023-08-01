      subroutine rfouro(data,m,nn)
      dimension data(*),ntot(3)
      ntot(1)=2*2**m
      xn=ntot(1)/2
      ntp2=ntot(1)+2
      if(nn.lt.0) go to 1
      nt=ntot(1)
      do 2 n=1,nt
    2 data(n)=data(n)/xn
      call four2(data,ntot,1,-1,0)
      return
    1 do 4 n=1,ntp2
    4 data(n)=data(n)/2.
      call four2(data,ntot,1,1,-1)
      return
      end
