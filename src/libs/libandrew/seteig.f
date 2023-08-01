c----------------------------------------------------------------------
      subroutine seteig(lu)
      save
      integer*2 jobuf
      real*4 rbuf
      dimension iobuf(1500),jobuf(3000),rbuf(1500)
      equivalence (iobuf(1),jobuf(1),rbuf(1))
      integer*2 indf(1000),nrec,maxn,maxl
      common/radii/r(222)
      common/get/nrec(500,2),maxn(500,2),maxl(2),llu
c
c
      maxl(1)=0
      maxl(2)=-1
      llu=lu
c     read(llu) nuf,irlast,(r(i),i=1,31),(indf(k),k=1,200)
c     read(llu) (indf(k),k=201,nuf)
c  must be replaced by call to bffi
c
      print *, "seteig"
      call bffi(llu,1,iobuf,1500*4,jstat,nread,1)
      print *,iobuf(1:10)
      call byswap4(iobuf(1),1500)
      print *, iobuf(1:10)
      nuf=iobuf(1)
      irlast=iobuf(2)
      no=iobuf(3)
      nslo=iobuf(4)
      ninf=iobuf(5)
      ind=5
      do 5533 i=1,no
      ind=1+ind
 5533 r(i)=rbuf(ind)

      numr=nuf/ninf
      if(numr*ninf.lt.nuf) numr=1+numr
      i2=0

      do 230 i=1,numr
      i1=i2+1
      i2=min0(i1+ninf-1,nuf)
      call bffi(llu,1,iobuf,(i2-i1+1)*2,jstat,nread,0)
      call byswap2(jobuf,(i2-i1+1))
      do 5511 k=i1,i2
 5511 indf(k)=jobuf(k-i1+1)
  230 continue

      do 10 i=1,nuf
      itors=2
      if(indf(i).lt.0) itors=1
      maxl(itors)=maxl(itors)+1
      ind=maxl(itors)+1
      iarg=indf(i)
      nrec(ind,itors)=iabs(iarg)+1
      maxn(ind,itors)=0
      if(i.eq.nuf) goto 10
      iarg1=indf(i+1)
      maxn(ind,itors)=iabs(iarg1)-nrec(ind,itors)
      nrec(ind,itors)=nrec(ind,itors)-1
   10 continue
      return
      end
