cprog wig0
      double precision function wig0(lone,ltwo,lthree)
      implicit double precision (a-h,o-z)
      dfloat(n)=n
      if((lone.ge.iabs(ltwo-lthree)).and.(lone.le.(ltwo+lthree)))goto 10
20    wig0=0.d0
      return
10    l=lone+ltwo+lthree
      lh=l/2
      if(lh*2.ne.l) goto 20
      l1=min0(lone,ltwo,lthree)
      l3=max0(lone,ltwo,lthree)
      l2=l-l1-l3
      sgn=1.d0
      if((lh/2)*2.ne.lh) sgn=-1.d0
      r1=1.d0/dfloat(l+1)
      if(l1.eq.0) goto 30
      do 40 j1=1,l1
      j12=2*j1
40    r1=r1*dfloat((j12-1)*(l-j12+2))/dfloat(j12*(l-j12+1))
30    numit=l3-lh+l1
      if(numit.eq.0) goto 50
      do 60 nit=1,numit
      j3=lh-l1+nit
      j2=lh-nit
60    r1=r1*dfloat((l-2*j2-1)*(l-2*j3+2))/dfloat((l-2*j2)*(l-2*j3+1))
50    wig0=dsqrt(r1)*sgn
      return
      end
