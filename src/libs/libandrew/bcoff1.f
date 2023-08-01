c-----------------------------------------------
      subroutine bcoff1(l1,is,l2,b)
      implicit double precision (a-h,o-z)
      double precision b
      dimension b(5)
      dfloat(n)=n
      ls=l1+l2+is
      l1c=l1*(l1+1)
      l2c=l2*(l2+1)
      isc=is*(is+1)
      if((ls/2)*2.ne.ls) goto 10
      b(1)=wig0(l1,is,l2)
      ld=l1c+l2c-isc
      b(2)=dfloat(ld)*b(1)*0.5d0
      b(3)=0.5d0*dfloat((ld-2)*ld-2*l1c*l2c)*b(1)
      b(4)=0.d0
      b(5)=0.d0
      return
10    b(1)=0.d0
      b(2)=0.d0
      b(3)=0.d0
      b(4)=wig0(l1+1,l2+1,is+1)
      if(is.lt.iabs(l1-l2).or.is.gt.(l1+l2))then
        t1=0.
        goto 11
      endif
      t1=dfloat((ls+2)*(ls+4))*dfloat((ls+1-2*l2)*(ls+1-2*l1))
     1     *dfloat(ls+1-2*is)/dfloat(ls+3)
11    b(4)=b(4)*0.5d0*dsqrt(t1)
      b(5)=dfloat(l1c+l2c-isc-2)*b(4)
      return
      end
