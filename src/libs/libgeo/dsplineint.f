      subroutine dsplineint(f,x,i1,i2,a,qwk1,qwk2)
c
c    splint does quadrature on a spline fitted integrand.
c
c     f = integrand array
c     x = abciscae
c     i1, i2 =indicies in x() of starting and ending points
c     a = answer
c     qwk1, qwk2 work arrays of dimension at least 3*(dimension of x)

      implicit double precision(a-h,p-z)
      dimension f(*),x(*),qwk1(3,*),qwk2(3,*)
      call drspln(i1,i2,x,f,qwk1,qwk2)
      a=0.d0
      do i=i1,i2-1
        b=x(i+1)-x(i)
        a=a+b*(f(i)+b*(.5d0*qwk1(1,i)+b*(t*qwk1(2,i)+.25d0*b*qwk1(3,i))))
      enddo
      return
      end

