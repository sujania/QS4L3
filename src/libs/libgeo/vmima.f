      subroutine vmima(x,inc,n,xmax,xmin)
      dimension x(inc,*)
      xmax=x(1,1)
      xmin=x(1,1)
      do i=2,n
        xmax=amax1(x(1,i),xmax)
        xmin=amin1(x(1,i),xmin)
      enddo
      return
      end
