      subroutine cvmima(xr,incr,xi,inci,n,xmax,xmin)
      dimension xr(incr,*),xi(inci,*)
      xmax=xr(1,1)**2+xi(1,1)**2
      xmin=xmax
      do i=2,n
        t=xr(1,i)**2+xi(1,i)**2
        xmax=amax1(t,xmax)
        xmin=amin1(t,xmin)
      enddo
      xmax=sqrt(xmax)
      xmin=sqrt(xmin)
      return
      end
