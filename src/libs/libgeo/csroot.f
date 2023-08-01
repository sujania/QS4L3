      subroutine csroot(xr,xi,yr,yi)
      real xr,xi,yr,yi
c
c     (yr,yi) = complex sqrt(xr,xi)
c     branch chosen so that yr .ge. 0.0 and sign(yi) .eq. sign(xi)
c
      real s,tr,ti,pythag
      tr = xr
      ti = xi
      s = sqrt(0.5e0*(pythag(tr,ti) + abs(tr)))
      if (tr .ge. 0.0e0) yr = s
      if (ti .lt. 0.0e0) s = -s
      if (tr .le. 0.0e0) yi = s
      if (tr .lt. 0.0e0) yr = 0.5e0*(ti/yi)
      if (tr .gt. 0.0e0) yi = 0.5e0*(ti/yr)
      return
      end
