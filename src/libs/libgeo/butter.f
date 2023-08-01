c---------------------------------------------------------------
      complex function butter(wbar,nord)
c     butterworth filter of order nord
c     acts as low pass filter if wbar=w/wcuthigh
c     acts as high pass filter if wbar=-wcutlow/w
c     acts as band pass if wbar=(w/wcuthigh-wcutlow/w)
c     f(t)=f(w)*exp(-iwt)    convention
      complex ci
      data ci/(0., 1.57079633)/
c     data ci/(0.,-1.57079633)/
      butter=cexp(-ci*nord)
      do 10 n=1,2*nord-1,2
10    butter=butter/(cmplx(wbar,0.)-cexp(ci*n/nord))
      return
      end
