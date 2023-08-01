      function grinta(afac,bfac,cfac,dfac,rr,fsi)

      if(fsi.eq.0.) then
        add1=afac*log(rr)
      else
        add1=afac/(-1.*fsi)
      endif
      if(fsi.eq.1) then
        add2=(bfac*rr)*log(rr)
      else
        add2=(bfac*rr)/(1.-fsi)
      endif
      if(fsi.eq.2.) then
        add3=(cfac*rr**2)*log(rr)
      else
        add3=(cfac*rr**2)/(2.-fsi)
      endif
      if(fsi.eq.3.) then
        add4=(dfac*rr**3)*log(rr)
      else
        add4=(dfac*rr**3)/(3.-fsi)
      endif
      grinta=add1+add2+add3+add4

      return
      end
