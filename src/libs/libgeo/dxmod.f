      double precision function dxmod(a,b)
      double precision a,b
      if(a.gt.0.d0) then
        dxmod=dmod(a,b)
        return
      else
        n=-a/b+2
        dxmod=dmod(a+n*b,b)
        return
      endif
      end
