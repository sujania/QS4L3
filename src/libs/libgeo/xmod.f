      function xmod(a,b)
      if(a.gt.0.) then
        xmod=amod(a,b)
        return
      else
        n=-a/b+2
        xmod=amod(a+n*b,b)
        return
      endif
      end
