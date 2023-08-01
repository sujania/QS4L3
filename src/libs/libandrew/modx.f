
      function modx(i,imod)
      if(i.gt.0) then
        modx=mod(i,imod)
      else
        ifac=(-i+imod-1)/imod
        modx=mod(i+ifac*imod,imod)
      endif
      return
      end
