      function jround(x)
      if(x.gt.0.) then
        jround=ifix(x+.5)
      else
        jround=-ifix(-x+.5)
      endif
      return
      end
