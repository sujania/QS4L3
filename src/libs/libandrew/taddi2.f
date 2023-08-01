      subroutine taddi2(isec1,ifrc1,isec2,ifrc2,add)
      double precision add
      if(add.gt.0.d0) then
        isadd=add
        iadd=10000.d0*(add-isadd)+.5d0
        ifrc2=ifrc1+iadd
        isec2=isec1+isadd+ifrc2/10000
        ifrc2=mod(ifrc2,10000)
      else
        isadd=-add
        isub=10000.d0*(-add-isadd)+.5
        ifrc2=ifrc1-isub
        isec2=isec1-isadd
        if(ifrc2.lt.0) then
          nm=(-ifrc2)/10000+1
          isec2=isec2-nm
          ifrc2=ifrc2+nm*10000
          isec2=isec2+ifrc2/10000
          ifrc2=mod(ifrc2,10000)
        else
          isec2=isec2+ifrc2/10000
          ifrc2=mod(ifrc2,10000)
        endif
      endif
      return
      end
