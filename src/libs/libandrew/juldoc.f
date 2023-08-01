
      subroutine juldoc(jy,jd,jdoc)
      iy=jy-1900
c !valid 1901-2099
      jdoc=365*iy+(iy-1)/4+jd
      return
      end
