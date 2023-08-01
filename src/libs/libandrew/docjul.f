
      subroutine docjul(jdoc,jy,jd)
      iy=jdoc/365
      j1=365*iy+(iy-1)/4
      if(jdoc.le.j1) then
        iy=iy-1
        j1=365*iy+(iy-1)/4
      endif
      jy=1900+iy
      jd=jdoc-j1
      end
