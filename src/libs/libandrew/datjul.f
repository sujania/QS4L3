
      subroutine datjul(im,id,jy,jd)
      integer*4 kday(12,0:1)
      data kday/0,31,60,91,121,152,182,213,244,274,305,335
     1         ,0,31,59,90,120,151,181,212,243,273,304,334/
      jd=id+kday(im,min0(1,mod(jy,4)))
      return
      end
