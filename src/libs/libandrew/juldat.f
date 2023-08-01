
      subroutine juldat(jy,jd,im,id)
      integer*4 kday(12,0:1)
      data kday/0,31,60,91,121,152,182,213,244,274,305,335
     1         ,0,31,59,90,120,151,181,212,243,273,304,334/
      leap=min0(1,mod(jy,4))
      im=min0(12,1+jd/30)
      id=jd-kday(im,leap)
      if(id.le.0) then
        im=im-1
        id=jd-kday(im,leap)
      endif
      return
      end
