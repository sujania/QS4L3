      subroutine timsec(jy,jd,ih,im,fsec,jsec)
      if(jy.lt.1910.or.jy.gt.2030) then
        jsec=z'7fffffff'
      else
        isec=fsec
        call juldoc(jy,jd,jj)
c !valid 1901-2036
        jsec=60*((jj-25221)*1440+(im+60*ih))+isec
      endif
      return
      end
