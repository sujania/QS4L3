
      subroutine sectim(jsec,jy,jd,ih,im,fsec)
      k=jsec
      isec=modx(k,60)
      k=(k-isec)/60
      mm=modx(k,1440)
      jj=25221+(k-mm)/1440
      im=modx(mm,60)
      ih=(mm-im)/60
      fsec=isec
      call docjul(jj,jy,jd)
      return
      end
