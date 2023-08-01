      subroutine tadd(iy1,jd1,ih1,im1,fs1,iy2,jd2,ih2,im2,fs2,secs)
      double precision secs,dfs2
      iy2=iy1
      jd2=jd1
      ih2=ih1
      im2=im1
      dfs2=fs1+secs
      if(dfs2.ge.60.d0) then
        is2=dfs2
        imadd=is2/60
        im2=im2+imadd
        fs2=dfs2-imadd*60.d0
        if(fs2.lt.0.) then
          fs2=fs2+60.
          im2=im2-1
        endif
        if(fs2.ge.60.) then
          fs2=fs2-60.
          im2=im2+1
        endif
        if(fs2.lt.0.or.fs2.ge.60.) pause 'tadd: error 0'
        if(im2.ge.60) then
          ihadd=im2/60
          ih2=ih2+ihadd
          im2=im2-60*ihadd
          if(ih2.ge.24) then
            jdadd=ih2/24
            jd2=jd2+jdadd
            ih2=ih2-24*jdadd
   10       lp=lpyr(iy2)
            if(jd2.gt.365+lp) then
              jd2=jd2-365-lp
              iy2=1+iy2
              goto 10
            endif
          endif
        endif
      else if(dfs2.lt.0.d0) then
        imadd=-dfs2/60.d0+1.d0
        im2=im2-imadd
        fs2=dfs2+imadd*60.d0
        if(fs2.ge.60.0) then
          im2=im2+1
          fs2=fs2-60.
        endif
        if(fs2.lt.0) then
          im2=im2-1
          fs2=fs2+60.
        endif
        if(fs2.lt.0.or.fs2.ge.60.d0) pause 'tadd: error 1'
        if(im2.lt.0) then
          ihadd=(-im2-1)/60+1
          ih2=ih2-ihadd
          im2=im2+60*ihadd
          if(im2.lt.0.or.im2.ge.60) pause 'tadd: error 2'
          if(ih2.lt.0) then
            jdadd=(-ih2-1)/24+1
            jd2=jd2-jdadd
            ih2=ih2+24*jdadd
            if(ih2.lt.0.or.ih2.ge.24) pause 'tadd: error 3'
   20       if(jd2.le.0) then
              iy2=iy2-1
              lp=lpyr(iy2)
              jd2=jd2+365+lp
              goto 20
            endif
          endif
        endif
      else
        fs2=dfs2
      endif
      return
      end
