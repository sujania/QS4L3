      subroutine tdif(iy1,jd1,ih1,im1,fs1,iy2,jd2,ih2,im2,fs2,secs)
      double precision secs
      if(iy1.eq.0.and.iy2.eq.0) pause 'difference of indeterminate times'
      imin=0
      if(iy1.eq.0) then
        secs=1.d12
      else if(iy2.eq.0) then
        secs=-1.d12
      else
        jdif=0
        if(iy1.ne.iy2) then
          do i=min0(iy1,iy2),max0(iy1,iy2)-1
            if(mod(i,4).eq.0.and.(mod(i,100).ne.0.or.mod(i,400).eq.0))then
              jdif=jdif+366
            else
              jdif=jdif+365
            endif
          enddo
          if(iy1.lt.iy2) jdif=-jdif
        endif
        secs=fs1-fs2+60.d0*(float(im1-im2)+60.d0*(float(ih1-ih2)
     1     +24.d0*float(jdif+jd1-jd2)))
      endif
      return
      end
