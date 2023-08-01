c-------------------------------------
      subroutine splhsetup()
      parameter (MXKNT=21)
      common/splhprm/spknt(MXKNT),qq0(MXKNT,MXKNT),qq(3,MXKNT,MXKNT)
      dimension qqwk(3,MXKNT)
      data spknt/
     1   -1.00000,-0.78631,-0.59207,-0.41550,-0.25499
     1  ,-0.10909, 0.02353, 0.14409, 0.25367, 0.35329
     1  , 0.44384, 0.52615, 0.60097, 0.66899, 0.73081
     1  , 0.78701, 0.83810, 0.88454, 0.92675, 0.96512
     1  , 1.00000
     1    /

      do i=1,MXKNT
        do j=1,MXKNT
          if(i.eq.j) then
            qq0(j,i)=1.
          else
            qq0(j,i)=0.
          endif
        enddo
      enddo
      do i=1,MXKNT
        call rspln(1,MXKNT,spknt(1),qq0(1,i),qq(1,1,i),qqwk(1,1))
      enddo
      return
      end
