
c-------------------------
      subroutine loadchr(iarr,nwords,string)
      integer*4 iarr(*)
      character*(*) string
      character*4 str4
      equivalence (str4,istr4)
      lstring=len(string)
      k=1
      do i=1,nwords
        ke=min0(k+3,lstring)
        if(k.le.ke) then
          str4=string(k:ke)
        else
          str4=' '
        endif
        iarr(i)=istr4
        k=k+4
      enddo
      return
      end
