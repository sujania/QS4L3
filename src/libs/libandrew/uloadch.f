c-------------------------
      subroutine uloadch(iarr,nwords,string)
      integer*4 iarr(*)
      character*(*) string
      character*4 str4
      equivalence (str4,istr4)
      lstring=len(string)
      string=' '
      k=1
      do i=1,nwords
        istr4=iarr(i)
        ke=min0(k+3,lstring)
        if(k.le.ke) then
          string(k:ke)=str4(1:ke-k+1)
        endif
        k=k+4
      enddo
      do i=1,lstring
        if(string(i:i).eq.'\0') then
          string(i:lstring)=' '
          return
        endif
      enddo
      return
      end
