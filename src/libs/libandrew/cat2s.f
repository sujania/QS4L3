      character*(*) function cat2s(s1,s2,lcat2s)
      character*(*) s1,s2
      if(len(s1)+len(s2).gt.len(cat2s) ) then
        write(6,*) 'lengths:',len(s1),len(s2),len(cat2s)
        stop 'cat2s: buffer would overflow'
      endif
      k=0
      do i=1,len(s1)
        k=k+1
        cat2s(k:k)=s1(i:i)
      enddo
      do i=1,len(s2)
        k=k+1
        cat2s(k:k)=s2(i:i)
      enddo
      lcat2s=k
      do i=k+1,len(cat2s)
        cat2s(i:i)=' '
      enddo
      return
      end

