c-----------------------------------------------------
      subroutine vlint(n,x,y,np,xarg,a)
      dimension x(*),y(n,*),a(*)
      if(xarg.lt.x(1)) then
        i1=1
        i2=2
      else if(xarg.ge.x(n)) then
        i1=n-1
        i2=n
      else
       do i=1,n-1
         if(xarg.ge.x(i).and.xarg.lt.x(i+1)) then
           i1=i
           i2=i+1
           goto 10
         endif
       enddo
      endif
   10 continue
      w1=(x(i2)-xarg)/(x(i2)-x(i1))
      w2=(xarg-x(i1))/(x(i2)-x(i1))
      do i=1,np
        a(i)=y(i1,i)*w1+y(i2,i)*w2
      enddo
      return
      end
