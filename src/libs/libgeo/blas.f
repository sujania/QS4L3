#if ( !defined(MachineT) )
      subroutine scopy(n,x,incx,y,incy)
      dimension x(incx,*),y(incy,*)
      do i=1,n
        y(1,i)=x(1,i)
      enddo
      return
      end

      subroutine dscopy(n,x,incx,y,incy)
      real*8 x(incx,*),y(incy,*)
      do i=1,n
        y(1,i)=x(1,i)
      enddo
      return
      end

      subroutine dcscopy(n,x,incx,y,incy)
      complex*16 x(incx,*),y(incy,*)
      do i=1,n
        y(1,i)=x(1,i)
      enddo
      return
      end

      subroutine saxpy(n,s,a,inca,d,incd)
      dimension a(inca,*),d(incd,*)
      do i=1,n
        d(1,i)=d(1,i)+s*a(1,i)
      enddo
      return
      end

      subroutine dsaxpy(n,s,a,inca,d,incd)
      real*8 a(inca,*),d(incd,*),s
      do i=1,n
        d(1,i)=d(1,i)+s*a(1,i)
      enddo
      return
      end

      subroutine dcsaxpy(n,s,a,inca,d,incd)
      complex*16 a(inca,*),d(incd,*),s
      do i=1,n
        d(1,i)=d(1,i)+s*a(1,i)
      enddo
      return
      end

      subroutine sscal(n,s,a,inca)
      dimension a(inca,*)
      do i=1,n
        a(1,i)=s*a(1,i)
      enddo
      return
      end

      subroutine dsscal(n,s,a,inca)
      real*8 a(inca,*),s
      do i=1,n
        a(1,i)=s*a(1,i)
      enddo
      return
      end

      subroutine dcsscal(n,s,a,inca)
      complex*16 a(inca,*),s
      do i=1,n
        a(1,i)=s*a(1,i)
      enddo
      return
      end

      function sdot(n,a,inca,b,incb)
      dimension a(inca,*),b(incb,*)
      sum=0.
      do i=1,n
        sum=sum+a(1,i)*b(1,i)
      enddo
      sdot=sum
      return
      end 

      real*8 function dsdot(n,a,inca,b,incb)
      real*8 a(inca,*),b(incb,*),sum
      sum=0.d0
      do i=1,n
        sum=sum+a(1,i)*b(1,i)
      enddo
      dsdot=sum
      return
      end 

c     complex*16 function dcsdot(n,a,inca,b,incb)
c     complex*16 a(inca,*),b(incb,*),sum
c     sum=(0.d0,0.d0)
c     do i=1,n
c       sum=sum+a(1,i)*b(1,i)
c     enddo
c     dcsdot=sum
c     return
c     end 
      complex*16 function dcsdot(n,a,inca,b,incb)
      real*8 a(*),b(*),dsdot
      dcsdot=dcmplx( dsdot(n,a(1),2*inca,b(1),2*incb)
     1              -dsdot(n,a(2),2*inca,b(2),2*incb)
     1             , dsdot(n,a(1),2*inca,b(2),2*incb)
     1              +dsdot(n,a(2),2*inca,b(1),2*incb) )
      return
      end 


      subroutine sclear(n,a,inca)
      dimension a(inca,*)
      do i=1,n
        a(1,i)=0.0
      enddo
      return
      end

      subroutine dsclear(n,a,inca)
      real*8 a(inca,*)
      do i=1,n
        a(1,i)=0.0d0
      enddo
      return
      end

      subroutine dcsclear(n,a,inca)
      complex*16 a(inca,*)
      do i=1,n
        a(1,i)=(0.0d0,0.d0)
      enddo
      return
      end

      
#else
      subroutine blasdummy
      end
#endif
