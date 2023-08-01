      subroutine rfour(data,m,nn)
      dimension data(*),work(4098)
cc#if ( defined(MachineS) )
cc      n=2**(m+1)
cc      if(n.le.1024) then
cc        if(nn.ge.0) then
cc          call R$RX1$FFT(data,n,0)
cc          data(n+1)=data(2)
cc          data(2)=0.
cc          call sscal(n+2,2.,data,1)
cc        else
cc          call sscal(n+2,.5,data,1)
cc          data(2)=data(n+1)
cc          data(n+1)=0.
cc          call R$RX1$FFT(data,n,-1)
cc        endif
cc      else if(n.le.4096) then
cc        if(nn.ge.0) then
cc          call R$RXII1$FFT(data,n,0)
cc          data(n+1)=data(2)
cc          data(2)=0.
cc          call sscal(n+2,2.,data,1)
cc        else
cc          call sscal(n+2,.5,data,1)
cc          data(2)=data(n+1)
cc          data(n+1)=0.
cc          call R$RXII1$FFT(data,n,-1)
cc        endif
cc      else
cc        call rfouro(data,m,nn)
cc      endif
#if ( defined(MachineT) )
      n=2**(m+1)
      fac=2./float(n)
      if(n.le.4096) then
        if(nn.gt.0) then
          call rfftn(nn,data,work,m+1)
          do i=1,n
            data(i)=fac*work(i)
          enddo
          data(n+1)=data(2)
          data(n+2)=0.
          data(2)=0.
        else
          sav1=data(n+1)
          sav2=data(n+2)
          data(2)=sav1
          call rfftn(nn,data,work,m+1)
          work(n+1)=sav1
          work(n+2)=sav2
          do i=1,n+2
            data(i)=0.5*work(i)
          enddo
        endif
      else
        write(6,*) 'rfour: Warning: old fft used -- insufficient work space'
        call rfouro(data,m,nn)
      endif
#else
      call rfouro(data,m,nn)
#endif
      return
      end
