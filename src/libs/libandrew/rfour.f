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
ca#if ( defined(MachineT) )
ca      n=2**(m+1)
ca      fac=2./float(n)
ca      if(n.le.4096) then
ca        if(nn.gt.0) then
ca          call rfftn(nn,data,work,m+1)
ca          do i=1,n
ca            data(i)=fac*work(i)
ca          enddo
ca          data(n+1)=data(2)
ca          data(n+2)=0.
ca          data(2)=0.
ca        else
ca          sav1=data(n+1)
ca          sav2=data(n+2)
ca          data(2)=sav1
ca          call rfftn(nn,data,work,m+1)
ca          work(n+1)=sav1
ca          work(n+2)=sav2
ca          do i=1,n+2
ca            data(i)=0.5*work(i)
ca          enddo
ca        endif
ca      else
ca        write(6,*) 'rfour: Warning: old fft used -- insufficient work space'
ca        call rfouro(data,m,nn)
ca      endif
ca#else
      call rfouro(data,m,nn)
ca#endif
      return
      end
