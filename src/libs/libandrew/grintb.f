      function grintb(afac,bfac,cfac,dfac,rr,fsi)

      add1=afac/(fsi+1.)
      add2=(bfac*rr)/(fsi+2.0)
      add3=(cfac*rr**2)/(fsi+3.)
      add4=(dfac*rr**3)/(fsi+4.)
      grintb=add1+add2+add3+add4

      return
      end
