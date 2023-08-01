      subroutine vabs(a,i,b,j,c,k,n)
      dimension a(i*n),b(j*n),c(k*n)
      ik=1
      jk=1
      kk=1
      do 1 l=1,n
      c(kk)=sqrt(b(jk)**2+a(ik)**2)
      ik=ik+i
      jk=jk+j
      kk=kk+k
    1 continue
      return
      end
