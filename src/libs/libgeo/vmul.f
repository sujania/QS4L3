      subroutine vmul(a,i,b,j,c,k,n)
      dimension a(i*n),b(j*n),c(k*n)
      ik=1
      jk=1
      kk=1
      do 1 l=1,n
      c(kk)=b(jk)*a(ik)
      ik=ik+i
      jk=jk+j
      kk=kk+k
    1 continue
      return
      end
