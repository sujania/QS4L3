c----------------------------------------------------------------
      subroutine cmatinv(a,n,m,det,ipvot,index,pivot)
      implicit complex(a-h,o-z)
      dimension a(m,m),index(m,2),pivot(m),ipvot(m)
      equivalence (irow,jrow),(icol,jcol)
      det=(1.0,0.0)
      do j=1,n
        ipvot(j)=0
      enddo
      do 135 i=1,n
      t=0.0
      do 9 j=1,n
      if(ipvot(j)-1)13,9,13
   13 do 23 k=1,n
      if(ipvot(k)-1)43,23,81
   43 if(cabs(t)-cabs(a(j,k))) 83,23,23
   83 irow=j
      icol=k
      t=a(j,k)
   23 continue
    9 continue
      ipvot(icol)=ipvot(icol)+1
      if(irow-icol) 73,109,73
   73 det=-det
      do 12 l=1,n
      t=a(irow,l)
      a(irow,l)=a(icol,l)
      a(icol,l)=t
   12 continue
  109 index(i,1)=irow
      index(i,2)=icol
      pivot(i)=a(icol,icol)
      det=det*pivot(i)
c     following 6 statements to divide pivot row by pivot element
      a(icol,icol)=(1.,0.)
      pivr=(1.,0.)/pivot(i)
      do l=1,n
        a(icol,l)=a(icol,l)*pivr
      enddo
      do 135 li=1,n
      if(li-icol) 21,135,21
   21 t=a(li,icol)
      a(li,icol)=(0.,0.)
      do l=1,n
        a(li,l)=a(li,l)-a(icol,l)*t
      enddo
  135 continue
      do 3 i=1,n
      l=n-i+1
      if(index(l,1)-index(l,2))19,3,19
   19 jrow=index(l,1)
      jcol=index(l,2)
      do 549 k=1,n
      t=a(k,jrow)
      a(k,jrow)=a(k,jcol)
      a(k,jcol)=t
  549 continue
    3 continue
   81 return
      end
