c-----------------------------------------------------------------
      subroutine curfit(x,y,sigmay,npts,nterms,mode,a,deltaa,
     1 sigmaa,flamda,yfit,chisqr,func)
      dimension x(*),y(*),sigmay(*),a(*),deltaa(*),sigmaa(*),yfit(*)
      parameter (MXYS=1000)
      parameter (MXAS=20)
      dimension weight(MXYS),alpha(MXAS,MXAS),beta(MXAS),deriv(MXYS,MXAS),
     1 array(MXAS,MXAS),b(MXAS)
      external func
      if(npts.gt.MXYS) stop 'curfit: npts too large'
      if(nterms.gt.MXAS) stop 'curfit: nterms too large'
      nterms1=0
      do j=1,nterms
        if(deltaa(j).ne.0.) nterms1=nterms1+1
      enddo
      nfree=npts-nterms1
      if(nfree.le.0) then
        chisqr=0.
        return
      endif
c
c     evaluate weights
c
      do i=1,npts
        if(mode.lt.0) then
          weight(i)=1.
          if(y(i).ne.0.) weight(i)=1./abs(y(i))
        else if(mode.eq.0) then
          weight(i)=1.
        else
          weight(i)=1./sigmay(i)**2
        endif
      enddo
c
c     evaluate alpha and beta matrices
c
      jj=0
      do j=1,nterms
        if(deltaa(j).ne.0.) then
          jj=jj+1
          beta(jj)=0.
          kk=0
          do k=1,j
            if(deltaa(k).ne.0.) then
              kk=kk+1
              alpha(jj,kk)=0.
            endif
          enddo
        endif
      enddo

      call func(x,yfit,a,npts)  
      call fderiv(x,a,deltaa,npts,nterms,deriv,MXYS,func)
      do i=1,npts
        jj=0
        do j=1,nterms
          if(deltaa(j).ne.0.) then
            jj=jj+1
            beta(jj)=beta(jj)+weight(i)*(y(i)-yfit(i))*deriv(i,j)
            kk=0
            do k=1,j
              if(deltaa(k).ne.0.) then
                kk=kk+1
                alpha(jj,kk)=alpha(jj,kk)+weight(i)*deriv(i,j)*deriv(i,k)
              endif
            enddo
          endif
        enddo
      enddo
      do jj=1,nterms1
        do kk=1,jj
          alpha(kk,jj)=alpha(jj,kk)
        enddo
      enddo

c
c     evaluate chi square at starting point
c
      chisq1=fchisq(y,sigmay,npts,nfree,mode,yfit)
c
c     invert modified curvature matrix to find new parameters
c
   71 continue
      do jj=1,nterms1
        do kk=1,nterms1
          array(jj,kk)=alpha(jj,kk)/sqrt(alpha(jj,jj)*alpha(kk,kk))
        enddo
        array(jj,jj)=1.+flamda
      enddo

      call cfmatinv(array,MXAS,nterms1,det)

      jj=0
      do j=1,nterms
        if(deltaa(j).ne.0.) then
          jj=jj+1
          b(j)=a(j)
          do kk=1,nterms1
            b(j)=b(j)+beta(kk)*array(jj,kk)/sqrt(alpha(jj,jj)*alpha(kk,kk))
          enddo
        else
          b(j)=a(j)
        endif
      enddo
c 
c     if chi squared increased, increase flamda and try again
c
c     write(6,*) 'b=',(b(i),i=1,nterms)
      call func(x,yfit,b,npts)
      chisqr=fchisq(y,sigmay,npts,nfree,mode,yfit)
      if(chisq1.le.0.) then
        test=0.
      else
        test=abs(chisqr-chisq1)/chisq1
      endif
      if(test.lt.1.e-5)  then
        write(6,'(''improvement saturates'')')
      else if(chisqr.gt.chisq1) then
        flamda=10.*flamda
        write(6,'(''chi sq increased this time:'',1p2e13.5)') chisq1,chisqr
        if(flamda.gt.1.)  then
          write(6,'(''flamda too large inversion ceases'')')
          do i=1,nterms
            a(i)=999
            sigmaa(i)=999
          enddo
          return
        else
          goto 71
        endif
      endif
      jj=0
      do j=1,nterms
        a(j)=b(j)
        if(deltaa(j).ne.0.) then
          jj=jj+1
          ratio=array(jj,jj)/alpha(jj,jj)
          ratio=ratio*chisqr
          sigmaa(j)=sqrt(abs(ratio))
        else
          sigmaa(j)=0.
        endif
      enddo
      flamda=flamda/10.
      return
      end

c-----------------------------------------------------------------
      subroutine curfitold(x,y,sigmay,npts,nterms,mode,a,deltaa,
     1 sigmaa,flamda,yfit,chisqr,func)
      dimension x(*),y(*),sigmay(*),a(*),deltaa(*),sigmaa(*),yfit(*)
      parameter (MXYS=1000)
      parameter (MXAS=20)
      dimension weight(MXYS),alpha(MXAS,MXAS),beta(MXAS),deriv(MXYS,MXAS),
     1 array(MXAS,MXAS),b(MXAS)
      external func
   11 nfree=npts-nterms
      if(nfree) 13,13,20
   13 chisqr=0.
      go to 110
c
c     evaluate weights
c
   20 do 30 i=1,npts
   21 if(mode) 22,27,29
   22 if(y(i)) 25,27,23
   23 weight(i)=1./y(i)
      go to 30
   25 weight(i)=1./(-y(i))
      go to 30
   27 weight(i)=1.
      go to 30
   29 weight(i)=1./sigmay(i)**2
   30 continue
c
c     evaluate alpha and beta matrices
c
   31 do j=1,nterms
        beta(j)=0.
        do k=1,j
          alpha(j,k)=0.
        enddo
      enddo
      call func(x,yfit,a,npts)  
      call fderiv(x,a,deltaa,npts,nterms,deriv,MXYS,func) 
c       write(6,1000)((deriv(i,j),j=1,nterms),i=1,npts)
1000  format('deriv', 2f10.4)
   41 do 50 i=1,npts
c      print*,'weight', weight(i)
      do 46 j=1,nterms
      beta(j)=beta(j)+weight(i)*(y(i)-yfit(i))*deriv(i,j)
      do 46 k=1,j
   46 alpha(j,k)=alpha(j,k)+weight(i)*deriv(i,j)*deriv(i,k)
   50 continue
   51 do 53 j=1,nterms
      do 53 k=1,j
   53 alpha(k,j)=alpha(j,k)
c
c     evaluate chi square at starting point
c
   63 chisq1=fchisq(y,sigmay,npts,nfree,mode,yfit)
c
c     invert modified curvature matrix to find new parameters
c
   71 do 74 j=1,nterms
      do 73 k=1,nterms
   73 array(j,k)=alpha(j,k)/sqrt(alpha(j,j)*alpha(k,k))
   74 array(j,j)=1.+flamda
c      write(6,3000) ((array(i,j),i=1,2),j=1,2)

   80 call cfmatinv(array,MXAS,nterms,det)
c      write(6,3000) ((array(i,j),i=1,2),j=1,2)
3000  format('array',2e10.4)
   81 do 84 j=1,nterms
      b(j)=a(j)
      do 84 k=1,nterms
   84 b(j)=b(j)+beta(k)*array(j,k)/sqrt(alpha(j,j)*alpha(k,k))
c 
c     if chi squared increased, increase flamda and try again
c
      call func(x,yfit,b,npts)
   93 chisqr=fchisq(y,sigmay,npts,nfree,mode,yfit)
      if(abs(chisqr-chisq1).lt.(1e-20))go to 8989
c      write(6,2000) b(1),b(2),chisqr
2000  format('b1 b2 and chisq', 3e15.7)
      if(chisq1-chisqr)95,101,101
   95 flamda=10.*flamda
c
      write(6,888)

  888 format(1x,'chi sq increased this time')
      if(flamda.gt.1.)go to 222
      go to 71
c
c     evaluate parameters and uncertainties
c
  101 do 103 j=1,nterms
      a(j)=b(j)
      ratio=array(j,j)/alpha(j,j)
c     maximum likelihood estimate
      ratio=ratio*chisqr
c  103 sigmaa(j)=sqrt(ratio)
c
  103 sigmaa(j)=sqrt(sqrt((ratio)**2))
c
      flamda=flamda/10.
  110 return
  222 write(6,7767)
 7767 format(1x,'flamda too large inversion ceases')
cccccccccccccccccccccccccc fix output if non convergent ccccccccccc
      do i=1,nterms
      a(i)=999
      sigmaa(i)=999
      enddo
 8989 continue
c
      write(6,6666)
 6666 format(1x,1x,'improvement saturates')
      return
      end
c_______________________________________________________________________________
      function fchisq(y,sigmay,npts,nfree,mode,yfit)
      dimension y(*),sigmay(*),yfit(*)
   11 chisq=0.
   12 if(nfree)13,13,20
   13 fchisq=0.
      go to 40
c
c      accumulate chi square
c
   20 do 30 i=1,npts
   21 if(mode)22,27,29
   22 if(y(i))25,27,23
   23 weight=1./y(i)
      go to 30
   25 weight=1./(-y(i))
      go to 30
   27 weight=1.
      go to 30
   29 weight=1./sigmay(i)**2
   30 chisq=chisq+weight*(y(i)-yfit(i))**2
c
c     divide by number degrees of freedom
c
   31 free=nfree
   32 fchisq=chisq/free
   40 return
      end
c_____________________________________________________________________________
      subroutine fderiv(x,a,deltaa,npts,nterms,deriv,ideriv,func)
      dimension x(*),a(*),deltaa(*),deriv(ideriv,*)
      parameter (MXYS=1000)
      dimension yfitp(MXYS),yfitm(MXYS)
      if(npts.gt.MXYS) stop 'fderiv: npts too large'
      if(npts.gt.ideriv) stop 'fderiv: npts too large'

      jj=0
      do j=1,nterms
        if(deltaa(j).ne.0.) then
          jj=jj+1
          aj=a(j)
          delta=deltaa(j)
          a(j)=aj+delta
          call func(x,yfitp,a,npts)
          a(j)=aj-delta
          call func(x,yfitm,a,npts)
          do k=1,npts
            deriv(k,j)=(yfitp(k)-yfitm(k))/(2.*delta)
          enddo
          a(j)=aj
        else
          do k=1,npts
            deriv(k,j)=0.
          enddo
        endif
      enddo

      return
      end
c_______________________________________________________________________________
c     subroutine cfmatinv
c     purpose invert a symetric matrix and calculate its inverse
c     use call cfmatinv(array,iarray,order,det)
c     description of parameters
c     array-input matrix which is replaced by its inverse
c     norder-degree of matrix
c     det-determinant
c     dimension statement valid for orders up to 10
      subroutine cfmatinv(array,iarray,norder,det)
c     double precision array,amax,save
      dimension array(iarray,*),ik(20),jk(20)
      if(norder.gt.20) stop 'cfmatinv: norder too large'
   10 det=1.
   11 do 100 k=1,norder
      amax=0.
   21 do 30 i=k,norder
      do 30 j=k,norder
   23 if(abs(amax)-abs(array(i,j)))24,24,30
   24 amax=array(i,j)
      ik(k)=i
      jk(k)=j
   30 continue
c
c
c
   31 if(amax)41,32,41
   32 det=0
      go to 140
   41 i=ik(k)
      if(i-k)21,51,43
   43 do 50 j=1,norder
      save=array(k,j)
      array(k,j)=array(i,j)
   50 array(i,j)=-save
   51 j=jk(k)
      if(j-k)21,61,53
   53 do 60 i=1,norder
      save=array(i,k)
      array(i,k)=array(i,j)
   60 array(i,j)=-save
c
c    accumulate elements of inverse matrix
c
   61 do 70 i=1,norder
      if(i-k)63,70,63
   63 array(i,k)=-array(i,k)/amax
   70 continue
   71 do 80 i=1,norder
      do 80 j=1,norder
      if(i-k)74,80,74
   74 if(j-k)75,80,75
   75 array(i,j)=array(i,j)+array(i,k)*array(k,j)
   80 continue
   81 do 90 j=1,norder
      if(j-k)83,90,83
   83 array(k,j)=array(k,j)/amax
   90 continue
      array(k,k)=1./amax
  100 det=det*amax
c
c     restore ordering of matrix
c
  101 do 130 l=1,norder
      k=norder-l+1
      j=ik(k)
      if(j-k)111,111,105
  105 do 110 i=1,norder
      save=array(i,k)
      array(i,k)=-array(i,j)
  110 array(i,j)=save
  111  i=jk(k)
       if(i-k)130,130,113
  113  do 120 j=1,norder
       save=array(k,j)
       array(k,j)=-array(i,j)
  120  array(i,j)=save
  130  continue
  140  return
       end
