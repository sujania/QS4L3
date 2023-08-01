c----------------------------------------------------------------
      subroutine cgrl(nm,n,a,w,zr,zl,fv,ind,ngroup,ierr)
      dimension a(2*nm*nm),w(2*nm),zr(nm*nm*2),zl(nm*nm*2),fv(5*nm),ind(0:nm)

      tol=1.e-5

      do i=1,nm
      do j=1,nm
      do k=1,2
        zl(j+nm*(i-1+nm*(k-1)))=a(k+2*(j-1+nm*(i-1)))
      enddo
      enddo
      enddo
      call cg(nm,n,zl(1),zl(1+nm*nm),w(1),w(1+nm),1
     1     ,zr(1),zr(1+nm*nm),fv(1),fv(1+nm),fv(1+nm*2),ierr)


c group and order the eigenvalues and eigenvectors
      ngroup=0
      do i=1,n
        fv(i)=0.0
      enddo
      mput=0
  100 imin=0
      do i=1,n
        if(fv(i).eq.0.) then
          if(imin.eq.0) then
            imin=i
            cmin=w(i)
          else
            if(w(i).lt.cmin) then
              imin=i
              cmin=w(i)
            endif
          endif
        endif
      enddo
      if(imin.ne.0) then ! new group
        ind(ngroup)=mput
        mput=mput+1
        ngroup=ngroup+1
        ind(ngroup)=mput
        fv(imin)=mput
        do i=1,n
          if( fv(i).eq.0..and.cabs(cmplx(w(i),w(i+nm))-cmplx(w(imin),w(imin+nm))).lt.tol) then
            mput=mput+1
            fv(i)=mput
            ind(ngroup)=mput
          endif
        enddo
        goto 100
      endif
      if(mput.ne.n) stop 'wrong mput value'
      do i=1,n
        mput=fv(i)
        do k=1,2
          fv(nm+k+2*(mput-1))=w(i+nm*(k-1))
          do j=1,n
            zl(k+2*(j-1+nm*(mput-1)))=zr(j+nm*(i-1+nm*(k-1)))
          enddo
        enddo
      enddo
      call scopy(2*nm,fv(nm+1),1,w,1)

c     ngroup=0
c     ind(ngroup)=0
c     do i=1,nm-1
c       if(cabs(cmplx(w(1+2*(i-1)),w(1+2*(i-1)))-cmplx(w(1+2*(i)),w(1+2*(i)))).gt.tol) then
c         ngroup=ngroup+1
c         ind(ngroup)=i
c       endif
c       ind(ngroup+1)=i+1
c     enddo
c     ngroup=ngroup+1


      call scopy(2*nm*nm,zl,1,zr,1)
      call cmatinv(zl,n,nm,det,fv(1),fv(1+nm),fv(1+3*nm))
      return
      end

