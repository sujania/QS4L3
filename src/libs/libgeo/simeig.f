c----------------------------------------------------------------
      subroutine simeig(nm,n,aa,ncom,w,zr,zl
     1      ,fv,ind,ngp,ind1,am,zrm,zlm)

c*********************************************
c     given ncom commuting matrices, find
c     simultaneous right and left eigenvectors
c*********************************************

c input and output

      implicit complex (a-h,o-z)
      dimension aa(nm,nm,*),w(nm,ncom),zr(nm,nm),zl(nm,nm)

c workspace

      real fv(nm,5)
      integer ind(0:nm,ncom),ngp(ncom),ind1(0:nm)
      complex am(nm,nm),zrm(nm,nm),zlm(nm,nm)
    
      ngroup=0
      igp=0
      ilev=0
      ngrp=0
      igp=0

   20 igp=igp+1
      if(igp.ge.ngrp) then ! all groups processed. go to next level
        if(ilev.gt.0) ngrp=ngp(ilev)
        ilev=ilev+1
        if(ilev.gt.ncom) goto 99
        ngp(ilev)=0
        igp=0
      endif
c     construct matrix to be decomposed at this level

      if(ilev.eq.1) then
        call cgrl(nm,n,aa(1,1,ilev),w(1,ilev),zr,zl,fv(1,1),ind(0,ilev),ngp(ilev),ierr)
      else
        ndim=ind(igp+1,ilev-1)-ind(igp,ilev-1)
        ibi=ind(igp,ilev-1)
        do i=1,ndim
        do j=1,ndim
          am(i,j)=(0.,0.)
          do k=1,n
          do l=1,n
             am(i,j)=am(i,j)+zl(ibi+i,k)*aa(k,l,ilev)*zr(l,ibi+j)
          enddo
          enddo
        enddo
        enddo
        call cgrl(nm,ndim,am,w(1+ibi,ilev),zrm,zlm,fv(1,1),ind1,ngroup1,ierr)

c       process the eigenvectors
        do j=1,ndim
          do i=1,n
            am(i,j)=(0.,0.)
            do k=1,ndim
              am(i,j)=am(i,j)+zr(i,k+ibi)*zrm(k,j)
            enddo
          enddo
        enddo
        do j=1,ndim
          do i=1,n
            zr(i,j+ibi)=am(i,j)
          enddo
        enddo

        do j=1,ndim
          do i=1,n
            am(j,i)=(0.,0.)
            do k=1,ndim
              am(j,i)=am(j,i)+zlm(j,k)*zl(k+ibi,i)
            enddo
          enddo
        enddo
        do j=1,ndim
          do i=1,n
            zl(j+ibi,i)=am(j,i)
          enddo
        enddo
c set up groups for the next level
        do i=1,ngroup1
          ind(ngp(ilev),ilev)=ind1(i-1)+ibi
          ngp(ilev)=ngp(ilev)+1
        enddo
        ind(ngp(ilev),ilev)=ind1(ngroup1)+ibi
      endif
      goto 20

   99 continue
       
      return
      end

