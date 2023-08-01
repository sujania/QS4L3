      parameter (MDIM=4)
      parameter (NGEN=15)
      complex gen(MDIM,MDIM,NGEN),c1(MDIM,MDIM)
      complex ci,dotmat,val,cwk(NGEN)
      dimension aar(NGEN,NGEN,MDIM),aai(NGEN,NGEN,MDIM)
      dimension wr(NGEN,MDIM),wi(NGEN,MDIM),zr(NGEN,NGEN),zi(NGEN,NGEN)
     1     ,fv1(NGEN),fv2(NGEN),fv3(NGEN),index(NGEN,2)
      ci=(0.,1.)
      do k=1,NGEN
        do j=1,MDIM
          do i=1,MDIM
            gen(i,j,k)=(0.d0,0.d0)
          enddo
        enddo
      enddo

c     define the diagonal generators
      n=1
      gen(1,1,n)=1.
      gen(2,2,n)=1.
      gen(3,3,n)=-1.
      gen(4,4,n)=-1.
      n=2
      gen(1,1,n)=1.
      gen(2,2,n)=-1.
      gen(3,3,n)=-1.
      gen(4,4,n)=1.
      n=3
      gen(1,1,n)=ci
      gen(2,2,n)=-ci
      gen(3,3,n)=ci
      gen(4,4,n)=-ci
      ncom=n
c     define the other generators
      nh=MDIM/2
      do i=1,nh
        do j=1,nh
          if(i.ne.j) then
            n=n+1
            gen(i,j,n)=1.
            gen(j+nh,i+nh,n)=-1.
            n=n+1
            gen(i,j,n)=ci
            gen(j+nh,i+nh,n)=ci
          endif
        enddo
      enddo
      do i=1,nh
        do j=i,nh
          if(i.eq.j) then
            n=n+1
            gen(i+nh,j,n)=1.
            n=n+1
            gen(i,j+nh,n)=1.
          else
            n=n+1
            gen(i+nh,j,n)=1.
            gen(j+nh,i,n)=1.
            n=n+1
            gen(i+nh,j,n)=ci
            gen(j+nh,i,n)=-ci
            n=n+1
            gen(j,i+nh,n)=1.
            gen(i,j+nh,n)=1.
            n=n+1
            gen(j,i+nh,n)=ci
            gen(i,j+nh,n)=-ci
          endif
        enddo
      enddo
      do i=1,n
        write(6,'(''Matrix '',i3,''.'')') i
        do k=1,MDIM
          write(6,'(10(''('',2f6.2,'')'',3x))') (gen(k,j,i),j=1,4)
        enddo
      enddo

      do icom=1,ncom
        do j=1,n
          call commut(gen(1,1,icom),MDIM,gen(1,1,j),MDIM,MDIM,c1,MDIM)
          do k=1,n
            val=dotmat(c1,MDIM,gen(1,1,k),MDIM,MDIM)
            aar(k,j,icom)=real(val)
            aai(k,j,icom)=aimag(val)
          enddo
        enddo
check it's real
        vmax=0.
        do j=1,n
          do k=1,n
            vmax=amax1(vmax,aai(k,j,icom))
          enddo
        enddo
        write(6,*) 'Reality check:',vmax
      enddo


c     do icom=1,ncom
c       write(6,'(''Matrix '',i3,''.'')') icom
c       do k=1,n
c         write(6,'(15(''('',2f5.1,'')'',3x))') (aar(k,j,icom),aai(k,j,icom),j=1,n)
c       enddo
c     enddo

      call simeig(NGEN,n,aar,aai,wr,wi,zr,zi,fv1,fv2,fv3,index,cwk,ncom)



      end
c----------------------------------------------------------------
      subroutine simeig(nm,n,aar,aai,wr,wi,zr,zi,fv1,fv2,fv3,index,cwk,ncom)

      parameter (NGEN=15)

      dimension savr(NGEN,NGEN),savi(NGEN,NGEN)
      dimension sa1r(NGEN,NGEN),sa1i(NGEN,NGEN)
      dimension sa2r(NGEN,NGEN),sa2i(NGEN,NGEN)
      dimension sa3r(NGEN,NGEN),sa3i(NGEN,NGEN)
      dimension zzr(NGEN,NGEN),zzi(NGEN,NGEN)
      dimension amr(NGEN,NGEN),ami(NGEN,NGEN)
      dimension ind(NGEN),w1r(NGEN),w1i(NGEN),w2r(NGEN),w2i(NGEN),istart(NGEN)

      dimension aar(nm,nm,*),aai(nm,nm,*),wr(nm,*),wi(nm,*)
     1     ,zr(nm,nm),zi(nm,nm),fv1(nm),fv2(nm),fv3(nm),index(nm,2)
      complex cwk(*)
      ilev=0
      igroup=0
  200 continue
      igroup=igroup+1
      if(igroup.gt.ngroup) then
        ilev=ilev+1
        igroup=1
      endif

      if(ilev.eq.1) then
        ndim=n
        do i=1,ndim
        do j=1,ndim
          amr(i,j)=aar(i,j,ilev)
          ami(i,j)=aai(i,j,ilev)
        enddo
        enddo
      else
        ndim=istart(igroup+1)-istart(igroup)
        do i=1,ndim
        ii=istart(igroup)+i-1
        do j=1,ndim
          jj=istart(igroup)+j-1
          amr(i,j)=0.
          ami(i,j)=0.
          do k=1,n
          do l=1,n
             tr=aar(k,l,ilev)*zr(l,jj)-aai(k,l,ilev)*zi(l,jj)
             ti=aar(k,l,ilev)*zi(l,jj)+aai(k,l,ilev)*zr(l,jj)
             amr(i,j)=amr(i,j)+zzr(ii,k)*tr-zzi(ii,k)*ti
             ami(i,j)=ami(i,j)+zzr(ii,k)*ti+zzi(ii,k)*tr
          enddo
          enddo
        enddo
        enddo
      endif


      do i=1,ndim
      do j=1,ndim
        sa1r(i,j)=aar(i,j,ilev)
        sa1i(i,j)=aai(i,j,ilev)
      enddo
      enddo
c find right eigenvectors (zr,zi)
      call cg(nm,ndim,sa1r(1,1),sa1i(1,1),w1r(1),w1i(1),1
     1    ,zr(1,1),zi(1,1),fv1,fv2,fv3,ierr)
c group and order the eigenvalues and eigenvectors
      tol=1.e-5
      do i=1,ndim
        w2r(i)=w1r(i)
        w2i(i)=w1i(i)
      enddo
      do i=1,ndim
        do j=1,ndim
          sa2r(i,j)=zr(i,j)
          sa2i(i,j)=zi(i,j)
        enddo
      enddo
      do i=1,ndim
        ind(i)=0
      enddo
      mput=0
  100 imin=0
      do i=1,ndim
        if(ind(i).eq.0) then
          if(imin.eq.0) then
            imin=i
            cmin=w2r(i)
          else
            if(w2r(i).lt.cmin) then
              imin=i
              cmin=w2r(i)
            endif
          endif
        endif
      enddo
      if(imin.ne.0) then
        mput=mput+1
        ind(imin)=mput
        ngroup=ngroup+1
        istart(ngroup)=mput
        istart(ngroup+1)=mput+1
        do i=1,ndim
          if( ind(i).eq.0.and.cabs(cmplx(w2r(i),w2i(i))-cmplx(w2r(imin),w2i(imin))).lt.tol) then
            mput=mput+1
            ind(i)=mput
            istart(ngroup+1)=mput+1
          endif
        enddo
        goto 100
      endif
      if(mput.ne.ndim) stop 'wrong mput value'
      do i=1,ndim
        mput=ind(i)
        w1r(mput)=w2r(i)
        w1i(mput)=w2i(i)
        do j=1,ndim
          zr(j,mput)=sa2r(j,i)
          zi(j,mput)=sa2i(j,i)
        enddo
      enddo
c find left eigenvectors (zzr,zzi)
      do i=1,ndim
      do j=1,ndim
        zzr(i,j)=zr(i,j)
        zzi(i,j)=zi(i,j)
      enddo
      enddo
      call cmatinv(zzr,zzi,ndim,nm,detr,deti,fv1,index,cwk)


      

       
      return
      end

c----------------------------------------------------------------
      subroutine cgrl(nm,n,a,w,zr,zl,fv,ierr)
      dimension a(2*nm*nm),w(2*nm),zr(nm*nm*2),zl(nm*nm*2),fv(4*nm)
c     a is assumed stored as a complex array, first rearrange it
c     as two real arrays.
      do i=1,nm
      do j=1,nm
      do k=1,2
        zr(i+nm*(j-1+nm*(k-1)))=a(k+2*(j-1+nm*(k-1)))
      enddo
      enddo
      enddo
      call scopy(nm*nm*2,zr,1,a,1)
      call cg(nm,n,a(1),a(1+nm*nm),w(1),w(1+nm),1
     1     ,zr(1),zr(1+nm*nm),fv(1),fv(1+nm),fv(1+nm*2),ierr)
      call scopy(nm*nm*2,zr,1,zl,1)
      call cmatinv(zl(1),zl(1+nm*nm),n,nm,detr,deti,fv(1),fv(1+nm),fv(1+2*nm))

      do i=1,nm
      do j=1,nm
      do k=1,2
        a(k+2*(j-1+nm*(k-1)))=zr(i+nm*(j-1+nm*(k-1)))
      enddo
      enddo
      enddo
      call scopy(nm*nm*2,a,1,zr,1)
      do i=1,nm
      do j=1,nm
      do k=1,2
        a(k+2*(j-1+nm*(k-1)))=zl(i+nm*(j-1+nm*(k-1)))
      enddo
      enddo
      enddo
      call scopy(nm*nm*2,a,1,zl,1)
      do i=1,nm
      do k=1,2
        fv(k+2*(i-1))=w(i+nm*(k-1))
      enddo
      enddo
      call scopy(nm*2,fv,1,w,1)
      return
      end



c----------------------------------------------------------------
      subroutine cmatinv(ar,ai,n,m,detr,deti,ipvot,index,pivot)
      implicit complex(a-h,o-z)
      real ar(m,m),ai(m,m),tr,ti,detr,deti,temp
      dimension index(m,2),pivot(m),ipvot(m)
      equivalence (irow,jrow),(icol,jcol)
      det=(1.0,0.0)
      do j=1,n
        ipvot(j)=0
      enddo
      do 135 i=1,n
c     following 12 statements for search for pivot element
      t=0.d0
      do 9 j=1,n
      if(ipvot(j)-1)13,9,13
   13 do 23 k=1,n
      if(ipvot(k)-1)43,23,81
   43 if(cabs(t)-cabs(cmplx(ar(j,k),ai(j,k))))83,23,23
   83 irow=j
      icol=k
      t=cmplx(ar(j,k),ai(j,k))
   23 continue
    9 continue
      ipvot(icol)=ipvot(icol)+1
c     following 15 statements to put pivot element on diagonal
      if(irow-icol)73,109,73
   73 det=-det
      do 12 l=1,n
      t=cmplx(ar(irow,l),ai(irow,l))
      ar(irow,l)=ar(icol,l)
      ai(irow,l)=ai(icol,l)
      ar(icol,l)=real(t)
   12 ai(icol,l)=aimag(t)
  109 index(i,1)=irow
      index(i,2)=icol
      pivot(i)=cmplx(ar(icol,icol),ai(icol,icol))
      det=det*pivot(i)
c     following 6 statements to divide pivot row by pivot element
      ar(icol,icol)=1.
      ai(icol,icol)=0.
      pivr=(1.,0.)/pivot(i)
      pivrr=real(pivr)
      pivri=aimag(pivr)
      do 205 l=1,n
        temp=ar(icol,l)*pivri+ai(icol,l)*pivrr
        ar(icol,l)=ar(icol,l)*pivrr-ai(icol,l)*pivri
  205   ai(icol,l)=temp
c      following 10 statements to reduce non-pivot rows
      do 135 li=1,n
      if(li-icol)21,135,21
   21 tr=ar(li,icol)
      ti=ai(li,icol)
      ar(li,icol)=0.
      ai(li,icol)=0.
      do 89 l=1,n
      temp=ai(li,l)-ar(icol,l)*ti-ai(icol,l)*tr
      ar(li,l)=ar(li,l)-ar(icol,l)*tr+ai(icol,l)*ti
   89 ai(li,l)=temp
  135 continue
c     following 11 statements to interchange columns
      do 3 i=1,n
      l=n-i+1
      if(index(l,1)-index(l,2))19,3,19
   19 jrow=index(l,1)
      jcol=index(l,2)
      do 549 k=1,n
      tr=ar(k,jrow)
      ti=ai(k,jrow)
      ar(k,jrow)=ar(k,jcol)
      ai(k,jrow)=ai(k,jcol)
      ai(k,jcol)=ti
      ar(k,jcol)=tr
  549 continue
    3 continue
   81 return
      end

c----------------------------------------------------------------
      subroutine commut(x,ix,y,iy,ndim,c,ic)
      complex x(ix,*),y(iy,*),c(ic,*)
      do i=1,ndim
      do j=1,ndim
        c(i,j)=(0.,0.)
        do k=1,ndim
          c(i,j)=c(i,j)+x(i,k)*y(k,j)-y(i,k)*x(k,j)
        enddo
      enddo
      enddo
      return
      end
c----------------------------------------------------------------
      complex function dotmat(a,ix,b,iy,ndim) ! coef for expansion of a
      complex a(ix,*),b(iy,*),sum,sum1
      sum=0.
      sum1=0.
      do j=1,ndim
      do i=1,ndim
        sum=sum+a(i,j)*conjg(b(i,j))
        sum1=sum1+b(i,j)*conjg(b(i,j))
      enddo
      enddo
      dotmat=sum/sum1
      return
      end

c----------------------------------------------------------------
      subroutine cg(nm,n,ar,ai,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)
c
      integer n,nm,is1,is2,ierr,matz
      real ar(nm,n),ai(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     x       fv1(n),fv2(n),fv3(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a complex general matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a=(ar,ai).
c
c        ar  and  ai  contain the real and imaginary parts,
c        respectively, of the complex general matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        wr  and  wi  contain the real and imaginary parts,
c        respectively, of the eigenvalues.
c
c        zr  and  zi  contain the real and imaginary parts,
c        respectively, of the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for comqr
c           and comqr2.  the normal completion code is zero.
c
c        fv1, fv2, and  fv3  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  cbal(nm,n,ar,ai,is1,is2,fv1)
      call  corth(nm,n,is1,is2,ar,ai,fv2,fv3)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  comqr(nm,n,is1,is2,ar,ai,wr,wi,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  comqr2(nm,n,is1,is2,fv2,fv3,ar,ai,wr,wi,zr,zi,ierr)
      if (ierr .ne. 0) go to 50
      call  cbabk2(nm,n,is1,is2,fv1,n,zr,zi)
   50 return
      end
c----------------------------------------------------------------
      subroutine cbal(nm,n,ar,ai,low,igh,scale)
c
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      real ar(nm,n),ai(nm,n),scale(n)
      real c,f,g,r,s,b2,radix
      logical noconv
c
c     this subroutine is a translation of the algol procedure
c     cbalance, which is a complex version of balance,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine balances a complex matrix and isolates
c     eigenvalues whenever possible.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex matrix to be balanced.
c
c     on output
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the balanced matrix.
c
c        low and igh are two integers such that ar(i,j) and ai(i,j)
c          are equal to zero if
c           (1) i is greater than j and
c           (2) j=1,...,low-1 or i=igh+1,...,n.
c
c        scale contains information determining the
c           permutations and scaling factors used.
c
c     suppose that the principal submatrix in rows low through igh
c     has been balanced, that p(j) denotes the index interchanged
c     with j during the permutation step, and that the elements
c     of the diagonal matrix used are denoted by d(i,j).  then
c        scale(j) = p(j),    for j = 1,...,low-1
c                 = d(j,j)       j = low,...,igh
c                 = p(j)         j = igh+1,...,n.
c     the order in which the interchanges are made is n to igh+1,
c     then 1 to low-1.
c
c     note that 1 is returned for igh if igh is zero formally.
c
c     the algol procedure exc contained in cbalance appears in
c     cbal  in line.  (note that the algol roles of identifiers
c     k,l have been reversed.)
c
c     arithmetic is real throughout.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      radix = 16.0e0
c
      b2 = radix * radix
      k = 1
      l = n
      go to 100
c     .......... in-line procedure for row and
c                column exchange ..........
   20 scale(m) = j
      if (j .eq. m) go to 50
c
      do 30 i = 1, l
         f = ar(i,j)
         ar(i,j) = ar(i,m)
         ar(i,m) = f
         f = ai(i,j)
         ai(i,j) = ai(i,m)
         ai(i,m) = f
   30 continue
c
      do 40 i = k, n
         f = ar(j,i)
         ar(j,i) = ar(m,i)
         ar(m,i) = f
         f = ai(j,i)
         ai(j,i) = ai(m,i)
         ai(m,i) = f
   40 continue
c
   50 go to (80,130), iexc
c     .......... search for rows isolating an eigenvalue
c                and push them down ..........
   80 if (l .eq. 1) go to 280
      l = l - 1
c     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj
c
         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (ar(j,i) .ne. 0.0e0 .or. ai(j,i) .ne. 0.0e0) go to 120
  110    continue
c
         m = l
         iexc = 1
         go to 20
  120 continue
c
      go to 140
c     .......... search for columns isolating an eigenvalue
c                and push them left ..........
  130 k = k + 1
c
  140 do 170 j = k, l
c
         do 150 i = k, l
            if (i .eq. j) go to 150
            if (ar(i,j) .ne. 0.0e0 .or. ai(i,j) .ne. 0.0e0) go to 170
  150    continue
c
         m = k
         iexc = 2
         go to 20
  170 continue
c     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1.0e0
c     .......... iterative loop for norm reduction ..........
  190 noconv = .false.
c
      do 270 i = k, l
         c = 0.0e0
         r = 0.0e0
c
         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + abs(ar(j,i)) + abs(ai(j,i))
            r = r + abs(ar(i,j)) + abs(ai(i,j))
  200    continue
c     .......... guard against zero c or r due to underflow ..........
         if (c .eq. 0.0e0 .or. r .eq. 0.0e0) go to 270
         g = r / radix
         f = 1.0e0
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
c     .......... now balance ..........
  240    if ((c + r) / f .ge. 0.95e0 * s) go to 270
         g = 1.0e0 / f
         scale(i) = scale(i) * f
         noconv = .true.
c
         do 250 j = k, n
            ar(i,j) = ar(i,j) * g
            ai(i,j) = ai(i,j) * g
  250    continue
c
         do 260 j = 1, l
            ar(j,i) = ar(j,i) * f
            ai(j,i) = ai(j,i) * f
  260    continue
c
  270 continue
c
      if (noconv) go to 190
c
  280 low = k
      igh = l
      return
      end
      subroutine cbabk2(nm,n,low,igh,scale,m,zr,zi)
c
      integer i,j,k,m,n,ii,nm,igh,low
      real scale(n),zr(nm,m),zi(nm,m)
      real s
c
c     this subroutine is a translation of the algol procedure
c     cbabk2, which is a complex version of balbak,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine forms the eigenvectors of a complex general
c     matrix by back transforming those of the corresponding
c     balanced matrix determined by  cbal.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by  cbal.
c
c        scale contains information determining the permutations
c          and scaling factors used by  cbal.
c
c        m is the number of eigenvectors to be back transformed.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors to be
c          back transformed in their first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
c
      do 110 i = low, igh
         s = scale(i)
c     .......... left hand eigenvectors are back transformed
c                if the foregoing statement is replaced by
c                s=1.0e0/scale(i). ..........
         do 100 j = 1, m
            zr(i,j) = zr(i,j) * s
            zi(i,j) = zi(i,j) * s
  100    continue
c
  110 continue
c     .......... for i=low-1 step -1 until 1,
c                igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
c
         do 130 j = 1, m
            s = zr(i,j)
            zr(i,j) = zr(k,j)
            zr(k,j) = s
            s = zi(i,j)
            zi(i,j) = zi(k,j)
            zi(k,j) = s
  130    continue
c
  140 continue
c
  200 return
      end
      subroutine comqr2(nm,n,low,igh,ortr,orti,hr,hi,wr,wi,zr,zi,ierr)
c
      integer i,j,k,l,m,n,en,ii,jj,ll,nm,nn,igh,ip1,
     x        itn,its,low,lp1,enm1,iend,ierr
      real hr(nm,n),hi(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     x       ortr(igh),orti(igh)
      real si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag
c
c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a complex upper hessenberg matrix by the qr
c     method.  the eigenvectors of a complex general matrix
c     can also be found if  corth  has been used to reduce
c     this general matrix to hessenberg form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        ortr and orti contain information about the unitary trans-
c          formations used in the reduction by  corth, if performed.
c          only elements low through igh are used.  if the eigenvectors
c          of the hessenberg matrix are desired, set ortr(j) and
c          orti(j) to 0.0e0 for these elements.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain further
c          information about the transformations which were used in the
c          reduction by  corth, if performed.  if the eigenvectors of
c          the hessenberg matrix are desired, these elements may be
c          arbitrary.
c
c     on output
c
c        ortr, orti, and the upper hessenberg portions of hr and hi
c          have been destroyed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors.  the eigenvectors
c          are unnormalized.  if an error exit is made, none of
c          the eigenvectors has been found.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  sqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
c     .......... initialize eigenvector matrix ..........
      do 101 j = 1, n
c
         do 100 i = 1, n
            zr(i,j) = 0.0e0
            zi(i,j) = 0.0e0
  100    continue
         zr(j,j) = 1.0e0
  101 continue
c     .......... form the matrix of accumulated transformations
c                from the information left by corth ..........
      iend = igh - low - 1
      if (iend) 180, 150, 105
c     .......... for i=igh-1 step -1 until low+1 do -- ..........
  105 do 140 ii = 1, iend
         i = igh - ii
         if (ortr(i) .eq. 0.0e0 .and. orti(i) .eq. 0.0e0) go to 140
         if (hr(i,i-1) .eq. 0.0e0 .and. hi(i,i-1) .eq. 0.0e0) go to 140
c     .......... norm below is negative of h formed in corth ..........
         norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)
         ip1 = i + 1
c
         do 110 k = ip1, igh
            ortr(k) = hr(k,i-1)
            orti(k) = hi(k,i-1)
  110    continue
c
         do 130 j = i, igh
            sr = 0.0e0
            si = 0.0e0
c
            do 115 k = i, igh
               sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
               si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
  115       continue
c
            sr = sr / norm
            si = si / norm
c
            do 120 k = i, igh
               zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
               zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
  120       continue
c
  130    continue
c
  140 continue
c     .......... create real subdiagonal elements ..........
  150 l = low + 1
c
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0e0) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0e0
c
         do 155 j = i, n
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
c
         do 160 j = 1, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
c
         do 165 j = low, igh
            si = yr * zi(j,i) + yi * zr(j,i)
            zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
            zi(j,i) = si
  165    continue
c
  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0e0
      ti = 0.0e0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 680
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1))
     x            + abs(hr(l,l)) + abs(hi(l,l))
         tst2 = tst1 + abs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0e0 .and. xi .eq. 0.0e0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0e0
      yi = (hi(enm1,enm1) - si) / 2.0e0
      call csroot(yr**2-yi**2+xr,2.0e0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0e0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
      si = 0.0e0
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
c
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0e0
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0e0
         hi(i,i-1) = sr / norm
c
         do 490 j = i, n
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
c
  500 continue
c
      si = hi(en,en)
      if (si .eq. 0.0e0) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0e0
      if (en .eq. n) go to 540
      ip1 = en + 1
c
      do 520 j = ip1, n
         yr = hr(en,j)
         yi = hi(en,j)
         hr(en,j) = sr * yr + si * yi
         hi(en,j) = sr * yi - si * yr
  520 continue
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
c
         do 580 i = 1, j
            yr = hr(i,j-1)
            yi = 0.0e0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
c
         do 590 i = low, igh
            yr = zr(i,j-1)
            yi = zi(i,j-1)
            zzr = zr(i,j)
            zzi = zi(i,j)
            zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  590    continue
c
  600 continue
c
      if (si .eq. 0.0e0) go to 240
c
      do 630 i = 1, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
c
      do 640 i = low, igh
         yr = zr(i,en)
         yi = zi(i,en)
         zr(i,en) = sr * yr - si * yi
         zi(i,en) = sr * yi + si * yr
  640 continue
c
      go to 240
c     .......... a root found ..........
  660 hr(en,en) = hr(en,en) + tr
      wr(en) = hr(en,en)
      hi(en,en) = hi(en,en) + ti
      wi(en) = hi(en,en)
      en = enm1
      go to 220
c     .......... all roots found.  backsubstitute to find
c                vectors of upper triangular form ..........
  680 norm = 0.0e0
c
      do 720 i = 1, n
c
         do 720 j = i, n
            tr = abs(hr(i,j)) + abs(hi(i,j))
            if (tr .gt. norm) norm = tr
  720 continue
c
      if (n .eq. 1 .or. norm .eq. 0.0e0) go to 1001
c     .......... for en=n step -1 until 2 do -- ..........
      do 800 nn = 2, n
         en = n + 2 - nn
         xr = wr(en)
         xi = wi(en)
         hr(en,en) = 1.0e0
         hi(en,en) = 0.0e0
         enm1 = en - 1
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 780 ii = 1, enm1
            i = en - ii
            zzr = 0.0e0
            zzi = 0.0e0
            ip1 = i + 1
c
            do 740 j = ip1, en
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
  740       continue
c
            yr = xr - wr(i)
            yi = xi - wi(i)
            if (yr .ne. 0.0e0 .or. yi .ne. 0.0e0) go to 765
               tst1 = norm
               yr = tst1
  760          yr = 0.01e0 * yr
               tst2 = norm + yr
               if (tst2 .gt. tst1) go to 760
  765       continue
            call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
c     .......... overflow control ..........
            tr = abs(hr(i,en)) + abs(hi(i,en))
            if (tr .eq. 0.0e0) go to 780
            tst1 = tr
            tst2 = tst1 + 1.0e0/tst1
            if (tst2 .gt. tst1) go to 780
            do 770 j = i, en
               hr(j,en) = hr(j,en)/tr
               hi(j,en) = hi(j,en)/tr
  770       continue
c
  780    continue
c
  800 continue
c     .......... end backsubstitution ..........
      enm1 = n - 1
c     .......... vectors of isolated roots ..........
      do  840 i = 1, enm1
         if (i .ge. low .and. i .le. igh) go to 840
         ip1 = i + 1
c
         do 820 j = ip1, n
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
  820    continue
c
  840 continue
c     .......... multiply by transformation matrix to give
c                vectors of original full matrix.
c                for j=n step -1 until low+1 do -- ..........
      do 880 jj = low, enm1
         j = n + low - jj
         m = min0(j,igh)
c
         do 880 i = low, igh
            zzr = 0.0e0
            zzi = 0.0e0
c
            do 860 k = low, m
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
  860       continue
c
            zr(i,j) = zzr
            zi(i,j) = zzi
  880 continue
c
      go to 1001
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine comqr(nm,n,low,igh,hr,hi,wr,wi,ierr)
c
      integer i,j,l,n,en,ll,nm,igh,itn,its,low,lp1,enm1,ierr
      real hr(nm,n),hi(nm,n),wr(n),wi(n)
      real si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag
c
c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr, num. math. 12, 369-376(1968) by martin
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
c
c     this subroutine finds the eigenvalues of a complex
c     upper hessenberg matrix by the qr method.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain
c          information about the unitary transformations used in
c          the reduction by  corth, if performed.
c
c     on output
c
c        the upper hessenberg portions of hr and hi have been
c          destroyed.  therefore, they must be saved before
c          calling  comqr  if subsequent calculation of
c          eigenvectors is to be performed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  sqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (low .eq. igh) go to 180
c     .......... create real subdiagonal elements ..........
      l = low + 1
c
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0e0) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0e0
c
         do 155 j = i, igh
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
c
         do 160 j = low, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
c
  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0e0
      ti = 0.0e0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 1001
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low e0 -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1))
     x            + abs(hr(l,l)) + abs(hi(l,l))
         tst2 = tst1 + abs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0e0 .and. xi .eq. 0.0e0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0e0
      yi = (hi(enm1,enm1) - si) / 2.0e0
      call csroot(yr**2-yi**2+xr,2.0e0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0e0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
      si = 0.0e0
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
c
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0e0
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0e0
         hi(i,i-1) = sr / norm
c
         do 490 j = i, en
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
c
  500 continue
c
      si = hi(en,en)
      if (si .eq. 0.0e0) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0e0
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
c
         do 580 i = l, j
            yr = hr(i,j-1)
            yi = 0.0e0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
c
  600 continue
c
      if (si .eq. 0.0e0) go to 240
c
      do 630 i = l, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
c
      go to 240
c     .......... a root found ..........
  660 wr(en) = hr(en,en) + tr
      wi(en) = hi(en,en) + ti
      en = enm1
      go to 220
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine corth(nm,n,low,igh,ar,ai,ortr,orti)
c
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      real ar(nm,n),ai(nm,n),ortr(igh),orti(igh)
      real f,g,h,fi,fr,scale,pythag
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure orthes, num. math. 12, 349-368(1968)
c     by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a complex general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     unitary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex input matrix.
c
c     on output
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the hessenberg matrix.  information
c          about the unitary transformations used in the reduction
c          is stored in the remaining triangles under the
c          hessenberg matrix.
c
c        ortr and orti contain further information about the
c          transformations.  only elements low through igh are used.
c
c     calls pythag for  sqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         h = 0.0e0
         ortr(m) = 0.0e0
         orti(m) = 0.0e0
         scale = 0.0e0
c     .......... scale column (algol tol then not needed) ..........
         do 90 i = m, igh
   90    scale = scale + abs(ar(i,m-1)) + abs(ai(i,m-1))
c
         if (scale .eq. 0.0e0) go to 180
         mp = m + igh
c     .......... for i=igh step -1 until m do -- ..........
         do 100 ii = m, igh
            i = mp - ii
            ortr(i) = ar(i,m-1) / scale
            orti(i) = ai(i,m-1) / scale
            h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
  100    continue
c
         g = sqrt(h)
         f = pythag(ortr(m),orti(m))
         if (f .eq. 0.0e0) go to 103
         h = h + f * g
         g = g / f
         ortr(m) = (1.0e0 + g) * ortr(m)
         orti(m) = (1.0e0 + g) * orti(m)
         go to 105
c
  103    ortr(m) = g
         ar(m,m-1) = scale
c     .......... form (i-(u*ut)/h) * a ..........
  105    do 130 j = m, n
            fr = 0.0e0
            fi = 0.0e0
c     .......... for i=igh step -1 until m do -- ..........
            do 110 ii = m, igh
               i = mp - ii
               fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
               fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
  110       continue
c
            fr = fr / h
            fi = fi / h
c
            do 120 i = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
               ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
  120       continue
c
  130    continue
c     .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........
         do 160 i = 1, igh
            fr = 0.0e0
            fi = 0.0e0
c     .......... for j=igh step -1 until m do -- ..........
            do 140 jj = m, igh
               j = mp - jj
               fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
               fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
  140       continue
c
            fr = fr / h
            fi = fi / h
c
            do 150 j = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
               ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
  150       continue
c
  160    continue
c
         ortr(m) = scale * ortr(m)
         orti(m) = scale * orti(m)
         ar(m,m-1) = -g * ar(m,m-1)
         ai(m,m-1) = -g * ai(m,m-1)
  180 continue
c
  200 return
      end
      subroutine cdiv(ar,ai,br,bi,cr,ci)
      real ar,ai,br,bi,cr,ci
c
c     complex division, (cr,ci) = (ar,ai)/(br,bi)
c
      real s,ars,ais,brs,bis
      s = abs(br) + abs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      end
      real function pythag(a,b)
      real a,b
c
c     finds sqrt(a**2+b**2) without overflow or destructive underflow
c
      real p,r,s,t,u
      p = amax1(abs(a),abs(b))
      if (p .eq. 0.0e0) go to 20
      r = (amin1(abs(a),abs(b))/p)**2
   10 continue
         t = 4.0e0 + r
         if (t .eq. 4.0e0) go to 20
         s = r/t
         u = 1.0e0 + 2.0e0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
      subroutine csroot(xr,xi,yr,yi)
      real xr,xi,yr,yi
c
c     (yr,yi) = complex sqrt(xr,xi)
c     branch chosen so that yr .ge. 0.0 and sign(yi) .eq. sign(xi)
c
      real s,tr,ti,pythag
      tr = xr
      ti = xi
      s = sqrt(0.5e0*(pythag(tr,ti) + abs(tr)))
      if (tr .ge. 0.0e0) yr = s
      if (ti .lt. 0.0e0) s = -s
      if (tr .le. 0.0e0) yi = s
      if (tr .lt. 0.0e0) yr = 0.5e0*(ti/yi)
      if (tr .gt. 0.0e0) yi = 0.5e0*(ti/yr)
      return
      end


