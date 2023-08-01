      program bspl2v

      parameter (MXLH=40)
      parameter (MXLENY=(MXLH+1)**2)
      parameter (MXDEP=21)
      parameter (MXSPL=MXDEP+3)
      parameter (MXMSZ=MXSPL*MXLENY)
      parameter (MAXP=1024)
      dimension wk1(MXLENY),wk2(MXLENY),wk3(MXLENY)
      dimension d0(MXLENY)
      dimension spl(MXDEP)
      dimension x(MXMSZ)
      real*8 rknots(3), rho_x
      character*100 mfl
      character*80 getunx

c     Read model
      Write(6,*) 'Give mantle model (i.e. S20RTS.sph)'
      read(5,111) mfl
111   format(a100)
      call rmod(mfl,x,lmx,nspl)
      natd=(lmx+1)**2

      write(6,*) 'Give latitude, longitude, depth (in km)'
10    read(5,*,end=99) xlt,xln,dep

c    calculate Y_k(xlt,xln) at random depth level
      call ylm(xlt,xln,lmx,d0,wk1,wk2,wk3)

c    calculate spline coefficients
      do i=1,nspl
       ind=(i-1)*natd+1
       spl(i)=sdot(natd,d0,1,x(ind),1)
      enddo

      rcmb=3480.
      rmoho=6346.
      rearth=6371.
      r=rearth-dep
      xd = r/rearth

      rknots(3)=24.0
      rknots(2)=1458.0
      rknots(1)=2892.0

      rknots(:) = (rearth - rknots(:)) / rearth ! normalize

      dv=0.
      do ip=1,nspl
        call bspl(ip-1,nspl,rknots(1:nspl),dble(xd),rho_x)
        dv=dv+rho_x*spl(ip)
      enddo

      write(6,112) xlt,xln,dep,dv*1000. !dqmu*1000
112   format(f9.2,f9.2,f9.2,f9.5)
      goto 10

99    continue

      end


      subroutine rmod(infl,x,smax,nspl)

      integer smax, i, j, k, nspl
      character*100 infl
      dimension x(*)
      dimension dum(2000)

      open(10,file=infl,status='old')
      read(10,*) smax,nspl
      natd=(smax+1)**2
      ind=1
      do k=1,nspl
       do i=0,smax
        ind1=ind+2*i
        if(i.eq.0) then
          x(ind)=0.000000
          read(10,*,end=100) (dum(j),j=ind,ind1)
        else 
          read(10,*,end=100) (x(j),j=ind,ind1)
        endif
        ind=ind1+1
       enddo
      enddo

      goto 200

 100  stop 'incompatible sph header'

 200  continue
c      smax=4
      end


c--------------------------------
c       B-SPLINES SHAPE FUNCTIONS
c--------------------------------
! Returns the ord'th spline's value defined at nknots points provided in
! rknots array at xi.

        subroutine bspl(ord, nknots, knot, xi, rho_x)

        integer, intent(in)            :: ord, nknots
        integer, parameter :: NKNOTS_MAX = 1000
        double precision, parameter    :: TOL = 1.5
        double precision, dimension(0:nknots-1), intent(in) :: knot
        double precision, intent(in)   :: xi
        double precision, intent(out)  :: rho_x

        integer             ::  ii, Nx
        double precision    ::  coefa, coefb, coefc, coefd, denom
        double precision    :: denomsum, dd, denom1, denom2
        double precision, dimension(0:NKNOTS_MAX) :: hh

        Nx = nknots - 1
        !! Compute vector hh of spacings */
        !! hh = farray1(0, Nx - 1); */
        call fill_hh(hh, knot, Nx)

        !! Consistency checks */
        if ((xi - TOL) .gt. knot(Nx)) then
            print*,"xi=",xi," / knot(",Nx,")=",knot(Nx)
            stop
            !return 0.0;
        else if ((xi + TOL) .lt. knot(0)) then
            print*,"xi=",xi," / knot(0)=", knot(0)
            stop
            !return 0.0;
        else if (ord .gt. Nx) then
            print*,"Warning: spl index ",ord," exceeds knot count", Nx
            stop
            !return 0.0;
        endif

        if (ord .eq. 0) then !          ! LHS */
            denom = 3. * hh(ord) * hh(ord) + 3. * hh(ord) * hh(ord + 1)
     1              + hh(ord + 1) * hh(ord + 1)
            if (xi .ge. knot(ord) .and. xi .le. knot(ord + 1)) then !!x0.le.x.le.x1 */
                coefa = 4. / (hh(ord) * (hh(ord) + hh(ord + 1)) * denom)
                coefb = 0.0
                coefc = -12 / denom
                coefd = 4 * (2 * hh(ord) + hh(ord + 1)) / denom
                rho_x = coefa * (xi - knot(ord)) * (xi - knot(ord)) *
     1                  (xi - knot(ord))
                rho_x = rho_x + coefb * (xi - knot(ord)) *
     1                  (xi - knot(ord))
                rho_x = rho_x + coefc * (xi - knot(ord))
                rho_x = rho_x + coefd
            else if (xi .gt. knot(ord + 1) .and. xi .le. knot(ord + 2)) then !  ! x1.le.x.le.x2 */
                coefa = -4. / (hh(ord + 1) * (hh(ord) + hh(ord + 1)) *
     1                  denom)
                coefb = 12 / ((hh(ord) + hh(ord + 1)) * denom)
                coefc = -12. * hh(ord + 1) / ((hh(ord) + hh(ord + 1)) *
     1                  denom)
                coefd = 4. * hh(ord + 1) * hh(ord + 1) / ((hh(ord) +
     1                  hh(ord + 1)) * denom)

                rho_x = coefa * (xi - knot(ord + 1)) * (xi - knot(ord +
     1                  1)) * (xi - knot(ord + 1))
                rho_x = rho_x + coefb * (xi - knot(ord + 1)) *
     1                  (xi - knot(ord + 1))
                rho_x = rho_x + coefc * (xi - knot(ord + 1))
                rho_x = rho_x + coefd
            else    ! x.gt.x2 */
                rho_x = 0.0
            endif
        else if (ord .eq. 1) then !     ! LHS+1 */
            denom = (3. * hh(ord - 1) * hh(ord - 1) + 4. * hh(ord - 1) *
     1               hh(ord) + hh(ord) * hh(ord) + 2. * hh(ord - 1) *
     2               hh(ord + 1) + hh(ord) * hh(ord + 1))
            denomsum = hh(ord - 1) + hh(ord) + hh(ord + 1)
            dd = denomsum * denom
            if (xi .ge. knot(ord - 1) .and. xi .le. knot(ord)) then !!x0.le.x.le.x1 */
                coefa = -4. * (3. * hh(ord - 1) + 2. * hh(ord) + hh(ord
     1                 +1))/(hh(ord - 1) * (hh(ord - 1) + hh(ord)) * dd)
                coefb = 0.
                coefc = 12. / denom
                coefd = 0.

                rho_x = coefa * (xi - knot(ord - 1)) * (xi - knot(ord -
     1                  1)) * (xi - knot(ord - 1))
                rho_x = rho_x + coefb * (xi - knot(ord - 1)) * (xi -
     1                  knot(ord - 1))
                rho_x = rho_x + coefc * (xi - knot(ord - 1))
                rho_x = rho_x + coefd
            else if (xi .ge. knot(ord) .and. xi .le. knot(ord + 1)) then!  ! x1.le.x.le.x2 */
                coefa = 4. * (2. * hh(ord - 1) * hh(ord - 1) + 6. *
     1                  hh(ord - 1) * hh(ord) + 3. * hh(ord) * hh(ord) +
     2                  3. * hh(ord - 1) * hh(ord + 1) + 3. * hh(ord) *
     3                  hh(ord + 1) + hh(ord + 1) * hh(ord + 1)) /
     4                  (hh(ord) * (hh(ord - 1) + hh(ord)) * (hh(ord) +
     5                  hh(ord + 1)) * dd)
                coefb = -12. * (3. * hh(ord - 1) + 2. * hh(ord) + hh(ord
     1                  + 1)) / ((hh(ord - 1) + hh(ord)) * dd)
                coefc = 12. * (-2. * hh(ord - 1) * hh(ord - 1) + hh(ord)
     1                  * hh(ord) + hh(ord) * hh(ord + 1)) /
     2                  ((hh(ord - 1) + hh(ord)) * dd)
                coefd = 4. * hh(ord - 1) * (4. * hh(ord - 1) * hh(ord) +
     1                  3. * hh(ord) * hh(ord) + 2. * hh(ord - 1) *
     2                  hh(ord + 1) + 3. * hh(ord) * hh(ord + 1)) /
     3                  ((hh(ord - 1) + hh(ord)) * dd)

                rho_x = coefa * (xi - knot(ord)) * (xi - knot(ord)) *
     1                  (xi - knot(ord))
                rho_x = rho_x +  coefb * (xi - knot(ord)) * (xi -
     1                  knot(ord))
                rho_x = rho_x +  coefc * (xi - knot(ord))
                rho_x = rho_x +  coefd
            else if (xi .ge. knot(ord + 1) .and. xi .le. knot(ord + 2)) then !  ! x2.le.x.le.x3 */
                dd = dd *  (hh(ord) + hh(ord + 1))
                coefa = -4. * (2. * hh(ord - 1) + hh(ord)) / (hh(ord +
     1                  1) * dd)
                coefb = 12. * (2. * hh(ord - 1) + hh(ord)) / dd
                coefc = -12. * (2. * hh(ord - 1) + hh(ord)) * hh(ord +
     1                  1) / dd
                coefd = 4. * (2. * hh(ord - 1) + hh(ord)) * hh(ord + 1)
     1                  * hh(ord + 1) / dd

                rho_x = coefa * (xi - knot(ord + 1)) * (xi - knot(ord +
     1                  1)) * (xi - knot(ord + 1))
                rho_x = rho_x +  coefb * (xi - knot(ord + 1)) * (xi -
     1                  knot(ord + 1))
                rho_x = rho_x +  coefc * (xi - knot(ord + 1))
                rho_x = rho_x +  coefd
            else ! x.gt.x3 */
                rho_x = 0.0
            endif
        else if (ord .eq. Nx - 1) then !        ! RHS-1 */
            denom = hh(ord - 2) * hh(ord - 1) + hh(ord - 1) * hh(ord -
     1              1) + 2. * hh(ord - 2) * hh(ord) + 4. * hh(ord - 1) *
     2              hh(ord) + 3. * hh(ord) * hh(ord)
            denomsum = hh(ord - 2) + hh(ord - 1) + hh(ord)
            dd = denomsum * denom
            if (xi .ge. knot(ord - 2) .and. xi .le. knot(ord - 1)) then !       ! x0.le.x.le.x1 */
                coefa = 4. * (hh(ord - 1) + 2. * hh(ord)) / (hh(ord - 2)
     1                  * (hh(ord - 2) + hh(ord - 1)) * dd)
                coefb = 0.0
                coefc = 0.0
                coefd = 0.0

                rho_x = coefa * (xi - knot(ord - 2)) * (xi - knot(ord -
     1                  2)) * (xi - knot(ord - 2))
                rho_x = rho_x +  coefb * (xi - knot(ord - 2)) * (xi -
     1                  knot(ord - 2))
                rho_x = rho_x +  coefc * (xi - knot(ord - 2))
                rho_x = rho_x +  coefd
            else if (xi .ge. knot(ord - 1) .and. xi .le. knot(ord)) then!       ! x1.le.x.le.x2 */
                coefa = -4. * (hh(ord - 2) * hh(ord - 2) + 3. * hh(ord -
     1                  2) * hh(ord - 1) + 3. * hh(ord - 1) * hh(ord -
     2                  1) + 3. * hh(ord - 2) * hh(ord) + 6. * hh(ord -
     3                  1) * hh(ord) + 2. * hh(ord) * hh(ord)) /
     4                  (hh(ord - 1) * (hh(ord - 2) + hh(ord - 1)) *
     5                  (hh(ord - 1) + hh(ord)) * dd)
                coefb = 12. * (hh(ord - 1) + 2. * hh(ord)) / ((hh(ord -
     1                  2) + hh(ord - 1)) * dd)
                coefc = 12. * hh(ord - 2) * (hh(ord - 1) + 2. * hh(ord))
     1                  / ((hh(ord - 2) + hh(ord - 1)) * dd)
                coefd = 4. * hh(ord - 2) * hh(ord - 2) * (hh(ord - 1) +
     1                  2. * hh(ord))/((hh(ord - 2) + hh(ord - 1)) * dd)

                rho_x = coefa * (xi - knot(ord - 1)) * (xi - knot(ord -
     1                  1)) * (xi - knot(ord - 1))
                rho_x = rho_x +  coefb * (xi - knot(ord - 1)) * (xi -
     1                  knot(ord - 1))
                rho_x = rho_x +  coefc * (xi - knot(ord - 1))
                rho_x = rho_x +  coefd
            else if (xi .ge. knot(ord) .and. xi .le. knot(ord + 1)) then!       ! x2.le.x.le.x3 */
                dd = dd *  (hh(ord - 1) + hh(ord))
                coefa = 4. * (hh(ord - 2) + 2. * hh(ord - 1) + 3. *
     1                  hh(ord)) / (hh(ord) * dd)
                coefb = -12. * (hh(ord - 2) + 2. * hh(ord - 1) + 3. *
     1                  hh(ord)) / dd
                coefc = 12. * (-hh(ord - 2) * hh(ord - 1) - hh(ord - 1)
     1                  * hh(ord - 1) + 2. * hh(ord) * hh(ord)) / dd
                coefd = 4. * hh(ord) * (3. * hh(ord - 2) * hh(ord - 1) +
     1                  3. * hh(ord - 1) * hh(ord - 1) + 2. * hh(ord -
     2                  2) * hh(ord) + 4. * hh(ord - 1) * hh(ord)) / dd

                rho_x = coefa * (xi - knot(ord)) * (xi - knot(ord)) *
     1                  (xi - knot(ord))
                rho_x = rho_x +  coefb * (xi - knot(ord)) * (xi -
     1                  knot(ord))
                rho_x = rho_x +  coefc * (xi - knot(ord))
                rho_x = rho_x +  coefd
            else! x.gt.x4 */
                rho_x = 0.0
            endif
        else if (ord .eq. Nx) then !    ! RHS */
            denom = (hh(ord - 2) + hh(ord - 1)) * (hh(ord - 2) * hh(ord
     1              - 2) + 3. * hh(ord - 2) * hh(ord - 1) + 3. * hh(ord
     2              - 1) * hh(ord - 1))
            if (xi .ge. knot(ord - 2) .and. xi .le. knot(ord - 1)) then !! x0.le.x.le.x1 */
                coefa = 4. / (hh(ord - 2) * denom)
                coefb = 0.0
                coefc = 0.0
                coefd = 0.0
                rho_x = coefa * (xi - knot(ord - 2)) * (xi -
     1                  knot(ord - 2)) * (xi - knot(ord - 2))
                rho_x = rho_x +  coefb * (xi - knot(ord - 2)) * (xi -
     1                  knot(ord - 2))
                rho_x = rho_x +  coefc * (xi - knot(ord - 2))
                rho_x = rho_x +  coefd
            else if (xi .ge. knot(ord - 1) .and. xi .le. knot(ord)) then!       ! x1.le.x.le.x2 */
                coefa = -4. / (hh(ord - 1) * denom)
                coefb = 12 / denom
                coefc = 12 * hh(ord - 2) / denom
                coefd = 4. * hh(ord - 2) * hh(ord - 2) / denom

                rho_x = coefa * (xi - knot(ord - 1)) * (xi -
     1                  knot(ord - 1)) * (xi - knot(ord - 1))
                rho_x = rho_x +  coefb * (xi - knot(ord - 1)) * (xi -
     1                  knot(ord - 1))
                rho_x = rho_x +  coefc * (xi - knot(ord - 1))
                rho_x = rho_x +  coefd
            else! x.gt.x2 */
                rho_x = 0.0
            endif
        else !          ! Away from borders */
            denom1 = hh(ord - 2) + hh(ord - 1) + hh(ord) + hh(ord + 1)
            if (xi .ge. knot(ord - 2) .and. xi .le. knot(ord - 1)) then !       ! x0.le.x.le.x1 */
                coefa = 4. / (hh(ord - 2) * (hh(ord - 2) + hh(ord - 1))*
     1                  (hh(ord - 2) + hh(ord - 1) + hh(ord)) * denom1)
                coefb = 0.0
                coefc = 0.0
                coefd = 0.0
                rho_x = coefa * (xi - knot(ord - 2)) * (xi -
     1                  knot(ord - 2)) * (xi - knot(ord - 2))
                rho_x = rho_x +  coefb * (xi - knot(ord - 2)) * (xi -
     1                  knot(ord - 2))
                rho_x = rho_x +  coefc * (xi - knot(ord - 2))
                rho_x = rho_x +  coefd
            else if (xi .ge. knot(ord - 1) .and. xi .le. knot(ord)) then!       ! x1.le.x.le.x2 */
                denom2 = (hh(ord - 2) + hh(ord - 1)) * (hh(ord - 2) +
     1                   hh(ord - 1) + hh(ord))
                denom = denom1 * denom2

                coefa = -4. * (hh(ord - 2) * hh(ord - 2) + 3. *
     1                  hh(ord - 2) * hh(ord - 1) + 3. * hh(ord - 1) *
     2                  hh(ord - 1) + 2. * hh(ord - 2) * hh(ord) +  4. *
     3                  hh(ord - 1) * hh(ord) + hh(ord) * hh(ord) +
     4                  hh(ord - 2) * hh(ord + 1) + 2. * hh(ord - 1) *
     5                  hh(ord + 1) + hh(ord) * hh(ord + 1)) /
     6                  (hh(ord - 1) * (hh(ord - 1) + hh(ord)) *
     7                  (hh(ord - 1) + hh(ord) + hh(ord + 1)) * denom)
                coefb = 12. / denom
                coefc = 12. * hh(ord - 2) / denom
                coefd = 4. * hh(ord - 2) * hh(ord - 2) / denom

                rho_x = coefa * (xi - knot(ord - 1)) * (xi -
     1                  knot(ord - 1)) * (xi - knot(ord - 1))
                rho_x = rho_x +  coefb * (xi - knot(ord - 1)) * (xi -
     1                  knot(ord - 1))
                rho_x = rho_x +  coefc * (xi - knot(ord - 1))
                rho_x = rho_x +  coefd
            else if (xi .ge. knot(ord) .and. xi .le. knot(ord + 1)) then!       ! x2.le.x.le.x3 */
                denom2 = (hh(ord - 1) + hh(ord)) * (hh(ord - 2) +
     1                   hh(ord - 1) + hh(ord)) * (hh(ord - 1) + hh(ord)
     2                   + hh(ord + 1))
                denom = denom1 * denom2
                coefa = 4. * (hh(ord - 2) * hh(ord - 1) + hh(ord - 1) *
     1                  hh(ord - 1) + 2. * hh(ord - 2) * hh(ord) + 4. *
     2                  hh(ord - 1) * hh(ord) + 3. * hh(ord) * hh(ord) +
     3                  hh(ord - 2) * hh(ord + 1) + 2. * hh(ord - 1) *
     4                  hh(ord + 1) + 3. * hh(ord) * hh(ord + 1) +
     5                  hh(ord + 1) * hh(ord + 1)) / (hh(ord) *
     6                  (hh(ord) + hh(ord + 1)) * denom)
                coefb = -12. * (hh(ord - 2) + 2. * hh(ord - 1) + 2. *
     1                  hh(ord) + hh(ord + 1)) / denom
                coefc = 12. * (-hh(ord - 2) * hh(ord - 1) - hh(ord - 1)
     1                  * hh(ord - 1) + hh(ord) * hh(ord) + hh(ord) *
     2                  hh(ord + 1)) / denom
                coefd = 4. * (2. * hh(ord - 2) * hh(ord - 1) * hh(ord) +
     1                  2. * hh(ord - 1) * hh(ord - 1) * hh(ord) +
     2                  hh(ord - 2) * hh(ord) * hh(ord) + 2. *
     3                  hh(ord - 1) * hh(ord) * hh(ord) + hh(ord - 2) *
     4                  hh(ord - 1) * hh(ord + 1) + hh(ord - 1) *
     5                  hh(ord - 1) * hh(ord + 1) + hh(ord - 2) *
     6                  hh(ord) * hh(ord + 1) + 2. * hh(ord - 1) *
     7                  hh(ord) * hh(ord + 1)) / denom
                rho_x = coefa * (xi - knot(ord)) * (xi - knot(ord)) *
     1                  (xi - knot(ord))
                rho_x = rho_x +  coefb * (xi - knot(ord)) * (xi -
     1                  knot(ord))
                rho_x = rho_x +  coefc * (xi - knot(ord))
                rho_x = rho_x +  coefd
            else if (xi .ge. knot(ord + 1) .and. xi .le. knot(ord + 2)) then !  ! x3.le.x.le.x4 */
                denom2 = (hh(ord) + hh(ord + 1)) * (hh(ord - 1) +
     1                   hh(ord) + hh(ord + 1))
                denom = denom1 * denom2

                coefa = -4. / (hh(ord + 1) * denom)
                coefb = 12 / denom
                coefc = -12 * hh(ord + 1) / denom
                coefd = 4. * hh(ord + 1) * hh(ord + 1) / denom

                rho_x = coefa * (xi - knot(ord + 1)) * (xi -
     1                  knot(ord + 1)) * (xi - knot(ord + 1))
                rho_x = rho_x +  coefb * (xi - knot(ord + 1)) * (xi -
     1                  knot(ord + 1))
                rho_x = rho_x +  coefc * (xi - knot(ord + 1))
                rho_x = rho_x +  coefd
            else! x.gt.x4 */
                rho_x = 0.0
            endif
        endif

        return
        end

c-----------------------------------
c       DISTANCES BETWEEN KNOT RADII
c-----------------------------------
!       Compute the distance between the b-spline knot radii

        subroutine fill_hh(hh, knot, Nx)

        integer, intent(in)                       :: Nx
        integer, parameter :: NKNOTS_MAX = 1000
        double precision, dimension(Nx), intent(in) :: knot
        double precision, dimension(1:NKNOTS_MAX), intent(out):: hh

        !LOCAL VARIABLES
        integer     :: ii

        do ii = 1, Nx
            hh(ii) = knot(ii + 1) - knot(ii)
        enddo

        return
        end


      subroutine ylm(xlat,xlon,lmax,y,wk1,wk2,wk3)
c
      complex temp,fac,dfac
      dimension wk1(1),wk2(1),wk3(1),y(1)
c
c     wk1,wk2,wk3 should be dimensioned at least (lmax+1)*4
c
      data radian/57.2957795/    ! 360./2pi
c
c     transform to spherical coordinates
      theta=(90.-xlat)/radian
      phi=xlon/radian
c
c    loop over l values
      ind=0
      lm1=lmax+1
      do 10 il1=1,lm1
      l=il1-1
      call legndr(theta,l,l,wk1,wk2,wk3)
c
      fac=(1.,0.)
      dfac=cexp(cmplx(0.,phi))
c
c    loop over m values
      do 20 im=1,il1
      temp=fac*cmplx(wk1(im),0.)
      ind=ind+1
      y(ind)=real(temp)
      if(im.eq.1) goto 20
      ind=ind+1
      y(ind)=aimag(temp)
   20 fac=fac*dfac   ! calculates exp(im phi)
c
   10 continue
      return
      end

      SUBROUTINE LEGNDR(THETA,L,M,X,XP,XCOSEC)
      DIMENSION X(*),XP(*),XCOSEC(*)
      DOUBLE PRECISION SMALL,SUM,COMPAR,CT,ST,FCT,COT,FPI,X1,X2,X3,
     1F1,F2,XM,TH,DFLOAT
      DATA FPI/12.56637062D0/
      DFLOAT(I)=FLOAT(I)
      SUM=0.D0
      LP1=L+1
      TH=THETA
      CT=DCOS(TH)
      ST=DSIN(TH)
      MP1=M+1
      FCT=DSQRT(DFLOAT(2*L+1)/FPI)
      SFL3=SQRT(FLOAT(L*(L+1)))
      COMPAR=DFLOAT(2*L+1)/FPI
      DSFL3=SFL3
      SMALL=1.D-16*COMPAR
      DO 1 I=1,MP1
      X(I)=0.
      XCOSEC(I)=0.
    1 XP(I)=0.
      IF(L.GT.1.AND.ABS(THETA).GT.1.E-5) GO TO 3
      X(1)=FCT
      IF(L.EQ.0) RETURN
      X(1)=CT*FCT
      X(2)=-ST*FCT/DSFL3
      XP(1)=-ST*FCT
      XP(2)=-.5D0*CT*FCT*DSFL3
      IF(ABS(THETA).LT.1.E-5) XCOSEC(2)=XP(2)
      IF(ABS(THETA).GE.1.E-5) XCOSEC(2)=X(2)/ST
      RETURN
    3 X1=1.D0
      X2=CT
      DO 4 I=2,L
      X3=(DFLOAT(2*I-1)*CT*X2-DFLOAT(I-1)*X1)/DFLOAT(I)
      X1=X2
    4 X2=X3
      COT=CT/ST
      COSEC=1./ST
      X3=X2*FCT
      X2=DFLOAT(L)*(X1-CT*X2)*FCT/ST
      X(1)=X3
      X(2)=X2
      SUM=X3*X3
      XP(1)=-X2
      XP(2)=DFLOAT(L*(L+1))*X3-COT*X2
      X(2)=-X(2)/SFL3
      XCOSEC(2)=X(2)*COSEC
      XP(2)=-XP(2)/SFL3
      SUM=SUM+2.D0*X(2)*X(2)
      IF(SUM-COMPAR.GT.SMALL) RETURN
      X1=X3
      X2=-X2/DSQRT(DFLOAT(L*(L+1)))
      DO 5 I=3,MP1
      K=I-1
      F1=DSQRT(DFLOAT((L+I-1)*(L-I+2)))
      F2=DSQRT(DFLOAT((L+I-2)*(L-I+3)))
      XM=K
      X3=-(2.D0*COT*(XM-1.D0)*X2+F2*X1)/F1
      SUM=SUM+2.D0*X3*X3
      IF(SUM-COMPAR.GT.SMALL.AND.I.NE.LP1) RETURN
      X(I)=X3
      XCOSEC(I)=X(I)*COSEC
      X1=X2
      XP(I)=-(F1*X2+XM*COT*X3)
    5 X2=X3
      RETURN
      END

      real function sdot(n,sx,incx,sy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      real sx(*),sy(*),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      stemp = 0.0e0
      sdot = 0.0e0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      sdot = stemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)*sy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
     *   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
   50 continue
   60 sdot = stemp
      return
      end
