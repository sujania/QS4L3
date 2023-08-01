      program mulAcst

c      double precision A(1000,5000,3),model(5000,3),cst(5000)
      double precision A(5000,5000),model(5000),cst(150)
      real dum1,dum2,dum41
      integer Ncst, Nmodel

      open(10,file='A.dat')
      open(20,file='../mnew.dat')
      open(30,file='cstnew.dat')

      write(6,*) 'Give number of cst coefficients:'
      read(5,*) Ncst !45

      write(6,*) 'Give number of model parameter:'
      read(5,*) Nmodel ! 1890

      do i=1,Ncst
         cst(i)=0.0
      enddo

      do j=1,Nmodel
         do i=1,Ncst
            read(10,*) A(i,j)
         enddo
      enddo

      do i=1,Nmodel
         read(20,*)model(i)
      enddo

      do i=1,Ncst
         do j=1,Nmodel
cs            write(*,*) 'mulAcst', cst(i), A(i,j), model(j)
            cst(i)=cst(i)+A(i,j)*model(j)
         enddo
         write(30,*) cst(i)
      enddo
      close(30)
      close(20)
      close(10)
      end
