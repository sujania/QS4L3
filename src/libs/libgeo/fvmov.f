
C$**********************************************************************
      SUBROUTINE FVMOV(A,I,C,K,N)
      DIMENSION A(*),C(*)
      II=1
      KK=1
      DO 10 NN=1,N
      C(KK)=A(II)
      KK=KK+K
   10 II=II+I
      RETURN
      END
