
C$**********************************************************************
      SUBROUTINE FVCLR(A,I,N)
      DIMENSION A(I*N)
      IK=1
      DO 1 K=1,N
      A(IK)=0.
      IK=IK+I
    1 CONTINUE
      RETURN
      END
