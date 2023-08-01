C$**********************************************************************
      SUBROUTINE FSVESQ(A,I,C,N)
      DIMENSION A(I*N)
      C=0.
      IK=1
      DO 1 K=1,N
      C=C+A(IK)*A(IK)
      IK=IK+I
    1 CONTINUE
      RETURN
      END
