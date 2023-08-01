
C$**********************************************************************
      SUBROUTINE FDOTPR(A,I,B,J,C,N)
      DIMENSION A(N*I),B(J*N)
      C=0.
      IK=1
      JK=1
      DO 1 K=1,N
      C=C+A(IK)*B(JK)
      IK=IK+I
      JK=JK+J
    1 CONTINUE
      RETURN
      END
