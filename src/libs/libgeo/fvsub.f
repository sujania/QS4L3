
C$**********************************************************************
      SUBROUTINE FVSUB(A,I,B,J,C,K,N)
      DIMENSION A(I*N),B(J*N),C(K*N)
      IK=1
      JK=1
      KK=1
      DO 1 L=1,N
      C(KK)=B(JK)-A(IK)
      IK=IK+I
      JK=JK+J
      KK=KK+K
    1 CONTINUE
      RETURN
      END
