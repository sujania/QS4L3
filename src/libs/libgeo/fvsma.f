C$**********************************************************************
      SUBROUTINE FVSMA(A,I,B,C,K,D,L,N)
      DIMENSION A(I*N),C(K*N),D(L*N)
#if ( defined(MachineS) || defined(MachineT) )
      CALL SCOPY(N,C,K,D,L)
      CALL SAXPY(N,B,A,I,D,L)
#else
      IK=1
      KK=1
      LK=1
      DO 1 M=1,N
      D(LK)=A(IK)*B+C(KK)
      IK=IK+I
      KK=KK+K
      LK=LK+L
    1 CONTINUE
#endif
      RETURN
      END

