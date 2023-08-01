C$**********************************************************************
      SUBROUTINE FVSMA2(A,I,B,D,L,N)
#if ( defined(MachineT)  )
      DIMENSION A(*),D(*)
      CALL SAXPY(N,B,A,I,D,L)
#elif ( defined(MachineT)  )
      DIMENSION A(I,*),D(L,*)
      bb=b
      do i=1,n
        d(1,i)=bb*a(1,i)+d(1,i)
      enddo
#else
      DIMENSION A(I*N),D(L*N)
      IK=1
      LK=1
      DO 1 M=1,N
      D(LK)=A(IK)*B+D(LK)
      IK=IK+I
      LK=LK+L
    1 CONTINUE
#endif
      RETURN
      END
