
      SUBROUTINE CSVECa(VS,SOUR)
      COMMON/PREMDATA/AMP(12),AKER(10),AR1(4),AR2(4),R1,R2,NO,IT,LO
      DIMENSION SOUR(9),X(150),XP(150),XCOSEC(150)
!     apv fix dimension 1->*
      COMPLEX SS,VS(*)
      
!      CALL SOPREM(SOUR,VS,X,XP,XCOSEC)
      CALL SOPREMa(SOUR,VS,X,XP,XCOSEC)
      
      IF(LO.EQ.0) RETURN
      
      DO 10 J=1,LO+1
   10 VS(2*LO+2-J)=VS(LO+2-J)
      DO 11 J=1,LO
      IF(MOD(J,2).EQ.0) SS=(1.,0.)
      IF(MOD(J,2).EQ.1) SS=(-1.,0.)
   11 VS(LO+1-J)=SS*CONJG(VS(LO+1+J))
  
      RETURN
      END
