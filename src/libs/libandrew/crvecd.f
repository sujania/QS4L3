      SUBROUTINE CRVECd(VR,station)
      !apv fix dimensions 1->*
      COMMON/PREMDATA/AMP(12),AKER(10),AR1(4),AR2(4),R1,R2,NO,IT,LO
      integer comp
      DIMENSION REC(5),X(150),XP(150),XCOSEC(150),station(5)
      COMPLEX SS,VR(*)
      COMMON/IDASTA/NIDA,IDAIND(30),IDANAM(30),SLAIDA(30),
     . SLOIDA(30),ELVIDA(30)
      COMMON/STDATA/NSTA,IND(60),STNA(60),STLAT(60),STLON(60),ELEV(60)

      CALL REPREMa(station,VR,X,XP,XCOSEC)



      IF(LO.EQ.0) RETURN

      DO 12 J=1,LO+1
   12 VR(2*LO+2-J)=VR(LO+2-J)
      DO 13 J=1,LO
      IF(MOD(J,2).EQ.0) SS=(1.,0.)
      IF(MOD(J,2).EQ.1) SS=(-1.,0.)
   13 VR(LO+1-J)=SS*CONJG(VR(LO+1+J))

      RETURN
      END


