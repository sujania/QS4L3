C     REPREM CALCULATES THE RECEIVER VECTOR REVEC
C     INPUT:  RECE(5) , LAT,LONG,VERT,N-S,E-W  (????)
C     OUTPUT: REVEC(LORD+1)  REVEC(1) - M=0 ,ETC.
C     WORKSPACES HAVE TO BE ALLOCATED TO (LORD+1)
C     REVEC(-M)=(-1)**M*(REVEC(M)

      SUBROUTINE REPREMa(REC,REVEC,X,XP,XCOSEC)
c     apv fix dimensions 1->*
      COMPLEX REVEC(*),FACT,EFAC
      DIMENSION REC(5),X(*),XCOSEC(*),XP(*)
      COMMON/PREMDATA/OM,QU,VAA,HAA,GVEL,ELLIP,DUM1(5),ROTSP
     + ,DUM2(10),U,UP,V,VP,U2,U2P,V2,V2P,R1,R2,NO,IQ,L
      DATA PI,PIB2/3.1415926     ,1.5707963     /
      THETA=PIB2-REC(1)*PI/180.E0
      PHI=REC(2)*PI/180.E0
      COSEC=1.E0/SIN(THETA)
      FACT=(1.,0.)
      EFAC=CEXP(CMPLX(0.E0,PHI))
      CALL LEGRG(THETA,L,L,X,XP,XCOSEC)
      LP1=L+1
      DO 10 I=1,LP1
      XM=I-1
      IF(IQ.EQ.2) GOTO 20
      REVEC(I)=(REC(3)*VAA*X(I)+REC(4)*HAA*XP(I)
     +       +REC(5)*CMPLX(0.E0,XM*COSEC*HAA)*X(I))*FACT
ca      WRITE(6,*) REC(3),VAA,X(I),REC(4),HAA,XP(I),REC(5),XM,COSEC,FACT
      GOTO 9
   20 REVEC(I)=(CMPLX(0.E0,XM*COSEC*X(I)*REC(4))
     1        -REC(5)*XP(I))*HAA*FACT
    9 FACT=FACT*EFAC
   10 CONTINUE
      RETURN
      END