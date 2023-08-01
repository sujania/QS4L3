C     SOPREM RETURNS THE EXITATION OF NORMAL MODES IN EVEC
C     INPUT: SOUR(I),I=1,9
C            LATITUDE,LONGITUDE,DEPTH,XM(1-6)
C     OUTPUT:COMPLEX EVEC(LORD+1) WHERE LORD IS ANGULAR ORDER
C     WORKSPACES X,XP,XCOSEC HAVE TO BE DIMENSIONED (LORD+1)
C     EVEC(1) - M=0   EVEC(2) - M=1   EVEC(LORD+1) - M=LORD
C     THIS IS A SINGLE PRECISION VERSION OF EXTAN
C     EVEC(-M)=(-1)**M*CONJG(EVEC(M)

      SUBROUTINE SOPREMa(SOUR,EVEC,X,XP,XCOSEC)
      ! APV evec was dimension 3
      COMPLEX EVEC(*),FACT,EFAC,E(6)
c     apv change dimension from 1 to *
      DIMENSION SOUR(9),X(*),XP(*),XCOSEC(*)
      COMMON/PREMDATA/OM,QU,VAA,HAA,GVEL,ELLIP,DUM1(5),ROTSP
     1 ,DUM2(10),U,UP,V,VP,U2,U2P,V2,V2P,R1,R2,NO,IQ,L
      DATA PI,PIB2/3.141592      ,1.570796      /

      RS=1.0-SOUR(3)/6371.0
      THETA=PIB2-SOUR(1)*PI/180.E0
      PHI=SOUR(2)*PI/180.E0
      H=RS-R1
      HN=R2-R1
      RHN=1.E0/HN
      HN2=RHN*RHN
      HN3=RHN*HN2
      A=HN3*(HN*(VP+V2P)+
     1         2.E0*(V-V2))
      B=HN2*(3.E0*(V2-V)-
     1         HN*(V2P+2.D0*VP))
      VV=V+H*(VP+H*(B+H*A))
      VVP=VP+H*(2.E0*B+3.E0*H*A)
ca ***
ca these two values have to be set to zero here
ca *** 
      uu=0.0
      uup=0.0
ca ***
      IF(IQ.EQ.2) GOTO 20
      A=HN3*(HN*(UP+U2P)+
     1         2.E0*(U-U2))
      B=HN2*(3.E0*(U2-U)-
     1         HN*(U2P+2.E0*UP))
      UU=U+H*(UP+H*(B+H*A))
      UUP=UP+H*(2.E0*B+3.E0*H*A)

   20 B0=SQRT(FLOAT(2*L+1)/(4.E0*PI))
      B1=-.5D0*B0*SQRT(FLOAT(L*(L+1)))
      B2=0.E0
      IF(L.GT.0) B2=-.25*B1*SQRT(FLOAT((L-1)*(L+2)))
      LP1=L+1
      DO 1 I=1,LP1
    1 EVEC(I)=(0.,0.)
      ES=VVP+(UU-VV)/RS
      IF(THETA.NE.0.E0.AND.THETA.NE.PI) GOTO 100
      FACT=(1.,0.)
      IF(THETA.EQ.PI) FACT=(-1.,0.)
      IF(IQ.EQ.2) GOTO 110
      EVEC(1)=SOUR(4)*B0*UUP+(SOUR(5)+SOUR(6))
     1     *(UU-.5E0*FLOAT(L*(L+1))*VV)*B0/RS
      IF(L.EQ.0) RETURN
      EVEC(2)=FACT*ES*B1*(SOUR(7)-FACT*CMPLX(0.E0,SOUR(8)))
      IF(L.EQ.1) RETURN
      EVEC(3)=(SOUR(5)-SOUR(6))*2.E0*B2*VV/RS+SOUR(9)*
     1     FACT*CMPLX(0.E0,-4.E0*B2*VV/RS)
      RETURN

  110 IF(L.EQ.0) RETURN
      EVEC(2)=-FACT*B1*ES*(FACT*CMPLX(0.E0,SOUR(7))-SOUR(8))
      IF(L.EQ.1) RETURN
      EVEC(3)=-FACT*B2*VV*2.E0/RS*CMPLX(0.E0,SOUR(5)-SOUR(6))
     1     -4.E0*B2*VV*SOUR(9)/RS
      RETURN
  100 CALL LEGRG(THETA,L,L,X,XP,XCOSEC)
      CT=COS(THETA)
      ST=SIN(THETA)
      COSEC=1.E0/ST
      COT=CT/ST
      FL3=L*(L+1)
      FACT=(1.,0.)
      EFAC=CEXP(CMPLX(0.E0,-PHI))
      DO 11 I=1,LP1
      XM=I-1
      Y=X(I)
      YP=XP(I)
      IF(IQ.EQ.2) GOTO 21
      E(1)=UUP*Y
      E(2)=(UU*Y-VV*(COT*YP-((XM*COSEC)**2-FL3)*Y))/RS
      E(3)=(UU*Y+VV*(COT*YP-((XM*COSEC)**2)*Y))/RS
      E(4)=ES*YP
      E(5)=CMPLX(0.E0,-XM*ES*COSEC*Y)
      E(6)=CMPLX(0.E0,-2.E0*XM*VV*COSEC*(YP-COT*Y)/RS)
      GOTO 22
   21 E(1)=(0.,0.)
      E(2)=CMPLX(0.E0,-VV*XM*COSEC*(YP-COT*Y)/RS)
      E(3)=-E(2)
      E(4)=CMPLX(0.E0,-XM*COSEC*Y*ES)
      E(5)=-YP*ES
      E(6)=VV*(2.E0*COT*YP-(2.E0*((XM*COSEC)**2)-FL3)*Y)/RS
   22 DO 23 J=1,6
   23 EVEC(I)=EVEC(I)+SOUR(J+3)*E(J)
      EVEC(I)=EVEC(I)*FACT
      FACT=FACT*EFAC
   11 CONTINUE
      RETURN
      END

c ------------------------------------------------------------------------
