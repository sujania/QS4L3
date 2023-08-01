      SUBROUTINE INSTR(T0,TG,H0,H1,SIGMA,PERIOD,IFGO,XMAG,PHASE)
      save
      PI=3.1415926
      TWOPI=2.*PI
      IF(IFGO.GT.0) GO TO 100
      XMU=TG/T0
      U1=25./T0
      IF(T0.EQ.15.)U1=15./T0
      GO TO 101
  100 U1=PERIOD/T0
  101  FT=1.-((1.+1./(XMU**2))+4.*H0*H1*(1.-SIGMA**2)/XMU)*(U1**2)
     1+(U1**4)/(XMU**2)
      ST=-2.*(H0+H1/XMU)*U1+2.*(H0/XMU+H1)*(U1**3)/XMU
      F=U1/SQRT(FT**2+ST**2)
      PHASE=ATAN2(FT,ST)+PI
      IF(IFGO.GT.0) GO TO 102
      Q=750.*TWOPI/(T0*F)
      XMAG=750.
      RETURN
  102 XMAG=Q*T0*F/TWOPI
      RETURN
      END
