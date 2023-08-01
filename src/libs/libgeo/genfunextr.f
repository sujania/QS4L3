c----------------------------------------------------------
      subroutine genfunextr(iopt,gxmin,gxmax,gymin,gymax)

      include 'genfun.h'

      call vmima(parms(1,indprm(iopt)+1),GFINC,numprm(iopt),gxmax,gxmin)
      call vmima(parms(2,indprm(iopt)+1),GFINC,numprm(iopt),gymax,gymin)

      return
      end

