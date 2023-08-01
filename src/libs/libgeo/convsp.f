      SUBROUTINE CONVSP(J)
      IF(J.GE.0) RETURN
      J=-J
      CALL BSET(J,0)
      RETURN
      END
