      complex function wwsp(omega)
      complex s
      s=cmplx(0.0,omega)
      wwsp=( 397.54767, 0.0 )*(s**3)/
     1     (  (s- ( -5.0136607, 6.4615109 ) )
     1       *(s- ( -5.0136607,-6.4615109 ) )
     1       *(s- ( -8.2981509, 0.0000000 ) )
     1       *(s- ( -8.6940765,-7.1968661 ) )
     1       *(s- ( -8.6940765, 7.1968661 ) )  )
      return
      end

      complex function wwspbn(omega)
      complex s
      s=cmplx(0.0,omega)
      wwspbn=( 432.83395, 0.0 )*(s**3)/
     1     (  (s- ( -4.04094, 6.47935 ) )
     1       *(s- ( -4.04094,-6.47935 ) )
     1       *(s- ( -9.25238, 0.00000 ) )
     1       *(s- ( -7.67430, 0.00000 ) )
     1       *(s- (-16.72981, 0.00000 ) )  )
      return
      end

      complex function wwlpbn(omega)
      complex s
      s=cmplx(0.0,omega)
      wwlpbn=( 0.5985275, 0.0 )*(s**3)/
     1     (  (s- ( -0.25700, 0.3376 ) )
     1       *(s- ( -0.25700,-0.3376 ) )
     1       *(s- ( -0.06283, 0.0000 ) )
     1       *(s- ( -0.06283, 0.0000 ) )  )
      return
      end

