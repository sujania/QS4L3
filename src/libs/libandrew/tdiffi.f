      double precision function tdiffi(itim1,itim2)
      dimension itim1(2),itim2(2)
      tdiffi=itim1(1)-itim2(1)+(ishft(itim1(2),-16)-ishft(itim2(2),-16))/10000.d0
      return
      end
