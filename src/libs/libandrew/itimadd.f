      subroutine itimadd(itim1,itim2,add)
      dimension itim1(*),itim2(*)
      double precision add
      isec1=itim1(1)
      ifrc1=ishft(itim1(2),-16)
      call taddi2(isec1,ifrc1,isec2,ifrc2,add)
      itim2(1)=isec2
      itim2(2)=or(ishft(ifrc2,16),and(z'0000ffff',itim1(2)))
      return
      end
