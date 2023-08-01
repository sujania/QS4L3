c----------------------------------------------------------
      function ginter(x1,x2,x,iop)
      dimension x1(*),x2(*)
      goto (1,2,3),iop
      stop 'ginter:invalid option'
    1 ginter=x1(2)+(x2(2)-x1(2))*(x-x1(1))/(x2(1)-x1(1))
      return
    2 h=x2(1)-x1(1)
      h2=h*h
      xx=(x-x1(1))/h
      ginter=x1(2)+xx*(x1(3)*h+.5*xx*x1(4)*h2)
     1       +(x2(2)-x1(2)-x1(3)*h-.5*x1(4)*h2)*xx**3*(6.*xx*xx-15.*xx+10.)
     1       +(x2(3)-x1(3)-x1(4)*h)*h*xx**3*(-3.*xx*xx+7.*xx-4.)
     1       +(x2(4)-x1(4))*h2*xx**3*(0.5*xx*xx-xx+.5)
      return
    3 ginter=exp(alog(x1(2))+alog(x2(2)/x1(2))*(x-x1(1))/(x2(1)-x1(1)))
      return
      end

