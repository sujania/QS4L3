      function krunge(n,y,f,x,h)
      dimension phi(6),savey(6),y(6),f(6)
      data m/0/
      save savey,phi
      m=m+1
      goto(1,2,3,4,5), m
1     krunge=1
      return
2     do 22 j=1,n
      savey(j)=y(j)
      phi(j)=f(j)
22    y(j)=savey(j)+0.5d0*h*f(j)
      x=x+0.5d0*h
      krunge=1
      return
3     do 33 j=1,n
      phi(j)=phi(j)+2.d0*f(j)
33    y(j)=savey(j)+0.5d0*h*f(j)
      krunge=1
      return
4     do 44 j=1,n
      phi(j)=phi(j)+2.0d0*f(j)
44    y(j)=savey(j)+h*f(j)
      x=x+0.5d0*h
      krunge=1
      return
5     do 55 j=1,n
55    y(j)=savey(j)+(phi(j)+f(j))*h/6.d0
      m=0
      krunge=0
      return
      end
