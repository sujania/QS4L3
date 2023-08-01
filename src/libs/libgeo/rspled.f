c----------------------------------
      function rspled(i1,i2,x,y,q,s)
      dimension x(*),y(*),q(3,*)
      data i/1/
      ii=i2-1
      i=max0(i,i1)
      i=min0(i,ii)
      if(x(i2)-x(i1))1,2,2
 1    if(s-x(i))3,3,4
 4    i=i-1
      if(i-i1)11,6,1
 3    if(s-x(i+1))5,6,6
 5    i=i+1
      if(i-ii)3,6,7
 2    if(s-x(i+1))8,8,9
 9    i=i+1
      if(i-ii)2,6,7
 8    if(s-x(i))10,6,6
 10   i=i-1
      if(i-i1)11,6,8
 7    i=ii
      go to 6
 11   i=i1
 6    h=s-x(i)
      rspled=q(1,i)+h*(2.*q(2,i)+h*3.*q(3,i))
      return
      end
