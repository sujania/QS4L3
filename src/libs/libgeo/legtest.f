      dimension x(100),xp(100),xcosec(100)
      dimension lord(10),xlat(10)

      data lord/3,7,9,12,15,21,0,0,0,0/
      data xlat/-20.,34.,67.,82.,-67.,-54.,0.,0.,0.,0./

      radian=45./atan(1.)

      do i=1,6
        theta=(90-xlat(i))/radian
        l=lord(i)
        mmax=lord(i)
        call legndr(theta,l,mmax,x,xp,xcosec)
        write(6,'(''Order:'',i6)') lord(i)
        write(6,'(''Latitude'',f8.4)') xlat(i)
        write(6,'(6f12.8)') (x(j),j=1,mmax+1)
      enddo

      end
