c----------------------------------------------------------
      function genfun(dep,ivop)

      include 'genfun.h'

      do i=indprm(ivop)+2,indprm(ivop)+numprm(ivop)-1
        if(dep.le.parms(1,i)) then
          ly=i
          goto 201
        endif
      enddo
      ly=indprm(ivop)+numprm(ivop)
  201 genfun=ginter(parms(1,ly-1),parms(1,ly),dep,ifix(parms(5,ly)))
      return
      end
