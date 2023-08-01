      subroutine tfour(data,lenh,idir)
      dimension data(*)
      parameter (MXAUX=100)
      dimension aux(MXAUX)
      if(idir.gt.0) then
        naux=MXAUX
        call tfg(data(1),data(2),lenh,-2,ier,aux,naux)
        if(ier.ne.0) write(6,'(''tfour: tfg - ier,naux ='',2i10)') ier,naux
        call tfr(data(1),data(2),lenh,-2,ier)
        do i=1,lenh
          data(2*i)=-data(2*i)
        enddo
      else
        do i=1,lenh
          data(2*i)=-data(2*i)
        enddo
        call tfr(data(1),data(2),lenh,2,ier)
        naux=MXAUX
        call tfg(data,data(2),lenh,2,ier,aux,naux)
        if(ier.ne.0) write(6,'(''tfour: tfg - ier,naux ='',2i10)') ier,naux
      endif
      return
      end
