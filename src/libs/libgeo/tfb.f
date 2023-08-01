      subroutine tfb(dat,nt,ni)
      real*4 dat(2)
      dimension aux(33)
      naux=33
      
c  forward transform - time to frequency

      if(ni.gt.0) go to 10
      call tfg(dat(1),dat(2),nt,ni,ier,aux,naux)
      if(ier.ne.0) write(5,501) ier,naux
 501  format('error return from direct tfg - ier,naux =',2i10)
      call tfr(dat(1),dat(2),nt,ni,ier)
      do 11 i=1,nt+1
   11 dat(i*2)=-dat(i*2)
      return

c  inverse transform - frequency to time

   10 continue
      do 12 i=1,nt+1
   12 dat(i*2)=-dat(i*2)
      call tfr(dat(1),dat(2),nt,ni,ier)
      call tfg(dat(1),dat(2),nt,ni,ier,aux,naux)
      if(ier.ne.0) write(5,502) ier,naux
 502  format('error return from inverse tfg - ier,naux =',2i10)
      return
      end

