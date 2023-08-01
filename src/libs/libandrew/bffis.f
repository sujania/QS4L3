c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine bffis(lufl,ifbin,ibuf,nbytes,istat,nread,irec)
      include 'openfile.h'
      dimension ibuf(*)
      krec=irec-1
      if(irec.eq.0) krec=jrec(lufl)
      if (jrecl(lufl).le.0) stop 'bffis: unsupported'
         nb=nbytes
         if(jfile(lufl).ne.200) nb=min0(nb,jrecl(lufl))
         call clseek(jchn(lufl),jrecl(lufl)*krec,0,ires,ierrno)
         if(ierrno.ne.0) call check('clseek in bffis')
         call cread(jchn(lufl),ibuf,nb,nread,ierrno)
         if(ierrno.ne.0) then
            call cperror('cread in bffis')
            write(6,'(''Error:'',i5)') ierrno
            istat=4
            call exit(2)
         else
            istat=2
            if(nread.eq.0) istat=3
            if(istat.ne.3) then
               jrec(lufl)=krec+(nread+jrecl(lufl)-1)/jrecl(lufl)
            else
               jrec(lufl)=krec
            endif
         endif
     
      return
      end
      
