c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine bffos(lufl,ifbin,ibuf,nbytes,istat,irec)
      include 'openfile.h'
      dimension ibuf(*)
      parameter (IOPT=0)
      krec=irec-1
      if(irec.eq.0) krec=jrec(lufl)
      kchn=jchn(lufl)
      if(jrecl(lufl).eq.0) then
        write(6,*) 'bffos: unsupported option'
      else if(jrecl(lufl).gt.0) then
        nb=nbytes
        if(jfile(lufl).ne.200) nb=min0(nbytes,jrecl(lufl))
        call clseek(jchn(lufl),jrecl(lufl)*krec,IOPT,ires,ierrno)
        if(ierrno.ne.0) call check('clseek in bffos')
        call cwrite(jchn(lufl),ibuf,nb,ires,ierrno)
        if(ierrno.ne.0.or.ires.ne.nbytes) call check('cwrite in bffos')
        istat=2
        jrec(lufl)=krec+(nb+jrecl(lufl)-1)/jrecl(lufl)
        lenglu(lufl)=max0(lenglu(lufl),jrec(lufl))
      else
       write(6,*) 'bffos: unsupported option'
      endif
      return
      end
