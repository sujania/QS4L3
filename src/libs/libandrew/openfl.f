c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine openfl(lufl,namein,iap,ifile,irec,istat,lrec)
      character*(*) namein
      istat=0
      call opnfil(lufl,namein,iap,ifile,irec,jstat,lrec)
      if(jstat.ne.0) then
	 write(6,'(''file='',a60)') namein
	 stop 'from openfl: file does not exist'
      endif
      return
      end
