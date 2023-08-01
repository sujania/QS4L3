c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine openfl(lufl,namein,iap,ifile,irec,istat,lrec)
      character*(*) namein
C$$$      istat=0
C$$$      call opnfil(lufl,namein,iap,ifile,irec,jstat,lrec)
C$$$      if(jstat.ne.0) then
C$$$	 write(6,'(''file='',a60)') namein
C$$$	 stop 'from openfl: file does not exist'
C$$$      endif
      return
      end
