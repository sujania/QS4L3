c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine bffi(lufl,ifbin,ibuf,nbytes,istat,nread,irec)
      include 'openfile.h'
      dimension ibuf(*)
      if(jrecl(lufl).eq.1.and.jfile(lufl).eq.200) then
         stop 'bffi: no longer supported'
C$$$        nread=0
C$$$        nget=nbytes
C$$$        irc=irec
C$$$        if(irc.eq.0) irc=jrec(lufl)+1
C$$$        iad=1
C$$$   10   nb=min0(nget,65532)
C$$$        call bffis(lufl,ifbin,ibuf(iad),nb,istat,nr,irc)
C$$$        if(istat.eq.2.and.nr.eq.nb) then
C$$$          if(nr.eq.nget) then
C$$$            nread=nread+nr
C$$$          else
C$$$            nr4=nr/4
C$$$            nr=4*nr4
C$$$            irc=irc+nr
C$$$            nread=nread+nr
C$$$            iad=iad+nr4
C$$$            nget=nget-nr
C$$$            goto 10
C$$$          endif
C$$$        else
C$$$          nread=nread+nr
C$$$        endif
      else
        call bffis(lufl,ifbin,ibuf,nbytes,istat,nread,irec)
      endif
      return
      end
