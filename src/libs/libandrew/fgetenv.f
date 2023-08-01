c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      character*80 function fgetenv(ident,nbyts)
      character*(*) ident
      call cgetenv(ident,fgetenv,nbyts)
      do i=nbyts+1,80
        fgetenv(i:i)=' '
      enddo
      return
      end
