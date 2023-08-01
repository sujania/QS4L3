c----------------------------------------------------------
      subroutine genfunstartup(filei)
      character*(*) filei
      character*100 file
      include 'genfun.h'

      character*80 line
      file=filei
      lfile=istlen(file)
      write(0,'(a)') 'genfunstartup from: '//file(1:lfile)
      irow=0
      nop=0
      do i=1,MXOPT
        numprm(i)=-1
      enddo
      open(87,file=file(1:lfile),status='old')
   10 read(87,'(a)',end=99) line
      if(line(1:1).eq.'#') goto 10
      if(line(1:6).eq.'option') then
        read(line(7:80),*) iop
        if(iop.le.0.or.iop.gt.MXOPT) stop 'genfunstartup: iop too big'
        indprm(iop)=irow
        if(numprm(iop).ne.-1) then
          write(6,'(''genfunstartup: Warning: option''
     1         ,i4,'' overwritten'')') iop
        endif
        numprm(iop)=0
      else
        irow=irow+1
        if(irow.gt.MXPRM) stop 'genfunstartup: irow too big'
        read(line,*) (parms(i,irow),i=1,GFINC)
        numprm(iop)=numprm(iop)+1
      endif
      goto 10
   99 close(87)
      return
      end
