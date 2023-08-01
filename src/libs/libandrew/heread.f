c-------------------------------------------------
      subroutine heread(lu,file,x,nstr,lmax,mask,lask,inorm)
      save
      character*(*) file
      dimension x(*),mask(*),lask(*)
      character*132 abuf
      character*80 form
      sqth=sqrt(.5)
      open(lu,file=file,status='old')
      read(lu,'(a132)') abuf
      llen=132
      do 30 i=1,132
      if(abuf(llen:llen).ne.' ') goto 40
   30 llen=llen-1
   40 read(abuf,'(10x,i5)') lmax
      write(6,*) 'lmax in heread is',lmax
      num=llen-15-1-(lmax+1)-1-3-1

      if(llen.gt.15) then
        idef=0
        write(form,'(''(16x,'',i3,''i1,1x,i3,1x,80i1)'')') lmax+1
!        write(6,'(a30,2i8)') form(1:30),lmax,num
        read(abuf,form) (lask(i+1),i=0,lmax),nstr,(mask(i),i=1,num)
        if(num.ne.nstr) pause 'error 1 in heread'
      else
        idef=1
        do 50 i=0,lmax
        lask(i+1)=1
   50   continue
      endif

      leny=0
      do 60 l=0,lmax
      if(lask(l+1).ne.0) leny=leny+2*l+1
   60 continue

      ip=0
  100 read(lu,'(a132)',end=99) abuf
      backspace lu
      ip=1+ip
      if(idef.ne.0) mask(ip)=1

      ind=1+(ip-1)*leny
      do 704 l=0,lmax
      if(lask(l+1).ne.0) then
        ind1=ind+2*l
        read(lu,'(11e12.4)')(x(i),i=ind,ind1)
        if(inorm.eq.1) x(ind)=x(ind)/sqth
        ind=ind1+1
      endif
  704 continue
      goto 100

   99 if(idef.ne.0) then
        nstr=ip
        nstruc=ip
      else
        nstruc=0
        do i=1,nstr
          if(mask(i).ne.0) nstruc=nstruc+1
        enddo
c       if(nstruc.ne.ip) pause 'error 3 in heread'
      endif

      ii=istlen(file)
!      write(6,'(''model from'',1x,a)') file(1:ii)
!      write(6,'(i4,'' parameters:  mask '',80i1)') nstr,(mask(i),i=1,nstr)
!      write(6,'(4x,'' lmax ='',i4,'':  mask '',80i1)') lmax,(lask(i+1),i=0,lmax)
      close(lu)
      return
      end
