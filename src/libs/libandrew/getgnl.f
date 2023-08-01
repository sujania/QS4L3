c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      character*80 function getgnl(ident,nni,nbyts
     1  ,swnam,cline,deflt
     1  ,nswt,lswnam,ldeflt,ifreq,nmin,nmax
     1  ,icnt,iopn,iptr,itable,lcom)
      character*(*) ident
      character*80 str80,fgetenv
      character*1 ket
      include 'getgnl.h'

      character*1 null
      null=char(0)
      nn=nni
      ll=len(ident)
      do while (ll.gt.0.and.(ident(ll:ll).eq.' '
     1          .or.ident(ll:ll).eq.null))
        ll=ll-1
      enddo
      imatch=0
      iswt=0
      ierr=0
      do i=1,nswt
        if((ll.eq.0.and.lswnam(i).eq.0)
     1   .or.(ll.gt.0.and.ident(1:ll).eq.swnam(i)(1:ll))) then
          imatch=1+imatch
          if(iswt.eq.0) iswt=i
          if(ll.eq.lswnam(i)) then
            imatch=1
            iswt=i
            goto 22
          endif
        endif
      enddo
   22 continue
      if(imatch.eq.0) then
        ierr=1
        goto 99
      else if(imatch.ne.1) then
        ierr=2
        goto 99
      else



        if(nn.eq.0) then
          if(ll.eq.0) then
            getgnl=cline(1:lcom)
            nbyts=lcom
            goto 100
          else
            if(icnt(iswt).lt.0) then
              getgnl=' '
              nbyts=-1
              goto 100
            else
              getgnl=ident(1:ll)
              nbyts=ll
              goto 100
            endif
          endif
        else if(icnt(iswt).lt.0) then
          nbyts=-1
          getgnl=' '
          if(ldeflt(iswt).ge.0) then
 202        continue    ! label for looping to replace multiple environment variables
            ldef=ldeflt(iswt)
            ipdol=0
            do i=1,ldef
              if(deflt(iswt)(i:i).eq.'$') ipdol=i
            enddo
            if(ipdol.ne.0) then
              str80=deflt(iswt)(1:ldef)
              iq1=ipdol+1
              if(str80(iq1:iq1).eq.'(') then
                 ket=')'
                 iq1=iq1+1
              else if(str80(iq1:iq1).eq.'{') then
                 ket='}'
                 iq1=iq1+1
              endif
              iq2=ldef
              iq3=iq2+1
              if(ket.ne.' ') then
                do i=iq1,ldef
                  if(str80(i:i).eq.ket) iq2=i-1
                enddo
                iq3=iq2+2
              endif
              ldeflt(iswt)=0
              if(ipdol.gt.1) then
                deflt(iswt)(1:ipdol-1)=str80(1:ipdol-1)
                ldeflt(iswt)=ipdol-1
              endif
              deflt(iswt)(ldeflt(iswt)+1:80)
     1            =fgetenv(str80(iq1:iq2),ladd)
              if(ladd.le.0) then
                write(0,*) 'chekgl: Environment variable '//str80(iq1:iq2)//' not set'
                call exit(2)
              endif

              ldeflt(iswt)=ldeflt(iswt)+ladd
              if(iq3.le.ldef) then
                 deflt(iswt)=deflt(iswt)(1:ldeflt(iswt))//str80(iq3:ldef)
                 ldeflt(iswt)=ldeflt(iswt)+ldef-iq3+1
              endif
              
              goto 202  ! iterate replacement of environment variables
              
            endif



            ii=1
            ik=0
            do while(ii.le.ldeflt(iswt).and.ik.lt.nn)
              do while(ii.le.ldeflt(iswt).and.deflt(iswt)(ii:ii).eq.' ')
                ii=ii+1
              enddo
              ik=ik+1
              ij=0
              do while(ii.le.ldeflt(iswt).and.deflt(iswt)(ii:ii).ne.' ')
                ij=ij+1
                if(ik.eq.nn) then
                  nbyts=ij
                  getgnl(ij:ij)=deflt(iswt)(ii:ii)
                endif
                ii=ii+1
              enddo
            enddo
          
          endif


          goto 100
        else if(nn.gt.icnt(iswt)) then
          nbyts=0
          getgnl=' '
          goto 100
        else




          ipt=iptr(iswt)
          ipt0=ipt
          do i=1,nn-1
            ipt=itable(ipt,3)
            if(ipt.eq.ipt0) then
              ierr=3
              goto 99
            endif
          enddo
          ip1=itable(ipt,1)
          ip2=itable(ipt,2)
          nbyts=ip2-ip1+1
          getgnl=cline(ip1:ip2)
          goto 100
        endif
      endif
   99 continue
      goto (1,2,3),ierr
    1 write(0,'(''getgnl: no match'')')
      goto 97
    2 write(0,'(''getgnl: non-unique abbreviation for identifier'')')
      goto 97
    3 write(0,'(''getgnl: required argument not found in table'')')
      goto 97
   97 call exit(2)
  100 continue
      return
      end
