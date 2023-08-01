c------------------------------------
      subroutine dofilt(data,ndata,dt,ahhdr,work,nwork)

c work() is used for workspace. Must be dimensioned in calling program.
c nwork is the number of words avaiable in work -- used for checking.
c The space needed will be 2 greater than the next power of 2 higher than ndata,
c OR 258, whichever is greater.

      dimension data(*),ahhdr(0:*)
      include 'filtpar.h'
      complex cv,respahx,butter,sroins,wwspbn,wwlpbn
      dimension cva(2)
      equivalence (cv,cva(1),cv1),(cva(2),cv2)
      dimension work(nwork)
      call remavl(ndata,data)
      if(nflstg.eq.0) return
      pi=4.*atan(1.)
      ntran=7
      do while(2**(ntran+1).lt.ndata)
        ntran=ntran+1
      enddo
      nfp=2**ntran
      la=2*nfp+2
      fp=nfp
      dom=pi/(fp*dt)
      if(nwork.lt.la) stop 'dofilt: workspace too small'
      ia=1


      i1=ia
      do i=1,ndata
        work(i1)=data(i)
        i1=i1+1
      enddo
      do i=ndata+1,la
        work(i1)=0.
        i1=i1+1
      enddo
      call rfour(work(ia),ntran,1)
      do istg=1,nflstg
        if     (ifltopt(istg).eq.JFOCHP) then
          om1=2.*pi*fltprm(1,istg)
          om2=2.*pi*fltprm(2,istg)
          xm=pi/(om2-om1)
          i1=min0(ifix(om1/dom+1.),nfp)
          i2=min0(ifix(om2/dom),nfp)
          ic=ia
          do i=0,i1-1
            work(ic)=0.
            work(ic+1)=0.
            ic=ic+2
          enddo
          do i=i1,i2
            fc=.5*(1.+cos( xm*(om2-i*dom) ) )
            work(ic)=work(ic)*fc
            work(ic+1)=work(ic+1)*fc
            ic=ic+2
          enddo
        else if(ifltopt(istg).eq.JFOCLP) then
          om1=2.*pi*fltprm(1,istg)
          om2=2.*pi*fltprm(2,istg)
          xm=pi/(om2-om1)
          i1=min0(ifix(om1/dom+1.),nfp)
          i2=min0(ifix(om2/dom),nfp)
          ic=ia+2*i1
          do i=i1,i2
            fc=.5*(1.+cos( xm*(i*dom-om1) ))
            work(ic)=work(ic)*fc
            work(ic+1)=work(ic+1)*fc
            ic=ic+2
          enddo
          do i=i2+1,nfp
            work(ic)=0.
            work(ic+1)=0.
            ic=ic+2
          enddo
        else if(ifltopt(istg).eq.JFOBHP) then
          om1=2.*pi*fltprm(1,istg)
          nord=fltprm(2,istg)
          ic=ia+2
          do i=1,nfp
            cv=butter(-om1/(i*dom),nord)
            a1=work(ic)
            a2=work(ic+1)
            work(ic)=a1*cv1-a2*cv2
            work(ic+1)=a1*cv2+a2*cv1
            ic=ic+2
          enddo
        else if(ifltopt(istg).eq.JFOBLP) then
          om1=2.*pi*fltprm(1,istg)
          nord=fltprm(2,istg)
          ic=ia+2
          do i=1,nfp
            cv=butter((i*dom)/om1,nord)
            a1=work(ic)
            a2=work(ic+1)
            work(ic)=a1*cv1-a2*cv2
            work(ic+1)=a1*cv2+a2*cv1
            ic=ic+2
          enddo
        else if(ifltopt(istg).eq.JFOBBP) then
          om1=2.*pi*fltprm(1,istg)
          om2=2.*pi*fltprm(2,istg)
          nord=fltprm(3,istg)
          ic=ia+2
          do i=1,nfp
            om=i*dom
            cv=butter((om/om2-om1/om),nord)
            a1=work(ic)
            a2=work(ic+1)
            work(ic)=a1*cv1-a2*cv2
            work(ic+1)=a1*cv2+a2*cv1
            ic=ic+2
          enddo
        else if(ifltopt(istg).eq.JFOPHP) then
          om1=2.*pi*fltprm(1,istg)
          eta=fltprm(2,istg)
          nord=fltprm(3,istg)
          oms=om1*om1
          tome=2.*om1*eta
          ic=ia+2
          do i=1,nfp
            om=i*dom
            cv=(cmplx(om*om,0.)/cmplx(om*om-oms,-tome*om))**nord
            a1=work(ic)
            a2=work(ic+1)
            work(ic)=a1*cv1-a2*cv2
            work(ic+1)=a1*cv2+a2*cv1
            ic=ic+2
          enddo
        else if(ifltopt(istg).eq.JFOPLP) then
          om1=2.*pi*fltprm(1,istg)
          eta=fltprm(2,istg)
          nord=fltprm(3,istg)
          oms=om1*om1
          tome=2.*om1*eta
          ic=ia+2
          do i=1,nfp
            om=i*dom
            cv=(cmplx(-oms,0.)/cmplx(om*om-oms,-tome*om))**nord
            a1=work(ic)
            a2=work(ic+1)
            work(ic)=a1*cv1-a2*cv2
            work(ic+1)=a1*cv2+a2*cv1
            ic=ic+2
          enddo
        else if(ifltopt(istg).eq.JFOSPW) then
          spow=fltprm(1,istg)
          cv=cexp( cmplx(0.,spow*.5*pi) )
          ic=ia+2
          do i=1,nfp
            amp=(i*dom)**spow
            a1=work(ic)
            a2=work(ic+1)
            work(ic)=amp*(a1*cv1-a2*cv2)
            work(ic+1)=amp*(a1*cv2+a2*cv1)
            ic=ic+2
          enddo
        else if(ifltopt(istg).eq.JFODTS) then
          dts=fltprm(1,istg)
          tcor=fltprm(2,istg)
          omfac1=tcor/(2.0*pi)
          omfac2=dts/pi
          ic=ia+2
          do i=1,nfp
            om=i*dom
            cv=cexp ( cmplx( 0.,om*omfac2) * clog( cmplx(0.,om*omfac1) ) )
            a1=work(ic)
            a2=work(ic+1)
            work(ic)=a1*cv1-a2*cv2
            work(ic+1)=a1*cv2+a2*cv1
            ic=ic+2
          enddo
        else if(ifltopt(istg).eq.JFOTSH) then
          tsh1=fltprm(1,istg)
          ic=ia+2
          do i=1,nfp
            om=i*dom
            cv=cexp ( cmplx( 0.,-om*tsh1) )
            a1=work(ic)
            a2=work(ic+1)
            work(ic)=a1*cv1-a2*cv2
            work(ic+1)=a1*cv2+a2*cv1
            ic=ic+2
          enddo
        else if(ifltopt(istg).eq.JFOSRO) then
          ic=ia+2
          do i=1,nfp
            cv=5000.0e+06*sroins(i*dom)
            a1=work(ic)
            a2=work(ic+1)
            work(ic)=a1*cv1-a2*cv2
            work(ic+1)=a1*cv2+a2*cv1
            ic=ic+2
          enddo
        else if(ifltopt(istg).eq.JFOWWS) then
          ic=ia+2
          do i=1,nfp
            cv=1.e06*wwspbn(i*dom)
            a1=work(ic)
            a2=work(ic+1)
            work(ic)=a1*cv1-a2*cv2
            work(ic+1)=a1*cv2+a2*cv1
            ic=ic+2
          enddo
        else if(ifltopt(istg).eq.JFOWWL) then
          ic=ia+2
          do i=1,nfp
            cv=1.e06*wwlpbn(i*dom)
            a1=work(ic)
            a2=work(ic+1)
            work(ic)=a1*cv1-a2*cv2
            work(ic+1)=a1*cv2+a2*cv1
            ic=ic+2
          enddo
        else if(ifltopt(istg).eq.JFOHTR) then
          ic=ia+2
          do i=1,nfp
            a1=work(ic)
            a2=work(ic+1)
            work(ic)=-a2
            work(ic+1)=a1
            ic=ic+2
          enddo


        else if(ifltopt(istg).eq.JFOGSS) then
          gsom0=2*pi*fltprm(1,istg)
          gsowd=2*pi*fltprm(2,istg)
          ic=ia+2
          do i=1,nfp
            arg=-.5*((i*dom-gsom0)/gsowd)**2
            if(arg.lt.-30) then
              rfc=0.
            else
              rfc=exp(arg)
            endif
            a1=work(ic)
            a2=work(ic+1)
            work(ic)=work(ic)*rfc
            work(ic+1)=work(ic+1)*rfc
            ic=ic+2
          enddo




        else if(ifltopt(istg).eq.JFOCVF) then
          ic=ia+2
          do i=1,nfp
            cv=(0.,0.)
            do j=1,lfilt
              cv=cv+filt(j)*cexp(cmplx(0.,-(tshff+smpinf*(j-1))*(i*dom)) )
            enddo
            a1=work(ic)
            a2=work(ic+1)
            work(ic)=a1*cv1-a2*cv2
            work(ic+1)=a1*cv2+a2*cv1
            ic=ic+2
          enddo
        else if(ifltopt(istg).eq.JFODCI) then
c         write(6,'(''Amplitude response at 1 rad/s'',1pe13.5)') cabs(respahx(1.,ahhdr))
          ic=ia+2
          do i=1,nfp
            cv=(1.,0.)/respahx(i*dom,ahhdr)
            a1=work(ic)
            a2=work(ic+1)
            work(ic)=a1*cv1-a2*cv2
            work(ic+1)=a1*cv2+a2*cv1
            ic=ic+2
          enddo
        endif
      enddo

      work(ia)=0.
      work(ia+1)=0.
      work(ia+la-2)=0.
      work(ia+la-1)=0.
      call rfour(work(ia),ntran,-1)
      i1=ia
      do i=1,ndata
        data(i)=work(i1)
        i1=i1+1
      enddo

      return
      end
