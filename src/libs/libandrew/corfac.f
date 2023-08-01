c----------------------------------------------------------------------
      subroutine corfac(iq,wcom,jcom,xac,xf,xln)
      save
      real lcon,ncon
      common/modl1/n,nic,noc,moho,nsl,ifanis,r(222)
     1            ,rho(222),qrho(3,222),g(222),ell(222),eta(222)
      common/modl2/acon(222),qacon(3,222),ccon(222),qccon(3,222)
     1            ,lcon(222),qlcon(3,222),ncon(222),qncon(3,222)
     2            ,fcon(222),qfcon(3,222)
      common/modl3/qshear(222),qkappa(222)
      common/nond/rn,wn,vn,gn,rhobar
c
      data pi/3.14159265358979/
      data fot/1.33333333333333333333/
c
      fct=2.*alog(.5*wcom/pi)/pi
      xmu=qshear(iq)*fct
      xln=1.+xmu
      if(jcom.eq.2) return
      xka=qkappa(iq)*fct
      rt=fot*(acon(iq)+ccon(iq)-2.*fcon(iq)+5.*ncon(iq)+6.*lcon(iq))
     1   /(8.*acon(iq)+3.*ccon(iq)+4.*fcon(iq)+8.*lcon(iq))
      xf=1.+((1.-rt)*xka-.5*rt*xmu)/(1.-1.5*rt)
      xac=1.+(1.-rt)*xka+rt*xmu


      ! turn off attenuation for the moment
!      xac = 1.
!      xf = 1.
!      xln = 1.

      return
      end
