      subroutine delaz(eplat,eplong,stlat,stlong,delta,azep,azst)
      data geoco,hpi,twopi,rad,reprad/.993277,1.5707963,
     16.2831853,.0174533,57.29578/
c     tan(x)=sin(x)/cos(x)
      arcos(x)=atan2(sqrt(1.-x*x),x)
      if(eplat.eq.stlat.and.eplong.eq.stlong) then
        delta=0.
        azep=0.
        azst=0.
        return
      endif
      el=atan(geoco*tan(eplat*rad))
      el=hpi-el
      stl=atan(geoco*tan(stlat*rad))
      stl=hpi-stl
      elon=eplong*rad
      slon=stlong*rad
      as=cos(stl)
      bs=sin(stl)
      cs=cos(slon)
      ds=sin(slon)
      a=cos(el)
      b=sin(el)
      c=cos(elon)
      d=sin(elon)
      cdel=a*as+b*bs*(c*cs+d*ds)
      if(abs(cdel).gt.1.) cdel=sign(1.,cdel)
      delt=arcos(cdel)
      delta=delt*reprad
      sdel=sin(delt)
      caze=(as-a*cdel)/(sdel*b)
      if(abs(caze).gt.1.) caze=sign(1.,caze)
      aze=arcos(caze)
      if(bs.gt.0.) cazs=(a-as*cdel)/(bs*sdel)
      if(bs.eq.0.) cazs=sign(1.,cazs)
      if(abs(cazs).gt.1.) cazs=sign(1.,cazs)
      azs=arcos(cazs)
      dif=ds*c-cs*d
      if(dif.lt.0.) aze=twopi-aze
      azep=reprad*aze
      if(dif.gt.0.) azs=twopi-azs
      azst=reprad*azs
      return
      end
