      integer AAHLCD,AAHACD,AAHLCH,AAHACH,AAHLTP,AAHATP,AAHSLT,AAHSLN
     1       ,AAHSEL,AAHDSN,AAHAZR,AAHNPL,AAHNL1,AAHNZR,AAHNL2,AAHPZS
     1       ,AAHELT,AAHELN,AAHEDP,AAHEYR,AAHEMO,AAHEDA,AAHEHR,AAHEMI
     1       ,AAHESC,AAHLEC,AAHAEC,AAHITP,AAHNDT,AAHDEL,AAHAMX,AAHSYR
     1       ,AAHSMO,AAHSDA,AAHSHR,AAHSMI,AAHSSC,AAHNL3,AAHLDC,AAHADC
     1       ,AAHLLG,AAHALG,AAHNXR,AAHXRS

      parameter (MAHPZS=29)
      parameter (MAHLOG=51)
      parameter (MAHCOM=20)
      parameter (MAHXRS=21)

      parameter (AAHLCD=0)                   ! Station code
      parameter (AAHACD=AAHLCD+1)
      parameter (AAHLCH=AAHACD+2)            ! Channel
      parameter (AAHACH=AAHLCH+1)
      parameter (AAHLTP=AAHACH+2)            ! Type
      parameter (AAHATP=AAHLTP+1)
      parameter (AAHSLT=AAHATP+2)            ! Station coordinates
      parameter (AAHSLN=AAHSLT+1)
      parameter (AAHSEL=AAHSLN+1)
      parameter (AAHDSN=AAHSEL+1)            ! Digital sensitivity
      parameter (AAHAZR=AAHDSN+1)            ! A0
      parameter (AAHNPL=AAHAZR+1)            ! Poles and zeros
      parameter (AAHNL1=AAHNPL+1)
      parameter (AAHNZR=AAHNL1+1)
      parameter (AAHNL2=AAHNZR+1)
      parameter (AAHPZS=AAHNL2+1)
      parameter (AAHELT=AAHPZS+4*MAHPZS)     ! Event parameters
      parameter (AAHELN=AAHELT+1)
      parameter (AAHEDP=AAHELN+1)
      parameter (AAHEYR=AAHEDP+1)
      parameter (AAHEMO=AAHEYR+1)
      parameter (AAHEDA=AAHEMO+1)
      parameter (AAHEHR=AAHEDA+1)
      parameter (AAHEMI=AAHEHR+1)
      parameter (AAHESC=AAHEMI+1)
      parameter (AAHLEC=AAHESC+1)             ! Event comment
      parameter (AAHAEC=AAHLEC+1)
      parameter (AAHITP=AAHAEC+MAHCOM)        ! Data parameters
      parameter (AAHNDT=AAHITP+1)
      parameter (AAHDEL=AAHNDT+1)
      parameter (AAHAMX=AAHDEL+1)
      parameter (AAHSYR=AAHAMX+1)
      parameter (AAHSMO=AAHSYR+1)
      parameter (AAHSDA=AAHSMO+1)
      parameter (AAHSHR=AAHSDA+1)
      parameter (AAHSMI=AAHSHR+1)
      parameter (AAHSSC=AAHSMI+1)
      parameter (AAHNL3=AAHSSC+1)
      parameter (AAHLDC=AAHNL3+1)             ! Data comment
      parameter (AAHADC=AAHLDC+1)
      parameter (AAHLLG=AAHADC+MAHCOM)        ! AH Log
      parameter (AAHALG=AAHLLG+1)
      parameter (AAHNXR=AAHALG+MAHLOG)        ! Extras
      parameter (AAHXRS=AAHNXR+1)
      parameter (LAHHDR=AAHXRS+MAHXRS)

      parameter (XDIPA=1818.18)
