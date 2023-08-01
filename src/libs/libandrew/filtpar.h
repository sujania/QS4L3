      parameter(JFOCHP=1)  ! cosine highpass        f1 f2
      parameter(JFOCLP=2)  ! cosine lowpass         f1 f2
      parameter(JFOBHP=3)  ! Butterworth highpass   f n
      parameter(JFOBLP=4)  ! Butterworth lowpass    f n
      parameter(JFOBBP=5)  ! Butterworth bandpass   f1 f2 n
      parameter(JFOPHP=6)  ! Pendulum highpass      f eta n
      parameter(JFOPLP=7)  ! Pendulum lowpass       f eta n
      parameter(JFOSPW=8)  ! Power of i omega       exponent
      parameter(JFODTS=9)  ! Differential t-star    dt* Tcorner
      parameter(JFODCI=10) ! Deconvolve instrument  none
      parameter(JFOHTR=11) ! Hilbert transform      none
      parameter(JFOSRO=12) ! Convolve SRO           none
      parameter(JFOWWS=13) ! Convolve WWSSN sp.     none
      parameter(JFOWWL=14) ! Convolve WWSSN lp.     none
      parameter(JFOCVF=15) ! Convolve with filter from file (none)
      parameter(JFOTSH=16) ! Time shift             dt
      parameter(JFOGSS=17) ! Gaussian               f0 fwidth

      parameter (MXFSTG=20)
      parameter (MXFPRM=16)
      parameter (MXFILT=200)
      common/filtpar/nflstg,ifltopt(MXFSTG),fltprm(MXFPRM,MXFSTG)
     1     ,lfilt,smpinf,tshff,filt(MXFILT)
