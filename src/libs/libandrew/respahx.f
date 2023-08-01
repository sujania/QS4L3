      complex function respahx(om,headr)
      real*4 headr(0:*)
      include '../libdbdo/ahheader.h'

      respahx=headr(AAHDSN)*headr(AAHAZR)
      k=AAHPZS+2
      do 10 i=1,ifix(headr(AAHNZR))
         respahx=respahx*cmplx(-headr(k),om-headr(k+1))
         k=k+4
 10   continue
      k=AAHPZS
      do 20 i=1,ifix(headr(AAHNPL))
         respahx=respahx/cmplx(-headr(k),om-headr(k+1))
         k=k+4
 20   continue
      return
      end
