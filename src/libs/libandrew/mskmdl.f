      subroutine mskmdl(xin,npin,lmxin,mskin,lskin
     1                 ,xou,npou,lmxou,mskou,lskou)
c     APV fix: these were all dimensioned (1) or (0:1), was breaking bounds check
      dimension
     1  xin(*),mskin(*),lskin(0:*)
     1 ,xou(*),mskou(*),lskou(0:*)

      lmax=max0(lmxin,lmxou)
      nmax=max0(npin,npou)
      kin=0
      kou=0
      do 10 n=1,nmax
      do 10 l=0,lmax
      do 10 m=0,2*l

      if(l.le.lmxin.and.n.le.npin
     1  .and.(mskin(1).eq.-1.or.mskin(n).ne.0)
     1  .and.(lskin(0).eq.-1.or.lskin(l).ne.0)) then
        kin=1+kin
        val=xin(kin)
      else
        val=0.
      endif

      if(l.le.lmxou.and.n.le.npou
     1  .and.(mskou(1).eq.-1.or.mskou(n).ne.0)
     1  .and.(lskou(0).eq.-1.or.lskou(l).ne.0)) then
        kou=1+kou
        xou(kou)=val
      endif

   10 continue

      return
      end
