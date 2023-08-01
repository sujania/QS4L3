      subroutine rdmdl(luhmdl)
      common/hetmdl/ifani,lmaxm,lmaxm1,leny,npmm,pertm(588),bmdl(588)
      rewind luhmdl
      read(luhmdl,702) ifani,ifmino,lmaxm
  702 format(3i5)
      lmaxm1=lmaxm+1
      leny=lmaxm1**2
      npmm=12+(ifani-1)*5
      do 703 ip=1,npmm
      ind=1+(ip-1)*leny
      do 704 l1=1,lmaxm1
      ind1=ind+2*(l1-1)
      read(luhmdl,705)(pertm(i),i=ind,ind1)
  705 format(11e12.5)
  704 ind=ind1+1
  703 continue
      return
      end
