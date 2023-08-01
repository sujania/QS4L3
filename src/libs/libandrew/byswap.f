      subroutine byswap8(iar,nel)
      integer*4 iar(*)
      integer*4 it
!#if ( defined(MachineA) || defined(Machinel) || defined(Machinec) || defined(Machinep) )
      call byswap4(iar,nel*2)
      k=0
      do i=1,nel
        k=k+1
        k1=k+1
        it=iar(k)
        iar(k)=iar(k1)
        iar(k1)=it
        k=k1
      enddo
!#endif
      return
      end

c------------------------------------
 
 
      subroutine byswap4(iar,nel)
      integer*4 iar(*)
      character*1 c(4),csave
      integer*4 i4
      equivalence (c,i4)
!#if ( defined(MachineA) || defined(Machinel) || defined(Machinec) || defined(Machinep) )
      do i=1,nel
        i4=iar(i)
        csave=c(1)
        c(1)=c(4)
        c(4)=csave
        csave=c(2)
        c(2)=c(3)
        c(3)=csave
        iar(i)=i4
      enddo
!#endif
      return
      end
c----------------------------------
      subroutine byswap2(iar,nel)
      integer*2 iar(*)
      character*1 c(2),csave
      integer*2 i2
      equivalence (c,i2)
!#if ( defined(MachineA) || defined(Machinel) || defined(Machinec) || defined(Machinep) )
      do i=1,nel
        i2=iar(i)
        csave=c(1)
        c(1)=c(2)
        c(2)=csave
        iar(i)=i2
      enddo
      return
!#endif
      end

        
