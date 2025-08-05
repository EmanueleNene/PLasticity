      subroutine BeamToBeam_2Fast (xa,ya,za,xb,yb,zb,
     &                             u0,v0,fifoa,fifob,dd,
     &                             ac,bc,dfa,dfb,kdf)
c --- primary function: beam-beam contact calculations (see parent function)
c
#include "impcom.inc"
c
      double precision xa(3),ya(3),za(3),xb(3),yb(3),zb(3),u0,v0
      integer kdf

      double precision HALF,ONE,TWO,dd,fifoa(3),fifob(3),
     &                 dfa(3),dfb(3),ac(3),bc(3),u,v
      
      parameter (HALF=0.5d0,ONE=1.0d0,TWO=2.0d0)
      integer i

      u = u0
      fifoa(1) = HALF*u*(u-ONE)
      fifoa(2) = HALF*u*(u+ONE)
      fifoa(3) = ONE-u*u
      if (kdf.eq.1) then
          dfa(1) = HALF*(TWO*u-ONE)
          dfa(2) = HALF*(TWO*u+ONE)
          dfa(3) = -TWO*u
      endif

      v = v0
      fifob(1) = HALF*v*(v-ONE)
      fifob(2) = HALF*v*(v+ONE)
      fifob(3) = ONE-v*v
      if (kdf.eq.1) then
          dfb(1) = HALF*(TWO*v-ONE)
          dfb(2) = HALF*(TWO*v+ONE)
          dfb(3) = -TWO*v
      endif

      ac(1) = fifoa(1)*xa(1)
      ac(2) = fifoa(1)*ya(1)
      ac(3) = fifoa(1)*za(1)
      bc(1) = fifob(1)*xb(1)
      bc(2) = fifob(1)*yb(1)
      bc(3) = fifob(1)*zb(1)
      do i = 2,3
         ac(1) = ac(1) + fifoa(i)*xa(i)
         ac(2) = ac(2) + fifoa(i)*ya(i)
         ac(3) = ac(3) + fifoa(i)*za(i)
         bc(1) = bc(1) + fifob(i)*xb(i)
         bc(2) = bc(2) + fifob(i)*yb(i)
         bc(3) = bc(3) + fifob(i)*zb(i)
      enddo
      dd = (ac(1)-bc(1))*(ac(1)-bc(1))+
     &     (ac(2)-bc(2))*(ac(2)-bc(2))+
     &     (ac(3)-bc(3))*(ac(3)-bc(3))
      dd = sqrt(dd)

      return
      end
