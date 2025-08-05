      subroutine BeamToBeam_2 (xa,ya,za,xb,yb,zb,alpha,locu,locv,
     &                         u0,v0,fifoa,fifob,dd,ac,bc,dfa,dfb,kdf)
c --- primary function: beam-beam contact calculations (see parent function)
c
#include "impcom.inc"

      double precision xa(7),ya(7),za(7),xb(7),yb(7),zb(7),alpha,u0,v0
      integer locu,locv,kdf

      double precision dd,fifoa(7),fifob(7),dfa(7),dfb(7),
     &                 ac(3),bc(3),alpha0,u,v

      integer i

      alpha0 = 1.d0-2.d0*alpha
      if (locu.ne.-1.and.locu.ne.2.and.u0.lt.-alpha0) then
c --- left side
          u = 0.5d0*(u0+1.d0)/alpha
          fifoa(1) = (8.d0*alpha*alpha-6.d0*alpha)*u*u-6.d0*alpha+4.d0
          fifoa(2) = alpha*(2.d0*alpha*u*u*u+
     &                     (4.d0*alpha-1.d0)*u*u+
     &                     (2.d0*alpha-2.d0)*u -1.d0)
          fifoa(3) = alpha*(-4.d0*alpha*u*u*u+
     &                     (-8.d0*alpha+4.d0)*u*u+
     &                     (-4.d0*alpha+8.d0)*u +4.d0)
          fifoa(4) = alpha*(-2.d0*alpha*u*u*u+
     &                     (4.d0*alpha-1.d0)*u*u+
     &                     (-2.d0*alpha+2.d0)*u -1.d0)
          fifoa(5) = 0.d0
          fifoa(6) = alpha*(4.d0*alpha*u*u*u+
     &                     (-8.d0*alpha+4.d0)*u*u+
     &                     (4.d0*alpha-8.d0)*u +4.d0)
          fifoa(7) = 0.d0
          do i = 1,7
             fifoa(i) = 0.25d0*fifoa(i)
          enddo
          if (kdf.eq.1) then
              dfa(1) = (16.d0*alpha*alpha-12.d0*alpha)*u
              dfa(2) = alpha*(6.d0*alpha*u*u+
     &                       (8.d0*alpha-2.d0)*u+
     &                       (2.d0*alpha-2.d0))
              dfa(3) = alpha*(-12.d0*alpha*u*u+
     &                       (-16.d0*alpha+8.d0)*u+
     &                       (-4.d0*alpha+8.d0))
              dfa(4) = alpha*(-6.d0*alpha*u*u+
     &                       (8.d0*alpha-2.d0)*u+
     &                       (-2.d0*alpha+2.d0))
              dfa(5) = 0.d0
              dfa(6) = alpha*(12.d0*alpha*u*u+
     &                       (-16.d0*alpha+8.d0)*u+
     &                       (4.d0*alpha-8.d0))
              dfa(7) = 0.d0
              do i = 1,7
                 dfa(i) = 0.25d0*dfa(i)
              enddo
          endif

      elseif (locu.ne.1.and.locu.ne.2.and.u0.gt.alpha0) then
c --- right side
          u = 0.5d0*(u0-1.d0)/alpha
          fifoa(1) = alpha*(-2.d0*alpha*u*u*u+
     &                     (4.d0*alpha-1.d0)*u*u+
     &                     (-2.d0*alpha+2.d0)*u -1.d0)
          fifoa(2) = (8.d0*alpha*alpha-6.d0*alpha)*u*u-6.d0*alpha+4.d0
          fifoa(3) = alpha*(4.d0*alpha*u*u*u+
     &                     (-8.d0*alpha+4.d0)*u*u+
     &                     (4.d0*alpha-8.d0)*u +4.d0)
          fifoa(4) = 0.d0
          fifoa(5) = alpha*(2.d0*alpha*u*u*u+
     &                     (4.d0*alpha-1.d0)*u*u+
     &                     (2.d0*alpha-2.d0)*u -1.d0)
          fifoa(6) = 0.d0
          fifoa(7) = alpha*(-4.d0*alpha*u*u*u+
     &                     (-8.d0*alpha+4.d0)*u*u+
     &                     (-4.d0*alpha+8.d0)*u +4.d0)

          do i = 1,7
             fifoa(i) = 0.25d0*fifoa(i)
          enddo

          if (kdf.eq.1) then
              dfa(1) = alpha*(-6.d0*alpha*u*u+
     &                       (8.d0*alpha-2.d0)*u+
     &                       (-2.d0*alpha+2.d0))
              dfa(2) = (16.d0*alpha*alpha-12.d0*alpha)*u
              dfa(3) = alpha*(12.d0*alpha*u*u+
     &                       (-16.d0*alpha+8.d0)*u+
     &                       (4.d0*alpha-8.d0))
              dfa(4) = 0.d0
              dfa(5) = alpha*(6.d0*alpha*u*u+
     &                       (8.d0*alpha-2.d0)*u+
     &                       (2.d0*alpha-2.d0))
              dfa(6) = 0.d0
              dfa(7) = alpha*(-12.d0*alpha*u*u+
     &                       (-16.d0*alpha+8.d0)*u+
     &                       (-4.d0*alpha+8.d0))
              do i = 1,7
                 dfa(i) = 0.25d0*dfa(i)
              enddo
          endif
      else
          u = u0
          fifoa(1) = 0.5d0*u*(u-1.d0)
          fifoa(2) = 0.5d0*u*(u+1.d0)
          fifoa(3) = 1.d0-u*u
          fifoa(4) = 0.d0
          fifoa(5) = 0.d0
          fifoa(6) = 0.d0
          fifoa(7) = 0.d0
          if (kdf.eq.1) then
              dfa(1) = 0.5d0*(2.d0*u-1.d0)
              dfa(2) = 0.5d0*(2.d0*u+1.d0)
              dfa(3) = -2.d0*u
              dfa(4) = 0.d0
              dfa(5) = 0.d0
              dfa(6) = 0.d0
              dfa(7) = 0.d0
          endif

      endif

      if (locv.ne.-1.and.locv.ne.2.and.v0.lt.-alpha0) then
c --- left side
          v = 0.5d0*(v0+1.d0)/alpha
          fifob(1) = (8.d0*alpha*alpha-6.d0*alpha)*v*v-6.d0*alpha+4.d0
          fifob(2) = alpha*(2.d0*alpha*v*v*v+
     &                     (4.d0*alpha-1.d0)*v*v+
     &                     (2.d0*alpha-2.d0)*v -1.d0)
          fifob(3) = alpha*(-4.d0*alpha*v*v*v+
     &                     (-8.d0*alpha+4.d0)*v*v+
     &                     (-4.d0*alpha+8.d0)*v +4.d0)
          fifob(4) = alpha*(-2.d0*alpha*v*v*v+
     &                     (4.d0*alpha-1.d0)*v*v+
     &                     (-2.d0*alpha+2.d0)*v -1.d0)
          fifob(5) = 0.d0
          fifob(6) = alpha*(4.d0*alpha*v*v*v+
     &                     (-8.d0*alpha+4.d0)*v*v+
     &                     (4.d0*alpha-8.d0)*v +4.d0)
          fifob(7) = 0.d0
          do i = 1,7
             fifob(i) = 0.25d0*fifob(i)
          enddo
          if (kdf.eq.1) then
              dfb(1) = (16.d0*alpha*alpha-12.d0*alpha)*v
              dfb(2) = alpha*(6.d0*alpha*v*v+
     &                       (8.d0*alpha-2.d0)*v+
     &                       (2.d0*alpha-2.d0))
              dfb(3) = alpha*(-12.d0*alpha*v*v+
     &                       (-16.d0*alpha+8.d0)*v+
     &                       (-4.d0*alpha+8.d0))
              dfb(4) = alpha*(-6.d0*alpha*v*v+
     &                       (8.d0*alpha-2.d0)*v+
     &                       (-2.d0*alpha+2.d0))
              dfb(5) = 0.d0
              dfb(6) = alpha*(12.d0*alpha*v*v+
     &                       (-16.d0*alpha+8.d0)*v+
     &                       (4.d0*alpha-8.d0))
              dfb(7) = 0.d0
              do i = 1,7
                 dfb(i) = 0.25d0*dfb(i)
              enddo
          endif
      elseif (locv.ne.1.and.locv.ne.2.and.v0.gt.alpha0) then
c --- right side
          v = 0.5d0*(v0-1.d0)/alpha
          fifob(1) = alpha*(-2.d0*alpha*v*v*v+
     &                     (4.d0*alpha-1.d0)*v*v+
     &                     (-2.d0*alpha+2.d0)*v -1.d0)
          fifob(2) = (8.d0*alpha*alpha-6.d0*alpha)*v*v-6.d0*alpha+4.d0
          fifob(3) = alpha*(4.d0*alpha*v*v*v+
     &                     (-8.d0*alpha+4.d0)*v*v+
     &                     (4.d0*alpha-8.d0)*v +4.d0)
          fifob(4) = 0.d0
          fifob(5) = alpha*(2.d0*alpha*v*v*v+
     &                     (4.d0*alpha-1.d0)*v*v+
     &                     (2.d0*alpha-2.d0)*v -1.d0)
          fifob(6) = 0.d0
          fifob(7) = alpha*(-4.d0*alpha*v*v*v+
     &                     (-8.d0*alpha+4.d0)*v*v+
     &                     (-4.d0*alpha+8.d0)*v +4.d0)
          do i = 1,7
             fifob(i) = 0.25d0*fifob(i)
          enddo
          if (kdf.eq.1) then
              dfb(1) = alpha*(-6.d0*alpha*v*v+
     &                       (8.d0*alpha-2.d0)*v+
     &                       (-2.d0*alpha+2.d0))
              dfb(2) = (16.d0*alpha*alpha-12.d0*alpha)*v
              dfb(3) = alpha*(12.d0*alpha*v*v+
     &                       (-16.d0*alpha+8.d0)*v+
     &                       (4.d0*alpha-8.d0))
              dfb(4) = 0.d0
              dfb(5) = alpha*(6.d0*alpha*v*v+
     &                       (8.d0*alpha-2.d0)*v+
     &                       (2.d0*alpha-2.d0))
              dfb(6) = 0.d0
              dfb(7) = alpha*(-12.d0*alpha*v*v+
     &                       (-16.d0*alpha+8.d0)*v+
     &                       (-4.d0*alpha+8.d0))
              do i = 1,7
                 dfb(i) = 0.25d0*dfb(i)
              enddo
          endif
      else
          v = v0
          fifob(1) = 0.5d0*v*(v-1.d0)
          fifob(2) = 0.5d0*v*(v+1.d0)
          fifob(3) = 1.d0-v*v
          fifob(4) = 0.d0
          fifob(5) = 0.d0
          fifob(6) = 0.d0
          fifob(7) = 0.d0
          if (kdf.eq.1) then
              dfb(1) = 0.5d0*(2.d0*v-1.d0)
              dfb(2) = 0.5d0*(2.d0*v+1.d0)
              dfb(3) = -2.d0*v
              dfb(4) = 0.d0
              dfb(5) = 0.d0
              dfb(6) = 0.d0
              dfb(7) = 0.d0
          endif
      endif

      ac(1) = fifoa(1)*xa(1)
      ac(2) = fifoa(1)*ya(1)
      ac(3) = fifoa(1)*za(1)
      bc(1) = fifob(1)*xb(1)
      bc(2) = fifob(1)*yb(1)
      bc(3) = fifob(1)*zb(1)
      do i = 2,7
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
