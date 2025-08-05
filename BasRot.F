      subroutine BasRot(vn_t,vn_tpdt,deltaPhi,numVectors)
c --- primary function: rotate a basis vector
c
#include "impcom.inc"
#include "constants.inc"
c
      double precision vn_t(3,*),vn_tpdt(3,*),deltaPhi(3)
      double precision Direction(3),t2,ONE,EPS,
     &                 Magnitude,cosPhi,sinPhi,workVar
      integer  kvec,numVectors
      parameter (ONE = 1.0d0, EPS = 1.0d-6)

      Magnitude = sqrt(deltaPhi(1)**2+deltaPhi(2)**2+deltaPhi(3)**2)
      if (Magnitude .le. EPS) then
          do kvec=1,numVectors
             vn_tpdt(1,kvec) = vn_t(1,kvec)
     &           - deltaPhi(3)*vn_t(2,kvec)+deltaPhi(2)*vn_t(3,kvec)
             vn_tpdt(2,kvec) = vn_t(2,kvec)
     &           + deltaPhi(3)*vn_t(1,kvec)-deltaPhi(1)*vn_t(3,kvec)
             vn_tpdt(3,kvec) = vn_t(3,kvec)
     &           - deltaPhi(2)*vn_t(1,kvec)+deltaPhi(1)*vn_t(2,kvec)
          enddo
      else
          cosPhi = cos(Magnitude)
          sinPhi = sin(Magnitude)
          workVar = ONE/Magnitude
          Direction(1) = workVar*deltaPhi(1)
          Direction(2) = workVar*deltaPhi(2)
          Direction(3) = workVar*deltaPhi(3)
          do kvec=1,numVectors
             t2 = (1.0d0-cosPhi)*(
     &            Direction(1)*vn_t(1,kvec) +
     &            Direction(2)*vn_t(2,kvec) +
     &            Direction(3)*vn_t(3,kvec))
             vn_tpdt(1,kvec) = t2*Direction(1) +
     &            cosPhi*vn_t(1,kvec) +
     &            sinPhi*(Direction(2)*vn_t(3,kvec) -
     &                    Direction(3)*vn_t(2,kvec))
             vn_tpdt(2,kvec) = t2*Direction(2) +
     &            cosPhi*vn_t(2,kvec) +
     &            sinPhi*(Direction(3)*vn_t(1,kvec) -
     &                    Direction(1)*vn_t(3,kvec))
             vn_tpdt(3,kvec) = t2*Direction(3) +
     &            cosPhi*vn_t(3,kvec) +
     &            sinPhi*(Direction(1)*vn_t(2,kvec)-
     &                    Direction(2)*vn_t(1,kvec))
          enddo
      endif

      return
      end
