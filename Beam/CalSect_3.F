      subroutine CalSect_3(nnodi , xy, phiNodes,
     &                     dA, phi, phix, phiy, xI, yI)
c --- primary function: compute section data
c
#include "impcom.inc"
c
      integer nnodi,MNODI
      parameter (MNODI = 9)
      double precision xy(2,*), phi(4), dA(4), phix(4), phiy(4),
     &     parInt(MNODI,3), eJacob(4), spaDer(MNODI,2), phiNodes(*)
      double precision ZERO, ONE, deltaV, detJ, work, xI(4), yI(4)
      integer kintk, knode, k1, k2, k3
      parameter (ZERO = 0.0d0, ONE = 1.0d0)

c --- a nine-node 2x2 element for torsion stress function
c

      do kintk = 1,4
         call CalSect_6(kintk, nnodi, parInt(1,1))
         eJacob(1) = ZERO
         eJacob(2) = ZERO
         eJacob(3) = ZERO
         eJacob(4) = ZERO
         k3=0
         do k1 = 2,3
            do k2 = 1,2
               k3 = k3+1
               do knode = 1,nnodi
                   eJacob(k3)=eJacob(k3)+parInt(knode,k1)*xy(k2,knode)
               enddo
            enddo
         enddo
         deltaV = eJacob(1)*eJacob(4)-eJacob(2)*eJacob(3)
         detJ = ONE/deltaV
         work = eJacob(1)
         eJacob(1) = eJacob(4)*detJ
         eJacob(2) =-eJacob(2)*detJ
         eJacob(3) =-eJacob(3)*detJ
         eJacob(4) = work*detJ
         xI(kintk) = ZERO
         yI(kintk) = ZERO
         do knode=1,nnodi
            spaDer(knode,1)=
     &           eJacob(1)*parInt(knode,2)+eJacob(2)*parInt(knode,3)
            spaDer(knode,2)=
     &           eJacob(3)*parInt(knode,2)+eJacob(4)*parInt(knode,3)
            xI(kintk) = xI(kintk) + parInt(knode,1)*xy(1,knode)
            yI(kintk) = yI(kintk) + parInt(knode,1)*xy(2,knode)
         enddo
         dA(kintk) = deltaV
         phi(kintk) = ZERO
         phix(kintk) = ZERO
         phiy(kintk) = ZERO
         do k1 = 1,nnodi
            phi(kintk) = phi(kintk)+parInt(k1,1)*phiNodes(k1)
            phix(kintk) = phix(kintk)+spaDer(k1,1)*phiNodes(k1)
            phiy(kintk) = phiy(kintk)+spaDer(k1,2)*phiNodes(k1)
         enddo
      enddo
      
      return
      end




