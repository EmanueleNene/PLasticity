      subroutine CalSect_9 (nnodi,gProc,numNodes, numElems,
     &                      elemConn, xNode, yNode, globalPhi,
     &                      ModWt, Results)
c --- primary function: calculate section data
c
#include "impcom.inc"
#include "sec_info.inc"
      integer MNODI, nnodi,gProc
      parameter (MNODI = 9)
      integer numNodes, numElems, elemConn(MNODI,numElems),
     &     kel, k1, i, j, numIntPts
      parameter (numIntPts = 4)
      double precision xNode(numNodes),
     &     yNode(numNodes),xyElem(2,MNODI),dA(4),phi(4),
     &     globalPhi(numNodes), elemPhi(9), phix(4),phiy(4),
     &     Results(B_NumResPts,*), xI(4), yI(4),ModWt(*)

      do kel = 1,numElems
         do k1 = 1,nnodi
            i = elemConn(k1,kel)
            xyElem(1,k1) = xNode(i)
            xyElem(2,k1) = yNode(i)
            elemPhi(k1) = globalPhi(elemConn(k1,kel))
         enddo
         call CalSect_3(nnodi, xyElem(1,1),elemPhi(1),dA(1),
     &                  phi(1), phix(1),phiy(1), xI(1), yI(1))
         j = (kel-1)*numIntPts
         do k1=1,numIntPts
            i = j+k1
            if(gProc.le.2) Results(2,i) = phi(k1)
            if (gProc.eq.1) then
                Results(1,i) = dA(k1)
                Results(3,i) = phix(k1)
                Results(4,i) = phiy(k1)
                Results(5,i) = xI(k1)
                Results(6,i) = yI(k1)
                Results(7,i) = ModWt(kel)*dA(k1)
            elseif(gProc.eq.2) then
                Results(3,i) = phix(k1)
                Results(4,i) = phiy(k1)
            elseif(gProc.eq.3) then
                Results(8,i) = phix(k1)
                Results(9,i) = phiy(k1)
                Results(7,i) = ModWt(kel)*dA(k1)
            elseif(gProc.eq.4) then
                Results(10,i) = phix(k1)
                Results(11,i) = phiy(k1)
            endif
         enddo
      enddo

      return
      end
