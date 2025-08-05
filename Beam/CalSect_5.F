      subroutine CalSect_5 (nnodi,len_stiff, numNodes, numElems,
     &                      address, elemConn, xNode, yNode,
     &                      Gxy, globalK, globalF, Summary, ErrorFlag, 
     &                      gProc )
c --- primary function: compute section data
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"

      integer MNODI, eSize,nnodi
      parameter (MNODI = 9, eSize = 45)
      integer numNodes, numElems, elemConn(MNODI,numElems),
     &     address(numNodes),
     &     kel, k1, idof, jdof, kloc, kloc0, i, j, len_stiff,
     &     ErrorFlag, gProc
      double precision globalK(len_stiff),xNode(numNodes),
     &     yNode(numNodes),xyElem(2,MNODI),elemK(eSize),elemF(MNODI),
     &     globalF(numNodes), Gxy(numElems), Summary(*)

      call v_zero(globalK(1),len_stiff)
      call v_zero(globalF(1),numNodes)

      do kel = 1,numElems
         do k1 = 1,nnodi
            xyElem(1,k1) = xNode(elemConn(k1,kel))
            xyElem(2,k1) = yNode(elemConn(k1,kel))
         enddo
        
         call CalSect_4(nnodi,xyElem(1,1),Gxy(kel),
     &                  elemK(1),elemF(1),gProc,Summary(1),ErrorFlag)
         if (ErrorFlag .eq. 1) goto 999
         do i = 1,nnodi
            idof = elemConn(i,kel)
            globalF(idof) = globalF(idof)+elemF(i)
            do j = 1,nnodi
               jdof = elemConn(j,kel)
               if (jdof.le.idof) then
                   kloc=address(idof)-idof+jdof
                   if (i.le.j) then
                       kloc0=((j-1)*j)/2 +i
                   else
                       kloc0=((i-1)*i)/2 +j
                   endif
                   globalK(kloc)=globalK(kloc)+elemK(kloc0)
               endif
            enddo
         enddo
      enddo
c
      kloc = address(1)
      globalK(kloc) = globalK(kloc) + HUGE
      globalF(1) = 0.0d0
c
      do j = 1,numNodes
         kloc = address(j)
         if (globalK(kloc).eq.0.0d0) then
             globalK(kloc) = HUGE
         endif
      enddo
 999  continue

      return
      end



