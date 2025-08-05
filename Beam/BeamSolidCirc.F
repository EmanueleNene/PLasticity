      subroutine BeamSolidCirc (numCirElems, numThruElems, rOuter,
     &                          ElemConn, rData, numNodes, numElems)
c --- primary function: prepare section data for circular cross section
c
#include "impcom.inc"
c
      double precision
     &     rData(*), rInner, rOuter, rMid, rHole, dtheta, theta,
     &     plusminus, sqCurRadius
      integer ElemConn(9,*),k2,k1,k0,ndiv,numCirElems,
     &        numElems, numNodes, numThruElems
      double precision EIGHT, ONE, TWO, ZERO
      parameter (EIGHT=8.0d0,ONE=1.0d0,TWO=2.0d0,ZERO=0.0d0)

#if defined(DEBUG)
      call SubBeg('BeamSolidCirc')
#endif

      if (numCirElems .lt. 8) numCirElems = 8
      ndiv = 2*numCirElems
      numNodes = 4*numThruElems*numCirElems+1
      dtheta = EIGHT*atan(ONE)/dble(ndiv)

c --- put node 1 at (0,0)
      rData(1) = ZERO
      rData(numNodes+1) = ZERO
      rHole = 0.0d0

c --- node numbering
      do k2 = 1,numThruElems
         theta = 0.0d0
         rInner = rHole
         sqCurRadius = dble(k2)/numThruElems
         rHole = rOuter*sqrt(sqCurRadius)
         rMid = (rInner+rHole)/TWO
         plusminus = 1.0d0
         do k1 = 1,ndiv
            k0 = 2*k1 + 2*(k2-1)*ndiv

c --- plusminus defines nodes alternatively around the circle
c --- in the +theta and -theta direction (performance issues)
c
            plusminus = -1.0d0*plusminus
            if (plusminus .gt. 0.0d0) then
                theta = theta+dtheta
            endif
            rData(k0) = rMid*cos(theta)
            rData(numNodes+k0) = plusminus*rMid*sin(theta)
            rData(k0+1) = rHole*cos(theta)
            rData(numNodes+k0+1) = plusminus*rHole*sin(theta)
         enddo
      enddo

c --- element connectivity
c --- define element 1 first, special case
c
      ElemConn(1,1) = 1
      ElemConn(2,1) = 3
      ElemConn(3,1) = 9
      ElemConn(4,1) = 1
      ElemConn(5,1) = 2
      ElemConn(6,1) = 5
      ElemConn(7,1) = 8
      ElemConn(8,1) = 1
      ElemConn(9,1) = 4
      rData(2*numNodes+1)    = ONE
c
      do k2 = 2,numThruElems
         k0 = (k2-1)*numCirElems + 1
         ElemConn(1,k0) = 2*(k2-2)*ndiv + 3
         ElemConn(4,k0) = ElemConn(1,k0) + 6
         ElemConn(8,k0) = ElemConn(1,k0) + 2
         ElemConn(2,k0) = 2*(k2-1)*ndiv + 3
         ElemConn(3,k0) = ElemConn(2,k0) + 6
         ElemConn(5,k0) = ElemConn(2,k0) - 1
         ElemConn(6,k0) = ElemConn(2,k0) + 2
         ElemConn(7,k0) = ElemConn(2,k0) + 5
         ElemConn(9,k0) = ElemConn(2,k0) + 1
         rData(2*numNodes+k0)    = ONE
      enddo

c --- plusminus is used here to keep track of elements
c --- being defined on the top or bottom of the section.
c --- at end of loop, plusminus=1 means even number of divisions.
c
      plusminus = -1.0d0
      do k1=2,numCirElems
         plusminus = -1.0d0*plusminus
c --- bottom
         if (plusminus .gt. 0.0d0) then
             ElemConn(1,k1) = 1
             ElemConn(2,k1) = 4*k1+3
             ElemConn(3,k1) = ElemConn(2,k1)-8
             ElemConn(4,k1) = 1
             ElemConn(5,k1) = ElemConn(2,k1)-1
             ElemConn(6,k1) = ElemConn(2,k1)-4
             ElemConn(7,k1) = ElemConn(2,k1)-9
             ElemConn(8,k1) = 1
             ElemConn(9,k1) = ElemConn(2,k1)-5
             rData(2*numNodes+k1)    = ONE
c --- top
         else
             ElemConn(1,k1) = 1
             ElemConn(2,k1) = 4*k1-3
             ElemConn(3,k1) = ElemConn(2,k1)+8
             ElemConn(4,k1) = 1
             ElemConn(5,k1) = ElemConn(2,k1)-1
             ElemConn(6,k1) = ElemConn(2,k1)+4
             ElemConn(7,k1) = ElemConn(2,k1)+7
             ElemConn(8,k1) = 1
             ElemConn(9,k1) = ElemConn(2,k1)+3
             rData(2*numNodes+k1)    = ONE
         endif
      enddo
      if (plusminus .gt. 0.0d0) then
          ElemConn(2,numCirElems) = ElemConn(2,numCirElems) - 2
          ElemConn(5,numCirElems) = ElemConn(5,numCirElems) - 2
      else
          ElemConn(3,numCirElems) = ElemConn(3,numCirElems) - 6
          ElemConn(7,numCirElems) = ElemConn(7,numCirElems) - 6
      endif
c
      do k2 = 2,numThruElems
         plusminus = -1.0d0
         do k1=2,numCirElems
            k0 = (k2-1)*numCirElems + k1
            plusminus = -1.0d0*plusminus
c --- bottom
            if (plusminus .gt. 0.0d0) then
                ElemConn(1,k0) = 2*(k2-2)*ndiv + 4*k1 + 3
                ElemConn(4,k0) = ElemConn(1,k0) - 8
                ElemConn(8,k0) = ElemConn(1,k0) - 4
                ElemConn(2,k0) = 2*(k2-1)*ndiv + 4*k1 + 3
                ElemConn(3,k0) = ElemConn(2,k0) - 8
                ElemConn(5,k0) = ElemConn(2,k0) - 1
                ElemConn(6,k0) = ElemConn(2,k0) - 4
                ElemConn(7,k0) = ElemConn(2,k0) - 9
                ElemConn(9,k0) = ElemConn(2,k0) - 5
                rData(2*numNodes+k0)    = ONE
c --- top
            else
                ElemConn(1,k0) = 2*(k2-2)*ndiv + 4*k1 - 3
                ElemConn(4,k0) = ElemConn(1,k0) + 8
                ElemConn(8,k0) = ElemConn(1,k0) + 4
                ElemConn(2,k0) = 2*(k2-1)*ndiv + 4*k1 - 3
                ElemConn(3,k0) = ElemConn(2,k0) + 8
                ElemConn(5,k0) = ElemConn(2,k0) - 1
                ElemConn(6,k0) = ElemConn(2,k0) + 4
                ElemConn(7,k0) = ElemConn(2,k0) + 7
                ElemConn(9,k0) = ElemConn(2,k0) + 3
                rData(2*numNodes+k0)    = ONE
            endif
         enddo
         if (plusminus .gt. 0.0d0) then
             ElemConn(1,numCirElems*k2) = ElemConn(1,numCirElems*k2) - 2
             ElemConn(2,numCirElems*k2) = ElemConn(2,numCirElems*k2) - 2
             ElemConn(5,numCirElems*k2) = ElemConn(5,numCirElems*k2) - 2
         else
             ElemConn(3,numCirElems*k2) = ElemConn(3,numCirElems*k2) - 6
             ElemConn(4,numCirElems*k2) = ElemConn(4,numCirElems*k2) - 6
             ElemConn(7,numCirElems*k2) = ElemConn(7,numCirElems*k2) - 6
         endif
      enddo
      numElems = numCirElems*numThruElems

#if defined(DEBUG)
      call SubEnd('BeamSolidCirc')
#endif

      return
      end
