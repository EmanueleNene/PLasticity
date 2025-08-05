      subroutine BeamPipe (numCirElems, rInner, rOuter, rData,
     &                     ElemConn, numNodes, numElems)
c --- primary function: calculate section data for pipe cross section
c
#include "impcom.inc"
c
      double precision
     &     rData(*), rInner, rOuter, rMid,  dtheta, theta, plusminus
      integer ElemConn(9,*),k1,k0,ndiv,numCirElems,
     &        numElems, numNodes
      double precision EIGHT, ONE, TWO
      parameter (EIGHT=8.0d0,ONE=1.0d0,TWO=2.0d0)

#if defined(DEBUG)
      call SubBeg('BeamPipe')
#endif

      if (numCirElems .lt. 8) numCirElems = 8
      ndiv = 2*numCirElems
      numNodes = ndiv*3
      theta = 0.0d0
      dtheta = EIGHT*atan(ONE)/dble(ndiv)
      rMid = (rInner+rOuter)/TWO
      plusminus = 1.0d0
      do k1 = 1,ndiv
c --- plusminus defines nodes alternatively around the circle
c --- in the +theta and -theta direction (performance issues)
         plusminus = -1.0d0*plusminus
         k0 = (k1-1)*3
         if (plusminus .gt. 0.0d0) then
             theta = theta+dtheta
         endif
         rData(k0+1) = rInner*cos(theta)
         rData(numNodes+k0+1) = plusminus*rInner*sin(theta)
         rData(k0+2) = rMid*cos(theta)
         rData(numNodes+k0+2) = plusminus*rMid*sin(theta)
         rData(k0+3) = rOuter*cos(theta)
         rData(numNodes+k0+3) = plusminus*rOuter*sin(theta)
      enddo
      
c --- element connectivity
c --- define element 1 first, special case
c
      ElemConn(1,1) = 1
      ElemConn(2,1) = 3
      ElemConn(3,1) = 12
      ElemConn(4,1) = 10
      ElemConn(5,1) = 2
      ElemConn(6,1) = 6
      ElemConn(7,1) = 11
      ElemConn(8,1) = 4
      ElemConn(9,1) = 5
      rData(2*numNodes+1) = ONE
c
      plusminus = -1.0d0
      do k1=2,numCirElems

c --- plusminus is used here to keep track of elements
c --- being defined on the top or bottom of the section.
c --- at end of loop, plusminus=1 means even number of divisions.
c
         plusminus = -1.0d0*plusminus
         if (plusminus .gt. 0.0d0) then
             k0 = 6*k1+1
             ElemConn(1,k1) = k0
             ElemConn(2,k1) = ElemConn(1,k1)+2
             ElemConn(3,k1) = ElemConn(1,k1)-10
             ElemConn(4,k1) = ElemConn(1,k1)-12
             ElemConn(5,k1) = ElemConn(1,k1)+1
             ElemConn(6,k1) = ElemConn(1,k1)-4
             ElemConn(7,k1) = ElemConn(1,k1)-11
             ElemConn(8,k1) = ElemConn(1,k1)-6
             ElemConn(9,k1) = ElemConn(1,k1)-5
         else
             k0 = 6*k1-8
             ElemConn(1,k1) = k0
             ElemConn(2,k1) = ElemConn(1,k1)+2
             ElemConn(3,k1) = ElemConn(1,k1)+14
             ElemConn(4,k1) = ElemConn(1,k1)+12
             ElemConn(5,k1) = ElemConn(1,k1)+1
             ElemConn(6,k1) = ElemConn(1,k1)+8
             ElemConn(7,k1) = ElemConn(1,k1)+13
             ElemConn(8,k1) = ElemConn(1,k1)+6
             ElemConn(9,k1) = ElemConn(1,k1)+7
         endif
         rData(2*numNodes+k1) = ONE
      enddo
      if (plusminus .gt. 0.0d0) then
          ElemConn(1,numCirElems) = ElemConn(1,numCirElems) - 3
          ElemConn(2,numCirElems) = ElemConn(2,numCirElems) - 3
          ElemConn(5,numCirElems) = ElemConn(5,numCirElems) - 3
      else
          ElemConn(3,numCirElems) = ElemConn(3,numCirElems) - 9
          ElemConn(4,numCirElems) = ElemConn(4,numCirElems) - 9
          ElemConn(7,numCirElems) = ElemConn(7,numCirElems) - 9
      endif
      numElems = numCirElems

#if defined(DEBUG)
      call SubEnd('BeamPipe')
#endif

      return
      end
