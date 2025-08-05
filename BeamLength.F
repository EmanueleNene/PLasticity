      function BeamLength (numElemNodes,x,Transpose,nr)
c --- primary function: calculate beam length
c
c  input arguments:
c     numElemNodes (int,sc,in)        - number of element nodes
c     x            (dp,ar(*),in)      - nodal coords
c     Transpose    (int,sc,in)        - need transpose? (1-yes)
c     nr           (int,sc,in)        - number of rows
c
c  output arguments:
c     BeamLength   (dp,sc,out)        - beam length
c
#include "impcom.inc"
c
      integer
     &     numElemNodes,Transpose,nr
      
      double precision 
     &     x(nr*3), BeamLength, CalDist2P
      
      integer
     &     domInt, i, nod

      double precision
     &     Gauss3(3), Gaussw3(3), shpe(3,2), dxdg(3), dL, Length,
     &     Coords(nr,numElemNodes)

      double precision ONE, ZERO, HALF, TWO
      parameter (ZERO = 0.0d0, HALF = 0.5D0, ONE = 1.0d0, TWO = 2.0d0)

c --- Gaussian points
      data Gauss3/-0.7745966692415d0,0.0d0,0.7745966692415d0/
c --- Gaussian weights
      data Gaussw3/0.555555555556d0,0.888888888889d0,0.555555555556d0/

      Length = ZERO

      if (numElemNodes.lt.2 .or. numElemNodes.gt.3) then 
c --- wrong data, return zero
          goto 990
      endif
      
      if (Transpose.eq.1) then
          do nod = 1,numElemNodes
             do i = 1,3
                Coords(i,nod) = x((i-1)*numElemNodes+nod)
             enddo
          enddo
      else
          call v_move(x(1),Coords(1,1),nr*numElemNodes)
      endif
      
      if (numElemNodes.eq.2) then
c --- for 2-node beam
          Length = CalDist2P(Coords(1,1), Coords(1,2))
          goto 990
      endif

c --- integration loop
      do domInt = 1,3
         shpe(1,1)=HALF*(Gauss3(domInt)-ONE)*Gauss3(domInt)
         shpe(2,1)=HALF*(Gauss3(domInt)+ONE)*Gauss3(domInt)
         shpe(3,1)=ONE-Gauss3(domInt)*Gauss3(domInt)
         shpe(1,2)=Gauss3(domInt)-HALF
         shpe(2,2)=Gauss3(domInt)+HALF
         shpe(3,2)=-TWO*Gauss3(domInt)
         dL = ZERO
         call v_zero (dxdg(1),3)
         do i = 1,3
            do nod = 1,numElemNodes
               dxdg(i) = dxdg(i)+shpe(nod,2)*Coords(i,nod)
            enddo
            dL = dL+dxdg(i)*dxdg(i)
         enddo
         Length = Length + sqrt(dL)*Gaussw3(domInt)
      enddo
 
 990  BeamLength = Length
      
      return
      end
