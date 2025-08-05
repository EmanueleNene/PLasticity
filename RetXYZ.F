      subroutine RetXYZ (nnod,nodes,xyzang)
c --- primary function: get the coordinates for an element
c
c  input arguments:
c        nnod     (int,sc,in)       - number of nodes
c        nodes    (int,ar(nnod),in) - the element nodes
c
c  output arguments:
c        xyzang   (dp,ar(6,*),out)  - the coordinates and rotation angles
c
#include "impcom.inc"
c
      integer  ndgall

      integer  nnod,nodes(nnod),node,i,j,kdef
      double precision   xyzang(6,*)

c --- get element nodal coordinates
      do node = 1,nnod
c --- retrieve coordinates for this node
         if (nodes(node).gt.0) then
             kdef = ndgall(nodes(node),xyzang(1,node))
             if (kdef.eq.0) call v_zero (xyzang(1,node),6)
         else
             call v_zero (xyzang(1,node),6)
         endif
      enddo

      return
      end
