      subroutine LineGeomTrans (xyz,nx,tr)
c --- primary function: geometric transform for 3-d line
c
c  input arguments:
c     xyz      (dp,ar(nx,*),in)  - coordinates (x,y,z down)
c     nx       (int,sc,in)       - row dimension of xyz array
c
c  output arguments:
c     tr       (dp,ar(3,3),out)  - transformation matrix
c
#include "impcom.inc"
c
      integer nx
c
      double precision
     & xyz(nx,*),tr(3,3),d12(3)

c --- compute current length and length components
      call v1mv2_v3 (xyz(1,2),xyz(1,1),d12(1),3)
      call TransBeam (d12(1),tr(1,1))
      return
      end
