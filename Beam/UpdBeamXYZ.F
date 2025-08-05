      subroutine UpdBeamXYZ (nnod,u,nr,xyz,nx,xyzup)
c --- primary function: update elem nodal coords for large deformation
c
c  input arguments:
c     nnod     (int,sc,in)          - number of nodes
c     u        (dp,ar(nr),in)       - displacement vector
c     nr       (int,sc,in)          - size of the u vector
c     xyz      (dp,ar(nx,nnod),in)  - coordinates to be updated
c     nx       (int,sc,in)          - row size of xy
c
c  output arguments:
c     xyzup    (dp,ar(3,nnod),out)  - updated coordinates
c
#include "impcom.inc"
c
      integer nnod,nr,nx, npr,i,j
      double precision u(nr),xyz(nx,nnod),xyzup(3,nnod)

#if defined(DEBUG)
      call SubBeg('UpdBeamXYZ')
#endif

c --- determine the number of dofs per node
      npr = nr/nnod
c --- assume u vector has at least ux,uy,uz in its DOF order
      j = 1
      do i = 1,nnod
         call v1pv2_v3 (xyz(1,i),u(j),xyzup(1,i),3)
         j = j + npr
      enddo

#if defined(DEBUG)
      call SubEnd('UpdBeamXYZ')
#endif

      return
      end
