      subroutine CalUdef (xyz,nx,nnod,trr,utot,udef,kerr)
c --- primary function: compute deformational displacements for beams
c
c     compute the deformational displacements
c     ref: rankin and brogan, 'an element independent corotational procedure
c          for treatment of large rotations', j.p.v. technology, v108 p165 may 86
c     note 1: define coords and translational disps wrt node 1
c     note 2: udef is in global coordinates
c
c  input arguments:
c     xyz    (dp,ar(nx,nnod),in)     - original element coordinates (in global cs) (x,y,z down)
c     nx     (int,sc,in)             - row dimension on xyz
c     nnod   (int,sc,in)             - number of nodes
c     trr    (dp,ar(3,3),in)         - rotation matrix
c     utot   (dp,ar(6*nnod),in)      - total displacements (with the rotations being Euler values)
c
c  output arguments:
c     udef   (dp,ar(6*nnod),out)     - deformational displacements (unrotated)
c                                      udef acts on current (rotated) orientation of the element
c     kerr   (int,sc,out)            - error flag, set if deformational rotation exceeds 30 degrees
c
#include "impcom.inc"
c
      integer nx,nnod,kerr,i,j
      double precision xyz(nx,nnod),trr(3,3),utot(6*nnod),udef(6*nnod),
     &                 urel(3),xrel(3),tx(3),trot(3,3),deltr(3,3),con
 
      kerr = 0
c --- assume the c.s. origin in the first node
      call v_zero (udef(1),3)
c --- translational DOF (eqn 2)
      call v1mv2_v3 (xyz(1,2),xyz(1,1),xrel(1),3)
      call v1mv2_v3 (utot(7),utot(1),urel(1),3)
      call v1pv2_v3 (xrel(1),urel(1),udef(7),3)
      call Atxv1_v2 (trr(1,1),xrel(1),tx(1), 3,3)
      call v1mv2_v1 (udef(7),tx(1),3)

      j = 1
      do i = 1,nnod
c --- rotational DOF (eqns 24 and 37)
         call CalRotTr (utot(j+3),trot(1,1))
         call AtxB_C (trot(1,1),trr(1,1),deltr(1,1), 3,3,3, 3,3,3)
         con = 1.0d0 + deltr(1,1) + deltr(2,2) + deltr(3,3)
         if (con.le.0.0d0) then
             con = 0.0d0
         else
             con = sqrt(1.0d0/con)
         endif
         udef(j+3) = con*(deltr(3,2) - deltr(2,3))
         udef(j+4) = con*(deltr(1,3) - deltr(3,1))
         udef(j+5) = con*(deltr(2,1) - deltr(1,2))
c --- check for element distortional rotation greater than 30 degrees
         con = udef(j+3)**2 + udef(j+4)**2 +  udef(j+5)**2
         if (con.gt.0.274d0) kerr = 1
         j = j + 6
      enddo
c
      return
      end
