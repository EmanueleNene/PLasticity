      subroutine CalUdef2D (xyz,nx,averot,utot,udef)
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
c     averot (dp,sc,in)              - average rotation of the element
c     utot   (dp,ar(6*nnod),in)      - total displacements (with the rotations being Euler values)
c
c  output arguments:
c     udef   (dp,ar(6*nnod),out)     - deformational displacements (unrotated)
c                                      udef acts on current (rotated) orientation of the element
c
      integer nx
      double precision
     & xyz(nx,2),averot,utot(6),udef(6),xrel(2),urel(2),cth,sth

      cth = cos(averot)
      sth = sin(averot)

c --- translational DOF (eqn 2)
      call v_zero (udef(1),2)
      call v1mv2_v3 (xyz(1,2),xyz(1,1),xrel(1),2)
      call v1mv2_v3 (utot(4),utot(1),urel(1),2)
      udef(4) = xrel(1) + urel(1) -  cth*xrel(1) + sth*xrel(2)
      udef(5) = xrel(2) + urel(2) -  sth*xrel(1) - cth*xrel(2)

c --- rotational DOF
      udef(3) = utot(3) - averot
      udef(6) = utot(6) - averot
c
      return
      end
