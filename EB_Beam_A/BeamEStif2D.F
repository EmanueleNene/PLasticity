      subroutine BeamEStif2D (ex,gxy,area,ai,shearf,aleng,nr,zs)
c --- primary function: compute elastic stiffness matrix for an euler beam
c
c     ref:  eqn. 5.116, przemieniecki, page 79
c     note: set zs matrix to zero before call to this routine
c
c  input arguments:
c     ex      (dp,sc,in)         - modulus of elasticity
c     gxy     (dp,sc,in)         - shear modulus
c     area    (dp,sc,in)         - cross-sectional area
c     ai      (dp,sc,in)         - area moment of inertia
c     shearf  (dp,sc,in)         - shear constant
c     aleng   (dp,sc,in)         - length of beam
c     nr      (int,sc,in)        - dimension of the zs matrix
c
c  output arguments:
c     zs      (dp,ar(nr,nr),out) - stiffness matrix
c
#include "impcom.inc"
c
      integer  nr
      double precision
     & ex,gxy,area,ai,shearf,aleng,zs(nr,nr),
     & eis,eil,eil2,eil3,con1
c
      con1 = ex*ai*shearf*12.0d0/(gxy*area*aleng**2)
      eis = ex*ai/(1.0d0 + con1)
      eil = eis/aleng
      eil2 = eil/aleng
      eil3 = eil2/aleng
      zs(1,1) = ex*area/aleng
      zs(4,4) = zs(1,1)
      zs(1,4) = -zs(1,1)
      zs(4,1) = zs(1,4)
      zs(2,2) = 12.0d0*eil3
      zs(5,2) = -zs(2,2)
      zs(3,3) = (4.0d0 + con1)*eil
      zs(3,6) = (2.0d0 - con1)*eil
      zs(2,3) = 6.0d0*eil2
      zs(2,6) = zs(2,3)
      zs(3,5) = -zs(2,3)
      zs(5,5) = zs(2,2)
      zs(2,5) = zs(5,2)
      zs(6,6) = zs(3,3)
      zs(6,3) = zs(3,6)
      zs(3,2) = zs(2,3)
      zs(5,6) = -zs(3,2)
      zs(6,5) = zs(5,6)
      zs(6,2) = zs(2,6)
      zs(5,3) = zs(3,5)
      return
      end
