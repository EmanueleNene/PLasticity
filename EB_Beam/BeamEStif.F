      subroutine BeamEStif (ex,gxy,area,aiz,aiy,aip,shearz,sheary,aleng,
     &                  nr,zs)
c --- primary function: compute elastic stiffness matrix for an euler beam
c
c     ref:  eqn. 5.116, przemieniecki, page 79
c     note: set zs matrix to zero before call to this routine
c
c  input arguments:
c     ex      (dp,sc,in)         - modulus of elasticity
c     gxy     (dp,sc,in)         - shear modulus
c     area    (dp,sc,in)         - cross-sectional area
c     aiz     (dp,sc,in)         - z-axis area moment of inertia
c     aiy     (dp,sc,in)         - y-axis area moment of inertia
c     aip     (dp,sc,in)         - polar moment of inertia
c     shearz  (dp,sc,in)         - z-axis shear constant
c     sheary  (dp,sc,in)         - y-axis shear constant
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
     & ex,gxy,area,aiz,aiy,aip,shearz,sheary,aleng,zs(nr,nr),
     & aleng2,eiz,eiy,con,con1,con2
c
      aleng2 = aleng*aleng
      eiz = ex*aiz
      eiy = ex*aiy
      con1 = 12.d0/(gxy*area*aleng2)
      con2 = con1*sheary*eiy
      con1 = con1*shearz*eiz
      zs(1,1) = ex*area/aleng
      zs(7,1) = -zs(1,1)
      zs(1,7) = zs(7,1)
      zs(7,7) = zs(1,1)
      con = 12.0d0/(aleng2*aleng)
      zs(2,2) = con*eiz/(1.0d0 + con1)
      zs(8,2) = -zs(2,2)
      zs(2,8) = zs(8,2)
      zs(8,8) = zs(2,2)
      zs(3,3) = con*eiy/(1.0d0 + con2)
      zs(9,3) = -zs(3,3)
      zs(3,9) = zs(9,3)
      zs(9,9) = zs(3,3)
      zs(4,4) = gxy*aip/aleng
      zs(10,4) = -zs(4,4)
      zs(4,10) = zs(10,4)
      zs(10,10) = zs(4,4)
      con = 6.0d0/aleng2
      zs(6,2) = con*eiz/(1.0d0 + con1)
      zs(2,6) = zs(6,2)
      zs(12,2) = zs(6,2)
      zs(2,12) = zs(12,2)
      zs(12,8) = -zs(6,2)
      zs(8,12) = zs(12,8)
      zs(8,6) = -zs(12,2)
      zs(6,8) = zs(8,6)
      zs(5,3) = -con*eiy/(1.0d0 + con2)
      zs(3,5) = zs(5,3)
      zs(11,3) = zs(5,3)
      zs(3,11) = zs(11,3)
      zs(11,9) = -zs(5,3)
      zs(9,11) = zs(11,9)
      zs(9,5) = -zs(11,3)
      zs(5,9) = zs(9,5)
      zs(5,5) = (4.0d0+con2)*eiy/(aleng*(1.0d0+con2))
      zs(11,11) = zs(5,5)
      zs(6,6) = (4.0d0+con1)*eiz/(aleng*(1.0d0+con1))
      zs(12,12) = zs(6,6)
      zs(11,5) = (2.0d0-con2)*eiy/(aleng*(1.0d0+con2))
      zs(5,11) = zs(11,5)
      zs(12,6) = (2.0d0-con1)*eiz/(aleng*(1.0d0+con1))
      zs(6,12) = zs(12,6)
c
      return
      end
