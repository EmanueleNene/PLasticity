      subroutine BeamSStif (force,aleng,nr,ioa,gstif)
c --- primary function: compute stress stiffness matrix for an euler beam
c
c     ref: przemieniecki, page 390
c     note - set gstif matrix to zero before call to this routine
c
c  input arguments:
c     force  (dp,sc,in)          - axial force
c     aleng  (dp,sc,in)          - length of beam
c     nr     (int,sc,in)         - dimension of the gstif matrix
c     ioa    (dp,sc,in)          - polar moment of inertia divided by the area
c                                  (polar moment of inertia = iyy + izz)
c  output arguments:    
c     gstif  (dp,ar(nr,nr),out)  - stress stiffness matrix
c
#include "impcom.inc"
c
      integer  nr
      double precision force,aleng,ioa,gstif(nr,nr),
     &                 aleng2,con,con1,DP6o5,THREE,DP7p5,TEN
      parameter (DP6o5=1.2d0,THREE=3.0d0,DP7p5=7.5d0,TEN=10.0d0)
c
      aleng2 = aleng*aleng
      con = force/aleng
      con1 = DP6o5*con
      gstif(2,2) = con1
      gstif(3,3) = con1
      gstif(8,8) = con1
      gstif(9,9) = con1
      gstif(8,2) = -gstif(9,9)
      gstif(2,8) = gstif(8,2)
      gstif(9,3) = -gstif(9,9)
      gstif(3,9) = gstif(9,3)
      con1 = con*aleng2/DP7p5
      gstif(5,5) = con1
      gstif(6,6) = con1
      gstif(11,11) = con1
      gstif(12,12) = con1
      con1 = con*aleng/TEN
      gstif(6,2) = con1
      gstif(5,3) = -gstif(6,2)
      gstif(3,5) = gstif(5,3)
      gstif(2,6) = gstif(6,2)
      gstif(12,2) = con1
      gstif(2,12) = gstif(12,2)
      gstif(11,3) = - gstif(12,2)
      gstif(3,11) = gstif(11,3)
      gstif(9,5) = con1
      gstif(5,9) = gstif(9,5)
      gstif(8,6) = -gstif(9,5)
      gstif(6,8) = gstif(8,6)
      gstif(12,8) = -gstif(9,5)
      gstif(8,12) = gstif(12,8)
      gstif(11,9) = con1
      gstif(9,11) = gstif(11,9)
      con1 = con*aleng2/(THREE*TEN)
      gstif(11,5) = -con1
      gstif(5,11) = gstif(11,5)
      gstif(12,6)= -con1
      gstif(6,12) = gstif(12,6)
c --- add torsional effects  - timoshenko + gere, section 5.5 (causes problems and therefore nulled)
      con1 = con*ioa
      gstif(4,4) = con1
      gstif(10,10) = con1
      gstif(4,10) = -con1
      gstif(10,4) = -con1
      return
      end
