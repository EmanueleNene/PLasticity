      subroutine BeamSStif2D (force,aleng,nr,gstif)
c --- primary function: compute stress stiffness matrix for an euler beam
c
c     ref: przemieniecki, page 390
c     note - set gstif matrix to zero before call to this routine
c
c  input arguments:
c     force  (dp,sc,in)          - axial force
c     aleng  (dp,sc,in)          - length of beam
c     nr     (int,sc,in)         - dimension of the gstif matrix
c
c  output arguments:    
c     gstif  (dp,ar(nr,nr),out)  - stress stiffness matrix
c
#include "impcom.inc"
c 
      integer  nr
      double precision
     & force,aleng,gstif(nr,nr),aleng2,con,con1,DP6o5,THREE,
     & DP7p5,TEN
      parameter (DP6o5=1.2d0,THREE=3.0d0,DP7p5=7.5d0,TEN=10.0d0)
c
      aleng2 = aleng*aleng
      con = force/aleng
      con1 = DP6o5*con
      gstif(2,2) = con1
      gstif(5,5) = con1
      gstif(5,2) = -gstif(2,2)
      gstif(2,5) = gstif(5,2)
      con1 = con*aleng2/DP7p5
      gstif(3,3) = con1
      gstif(6,6) = con1
      con1 = con*aleng/TEN
      gstif(3,2) = con1
      gstif(2,3) = gstif(3,2)
      gstif(6,2) = con1
      gstif(2,6) = gstif(6,2)
      gstif(5,3) = - gstif(6,2)
      gstif(3,5) = gstif(5,3)
      gstif(6,5) = - gstif(6,2)
      gstif(5,6) = gstif(6,5)
      con1 = con*aleng2/(TEN*THREE)
      gstif(6,3) = - con1
      gstif(3,6) = gstif(6,3)
c
      return
      end
