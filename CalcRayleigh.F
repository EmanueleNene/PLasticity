      subroutine CalcRayleigh (kelfil,nr,Kmtrx,uelm,kEnergy)
c --- primary function:   accumulate terms for Rayleigh frequency estimate
c     secondary function: calculate nonlinear transient kinetic energy
c
c  input arguments:
c     kelfil        (int,ar(10),in)     - element matrix availablity keys
c     nr            (int,sc,in)         - number of rows in element matrices
c     Kmtrx         (dp,ar(nr,nr,3),in) - the element matrices
c     uelm          (dp,ar(nr),in)      - element nodal displacements
c
c  output arguments:
c     kEnergy       (dp,sc,out)         - element kinetic energy
c
#include "impcom.inc"
#include "constants.inc"
#include "smp_info.inc"
#include "solu_info.inc"
#include "dyn_info.inc"
c
      double precision v_dot

      integer nr,kelfil(10)
      double precision   Kmtrx(nr,nr,*),uelm(nr,*),kEnergy,con1,con2

      double precision work(W_NrMax)

      con1 = 0.0d0
      con2 = 0.0d0

c --- numerator (stiffness)
      if (kelfil(1).eq.1) then
          call Axv1_v2(Kmtrx(1,1,1),uelm(1,A_RowIncS),work(1),nr,nr)
          con1 = v_dot(uelm(1,A_RowIncS),work(1),nr)
      endif

c --- denominator (mass)
      if (kelfil(2).eq.1) then
          call Axv1_v2(Kmtrx(1,1,2),uelm(1,A_RowIncS),work(1),nr,nr)
          con2 = v_dot(uelm(1,A_RowIncS),work(1),nr)
c --- calculate kinetic energy as 1/2*v*M*v
          call Axv1_v2(Kmtrx(1,1,2),uelm(1,A_RowVelo),work(1),nr,nr)
          kEnergy = 0.5d0*v_dot(uelm(1,A_RowVelo),work(1),nr)
      endif

      call lockSMP (W_LockTrans)
      I_RespStiff = I_RespStiff + con1
      I_RespMass = I_RespMass + con2
      call unlockSMP (W_LockTrans)

      return
      end

