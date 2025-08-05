      subroutine PlCrDriver  (keycut, moist,
     &                        nTens, elemId, matId, prop,
     &                        nlkword, matMW, matMWIndex,
     &                        matWork, matWorkIndex,
     &                        stress, dsdePl,
     &                        sedEl, sedPl, sedCr, sdv)
c --- primary function: visco-plasticity model driver
c
c  input arguments:
c       keycut      (int,sc,loc)       -  bisec/cut flag
c       matMW       (dp,ar(*),in)      -  mat array (modified copy of matrec)
c       matMWIndex  (int,ar(*),in)     -  array indexing into matMW array
c       matWork     (dp,ar(*),in)      -  mat work array
c       matWorkIndex(int,ar(*),in)     -  array indexing into matWork array
c       moist       (dp,sc,in)         -  moisture at the materila point
c       nTens       (int,sc,in)        -  nDirect + nShear
c       elemId      (int,sc,in)        -  element #
c       matId       (int,sc,in)        -  material #
c
c  output arguments:
c       stress      (dp,ar(nTens),out) -  Cauchy stress at t+Dt
c       strain      (dp,ar(nTens),out) -  total strain at t+Dt
c       siglv       (dp,ar(4),out)     -  hygro load vector
c       siglvi      (dp,ar(4),out)     -  initial stress load vector
c       eEnerPt     (dp,ar(7),out)     -  energy at this point
c
#include "impcom.inc"
#include "constants.inc"
c
      double precision ZERO
      parameter       (ZERO = 0.d0)

      integer          keycut, nTens,
     &                 matMWIndex(*), matWorkIndex(*),
     &                 elemId, matId
      integer*8        nlkword(W_nMatCol,W_nMatRow)
      double precision matMW(*), matWork(*),
     &                 dsdePl(nTens,*), stress(*),
     &                 sedEl, sedPl, sedCr, sdv(*)
      double precision moist, prop(13)

      integer          i1 , i2 , i3 , i4 , i5 , i6 , i7 ,
     &                 i15, w1 , w2 , w3 , w4 , w5 , w21,
     &                 option
      integer*8        keycrL

      integer          npropcr
      double precision propcr(1000)


      keycrL   = nlkword(1,2)
      sedEl = ZERO
      sedPl = ZERO
      sedCr = ZERO
c --- define indices
      i1  = matMWIndex(W_iStrain)
      i2  = max(matMWIndex(W_iStrainPl),1)
      i3  = max(matMWIndex(W_iStrainCr),1)
      i4  = max(matMWIndex(W_iEqvPl),1)
      i5  = max(matMWIndex(W_iEqvCr),1)
      i6  = matMWIndex(W_iMatMatrix)
      i7  = matMWIndex(W_iStress)
      i15 = matMWIndex(W_iMisc)

c --- retrieve creep constants
      if (keycrL .gt. 0)
     &    call GetTabParam(keycrL,moist,npropcr,propcr(1),option)
      
      w1  = matWorkIndex(W_WrkIndx1)
      w2  = matWorkIndex(W_WrkIndx2)
      w3  = matWorkIndex(W_WrkIndx3)
      w4  = matWorkIndex(W_WrkIndx4)
      w5  = matWorkIndex(W_WrkIndx5)
      w21 = matWorkIndex(W_WrkIndx6)
      
      call CreepPlastDriver (keycut, moist, nTens, 
     &                       elemId, matId, prop(1),
     &                       propcr(1), matMW(i1), matMW(i2), matMW(i3),
     &                       matMW(i4), matMW(i5), matMW(i15),
     &                       matWork(w1), matWork(w2), matWork(w3),
     &                       matWork(w4), matWork(w5), matWork(w21),
     &                       matMW(i7), stress(1), dsdePl(1,1),
     &                       nlkword(1,1), sedEl, sedPl, sedCr, 
     &                       sdv(1))


      return
      end
