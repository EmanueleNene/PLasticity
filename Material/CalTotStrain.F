      subroutine CalTotStrain (nTens, nlkword, khygro,
     &                         matMW, matMWIndex,
     &                         matRec, matRecIndex,
     &                         strain)
c --- primary function: compute the total strain based on various options
c
c  input arguments:
c      nTens   (int,sc,in)           - number of tensor components
c      matMW   (dp,ar(*),in)         - work array for mat record
c      matMWIndex(int,ar(*),in)      - array indexing into matMW
c      matRec  (dp,ar(*),in)         - mat record
c      matRecIndex(int,ar(*),in)     - mat record index
c      nlkword (int*8,ar(2,in),in)   - nonlinear material index
c
c  output arguments:
c      strain  (dp,ar(*),out)        - total strain
c
#include "impcom.inc"
#include "constants.inc"
c
      integer          nTens, khygro,
     &                 matRecIndex(*), matMWIndex(*)
      integer*8        nlkword(W_nMatCol,W_nMatRow)

      double precision matMW(*), strain(*), matRec(*)


      integer*8        keyplL, keycrL


      keyplL  = nlkword(1,1)
      keycrL  = nlkword(1,2)

      call v_move (matMW(matMWIndex(W_iStrain)),
     &             strain(1), nTens)

      if (keyplL.gt.0) then
c --- update total strain  E = E + E_pl
          call v1pv2_v1 (strain(1),
     &                   matMW(matMWIndex(W_iStrainPl)), nTens)
      endif

      if (keycrL.gt.0) then
c --- update total strain  E = E + E_cr
          call v1pv2_v1 (strain(1),
     &                   matMW(matMWIndex(W_iStrainCr)), nTens)
      endif


      if (khygro.gt.0) then
c --- update total strain  E = E + E_th
          call v1pv2_v1 (strain(1),
     &                   matMW(matMWIndex(W_iStrainHyg)), nTens)

          call v1pv2_v1 (strain(1),
     &                   matRec(matRecIndex(W_iStrainHyg)), nTens)
      endif
    
      return
      end
