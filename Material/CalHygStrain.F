      subroutine CalHygStrain(elemId, matId, nDirect, nTens,
     &                        moist, mref, matMW,  matWrkIndex, khygro)
c --- primary function: calculate linear hygro strains
c
c  input arguments:
c      elemId      (int,sc,in)         - element number
c      matId       (int,sc,in)         - material number
c      nDirect     (int,sc,in)         - number of direct components
c      nTens       (int,sc,in)         - number of tensor components
c      moist       (dp,sc,in)          - current moisture
c      mref        (dp,sc,in)          - reference moisture
c      matMW       (dp,ar(*),in)       - hygrostrain (material work array)
c      matWrkIndex (int,ar(*),in)      - material work array index
c
c  output arguments:
c      khygro      (int,sc,out)        - 1: hygro strains, 0: otherwise
c
#include "impcom.inc"
#include "constants.inc"
c
      double precision ZERO
      parameter       (ZERO  = 0.0d0)

      integer          v_last

      integer          elemId, matId, khygro, nDirect, nTens,
     &                 matWrkIndex(*)
      double precision moist, mref, matMW(*)
c

c
      integer          i, kstart
      double precision dmoistRef, beta(3)

c --- hygro expansion index on MP
      integer lp(3)
      data    lp /10,11,12/

c --- evaluate hygro coefficients

      call GetMatProps (elemId, matId, lp(1), moist, beta(1), nDirect)
c
      dmoistRef = moist - mref
      kstart = matWrkIndex(W_iStrainHyg)
      khygro = 0
      call v_zero(matMW(kstart),nTens)

      if (v_last(beta(1),nDirect).gt.0) then
c --- when hygro expansion is defined
          khygro = 1
          call cv1_v2 (beta(1),matMW(kstart),nDirect,dmoistRef)
      endif
c
      return
      end
