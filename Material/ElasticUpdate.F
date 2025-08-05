      subroutine ElasticUpdate (nTens, matMW, matMWIndex, stress, sedEl)
c --- primary function: performs elastic constitutive calc
c
c  input arguments:
c      nTens   (int,sc,in)           -  # of tensor components
c      matMWIndex(int,ar(*),in)      -  array indexing into matMW
c
c  output arguments:
c      matMW   (dp,ar(*),inout)      -  mat record in work array
c      stress  (dp,ar(*),out)        -  stress based on el strain
c      siglvi  (dp,ar(4),out)        -  initial stress load vector
c      sedEl   (dp,sc,out)           -  elastic strain energy
c
c
#include "impcom.inc"
#include "constants.inc"
c
      double precision ZERO, HALF
      parameter       (ZERO = 0.d0,
     &                 HALF = 0.5d0)
c

c
      integer          nTens, matMWIndex(*)
      double precision matMW(*), stress(nTens),sedEl(7)
c

c
      integer          i1, i
      double precision prs, evol


      call v_move (matMW(matMWIndex(W_iStress)), stress, nTens)

c --- elastic strain energy density
c
      i1 = matMWIndex(W_iStrain)
      sedEl(1) = ZERO
      do i=0,nTens-1
         sedEl(1) = sedEl(1) + stress(i+1)*matMW(i1+i)
      enddo
      sedEl(1) = sedEl(1) * HALF
      evol= matMW(i1)
      prs = stress(1)/3.d0
      sedEl(6) = prs*evol*0.5d0
      sedEl(7) = sedEl(1)-sedEl(6)

      return
      end

