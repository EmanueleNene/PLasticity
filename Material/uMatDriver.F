      subroutine uMatDriver (keycut, kMatPoint, 
     &                       moist, dmoist,
     &                       nTens, elemId, matId,
     &                       nlkword, matMW, matMWIndex,
     &                       stress, dsdePl,
     &                       sedEl, sedPl, sdv)
c --- primary function: interface to ANSYS usermat format
c
#include "impcom.inc"
#include "constants.inc"

      integer          userpropMax
      parameter        (userpropMax=1000)
      integer          keycut, nTens, kMatPoint,
     &                 matMWIndex(*),
     &                 elemId, matId, nStatev
      integer*8        nlkword(W_nMatCol,W_nMatRow)
      double precision matMW(*), dsdePl(nTens,nTens), 
     &                 stress(nTens), sedEl, sedPl,
     &                 usrprop(userpropMax), sdv(*)
      double precision moist, dmoist

      integer          i1 , i2 , i3 ,i4 , i6, i5,
     &                 nusrprop, option
      integer*8        keyplL
      double precision strain(6), dstrain(6), eppl(6),
     &                 strain0(6), epseq_t, epseq

      keycut   = 0
      i1 = matMWIndex (W_iStrain)
      i2 = matMWIndex (W_iStrainPl)
      i3 = matMWIndex (W_iEqvPl)
      i4 = matMWIndex (W_iStrainInc)
      i5 = matMWIndex (W_iMisc)
      i6 = matMWIndex (W_iStateVar)

c --- retrive state variables
      nstatev = nlkword(2,4)
      call RetMatRec (W_iStateVar,kMatPoint,nstatev,matMW(i6))

      keyplL = nlkword(1,1)
c --- retrieve material constants
      call v_zero (usrprop(1),10)
      call GetTabParam(keyplL,moist,nusrprop,usrprop(1),option)
      
      sedEl = matMW(i5)
      sedPl = matMW(i3+1)
      call v_move (matMW(i1), strain0(1), nTens)
      call v_move (matMW(i4), dstrain(1), nTens)
      call v_move (matMW(i2), eppl(1), nTens)

      call v_move (strain0(1), strain(1), nTens)
      call v1mv2_v1 (strain0(1), dstrain(1), nTens)

      epseq = abs(matMW(i3))
      epseq_t = epseq

      call ChabocheSV (matId, elemId, keycut,
     &                 nTens, nStatev, nusrprop,
     &                 stress(1), matMW(i6), dsdePl(1,1), sedEl, sedPl,
     &                 strain0(1), dstrain(1), eppl(1), usrprop(1))
     
c --- update maximum plastic strain and incremental plastic strain
      sdv(1) = max(sdv(1),epseq-epseq_t)
      sdv(2) = max(sdv(2),epseq)

      matMW(i3) = -epseq
      matMW(i3+1) = sedPl
      matMW(i5) = sedEl
      call v_move (eppl(1), matMW(i2),nTens)
      call v1mv2_v3 (strain(1), matMW(i2), matMW(i1), nTens)

c --- update state variables
      call PutMatRec (W_iStateVar, kMatPoint, nstatev, matMW(i6))

      return
      end
