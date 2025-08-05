      subroutine StresIni(elemId, matId, nTens,kDomIntPt, iSectPt,
     &                    kMatPoint, coords, basis, siglvi)
c --- primary function: provides initial stress state
c
c  input arguments:
c       nTens       (int,sc,in)       - nDirect + nShear
c       elemId      (int,sc,in)       - element Id
c       matId       (int,sc,in)       - material Id
c       kMatPoint   (int,sc,in)       - "k"th domain integration point
c       iSectPt     (int,sc,in)       - "k"th section point
c       coords      (dp,ar(3),in)     - coordinates of this pt
c       basis       (dp,ar(3,3),in)   - basis directions at this pt
c
c  output arguments:
c       tMatProp    (dp,ar(*),out)    - temporary array for material properties
c       siglvi      (dp,ar(4),out)    - initial stress load vector
c       eEnerPt     (dp,ar(7),out)    - energy at this point
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
c
      integer          elemId, nTens, kDomIntPt, iSectPt, kMatPoint,
     &                 matId
      double precision siglvi(6), coords(3), basis(3,3)

c --- get initial stress
      if (A_LoadStepNum.eq.1 .and. A_SubStepNum.eq.1) then
c --- get initial stress
          call GetIniStress (elemId, matId, nTens, basis(1,1),
     &                       kDomIntPt, iSectPt, coords(1), siglvi(1))
          if (A_SubStepNum.eq.1) then
              call PutMatRec(W_iStrainIni, kMatPoint, nTens, siglvi(1))
          endif
      else
          call RetMatRec (W_iStrainIni, kMatPoint, nTens, siglvi(1))
      endif

      return
      end
