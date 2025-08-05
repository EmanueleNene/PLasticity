      subroutine GetIniStress (elemId, matId, nTens, basis,
     &                         curIntPt, curSecPt, curCoords, iniData)
c --- primary function: get the user supplied stress state for
c                       initial stress analysis
c
c  input arguments:
c      elemId         (int,sc,in)     - element number
c      matId          (int,sc,in)     - material number
c      basis          (dp,ar(3,3),in) - basis directions (for the element)
c      curIntPt       (int,sc,in)     - current Integration point
c      curSecPt       (int,sc,in)     - current point through the thickness in the
c                                       current layer
c      curCoords      (dp,ar(3),in)   - current coordinates
c
c  output arguments:
c      iniData        (dp,ar(nDirect+nShear),out)
c
#include "impcom.inc"
c
      integer          elemId, elemType, matId, curIntPt,
     &                 curSecPt, nTens, lvl, kerr
      double precision curCoords(3), basis(3,3), ZERO,
     &                 iniData(nTens), dpValue, subc(1)
      character*(32)   string, pname
      data string /'    '/
      data pname  /'_INISTRESS'/
      parameter (ZERO = 0.0d0)

      dpValue = ZERO
      call v_zero (iniData(1),nTens)

      subc(1) = dble(matId)
      call parevl (pname,1,subc(1),lvl,dpValue,string,kerr)
      
      if (kerr.eq.0) then
          iniData(1) = dpValue
      endif

      return
      end
