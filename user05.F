      function user05()
c --- primary function: define real constant
c
#include "impcom.inc"
#include "constants.inc"
#include "file_sys.inc"
c
      integer user05,intinfun,elmiqr,elmget,rlgetn,etyiqr
      double precision dpinfun
      integer i,j,elData(EL_DIM),nnod,elmax,rch,ireal,etyp
      double precision nodes(W_NdMax),rvr(W_NumContParam),rnew

      user05 = 0
c --- real constant number
      rch = intinfun(2)
      if (rch.le.0 .or. rch.gt.W_NumContParam) goto 999
c --- real constant value
      rnew = dpinfun(3)
      elmax = elmiqr(0,DB_MAXDEFINED)
      do i = 1,elmax
         nnod = elmget (i,elData(1),nodes(1))
c --- only master elements are affected by this command!
         etyp = etyiqr(elData(EL_TYPE),-JETYP)
         if (etyp.ne.W_eMaster) cycle
         ireal = elData(EL_REAL)
         if (ireal.le.0) cycle
         j = rlgetn(ireal,1,W_NumContParam,rvr(1))
         if (j.le.0) cycle ! not initialized apperently
         if (j.lt.(rch - 1)) call v_zero(rvr(j+1),rch - j - 1)
         rvr(rch) = rnew
         call rlput (ireal,max(j,rch),rvr(1)) 
      enddo
  999 return
      end