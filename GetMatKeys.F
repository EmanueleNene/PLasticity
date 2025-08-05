      subroutine GetMatKeys (matId,nlkword)
c --- primary function:  obtain nonlinear material property information
c
c input arguments:
c    matId    (int,sc,in)    - element material number
c
c output arguments:
c    nlkword (int,ar(2,*),out)
c             nlkword(1,1)   - 1: rate-independent plasticity exists
c             nlkword(2,1)   - plasticity option
c
c             nlkword(1,2)   - 1: creep exists
c             nlkword(2,2)   - creep option
c
#include "impcom.inc"
#include "constants.inc"
c
      integer  nlgetPtrs

      integer*8 nlkword(W_nMatCol,W_nMatRow), MakeLongInt
      integer   matId, i, j
      double precision propnl(*), Header(*)
      pointer  (dPtrHeader,Header)
      pointer  (dPtrTable,propnl)

#if defined(DEBUG)
      call SubBeg('GetMatKeys')
#endif

c --- initialize output data
      do j = 1,W_nMatRow
         do i = 1,W_nMatCol
            nlkword(i,j) = 0
         enddo
      enddo

c --- get the pointer to the start of the table of property data
      i = nlgetPtrs(matId,NLMP_USERMAT,dPtrHeader,dPtrTable)
      if (i.gt.0) then
          nlkword(1,1) = MakeLongInt(matId,NLMP_USERMAT)
          nlkword(2,1) = nint(Header(1))
      endif

      i = nlgetPtrs(matId,NLMP_CREEP,dPtrHeader,dPtrTable)
      if (i.gt.0) then
          nlkword(1,2) = MakeLongInt(matId,NLMP_CREEP)
          nlkword(2,2) = nint(Header(13))
      endif

      i = nlgetPtrs(matId,NLMP_CYCPLAST,dPtrHeader,dPtrTable)
      if (i.gt.0) then
          nlkword(1,1) = MakeLongInt(matId,NLMP_CYCPLAST)
          nlkword(2,1) = nint(Header(1))
          nlkword(1,3) = nlkword(1,1)
      endif

      i = nlgetPtrs(matId,NLMP_PLAST,dPtrHeader,dPtrTable)
      if (i.gt.0) then
          nlkword(1,1) = MakeLongInt(matId,NLMP_PLAST)
          nlkword(2,1) = nint(Header(1))
      endif

      i = nlgetPtrs(matId,NLMP_USTATE,dPtrHeader,dPtrTable)
      if (i.gt.0) then
          nlkword(1,4) = MakeLongInt(matId,NLMP_USTATE)
          nlkword(2,4) = nint(Header(13))              
      endif       

#if defined(DEBUG)
      call SubEnd('GetMatKeys')
#endif

      return
      end


