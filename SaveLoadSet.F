      function SaveLoadSet (gmask,timfrq,itter,itime,ncumit,
     &                      u,f,locrfL,lochedL)
c --- primary function: save a load set on the result file
c
c  input arguments:
c     gmask    (int,sc,in)       - global mask
c     timfrq   (dp,sc,in)        - the time for this load set
c     itter    (int,sc,in)       - the iteration for this load set
c     itime    (int,sc,in)       - the load set for this set
c     ncumit   (int,sc,in)       - the cumulative load set number
c     u        (dp,ar(*),in)     - nodal displacement vector
c     f        (dp,ar(*),in)     - nodal force vector
c
c  output arguments:
c     locrfL       (intL,sc,out) - location of reaction forces
c     SaveLoadSet  (int,sc,out)  - 0: successful
c     lochedL      (intL,sc,out) - location of header
c
#include "impcom.inc"
#include "constants.inc"
#include "RST_file.inc"
c
      integer*8   Bin_Get
      PTRFTN    MemAlloc

      integer*8  loc8, curloc8, leng8, lochedL, locrfL, iL
      integer  SaveLoadSet,gmask,itter,itime,ncumit,
     &         kbf,i,j,ivect(100),nelm,mask
      double precision  u(*),f(*),timfrq
      logical nooutput

      pointer (piOrderLocalL,OrderLocal)
      pointer (piWorkL,iWork)
      pointer (pdWorkL,dWork)
      integer  OrderLocal(*), iWork(*)
      double precision  dWork(*)

#if defined(DEBUG)
      call SubBeg('SaveLoadSet')
#endif

c --- see if file is being written
      SaveLoadSet = 0
      if (R_BlockNum .le. 0) goto 999

c --- set function value for unsuccessful operation
      SaveLoadSet = 1
      lochedL = -1

c --- check if anything to write (first block out print control bits)
      mask = gmask
      i = iand(mask,W_AllResMask)

      nooutput = i.eq.0
      if (nooutput) goto 999
      SaveLoadSet = -1
      if (R_Unit .le. 0) goto 999

      R_CurSet = R_CurSet + 1

c --- update index vectors
      curloc8 = Bin_Get(R_BlockNum) + R_FilePos
      loc8 = 0
      call Bin_Put (R_BlockNum,loc8)
      piWorkL = MemAlloc(2*R_MaxSet,W_LenInt,'piWorkL         ')
      iL = R_ResIndxLoc
      j = R_MaxSet
      call Buf_Read (R_BlockNum,iL,j,iWork(1),kbf)

      call MakeShortInt (curloc8,iWork(R_CurSet),
     &                   iWork(R_CurSet+R_MaxSet))
      if (R_CurSet .lt. R_MaxSet) then
          iWork(R_CurSet+1) = 0
          iWork(R_MaxSet+R_CurSet+1) = 0
      endif
      iL = R_ResIndxLoc

      call Buf_Write (R_BlockNum,iL,2*R_MaxSet,iWork(1),2*R_MaxSet)

      call MemFree (piWorkL)

      pdWorkL = MemAlloc(R_MaxSet,W_LenDouble,'pdWorkL         ')
      iL = R_TimIndxLoc
      j = R_MaxSet*W_IntPerDp
      call Buf_Read (R_BlockNum,iL,j,dWork(1),kbf)
      dWork(R_CurSet) = timfrq
      iL = R_TimIndxLoc
      call Buf_Write (R_BlockNum,iL,R_MaxSet*W_IntPerDp,dWork(1),0)
      call MemFree (pdWorkL)

      piWorkL = MemAlloc(3*R_MaxSet,W_LenInt,'piWorkL         ')
      iL = R_LdStepLoc
      j = 3*R_MaxSet
      call Buf_Read (R_BlockNum,iL,j,iWork(1),kbf)
      iWork(3*R_CurSet-2) = itime
      iWork(3*R_CurSet-1) = itter
      iWork(3*R_CurSet)   = ncumit
      iL = R_LdStepLoc
      call Buf_Write (R_BlockNum,iL,R_MaxSet*3,iWork(1),R_MaxSet*3)
      call MemFree (piWorkL)

c --- get element order vector
      iL = W_HedLen
      j = 40
      call Buf_Read (R_BlockNum,iL,j,ivect(1),kbf)
      nelm = ivect(7)
      piOrderLocalL = MemAlloc(nelm,W_LenInt,'piOrderLocalL   ')
      iL = ivect(14)
      j = nelm
      call Buf_Read (R_BlockNum,iL,j,OrderLocal(1),kbf)

c --- write load step records
      call Bin_Put (R_BlockNum,curloc8)
      R_FilePos= 0
      lochedL= R_FilePos
      SaveLoadSet = 0

      call SaveDisp(R_BlockNum,R_FilePos,nelm,mask,itter,
     &              itime,ncumit,u(1),f(1),locrfL)

      call MemFree (piOrderLocalL)

c --- update header record
      iL = W_HedLen
      curloc8 = Bin_Get(R_BlockNum)
      loc8 = 0
      call Bin_Put(R_BlockNum,loc8)
      i = 40
      call Buf_Read (R_BlockNum,iL,i,ivect(1),kbf)
      ivect(9) = R_CurSet
      leng8 = curloc8 + R_FilePos
      ivect(10) = leng8
      call MakeShortInt (leng8,ivect(23),ivect(24))
      iL = W_HedLen
      i = 40
      call Buf_Write (R_BlockNum,iL,i,ivect(1),i)
      call Bin_Put (R_BlockNum,curloc8)

  999 continue

#if defined(DEBUG)
      call SubEnd('SaveLoadSet')
#endif

      return
      end
