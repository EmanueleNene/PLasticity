      subroutine IndxEmat (nRow,Emat,LV,Idx)
c --- primary function: sort element matrix into index order
c
c  input arguments:
c     nRow     (int,sc,in)             - size of element matrix
c
c  output arguments:
c     Emat     (dp,ar(nRow,nRow,inout) - sorted element matrix
c     LV       (dp,ar(nRow),inout)     - sorted element load vector
c     Idx      (int,ar(nRow),inout)    - sorted element index vector
c
#include "impcom.inc"
c
      integer   nRow, Idx(nRow), iRow, jRow, imin, vmin
      double precision   Emat(nRow,nRow), LV(nRow), con

#if defined(DEBUG)
      call SubBeg('IndxEmat')
#endif

      do iRow = 1,nRow
c --- find minimum
         if (Idx(iRow).eq.0) cycle
         imin = iRow
         vmin = Idx(iRow)
         do jRow = iRow+1,nRow
            if (Idx(jRow).gt.0 .and. Idx(jRow).lt.vmin) then
                imin = jRow
                vmin = Idx(jRow)
            endif
         enddo

c --- reorder
         if (imin.ne.iRow) then

             call m_swapCol (Emat(1,iRow),Emat(1,imin),nRow)
             call m_swapRow (Emat(iRow,1),Emat(imin,1),nRow)

             jRow = Idx(iRow)
             Idx(iRow) = Idx(imin)
             Idx(imin) = jRow

             con = LV(iRow)
             LV(iRow) = LV(imin)
             LV(imin) = con

         endif

      enddo
      
#if defined(DEBUG)
      call SubEnd('IndxEmat')
#endif

      return
      end
