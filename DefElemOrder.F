      subroutine DefElemOrder (n,eorder,NumEl)
c --- primary function: define the element order vector
c
c  input arguments:
c     n        (int,sc,in)       - highest element number defined
c                                  (size of eorder vector)
c
c  output arguments:
c     eorder   (int,ar(n),out)   - element order vector
c     NumEl     (int,sc,out)     - number of elements to process in solver
c
#include "impcom.inc"
#include "constants.inc"
c
      integer  etyiqr, elmget

      integer  n, eorder(n), NumEl, i, nnod,
     &         ityp, itld, ndst, nodes(W_NdMax), elData(EL_DIM)
      double precision dperr(3)

#if defined(DEBUG)
      call SubBeg('DefElemOrder')
#endif

c --- set initial order
      NumEl = 0
      itld = 0
      ndst = 0
      
      do i = 1,n
c --- get the number of stiffness nodes for this element
         nnod = elmget(i,elData(1),nodes(1))
         if (nnod .le. 0) cycle
         if (nnod .gt. W_NdMax) then
             dperr(1) = dble(nnod)
             dperr(2) = dble(i)
             dperr(3) = dble(W_NdMax)
             call erhandler('DefElemOrder',5000,3,'Number of nodes %I in
     &                       element %I exceeds maximum of %I.',
     &                       dperr(1),' ')
             call StopAll(-1)
         endif
         ityp = elData(EL_TYPE)
         if (ityp .le. 0) then
             dperr(1) = dble(i)
             dperr(1) = dble(ityp)
             call erhandler('DefElemOrder',5000,3,'Problems with element
     &                       %I having a wrong type of %I.',
     &                       dperr(1),' ')
             call StopAll(-1)
         endif
         if (ityp .ne. itld) then
             itld = ityp
             ndst = etyiqr(ityp,-NMNDST)
         endif
         
c --- only count element if it has stiffness nodes
c         if (ndst .gt. 0) then
             NumEl = NumEl + 1
             eorder(NumEl) = i
c         endif
      enddo


#if defined(DEBUG)
      call SubEnd('DefElemOrder')
#endif

      return
      end
