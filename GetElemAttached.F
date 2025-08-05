      subroutine GetElemAttached (node,numElems,ptr_elems)
c --- primary function: get the elements attached to the node
c
c  input arguments:
c     node      (int,sc,in)         - node number
c
c  output arguments:
c     numElems  (int,sc,out)        - number of elements
c     ptr_elems (int,sc,out)        - pointer to element array
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
#include "in_mem.inc"
c
      integer  elmget
      pointer (ptr_elems,Elems)

      PTRFTN  MemAlloc

      integer  iStart, iEnd, elData(EL_DIM), nodes(W_NdMax), nnod, node,
     &         numElems, elord, elem, i, AllElem(D_NumEl), Elems(*)


#if defined(DEBUG)
      call SubBeg('GetElemAttached')
#endif

      numElems = 0

c --- loop over all elements and get general contact elements
      call RangeSMP (D_NumEl,iStart,iEnd)
      do 500 elord = iStart,iEnd

c --- get the actual element number
         elem = M_ElemOrd(elord)
         AllElem(elord) = 0

         if (elem .le. 0) goto 500

c --- get the element data and nodes
         call vi_zero (nodes(1),W_NdMax)
         nnod = elmget(elem,elData(1),nodes(1))

         if (nnod .le. 0) goto 500

         do i = 1, min(3,nnod)
            if (node.eq.nodes(i)) then
                numElems = numElems + 1
                AllElem(elord) = elem
                goto 500
            endif
         enddo

 500  continue

      ptr_elems = MemAlloc (numElems,W_LenInt,'ptr_elems       ')
      i = 1
      do elord = 1, D_NumEl
         if (AllElem(elord).gt.0) then
             Elems(i) = AllElem(elord)
             i = i + 1
         endif
      enddo

      call vi_sort(Elems(1),NumElems)


#if defined(DEBUG)
      call SubEnd('GetElemAttached')
#endif

      return
      end
