      function ElemAttachedSet (nNodes,Nodes,ElemList)
c --- primary function: get list of elements attached to node set
c
c  input arguments:
c        nNodes          (int,sc,in)         - number of nodes
c        Nodes           (int,ar(nNodes),in) - external node numbers for operation
c
c  output arguments:
c        ElemAttachedSet (int,sc,out)        - number of elements in list
c        ElemList        (int,ar(*),out)     - attached external element numbers
c
#include "impcom.inc"
c
      integer  ElemAttachedSet, nNodes, Nodes(nNodes), ElemList(*),
     &         jNode, i, j, NxrefElSet,
     &         numElems, Elems(*), work
      pointer (ptr_elems,Elems)


      ElemAttachedSet = NxrefElSet(nNodes,Nodes(1),ElemList(1))
      goto 999  ! change that

      ElemAttachedSet = 0
c --- get data from node-elem xref
      ptr_elems = PTRFTNNULL

c --- NOTE: performace tuning needed

      do i = 1,nNodes
         jNode = Nodes(i)
         if (jNode .eq. 0) cycle

         call GetElemAttached(jNode,numElems,ptr_elems)

         if (i .eq. 1) then
c --- for first node simply copy its elements into final list
             do j = 1,numElems
                ElemAttachedSet = ElemAttachedSet + 1
                ElemList(ElemAttachedSet) = Elems(j)
             enddo
         else
c --- for other nodes vi_intersect final list with this node's elements
             work = ElemAttachedSet
	       call vi_intersect(work,ElemList(1),numElems,Elems(1),
     &                         ElemAttachedSet,ElemList(1))
         endif

         if (ptr_elems .ne. PTRFTNNULL) call MemFree(ptr_elems)

      enddo

  999 continue

      return
      end

