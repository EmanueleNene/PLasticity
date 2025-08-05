      subroutine RemoveEnds (Elptch,ndel)
c --- primary function: remove free ends of the fibers
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "dof_info.inc"
c
      logical done
      integer ElemAttachedSet,elmget
      integer Elptch(*),nodes(W_NdMax),elData(EL_DIM),
     &        nElem,node,i,j,k,ka,kb,nnod,ndel

      ndel = 0
c --- re-entry point
10    continue  
      done = .true.
c --- go through all the nodes
      do 100 i = 1,D_FwdLen
         node = M_FrwdList(i)
         if (node.le.0) cycle
c --- check if node is totally free of any contraint or force
         ka = (node-1)*D_NumDof
         do j = 1,2 ! D_NumDof count only x and y
            kb = ka+j
            if (M_DofBits(kb).ne.0) goto 100
            if (D_NumCp.gt.0) then
                if (M_SlaveCpDofs(kb).ne.kb) goto 100
            endif
         enddo
         nElem = ElemAttachedSet(1,i,Elptch(1))
         k = 0
         do 200 j = 1,nElem
            nnod = elmget(Elptch(j),elData(1),nodes(1))
            if (nnod.le.0) goto 200
c --- check if it is a mid-node
            if (nnod.gt.2 .and. nodes(3).eq.i) goto 100
c --- found an alive element
            if (elData(EL_DEAD).eq.0) then
c --- if it is the second alive element attached to this node then skip it
                if (k.ne.0) goto 100
c --- store the element index number
                k = j
            endif
  200    continue
c --- if this node is empty, go to the next one
         if (k.eq.0) goto 100
c --- mark this node as inactive
         M_FrwdList(i) = 0
         M_BackList(i) = 0
         ndel = ndel + 1
c --- delete the element
         nnod = elmget(Elptch(k),elData(1),nodes(1))
         elData(EL_DEAD) = 1
         call elmput(Elptch(k),elData(1),nnod,nodes(1))
         done = .false.
  100  continue
                  
      if (.not.done) goto 10

      return
      end