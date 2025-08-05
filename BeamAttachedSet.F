      function BeamAttachedSet (nNodes,Nodes,ElemList)
c --- primary function: get list of alive beam elements attached to node set
c                        if nNodes.lt.0 returns a number of beam elements connected
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
#include "constants.inc"
c
      integer  BeamAttachedSet      
      integer  ElemAttachedSet, nNodes, Nodes(nNodes), ElemList(*), 
     &         elmget, jelm, i, j, k, elData(EL_DIM), tmp(W_NdMax),
     &         etyiqr, vi_compress, nelem

      nelem = ElemAttachedSet(abs(nNodes),Nodes(1),ElemList(1))
      
      if (nNodes.gt.0) then
          do i = 1,nelem
             jelm = ElemList(i)
             j = elmget(jelm,elData(1),tmp(1))
             if (elData(EL_DEAD).eq.1 .or. 
     &           etyiqr(elData(EL_TYPE),-JBEAM).eq.0) ElemList(i) = 0
          enddo
          BeamAttachedSet = vi_compress(ElemList(1),nelem)
      else
          k = nelem
          do i = 1,k
             jelm = ElemList(i)
             j = elmget(jelm,elData(1),tmp(1))
             if (elData(EL_DEAD).eq.1 .or. 
     &           etyiqr(elData(EL_TYPE),-JBEAM).eq.0) nelem = nelem - 1
          enddo
          BeamAttachedSet = nelem
      endif
          
      
      return
      end

