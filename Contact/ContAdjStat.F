      function ContAdjStat(nodes,nodec,isegda,iseg0,nodbac,nPerSeg,nSeg)
c --- primary function: compute the relations between an integration point
c                       of the contact element and the corresponding foundation
c
c  input arguments:
c          isegda(1,iseg)        = code of segment type (icode)
c                                = 1   : straight line
c                                = 2   : parabolla
c                                = 7   : quadrangle
c          isegda(2,iseg)        = fondation number
c          isegda(j,iseg)        = 4 nodes (j=3-10) on seg. iseg
c          isegda(j,iseg)        = 4 adjacent segments (j=11-14)
c                                = 2 adjacent segments (j=11,12)
c          isegda(15,iseg)       = segment element number
c                                  < slave element death
c          isegda(16,iseg)       = attached beam element number, if it is
c          isegda(17,iseg)       = attached contact element number, if it is
c
c  output arguments:
c          ContAdjStat   (int,sc,out)
c                                = 0  : none connection
c                                = 1  : connection
c
#include "impcom.inc"
#include "dof_info.inc"
c
      integer nodec,nnode,nodes(nodec),iseg0,iseg,nodbac(*),i,
     &        ContAdjStat,nloop,nPerSeg,nSeg,isegda(nPerSeg,nSeg),
     &        inode,i1,icode

#if defined(DEBUG)
      call SubBeg('ContAdjStat')
#endif

      ContAdjStat = 0
      if (iseg0.eq.0) goto 999
      do 900 nloop = 1,5
         if (nloop.eq.1) iseg = iseg0
         if (nloop.eq.2) iseg = isegda(11,iseg0)
         if (nloop.eq.3) iseg = isegda(12,iseg0)
         if (nloop.eq.4) iseg = isegda(13,iseg0)
         if (nloop.eq.5) iseg = isegda(14,iseg0)
         if (iseg.eq.0) goto 900
         icode = isegda(1,iseg)
         nnode = 3
         if (icode.eq.7) nnode = 4
         do 120 inode = 1,nnode
            i1 = isegda(2+inode,iseg)
            if (i1.le.0) goto 120
            i1 = nodbac((i1-1)/D_NumDof+1)
            do 100 i = 1,nodec
               if (i1.eq.nodes(i)) then
                   ContAdjStat = 1
                   if (nloop.eq.1) ContAdjStat = 2
                   goto 999
               endif
 100        continue
 120     continue
 900  continue


 999  continue

#if defined(DEBUG)
      call SubEnd('ContAdjStat')
#endif

      return
c
      end

