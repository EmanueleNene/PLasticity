      subroutine AdjForCoup
c --- primary function: modify vector for couplings
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
#include "in_mem.inc"
c
      integer  vi_find,cpget,forget,disget
      PTRFTN   MemAlloc

      integer  maxe,nval,nfor,mloc,knode,maxterm,
     &         i,j,idf,kdf,ldf,ii,node,mast,node1,nodrl1,ndis,nc
      double precision  dvalue(4),Fmaster(4),Dmaster(4)

      pointer (ptr_iLoc,iLoc)
      integer iLoc(*)

#if defined(DEBUG)
      call SubBeg('AdjForCoup')
#endif

c --- get maximum number of terms
      maxterm = D_LenCp + 1
      ptr_iLoc = MemAlloc(maxterm,W_LenInt,'ptr_iLoc        ')

c --- loop on the coupled sets
      do i = 1,D_NumCp
         nval = cpget(i,iLoc(1))
         if (nval .le. 0) cycle

c --- determine the active dof
         j = iLoc(nval+1)
         idf = vi_find(D_DofList(1),j,D_NumDof)
         if (idf .le. 0) cycle

         ldf = D_DofList(idf)

c --- determine the last element on this set
         nc = 0
         node1 = 0
         maxe = 0
         mloc = 0

         do j = 1,nval
            node = iLoc(j)
            knode = M_FrwdList(node)
            if (knode .le. 0) cycle

            if (node1 .le. 0) then
                node1  = knode
                nodrl1 = node
                mloc = j
            endif
            kdf = (knode-1)*D_NumDof + idf
            nc = nc + 1
         enddo

         if (nc.le.1) cycle

c --- get any forces or constraints on the master
         mast = (node1-1)*D_NumDof + idf
         call v_zero (Fmaster(1),4)
         if (btest(M_DofBits(mast),W_ForcBit)) then
             nfor = forget(nodrl1,ldf,Fmaster(1))
c --- in case where F is on slave only we will leave F=0 on master
             if (Fmaster(1).eq.0.0d0) nfor = 0
         else
             nfor = 0
         endif
         call v_zero (Dmaster(1),4)
         if (btest(M_DofBits(mast),W_DispBit)) then
             ndis = disget(nodrl1,ldf,Dmaster(1))
         else
             ndis = 0
         endif

c --- modify M_DofBits for all slave DOFs
         do j = mloc+1,nval
            node  = iLoc(j)
            knode = M_FrwdList(node)
            if (knode.le.0) cycle
            kdf = (knode-1)*D_NumDof + idf
c --- identify the master for this slave
            if (M_PtrMasterCp.ne.PTRFTNNULL) M_MasterCp(kdf) = mast

c --- set the coupling bit
            M_DofBits(kdf) = ibset(M_DofBits(kdf),W_CpBit)

c --- sum forces on master node
            if (btest(M_DofBits(kdf),W_ForcBit)) then
                ii = forget(node,ldf,dvalue(1))
                if (ii.gt.0) then
                    nfor = 1
                endif
            endif

c --- move constraints to master node
            if (ndis.le.0) then
                if (btest(M_DofBits(kdf),W_DispBit)) then
                    ndis = 1
                    M_DofBits(kdf) = ibclr(M_DofBits(kdf),W_DispBit)
                endif
            endif
         enddo

c --- put total force on master node
         if (nfor .gt. 0) then
             M_DofBits(mast) = ibset(M_DofBits(mast),W_ForcBit)
         endif

c --- put constraint (if any) on master node
         if (ndis .gt. 0) then
             M_DofBits(mast) = ibset(M_DofBits(mast),W_DispBit)
         endif

      enddo

      call MemFree (ptr_iLoc)

#if defined(DEBUG)
      call SubEnd('AdjForCoup')
#endif

      return
      end
