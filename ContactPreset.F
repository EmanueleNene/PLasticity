      subroutine ContactPreset (eOrder,maxNodes)
c --- primary function: allocate space for contact
c
c  input arguments:
c     eOrder     (int,ar(D_NumEl),in) - element order vector
c     maxNodes   (int,sc,in)          - maximum number of nodes (compressed)
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "cont_info.inc"
#include "dof_info.inc"
#include "solu_info.inc"
c
      integer  elmget, etyget, rlinqr, InqSMP
      PTRFTN   MemAlloc

      integer eOrder(D_NumEl),maxNodes,i,j,ielm,iel,ival(3),
     & nodes(W_NdMax),elData(EL_DIM),nnod,ityp,jtyp,elChar(IELCSZ),
     & ifoun,ireal,NumR,ityold,imast,islav,icont0,itarg0,nproc


      pointer (piRlsetL,Rlset)
      integer  Rlset(*)

#if defined(DEBUG)
      call SubBeg('ContactPreset')
#endif

c --- determine number of foundations and initialize
      K_nSeg = 0
      K_nFoun = 0
      K_nMast = 0
      ityold = 0
c
      j = 0
      do iel = 1,D_NumEl
         ielm = eOrder(iel)
         if (ielm .le. 0) cycle
         nodes(3) = 0    !  fix for line slave
         nnod = elmget(ielm,elData(1),nodes(1))
         if (nnod .le. 0) cycle
         ityp = elData(EL_TYPE)
         if (ityp .ne. ityold) then
             i = etyget(ityp,elChar(1))
             ityold = ityp
         endif
         jtyp = elChar(JETYP)
         if (jtyp.eq.W_eSlave .or. jtyp.eq.W_eMaster) j = j + 1
         if (jtyp.eq.W_eSlave) then
             if (elChar(KYOP1).eq.1) then
                 elData(EL_SHAPE) = 7
             elseif (elChar(KYOP1).eq.2) then
                 elData(EL_SHAPE) = 99
             else
                 if (nodes(3).eq.0 .or. nnod.eq.2) then
                     elData(EL_SHAPE) = 1
                 else
                     elData(EL_SHAPE) = 2
                 endif
             endif
             call elmput (ielm,elData(1),nnod,nodes(1))
         endif
      enddo

      if (j .eq. 0) goto 999

      ityold = 0
      NumR = rlinqr(0,DB_MAXDEFINED)
      NumR = max(1,NumR)
      piRlsetL = MemAlloc(3*NumR,W_LenInt,'piRlsetL        ')
      call vi_zero (Rlset(1),3*NumR)

      M_PtrElSet = MemAlloc(j,W_LenInt,'M_PtrElSet      ')
      call vi_zero (M_ElSet(1),j)

      do iel = 1,D_NumEl
         ielm = eOrder(iel)
         if (ielm .le. 0) cycle
         call vi_zero (nodes(1),W_NdMax)
         nnod = elmget(ielm,elData(1),nodes(1))
         if (nnod .le. 0) cycle
         ityp = elData(EL_TYPE)
         if (ityp .ne. ityold) then
             i = etyget(ityp,elChar(1))
             ityold = ityp
         endif
         jtyp = elChar(JETYP)
c --- count the number of contact elements, foundations, and segments
         if (jtyp.eq.W_eMaster) then

             K_nMast = K_nMast + 1
             ireal = elData(EL_REAL)
             if (ireal.gt.NumR) then
                 call erhandler ('ContactPreset',5000,3,
     &                           'Element %i uses undefined real
     &                            constant set.',
     &                            dble(ielm),' ') 
                 call StopAll(-1)
             endif
             
             call conget (ielm,ival(1))

             if (A_KeyFirstPass.eq.1) then
                 ifoun = Rlset(3*ireal-2)
                 if (ifoun.le.0) then
                     ifoun = K_nFoun + 1
                     Rlset(3*ireal-2) = ifoun
                 endif
                 Rlset(3*ireal-1) = Rlset(3*ireal-1) + 1
                 if (ival(3).lt.0) then
                     ival(3) = -ifoun
                 else
                     ival(3) = ifoun
                 endif
                 call conput (ielm,ival(1))
             else
                 ifoun = abs(ival(3))
                 if (ifoun.le.0) ifoun = Rlset(3*ireal-2)
                 Rlset(3*ireal-2) = ifoun
                 Rlset(3*ireal-1) = Rlset(3*ireal-1) + 1
                 if (ival(3).eq.0) then
                     ival(3) = ifoun
                     call conput (ielm,ival(1))
                 endif
             endif

             if (ifoun .gt. K_nFoun) K_nFoun = ifoun

c --- count the number of contact elements, foundations, and segments
         elseif (jtyp.eq.W_eSlave) then
             K_nSeg = K_nSeg + 1
             ireal = elData(EL_REAL)
             call conget (ielm,ival(1))
             if (A_KeyFirstPass.eq.1) then
                 ifoun = Rlset(3*ireal-2)
                 if (ifoun.le.0) then
                     ifoun = K_nFoun + 1
                     Rlset(3*ireal-2) = ifoun
                 endif
                 Rlset(3*ireal) = Rlset(3*ireal) + 1
                 ival(1) = ifoun
                 call conput (ielm,ival(1))
             else
                 ifoun = abs(ival(1))
                 if (ifoun.le.0) ifoun = Rlset(3*ireal-2)
                 Rlset(3*ireal-2) = ifoun
                 Rlset(3*ireal) = Rlset(3*ireal) + 1
                 if (ival(1).eq.0) then
                     ival(1) = ifoun
                     call conput (ielm,ival(1))
                 endif
             endif
             if (ifoun .gt. K_nFoun) K_nFoun = ifoun
         endif
      enddo


      if ((K_nMast+K_nSeg).le.0.and.M_PtrElSet.ne.PTRFTNNULL) then
         call MemFree (M_PtrElSet)
      else
c --- correct the contact and slave element order
         imast = 0
         islav = K_nMast
         do ireal = 1,NumR
            if (Rlset(3*ireal-2).gt.0) then
               imast = imast + Rlset(3*ireal-1)
               islav = islav + Rlset(3*ireal)
               Rlset(3*ireal-1) = imast - Rlset(3*ireal-1)
               Rlset(3*ireal) = islav - Rlset(3*ireal)
            endif
         enddo
         ityold = 0

         do iel = 1,D_NumEl
            ielm = eOrder(iel)
            if (ielm .le. 0) cycle
            nnod = elmget(ielm,elData(1),nodes(1))
            if (nnod .le. 0) cycle
            ireal = elData(EL_REAL)
            ifoun = Rlset(3*ireal-2)
            if (ifoun.le.0) cycle
            ityp = elData(EL_TYPE)
            if (ityp .ne. ityold) then
               i = etyget(ityp,elChar(1))
               ityold = ityp
            endif
            jtyp = elChar(JETYP)

c --- count the number of contact elements, foundations, and segments
            if (jtyp.eq.W_eMaster) then
               Rlset(3*ireal-1) = Rlset(3*ireal-1) + 1
               imast = Rlset(3*ireal-1)
               M_ElSet(imast) = ielm
            elseif (jtyp.eq.W_eSlave) then
               Rlset(3*ireal) = Rlset(3*ireal) + 1
               islav = Rlset(3*ireal)
               M_ElSet(islav) = ielm
            endif
c
         enddo
         icont0 = 0
         itarg0 = K_nMast
         do ireal = 1,NumR
            if (Rlset(3*ireal-2).gt.0) then
               imast = Rlset(3*ireal-1) - icont0
               islav = Rlset(3*ireal) - itarg0
c --- correct the contact order (sort array into ascending order)
               call vi_sort(M_ElSet(icont0+1),imast)
c --- correct the slave order (sort array into ascending order)
               call vi_sort(M_ElSet(itarg0+1),islav)
               icont0 = Rlset(3*ireal-1)
               itarg0 = Rlset(3*ireal)
            endif
         enddo
      endif
      call MemFree (piRlsetL)

c --- space for contact data
      if (K_nMast.gt.0 .or. K_nFoun.gt.0) then
c --- set number of items per database entry
          K_nPerFounI = 26
          K_nPerFounDp = 41
          K_nPerSeg = 19
          K_nPerMast = 9
          nproc = InqSMP(0)
c --- get dp-memory space for foundation and contact databases
          j = K_nFoun*K_nPerFounDp + maxNodes*D_NumDof

          K_PtrFounDp = MemAlloc(j,W_LenDouble,'K_PtrFounDp     ')
          K_PtrNdCoord = K_PtrFounDp + 
     &                   K_nFoun*K_nPerFounDp*W_AddrStrideDp
         
c --- get in-memory space for foundation and contact databases
          j = K_nSeg*K_nPerSeg + K_nFoun*K_nPerFounI

          K_PtrSeg= MemAlloc(j,W_LenInt,'K_PtrSeg        ')
          K_PtrFounI  = K_PtrSeg + K_nSeg*K_nPerSeg*W_AddrStrideI

          j = K_nMast*K_nPerMast
          K_PtrContInfo = MemAlloc(j,W_LenInt,'K_PtrContInfo   ')
         
          j = maxNodes
          K_PtrPiloNd  = MemAlloc(j,W_LenInt,'K_PtrPiloNd     ')
         
          j = K_nSeg*nproc
          K_PtrConSrch  = MemAlloc(j,W_LenByte,'K_PtrConSrch    ')
                  
      endif

 999  continue

#if defined(DEBUG)
      call SubEnd('ContactPreset')
#endif

      return
      end
