      subroutine Cont_4 (NumEl,BacLen,NumDof,nodbac,nPerSeg,
     &                   nSeg,nPerFounI,nPerFounDp,nFoun,isegda,ififo,
     &                   founda,icnda,nopilo,ncon,nPerMast,nMast)
c --- primary function: verify whether slave element belongs to two pass
c                       find the adjacent segments for slave element
c
c  input arguments:
c     BacLen    (int,sc,in)         - number of active nodes
c     NumDof    (int,sc,in)         - number of active dof per node
c     nodbac    (int,ar(BacLen),out)- internal to exter. list
c     nPerSeg   (int,sc,in)         - number of items per segment database entry
c     nSeg      (int,sc,in)         - number of segments on foundation
c     nPerFounI (int,sc,in)         - number of integer items per foundation
c                                     database entry
c     nFoun    (int,sc,in)          - number of foundations
c
c  output arguments:
c     isegda   (int,ar(nPerSeg,nSeg),out)
c                             - segment information
c          isegda(1,iseg)   = code of segment type (icode)
c                           = 1   : straight line
c                           = 2   : parabolic arc
c          isegda(2,iseg)   = fondation number
c          isegda(j,iseg)   odes (j=3-10) on seg. iseg
c          isegda(j,iseg)   = 4 adjacent segments (j=11-14)
c                           = 2 adjacent segments (j=11,12
c                                for seg. iseg
c                           = 0  (icode > 10)
c                             --> none adjacent segment
c          isegda(15,iseg)  = segment element number
c                             < slave element death
c          isegda(16,iseg)  = attached beam element number, if it is
c          isegda(17,iseg)   attached contact element number, if it is
c          isegda(18,iseg)        = remark
c
c     ififo    (dp,ar(nPerFounI,nFoun),out)
c                              - foundation information (integer)
c          ififo(1,ifoun)   = number of segments
c          ififo(4,ifoun)   = constrain condition in x-dir.
c          ififo(5,ifoun)   = constrain condition in y-dir.
c          ififo(6,ifoun)   = constrain condition in z-dir.
c          ififo(7,ifoun)   = dimensions
c          ififo(8,ifoun)   = 1st segment numbering on this foundation
c          ififo(9,ifoun)   = 2nd segment numbering on this foundation
c          ififo(10,ifoun)  = lowest nodal numbering on this foundation
c          ififo(12,ifoun)  = real constant number (ireal)
c          ififo(13,ifoun)  = contact pair key
c                             0 - no contact pair
c                             1 - only slave element is defined
c                             2 - only contact element is defined
c                             3 - contact pair is defined
c          ififo(14,ifoun)  = type of associated contant element type (ityp)
c          ififo(15,ifoun)  = type of deformable slave surface
c                            -3 - rigid slave surface no node fix
c                            -2 - mixed type surface (may wrong)
c                            -1 - rigid slave surface
c                             0 - none
c                             1 - deformable slave surface
c                             2 - deformable slave surface (two pass)
c                             3 - self contact
c          ififo(16,ifoun)  = contact element number who is closest to
c                                current slave
c          ififo(17,ifoun)  = slave element number who is pair of closest
c                             contact element
c          ififo(18,ifoun)  = iteration number of surface adjust
c
#include "impcom.inc"
#include "constants.inc"
#include "conv_stat.inc"
c
      PTRFTN   MemAlloc
      integer BacLen,NumDof,nPerSeg,nSeg,nPerFounI,
     &  nFoun,elmget,etyiqr,rlget,ElemAttachedSet,elmgct,
     &  nPerMast,nMast,nPerFounDp,elmgcn
      integer isegda(nPerSeg,nSeg),ififo(nPerFounI,nFoun),
     &  nodbac(BacLen),icnda(nPerMast,nMast),nopilo(BacLen),ncon
      integer elData(EL_DIM),nodes(W_NdMax),ival(3),ifoun,
     &  iseg,jseg,ielm,nnod,ityp,inode,iseg1,iseg2,iseg3,
     &  nvois,ln1,ln3,ln1v,ln3v,ivois1,ivois3,ietypz,
     &  icon,ibeam,j1,knode,lnode,iddl,n1,n2,
     &  coninf(3),nlist(8),n,nElem,jnode,NumEl,jcon,
     &  jelm,ifoun1,ielas,ietypt,kfoun,
     &  inopil,inopil0,icode,jcode,n3,n4,ino,inoexc,iropil,
     &  nside,nnode,iielm,k,k2,jetypz,jjelm,
     &  elData1(EL_DIM),nodes1(W_NdMax),maxn

      double precision founda(nPerFounDp,nFoun),rvr(W_NumContParam),
     &  dvet(11)

      pointer (piElptch,Elptch)
      integer Elptch(*)

#if defined(DEBUG)
      call SubBeg('Cont_4')
#endif

      piElptch = PTRFTNNULL

      if (NumEl.gt.0)
     &    piElptch = MemAlloc(NumEl,W_LenInt,'piElptch        ')

c --- foundation information
      do 500 kfoun = 1,nFoun
         ifoun = ififo(24,kfoun)
         if (ififo(12,ifoun).le.0) goto 500
         call v_zero (rvr(1),W_NumContParam)
         iseg1 = rlget(ififo(12,ifoun),rvr(1))
         if (ififo(1,ifoun).le.0) goto 490
         iseg1 = ififo(8,ifoun)
         iseg2 = ififo(9,ifoun)
         iseg3 = iseg1 + ififo(1,ifoun) - 1
         if (iseg2 .ne. iseg3)
     &       call erhandler ('Cont_4',5000,3,
     &                       'Contact specs are corrupted.', 0.0d0,' ')
         inopil = (ififo(2,ifoun)-1)/NumDof + 1
         inopil0 = inopil
         ififo(4,ifoun) = 1
         ififo(5,ifoun) = 1
         ififo(6,ifoun) = 1
         ififo(18,ifoun) = 0
         ififo(19,ifoun) = 0
         founda(24,ifoun) = 0.d0
         ielas = 0
         do 200 iseg = iseg1,iseg2
            icode = isegda(1,iseg)
            nnod = 3
            if (icode.eq.7)  nnod = 4
            if (icode.eq.99) nnod = 2 
c --- check boundary condition
            if (inopil0.gt.0) goto 35
            do 30 inode = 1,nnod
               knode = (isegda(2+inode,iseg)-1)/NumDof + 1
               if (knode .le. 0) goto 30
               knode = nodbac(knode)
               call conchk1 (knode,coninf(1),ififo(7,ifoun))
               do 20 iddl =1,ififo(7,ifoun)
                  if (ififo(3+iddl,ifoun) .eq. 0) goto 20
                  ififo(3+iddl,ifoun) = min(ififo(3+iddl,ifoun),
     &                                      coninf(iddl))
   20          continue
   30       continue

   35       continue
   
            if (ififo(15,ifoun).lt.1) goto 100
            if (icode.ge.10) goto 100
            ibeam = abs(isegda(16,iseg))

c --- verify whether slave belongs to two pass
            isegda(17,iseg) = 0
            if (abs(isegda(15,iseg)).gt.0) then
c --- get stiffness
                ifoun1 = elmgct (abs(isegda(15,iseg)),dvet(1))
                ielas = ielas + 1
                founda(24,ifoun) = founda(24,ifoun) + dvet(4)
            endif

            if (icode.eq.1.or.icode.eq.2) then
                n1 = (isegda(5,iseg)-1)/NumDof + 1
                n1 = nodbac(n1)
                nlist(1) = n1
                n2 = (isegda(3,iseg)-1)/NumDof + 1
                n2 = nodbac(n2)
                nlist(2) = n2
                n = 2
                n3 = 0
                n4 = 0
            elseif (icode.eq.7) then
                n1 = (isegda(3,iseg)-1)/NumDof + 1
                n1 = nodbac(n1)
                nlist(1) = n1
                n2 = (isegda(4,iseg)-1)/NumDof + 1
                n2 = nodbac(n2)
                nlist(2) = n2
                n3 = (isegda(5,iseg)-1)/NumDof + 1
                n3 = nodbac(n3)
                nlist(3) = n3
                n4 = (isegda(6,iseg)-1)/NumDof + 1
                n4 = nodbac(n4)
                nlist(4) = n4
                n = 4
            endif

            do inode = 1,n
               nlist(inode) = abs(nlist(inode))
            enddo

            nElem = ElemAttachedSet(n,nlist(1),Elptch(1))

            do 50 icon = 1,nElem
               ielm = Elptch(icon)
               if (ielm.eq.isegda(15,iseg)) goto 50
               nnod = elmget(ielm,elData(1),nodes(1))
               if (nnod .le. 0) goto 50
               ityp = elData(EL_TYPE)
               ietypz = etyiqr(ityp,-JETYP)
               if (ietypz.ne.W_eMaster) goto 50
               call conget (ielm,ival(1))
               if (icnda(2,abs(ival(2))).ne.ibeam) goto 50
               nnod = 2
               do inode = 1,nnod
                  if (nodes(inode).eq.n1) then
                      do jnode = 1,nnod
                         if (nodes(jnode).eq.n2) then
                             if (n3.eq.0) goto 40
                             do knode = 1,nnod
                                if (nodes(knode).eq.n3) then
                                    if (n4.eq.0) goto 40
                                    do lnode = 1,nnod
                                       if (nodes(lnode).eq.n4) goto 40
                                    enddo
                                endif
                             enddo
                         endif
                      enddo
                  endif
               enddo
               goto 50
  40           continue 
               ifoun1 = abs(ival(3))
               if (ifoun1.le.0.or.ifoun1.gt.nFoun) goto 50
               if (ififo(15,ifoun1).le.0) goto 50
               if (ififo(13,ifoun1).ne.3) goto 50
c --- we do not handle symmetric contact, but only self-contact
               if (ifoun1.ne.ifoun) goto 50
               isegda(17,iseg) = ielm
               call conget (abs(isegda(15,iseg)),ival(1))
               ival(3) = ielm
               call conput (abs(isegda(15,iseg)),ival(1))
               if (ififo(15,ifoun).eq.1) ififo(15,ifoun) = 2
               if (ififo(19,ifoun).eq.0) then
                   ififo(19,ifoun) = elData(EL_REAL)
                   ififo(20,ifoun) = ifoun1
               endif
               if (ififo(19,ifoun).ne.elData(EL_REAL)) then
                   if (ififo(20,ifoun1).ne.0.and.ififo(20,ifoun1).ne.
     &                 ifoun) goto 50
                   ififo(19,ifoun) = elData(EL_REAL)
                   ififo(20,ifoun) = ifoun1
               endif
               if (ififo(20,ifoun1).ne.0.and.ififo(20,ifoun1).ne.ifoun)
     &             goto 50
               goto 90
   50       continue
 
   90       continue

            if (isegda(17,iseg).ne.0) then
                if (ififo(19,ifoun).ne.elData(EL_REAL)) then
                    ififo(15,ifoun) = 1
                endif
            else
                if (ififo(19,ifoun).eq.0) ififo(19,ifoun) =
     &                                   -ififo(12,ifoun)
                if (ififo(19,ifoun).ne.-ififo(12,ifoun)) then
                    ififo(15,ifoun) = 1
                endif
            endif

  100       if (icode.eq.99) goto 200
            nnod = 3
            if (icode.eq.7) nnod = 4
            do 150 ino = 1,nnod
               inoexc = (isegda(ino+2,iseg)-1)/NumDof + 1
               if (inoexc.le.0 .or. inoexc.eq.inopil) goto 150
               if (ififo(3,ifoun).ne.0) then
                   iropil = (abs(ififo(3,ifoun))-1)/NumDof + 1
                   if (inoexc.eq.iropil) goto 150
               endif
c --- form the pilot node list
               if (inopil0.gt.0) then
                   nopilo(inoexc) = ifoun
               else
                   nopilo(inoexc) = -ifoun
               endif
  150       continue
  200    continue
         if (ielas.gt.0) then
             founda(24,ifoun) = founda(24,ifoun)/ielas
         endif

c --- check the adjacent segments and write on isegda (7-->9).
         if (ififo(1,ifoun).eq.1) goto 490
         iseg3 = iseg2 - 1

c --- for beam-beam
         do 480 iseg = iseg1,iseg3
            icode = isegda(1,iseg)
            if (icode.gt.2) goto 480
            nvois = 0 
            ivois1 = isegda(11,iseg)
            ivois3 = isegda(12,iseg)
            if (ivois1.ne.0) nvois=nvois+1
            if (ivois3.ne.0) nvois=nvois+1
            if (nvois.eq.2) goto 480
            ln1 = isegda(3,iseg)
            ln3 = isegda(5,iseg)
            j1 = iseg + 1
            do 470 jseg = j1,iseg2
               jcode = isegda(1,jseg)
               if (isegda(1,jseg).gt.4) goto 470
               ln1v = isegda(3,jseg)
               ln3v = isegda(5,jseg)

               if (ivois1.eq.0 .and. ln1.eq.ln1v) then
                   isegda(11,iseg) = jseg
                   isegda(11,jseg) = iseg
                   nvois = nvois+1
               endif

               if (ivois3.eq.0 .and. ln3.eq.ln3v) then
                   isegda(12,iseg) = jseg
                   isegda(12,jseg) = iseg
                   nvois = nvois+1
               endif

               if (ivois1.eq.0 .and. ln1.eq.ln3v) then
                   isegda(11,iseg) = jseg
                   isegda(12,jseg) = iseg
                   isegda(6,iseg) = ln1v
                   isegda(7,jseg) = ln3
                   if (icode.eq.2) then
                       isegda(9,jseg) = isegda(4,iseg)
                   else
                       isegda(9,jseg) = 0
                   endif
                   if (jcode.eq.2) then
                       isegda(8,iseg) = isegda(4,jseg)
                   else
                       isegda(8,iseg) = 0
                   endif
                   nvois = nvois + 1
               endif
               if (ivois3.eq.0 .and. ln3.eq.ln1v) then
                   isegda(12,iseg) = jseg
                   isegda(11,jseg) = iseg
                   isegda(7,iseg) = ln3v
                   isegda(6,jseg) = ln1
                   if (icode.eq.2) then
                       isegda(8,jseg) = isegda(4,iseg)
                   else
                       isegda(8,jseg) = 0
                   endif
                   if (jcode.eq.2) then
                       isegda(9,iseg) = isegda(4,jseg)
                   else
                       isegda(9,iseg) = 0
                   endif
                   nvois = nvois + 1
               endif
               if (nvois.eq.2) goto 480
  470       continue
  480    continue
  490    continue
  500 continue
      do 700 icon = 1,ncon
         ielm = icnda(8,icon)
         nnod = elmget(ielm,elData(1),nodes(1))
         if (nnod .le. 0) goto 700
         ityp = elData(EL_TYPE)
         ietypz = etyiqr(ityp,-JETYP)
         if (ietypz.ne.W_eMaster) goto 700
         ifoun = icnda(9,icon)
         if (ifoun.eq.0) goto 700
         n = 2
         n1 = nodes(1)
         n2 = nodes(2)
         nlist(1) = n1
         nlist(2) = n2

         if (ififo(19,ifoun).le.0) goto 700
         if (ififo(15,ifoun).ne.2) goto 700
 600     continue

         nElem = ElemAttachedSet(n,nlist(1),Elptch(1))

         do 650 jcon = 1,nElem
            jelm = Elptch(jcon)
            if (ielm.eq.jelm) goto 650
            nnod = elmget(jelm,elData(1),nodes(1))
            if (nnod .le. 0) goto 650
            if (ififo(19,ifoun).ne.elData(EL_REAL)) goto 650
            ityp = elData(EL_TYPE)
            ietypt = etyiqr(ityp,-JETYP)
            if (ietypt.ne.W_eSlave) goto 650
            call conget (jelm,ival(1))
            ifoun1 = abs(ival(1))
            if (ifoun1.le.0.or.ifoun1.gt.nFoun) goto 650
            if (ififo(15,ifoun1).le.0) goto 650
            if (ififo(13,ifoun1).ne.3) goto 650
            nnod = 4
            if (nodes(3).eq.nodes(4).or.nodes(4).eq.0) then
                nnod = 3
            endif            
            do 640 inode = 1,nnod
               if (nodes(inode).eq.0) goto 640
               if (nodes(inode).eq.n1) then
                   if (n2.eq.0) goto 700
                   do jnode = 1,nnod
                      if (nodes(jnode).eq.n2) then
                          if (n3.eq.0) goto 700
                          do knode = 1,nnod
                              if (nodes(knode).eq.n3) then
                                  if (nnod.eq.3) goto 700
                                  if (n4.eq.0) goto 700
                                  do lnode = 1,nnod
                                     if (nodes(lnode).eq.n4) goto 700
                                  enddo
                              endif
                              if (nodes(knode).eq.n4) then
                                  if (nnod.eq.3) goto 700
                              endif
                          enddo
                      endif
                  enddo
               endif
 640        continue
 650     continue
         ififo(19,ifoun) = -abs(ififo(19,ifoun))
         ififo(15,ifoun) = 1
 700  continue

c --- assign elements responsible for contact nodes
      nside = 2
      do 800 icon = 1,nMast
         ielm = icnda(8,icon)
         nnode = elmget(ielm,elData(1),nodes(1))
         if (nnode .le. 0) goto 800
         ietypz = etyiqr(elData(EL_TYPE),-JETYP)
         if (ietypz.ne.W_eMaster) goto 800
         ifoun = icnda(9,icon)
         if (ifoun.eq.0) goto 800
         if (ififo(12,ifoun).le.0) goto 800
         iielm = 0
         do 880 k = 1,nside
            iielm = iielm + 1
            k2 = elmgcn (ielm,dvet(1))
            if (k2.le.0) call v_zero (dvet(1),2)
            if (abs(dvet(k)).gt.0) goto 880
            dvet(k) = ielm
            n1 = abs(nodes(iielm))
c --- find elements attached to contact node
            nElem = ElemAttachedSet(1,n1,Elptch(1))
            do 870 jcon = 1,nElem
               jelm = Elptch(jcon)
               nnode = elmget(jelm,elData1(1),nodes1(1))
               if (nnode .le. 0) goto 870
               if (ififo(12,ifoun).ne.elData1(EL_REAL)) goto 870
               jetypz = etyiqr(elData1(EL_TYPE),-JETYP)
               if (jetypz.ne.W_eMaster) goto 870
               k2 = elmgcn (jelm,dvet(1))
               if (k2.le.0) call v_zero (dvet(1),2)
               jjelm = 0
               do 860 inode = 1,nside
                  jjelm = jjelm + 1
                  if (abs(dvet(inode)).gt.0) goto 860
                  if (nodes1(jjelm).eq.n1) then
                      dvet(inode) = ielm
                      call elmpcn (jelm,dvet(1),2)
                      goto 870
                  endif
 860           continue
 870        continue
 880     continue
 800  continue

c --- add boundary nodes
      do 990 icon = 1,ncon
         ielm = icnda(8,icon)
         nnod = elmget(ielm,elData(1),nodes(1))
         if (nnod .le. 0) goto 990
         ityp = elData(EL_TYPE)
         ietypz = etyiqr(ityp,-JETYP)
         if (ietypz.ne.W_eMaster) goto 990
c --- skip beam-surface elements  
         if (etyiqr(ityp,-KYO10).eq.1) goto 990
         ifoun = icnda(9,icon)
         if (ifoun.eq.0) goto 990
         if (ififo(13,ifoun).ne.3) goto 990
         nlist(1) = nodes(1)

         nElem = ElemAttachedSet(1,nlist(1),Elptch(1))

         do 950 jcon = 1,nElem
            jelm = Elptch(jcon)
            if (ielm.eq.jelm) goto 950
            nnod = elmget(jelm,elData1(1),nodes1(1))
            if (nnod .le. 0) goto 950
            ityp = elData1(EL_TYPE)
            ietypz = etyiqr(ityp,-JETYP)
            if (ietypz.ne.W_eMaster) goto 950
c --- skip beam-surface elements            
            if (etyiqr(ityp,-KYO10).eq.1) goto 950
            call conget (jelm,ival(1))
            ifoun1 = abs(ival(3))
            if (ifoun1.ne.ifoun) goto 950
            if (nodes1(2).ne.nodes(1)) goto 950
            nodes(4) = nodes1(1)
            if (nodes1(3).ne.0.and.nodes1(3).ne.nodes1(2)) then
                nodes(6) = nodes1(3)
            else
                nodes(6) = 0
            endif
            goto 960
 950     continue
 960     continue
         nlist(1) = nodes(2)

         nElem = ElemAttachedSet(1,nlist(1),Elptch(1))

         do 970 jcon = 1,nElem
            jelm = Elptch(jcon)
            if (ielm.eq.jelm) goto 970
            nnod = elmget(jelm,elData1(1),nodes1(1))
            if (nnod .le. 0) goto 970
            ityp = elData1(EL_TYPE)
            ietypz = etyiqr(ityp,-JETYP)
            if (ietypz.ne.W_eMaster) goto 970
c --- skip beam-surface elements
            if (etyiqr(ityp,-KYO10).eq.1) goto 970
            call conget (jelm,ival(1))
            ifoun1 = abs(ival(3))
            if (ifoun1.ne.ifoun) goto 970
            if (nodes1(1).ne.nodes(2)) goto 970
            nodes(5) = nodes1(2)
            if (nodes1(3).ne.0.and.nodes1(3).ne.nodes1(2)) then
                nodes(7) = nodes1(3)
            else
                nodes(7) = 0
            endif
            goto 980
 970     continue
 980     continue
         nnod = 7
         maxn = etyiqr(ityp,-NMNDST) 
         call vi_zero (nodes(nnod+1),maxn-nnod)
         call elmpnd (ielm,nnod,nodes(1))
 990  continue
 999  continue

      call MemFree(piElptch)
      C_nPen = 0
      C_nConA = 0

#if defined(DEBUG)
      call SubEnd('Cont_4')
#endif
      return
      end