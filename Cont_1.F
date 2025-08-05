      subroutine Cont_1 (NumEl,BacLen,NumDof,nodfwd,nodbac,
     &                   nSeg,nPerFounI,nMast,nPerMast,nFoun,
     &                   ififo,icnda,icon,ipara)
c --- primary function: find underlying elements for contact and slave elements
c
c  input arguments:
c     NumEl    (int,sc,in)         - number of elements
c     BacLen   (int,sc,in)         - number of active nodes
c     NumDof   (int,sc,in)         - number of active dof per node
c     nodfwd   (int,ar(BacLen),in) - forward list for compressed nodes
c     nSeg     (int,sc,in)         - number of segments on foundation
c     nPerMast (int,sc,in)         - number of items per contact element
c     nPerFounI(int,sc,in)         - number of integer items per foundation
c                                    database entry
c     nMast   (int,sc,in)          - number of contact elements
c     nFoun    (int,sc,in)         - number of foundations
c
c  output arguments:
c     ififo    (int,ar(nPerFounI,nFoun),out)
c                                  - foundation information (integer)
c          ififo(1,ifoun)       = number of segments
c          ififo(4,ifoun)       = constrain condition in x-dir.
c          ififo(5,ifoun)       = constrain condition in y-dir.
c          ififo(6,ifoun)       = constrain condition in z-dir.
c          ififo(7,ifoun)       = dimensions
c          ififo(8,ifoun)       = 1st segment numbering on this foundation
c          ififo(9,ifoun)       = 2nd segment numbering on this foundation
c          ififo(10,ifoun)      = lowest nodal numbering on this foundation
c          ififo(12,ifoun)      = real constant number (ireal)
c          ififo(13,ifoun)      = contact pair key
c                                 0 - no contact pair
c                                 1 - only slave element is defined
c                                 2 - only contact element is defined
c                                 3 - contact pair is defined
c          ififo(14,ifoun)      = type of associated contant element type (ityp)
c          ififo(15,ifoun)      = type of deformable slave surface
c                                -3 - rigid slave surface no node fix
c                                -2 - mixed type surface (may wrong)
c                                -1 - rigid slave surface
c                                 0 - none
c                                 1 - deformable slave surface
c                                 2 - deformable slave surface (two pass)
c          ififo(16,ifoun)      = contact element number who is closest to
c                                 current slave
c          ififo(17,ifoun)      = slave element number who is pair of
c                                 closest contact element
c          ififo(18,ifoun)      = iteration number of surface adjust
c
c
c     icnda    (int,ar(nPerMast,nMast),out) - contact element information
c                               = 1 element type (ktyp)
c                               = 2 attached beam element number (ibeam)
c                               = 3 element in contact
c                               = 4 - 7 nodes number (nodes)
c                               = 8 contact element number (ielm)
c                               = 9 slave surface number (ifoun)
c     icon     (int,sc,in)         - active number of contact elements
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
c
      integer  elmgct,rlget,ElemAttachedSet,etyget,elmget,
     &         GetWriteUnit,etyiqr
      PTRFTN   MemAlloc
      integer  NumEl,BacLen,NumDof,nSeg,
     &         nPerFounI,nFoun,nPerMast,nMast
      integer  nodfwd(BacLen),nodbac(BacLen),
     &         icnda(nPerMast,nMast),ival(3),
     &         elData(EL_DIM),nodes(W_NdMax),
     &         ifoun,iel,ielm,nnod,ityp,icode,inode,ietypz,icon,
     &         ibeam,jelm,ireal,nrvr,mattab,
     &         elChar(IELCSZ),ififo(nPerFounI,nFoun),
     &         nnode,matsld,i,iott,
     &         n,nlist(2),nElem,itypold,ipara,
     &         kmin,kfoun,iwork
      double precision rvr(W_NumContParam),xyzc(11),fsglb,ebeam,fsmat

      pointer (piElptchL,Elptch)
      integer  Elptch(*)
      logical  Reod_needed
      double precision  time_e,time_b,time_we,time_wb

#if defined(DEBUG)
      call SubBeg('Cont_1')
#endif

      icon = 0
      fsglb = 1.0d0
      itypold = 0

      piElptchL = PTRFTNNULL

c --- set up bit map for patch of elements

      if (NumEl .gt. 0) 
     &    piElptchL = MemAlloc(NumEl,W_LenInt,'piElptchL       ')
      do 700 iel = 1,nMast+nSeg

         ielm = M_ElSet(iel)

         nnod = elmget(ielm,elData(1),nodes(1))
         ityp = elData(EL_TYPE)
         if (ityp .ne. itypold) then
             i = etyget(ityp,elChar(1))
             itypold = ityp
         endif
         ietypz = elChar(JETYP)
c --- foundation and segment information
         if (ietypz.eq.W_eSlave) then
             ireal = elData(EL_REAL)
             icode = elData(EL_SHAPE)
             call conget(ielm,ival(1))
c --- get the element real constant data
             call v_zero (rvr(1),W_NumContParam)
             nrvr = rlget (ireal,rvr(1))
             ifoun = abs(ival(1))
             if (ival(3).lt.0) ififo(11,ifoun) = ival(3)
             
             nnod = 3
             if (icode.eq.99) nnod = 1
             if (icode.eq.7) nnod = 4
            
c --- get max/min node numbering
             do inode = 1,nnod
                if (nodes(inode).gt.0) then
                    if (ififo(10,ifoun).eq.0) then
                        ififo(10,ifoun) = nodfwd(nodes(1))
                    endif
                    if (nodbac(ififo(10,ifoun)).gt.nodes(inode)) then
                        ififo(10,ifoun) = nodfwd(nodes(inode))
                    endif
                endif
             enddo

             ififo(1,ifoun) = ififo(1,ifoun) + 1
             if (icode.eq.99) then
                 ififo(2,ifoun) = (nodfwd(nodes(1))-1)*numdof + 1
             endif
             if (ififo(13,ifoun).eq.0) ififo(13,ifoun) = 1
             if (ififo(13,ifoun).eq.2) ififo(13,ifoun) = 3

             if (ififo(12,ifoun).eq.0) then
                 ififo(12,ifoun) = ireal
             endif

             if (ipara.eq.1) then
                 if (ififo(15,ifoun).eq.0) then
                     if (ival(2).gt.0) then
                         ififo(15,ifoun) = 1
                     else
                         ififo(15,ifoun) = -1
                         if (elChar(KYOP2).eq.0) ififo(15,ifoun)= -3
                     endif
                 endif
                 goto 700
             endif
c --- check whether the slave element attaches to others
             nnode = 2
             nlist(1) = nodes(1)
             nlist(2) = nodes(2)
c --- find attached beam element

             nElem = ElemAttachedSet(nnode,nlist(1),Elptch(1))
            
             do 50 nnode = 1,nElem
                jelm = Elptch(nnode)
                if (jelm.eq.ielm) goto 50
                i = elmget(jelm,elData(1),nodes(1))
                if (etyiqr(elData(EL_TYPE),-JBEAM).eq.0) goto 50 

                ibeam = -abs(jelm)
                if (ififo(15,ifoun).eq.0) ififo(15,ifoun) = 1
                ival(1) = ifoun
                ival(2) = ibeam
                ival(3) = 0
                call conput (ielm,ival(1))

c --- put stiffness into slave element data base
                matsld = elData(EL_MAT)
c --- find Young's modulus and the volume from attached element
c --- get material table
                call GetEx (abs(jelm),ebeam,matsld,mattab)

                fsmat = 2.d0
                if (mattab.eq.1) fsmat = fsmat*0.01d0
                if (mattab.eq.2) fsmat = fsmat/2.0d0

                if (ibeam.gt.0) then
                    fsmat=fsmat*10.d0
                endif

                fsmat = fsmat/fsglb
                ebeam = abs(ebeam)*fsmat

                call v_zero (xyzc(1),6)
                nrvr = elmgct(ielm,xyzc(1))
                xyzc(4) = ebeam
                call elmpct (ielm,xyzc(1))

                goto 700

  50         continue

             if (ififo(15,ifoun).eq.0) then
                 ififo(15,ifoun)= -1
                 if (elChar(KYOP2).eq.0) ififo(15,ifoun)= -3
             endif

             goto 700
         endif


c --- contact information
         if (ietypz.eq.W_eMaster) then

            ireal = elData(EL_REAL)
c --- get the element real constant data
            call v_zero(rvr(1),W_NumContParam)
            i = rlget(ireal,rvr(1))
            call conget (ielm,ival(1))

            ifoun = abs(ival(3))
            icon = icon + 1
            if (ififo(21,ifoun).eq.0) then
                ififo(21,ifoun) = icon
            else
                ififo(21,ifoun) = -icon
            endif
            ival(2) = icon
            call conput (ielm,ival(1))
            icnda(1,icon) = 3
            if (ipara.eq.1) then
                if (ival(3).lt.0) icnda(1,icon) = 1
                if (ival(1).lt.0) icnda(1,icon) = 3
            endif
            icnda(2,icon) = 0
            icnda(8,icon) = ielm
            icnda(9,icon) = ifoun
            if (ififo(13,ifoun).eq.0) ififo(13,ifoun) = 2
            if (ififo(13,ifoun).eq.1) ififo(13,ifoun) = 3
            if (ififo(12,ifoun).eq.0) then
                ififo(12,ifoun) = ireal
            endif
            if (ififo(14,ifoun).le.0) then
                ififo(14,ifoun) = ityp
            endif

            if (ipara.eq.1) goto 700

            nnode = 2
            nlist(1) = nodes(1)
            nlist(2) = nodes(2)

            do inode = 1,nnode
               if (nodes(inode).gt.0) then
                   n = nodfwd(nodes(inode))
                   nodes(inode) = (n-1)*NumDof + 1
               endif
               icnda(3+inode,icon) = nodes(inode)
            enddo

c --- find attached beam element

            nElem = ElemAttachedSet(nnode,nlist(1),Elptch(1))

            do 150 n = 1,nElem
               jelm = Elptch(n)
               if (jelm.eq.ielm) goto 150
               i = elmget(jelm,elData(1),nodes(1))
               if (etyiqr(elData(EL_TYPE),-JBEAM).eq.0) goto 150
               ibeam = abs(jelm)
               icnda(2,icon) = ibeam
               ival(2) = icon
               call conput (ielm,ival(1))

c --- put stiffness into element data base
               matsld = elData(EL_MAT)
c --- retrieve Young's modulus and the volume from attached element
               call GetEx (ibeam,ebeam,matsld,mattab)

               fsmat = 2.d0
               if (mattab.eq.1) fsmat = fsmat*0.01d0
               if (mattab.eq.2) fsmat = fsmat/2.0d0
               fsmat = fsmat/fsglb
               ebeam = abs(ebeam)*fsmat
               call v_zero (xyzc(1),6)
               nrvr = elmgct(ielm,xyzc(1))
               xyzc(5) = ebeam
               call elmpct (ielm,xyzc(1))
               goto 700
 150        continue
         endif
 700  continue

      call MemFree (piElptchL)

      do 710 kfoun = 1,nFoun
         ififo(24,kfoun) = kfoun
 710  continue

c --- reorder contact pair
c
      iott = GetWriteUnit()
      write (iott,3000)
 3000 format (/'--> CONTACT ELEMENT REODERING...')

c --- check if reodering is needed
      Reod_needed = .False.

      do 720 i = 1,nFoun-1
         if (ififo(12,i).gt.ififo(12,i+1))then
             Reod_needed = .True.
             goto 730
         endif
 720  continue

 730  if (.not.Reod_needed) then
          write (iott,3001)
 3001     format ('--> REODERING NOT NEEDED'/)
          goto 999
      endif

c --- the code below is very akward for optimizer
c     probably because of recursive column addressing
c
      call KCL_time(time_b,time_wb)
      do 750 kfoun = i,nFoun-1
         kmin = kfoun
         do 740 ifoun = kfoun+1,nFoun
            if (ififo(12,ififo(24,ifoun)).lt.ififo(12,ififo(24,kmin)))
     &      then
                kmin = ifoun
            endif
 740     continue
         if (kmin.ne.kfoun) then
             iwork = ififo(24,kmin)
             ififo(24,kmin) = ififo(24,kfoun)
             ififo(24,kfoun) = iwork
         endif
 750  continue
      call KCL_time(time_e,time_we)
c
      write (iott,3002) time_we - time_wb
 3002 format ('--> TIME ELAPSED FOR REODERING  ',g10.3,' SECONDS'/)


  999 continue

#if defined(DEBUG)
      call SubEnd('Cont_1')
#endif

      return
      end
