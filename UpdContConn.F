      subroutine UpdContConn (NumEl,disp,SVRIndex,pdSaveFnr,
     &                        ififo,icnda,isegda,founda,nprocs,
     &                        FirstRound)
c --- primary function: update connectivity of contact elems
c
c  input arguments:
c     NumEl        (int,sc,in)             - number of elements (D_NumEl)
c     disp         (dp,ar(D_ULen,*),in)    - global disp vectors
c     SVRIndex     (int,ar(D_NumEl,2),out) - SVR file index vector
c     FirstRound   (log,sc,in)             - indicator of elimination round
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "cont_info.inc"
#include "in_mem.inc"
#include "conv_stat.inc"
#include "dof_info.inc"
#include "SVR_file.inc"
#include "err_info.inc"
c
      external UpdElemConn
      integer elmget, GetWriteUnit, ElemAttachedSet, BeamAttachedSet, 
     &        etyget,vi_compress, etyiqr, rlinqr
      PTRFTN  MemAlloc, pdSaveFnr
      integer NumEl,elem,SVRIndex(NumEl,2),MaxChatter,ContTot,
     &        ityp,nnod,kupsvr_save,j,jelm,i,ibeam,k,ndel,
     &        ififo(K_nPerFounI,K_nFoun),icnda(K_nPerMast,K_nMast),
     &        isegda(K_nPerSeg,K_nSeg),nodes(W_NdMax),elData(EL_DIM),
     &        elData2(EL_DIM),elChar(IELCSZ),keyopt_1,n1,etyp2,
     &        ival(3),ifoun,iott,nElem,itypold,ncdele,etyp,
     &        idf,iconn,nconn,maxconn,iseg,abeam,nodes2(W_NdMax),
     &        maxi,icon,NumR,ireal,areal,nElemS,nldele,InqSMP,nprocs
      double precision MaxForceNR, disp(*), MaxContForc, con, 
     &                 MaxChatterD
      pointer (prtBeamStat,BeamStat)
      pointer (prtBeamConn,BeamConn)
      pointer (piElptch,Elptch)
      pointer (piInterE,InterE)
      pointer (pdFnr,Fnr)
      pointer (pbcBeams,cBeams)
      integer  Elptch(*),BeamStat(*),InterE(*),BeamConn(*),cBeams(*)
      double precision Fnr(D_ULen), founda(K_nPerFounDp,K_nFoun)
      logical  IfEqvDP, CoreExists, IfMinTime, FirstRound, parallelOff, 
     &         IsLinks
      integer*8 SVRfp_saveL
      double precision ZERO
      parameter (ZERO=0.0d0)
      
      
#if defined(DEBUG)
      call SubBeg('UpdContConn')
#endif
      
      if (K_nFoun.le.0) goto 999
      iott = GetWriteUnit()      
      
      if (A_KeySvrUpd .eq.1 .and. A_InLineSrch.eq.0) then
          if (C_nFraBo.gt.0) then
              write (iott,1000) C_nFraBo
 1000         format('    -> Number of fractured joints so far ',i6)
          endif
          if (C_nDeBo.gt.0) then
              write (iott,1001) C_nDeBo
 1001         format('    -> Number of delaminated joints at '
     &                      'current iteration=',i6)
          endif
      endif
      ContTot = C_nDeBo
      C_nDeBo = 0
      C_nFraBo = 0
      K_ContStat = 0
      
c --- only before first interation)
      if (A_KeyFirstIter.eq.0) goto 300
c --- this must be the test case
      if (NumEl.le.4) goto 300
c --- do not delete anything unless it is a first pass            
      if (A_SpecVer.eq.1 .and. ContTot.eq.0) then
c          if (.not. FirstRound) goto 300
      elseif (A_SpecVer.eq.2) then
          goto 300
      elseif (E_Key_Prev.eq.1) then
          E_Key_Prev = 0
          goto 300
      else
          if (A_BisecNum.ne.0) goto 300
      endif
      
      if (ContTot.lt.0) goto 300
      
      piElptch = PTRFTNNULL
      prtBeamStat = PTRFTNNULL
      prtBeamConn = PTRFTNNULL

c --- delete all inactive contact pairs
      ncdele = 0
      nldele = 0
      itypold = 0
      keyopt_1 = 0

c          i = 0
c          j = 0
c          do ifoun = 1,K_nFoun
c             if (ififo(22,ifoun).gt.i) then
c                 i = ififo(22,ifoun)
c                 j = ifoun 
c             endif
c          enddo
c      write (6,*) 'OGO', i,j


      if (A_BisecNum.ne.0) goto 111
c      write (8,*) 'Enter Upd 1'
    
      do i = 1,K_nMast+K_nSeg
         elem = M_ElSet(i)
         nnod = elmget(elem,elData(1),nodes(1))
         if (elData(EL_DEAD).eq.1) cycle
         call conget (elem,ival(1))
c --- get the foundation number
         if (i.le.K_nMast) then
             ifoun = abs(ival(3))
         else
             ifoun = abs(ival(1))
         endif 
c --- there is a contact
         if (ififo(26,ifoun).eq.0) cycle
c --- there is a near contact
c         if (ififo(26,ifoun).eq.-1) cycle
c --- selection depending on contact type
         ityp = elData(EL_TYPE)
         if (ityp .ne. itypold) then
             itypold = ityp
             keyopt_1 = etyiqr(ityp,-KYOP1)
         endif
c --- skip if near contact for sliding case
         if (keyopt_1.le.2 .and. ififo(26,ifoun).eq.-1) cycle
c --- skip the rigid contact
         if (ififo(15,ifoun).le.0) cycle
c --- mark as killed
         elData(EL_DEAD) = 1
         call elmput(elem,elData(1),nnod,nodes(1))

         if (i.le.K_nMast) then
c --- zero out all primary contact data base entries
c --- removal of elements will be done later
             ncdele = ncdele + 1
             ififo(1,ifoun) = 0
             ififo(13,ifoun) = 0
             ififo(22,ifoun) = 0
c --- fix for the slave
             icon = abs(ival(2))
             if (icon.gt.0) then
                 iseg = icnda(3,icon)
c --- no element - no contact             
                 icnda(3,icon) = 0
c --- since the entire contact is removed, it is OK to zero out without checking
                 if (iseg.gt.0) isegda(19,iseg) = 0
             endif    
          endif
      enddo

c      write (6,*) 'Enter Upd 2', ncdele
      
c --- if this is not the first ipass, no elements removed and everything is OK
      if (A_KeyFirstPass.eq.0 .and. C_ContChang+ncdele+A_BisecNum.eq.0)
     &    goto 300
      
c      write (8,*) 'Enter Upd 3'
      
c --- set up bit map for patch of elements
 111  piElptch = MemAlloc(2*NumEl,W_LenInt,'piElptch        ')

      IfMinTime = IfEqvDP(A_TimeInc,A_MinTimeStep) .or. 
     &                   (A_TimeInc.lt.A_MinTimeStep)

c --- maximum real number defined
      NumR = rlinqr(0,DB_MAXDEFINED)
c --- reserve space for the connectivity array according to real number
c     beams will be distinguished with real ID
      prtBeamStat = MemAlloc(NumR,W_LenInt,'prtBeamStat     ')
      pbcBeams = MemAlloc(NumEl,W_LenInt,'pbcBeams        ')
      
      call vi_zero (BeamStat(1),NumR)
      call vi_zero (cBeams(1),NumEl)

c --- only beams will be considered, mark the rest with -1
      do i = 1,K_nMast+K_nSeg
         elem = M_ElSet(i)
         nnod = elmget(elem,elData(1),nodes(1))         
c --- get real number
         ireal = elData(EL_REAL)
c --- mark as not available
         if (ireal.gt.0) BeamStat(ireal) = -1
      enddo
c      write (8,*) 'Enter Upd 8'
c --- fill in "if attached" array
c --- go through all contact elements


      do i = 1,K_nMast
         elem = M_ElSet(i)
         nnod = elmget(elem,elData(1),nodes(1))
         if (elData(EL_DEAD).eq.1) cycle
         call conget (elem,ival(1))
c --- contact number
         icon = abs(ival(2))
c --- foundation number
         ifoun = abs(ival(3))
c --- attached beam
         ibeam = icnda(2,icon)
c --- if no beam attached - cycle
         if (ibeam.eq.0 .or. ifoun.eq.0) cycle
         if (ififo(1,ifoun).le.0) cycle
c --- if not in contact
         if (ififo(26,ifoun).ne.0) cycle
c --- retrive segment number
         iseg = icnda(3,icon)
c --- this should also take care of unconnected surface-beam element
         if (iseg.eq.0) cycle
c --- beam in contact
         abeam = abs(isegda(16,iseg))
c --- get real number for ibeam (CAUTION, elData will be changed)
         nnod = elmget(ibeam,elData(1),nodes(1))
         ireal = elData(EL_REAL)

c --- if beam to beam contact
         if (abeam.gt.0) then
c --- mark that beams are connected
             cBeams(ibeam) = -1
             cBeams(abeam) = -1
c --- get real number for abeam (CAUTION, elData will be changed again)
             nnod = elmget(abeam,elData(1),nodes(1))
             areal = elData(EL_REAL)
c --- this is important, all beams with only one contact are marked for deletion here
             BeamStat(ireal) = BeamStat(ireal) + 1
             BeamStat(areal) = BeamStat(areal) + 1
         elseif (ififo(15,ifoun).eq.-1.or.ififo(15,ifoun).le.-3) then
c --- if beam to surf contact
c --- select if fixed (surface is always fixed)
             BeamStat(ireal) = BeamStat(ireal) + 1
         endif
      enddo
c      write (8,*) 'Enter Upd 9'      
c --- calculate number of connected fibers
c --- find maximum number of connections per fiber
c --- NOTE, the function of BeamStat will be changed here
      nconn = 0
      iconn = 1
      do ireal = 1,NumR
c --- if there is no contact found
         if (BeamStat(ireal).le.0) cycle
c --- increment total number of contact
         nconn = nconn + 1
c --- maximum number of contacts found for this real ID
         maxconn = BeamStat(ireal)
c --- pointer to the array
         BeamStat(ireal) = iconn
c --- increment pointer (beam number + number of contacts + beams in contact)
         iconn = iconn + 2 + maxconn
      enddo
c --- length of the new array
      maxconn = iconn - 1
      CoreExists = .true.
      if (maxconn.le.0 .or. nconn.eq.0) then
c --- if there no contacts at all
          CoreExists = .false.
          goto 444
      endif
      
c      write (8,*) 'Enter Upd 10'     
c --- we will form a new advanced connectivity array
      prtBeamConn = MemAlloc(maxconn,W_LenInt,'ptrBeamConn     ')
c --- zero out all entries
      call vi_zero (BeamConn(1),maxconn)

c --- fill in all real numbers to the first place
      do i = 1,NumR
c --- any contacts?
         j = BeamStat(i)
c --- yes, then make the stamp
         if (j.gt.0) BeamConn(j) = i
      enddo

c      write (8,*) 'Enter Upd 11'     

      call FillConn (BeamStat(1),NumR,BeamConn(1),maxconn,
     &               ififo(1,1),icnda(1,1),isegda(1,1))

c      write (8,*) '2. BeamStat', BeamStat (1:NumR)

c      write (8,*) 'Enter Upd 12'     
      
      call CoreSelect (BeamStat(1),BeamConn(1),cBeams(1),
     &                 NumR,NumEl,nconn,maxconn)

      
c      write (9,*) '3. BeamStat', BeamStat (1:NumR)
            
  444 continue  
c      write (8,*) 'Enter Upd 13'     
      
c --- kill free fibers (only if one beam per fiber)
      ndel = 0
c      write (6,*) '3. BeamStat final', BeamStat (1:NumR)
      if (FirstRound .and. A_SpecVer.eq.3) then
c --- remove tangling ends
          elem = 1
          itypold = 0      
          do ireal = 1,NumR
             if (BeamStat(ireal).le.0) cycle
             do while (elem.le.NumEl .and. abs(cBeams(elem)).ne.ireal)
                elem = elem + 1
             enddo
c --- go forward until the first contact
             do while (elem.le.NumEl .and. cBeams(elem).eq.ireal)
                cBeams(elem) = 0
                elem = elem + 1
             enddo
c --- go forward until the end
             do while (elem.le.NumEl .and. abs(cBeams(elem)).eq.ireal)
                elem = elem + 1
             enddo
c --- go backwards untill first contact
             j = elem - 1
             do while (j.ge.1 .and. cBeams(j).eq.ireal)
                cBeams(j) = 0
                j = j - 1
             enddo
          enddo
      endif
      
      IsLinks = .false.
      do 200 elem = 1,NumEl
c --- again, the entry can be negative if the element touches the surface 
c     in this case, the element was already selected in FillConn.F
         nnod = elmget(elem,elData(1),nodes(1))
         if (nnod.le.0) cycle
         if (elData(EL_DEAD).eq.1) cycle
         ityp = elData(EL_TYPE)
         if (ityp .ne. itypold) then
             itypold = ityp
             ityp = etyget(ityp,elChar(1))
             etyp = elChar(EL_TYPE)
         endif
         if (etyp.eq.W_eLink) IsLinks = .true.
c --- only beams are processed
         if (elChar(JBEAM).eq.0) cycle
         ireal = elData(EL_REAL)
c --- bypass the check in the minimum step 
         if (abs(ireal).eq.A_DelRealBeam) goto 123
         
         if (ireal.le.0) cycle
         i = BeamStat(ireal)
c --- make a choice about deletion
         if (i.lt.0) then
c --- definitely not a beam
             cycle
         elseif (FirstRound .and. cBeams(elem).eq.0) then
c --- a tangling end
         elseif (i.gt.0) then
c --- negative means connected to the fixed boundary
c --- treshold value for number of contacts CANNOT be varied via BeamConn(i+1) here
c --- it should be changed in CoreSelect.F
             if (BeamConn(i).eq.0 .or. (.not.CoreExists)) cycle
         else
c --- NOTE: no contacts at all means i.eq.0 and it should go through too         
         endif
 123     continue         
c --- form the set of elements attached to the beam (including links)
         nElemS  = ElemAttachedSet(nnod,nodes(1),Elptch(1))
         
c --- proceed only with unselected fibres disconnected from the boundaries
c --- kill all attached elements
         
         do 210 k = 1,nElemS
            jelm = Elptch(k)
            j = elmget(jelm,elData2(1),nodes2(1))

            if (elData2(EL_DEAD).eq.1) cycle
            etyp2 = etyiqr(elData2(EL_TYPE),-JETYP)

            if (etyp2.eq.W_eMaster) then
                call conget (jelm,ival(1))
                ifoun = abs(ival(3))
c --- do not touch rigid surface
                if (ififo(15,ifoun).gt.0) then
                    ififo(1,ifoun) = 0
                    ififo(13,ifoun) = 0
                    ififo(22,ifoun) = 0
                endif
                ncdele = ncdele + 1
            elseif (etyp2.eq.W_eSlave) then
                call conget (jelm,ival(1))
                ifoun = abs(ival(1))
                ififo(1,ifoun) = 0
                ififo(13,ifoun) = 0
                ififo(22,ifoun) = 0
                ncdele = ncdele + 1
            elseif (etyp2.eq.W_eLink) then
                nldele = nldele + 1
            endif
            elData2(EL_DEAD) = 1
            call elmput (jelm,elData2(1),j,nodes2(1))
  210    continue
         ndel = ndel + 1
  200 continue
  
      call MemFree(prtBeamConn)
c      write (6,*) 'Enter Upd 14', ncdele,nldele,ndel

      if (ndel.gt.0) then
c --- reform the global DOF vector and trigger matrix reordering
          K_ContStat = 1         
c --- sweep across database and make sure all pairs are killed
          do 230 i = 1,K_nMast+K_nSeg
             elem = M_ElSet(i)
             nnod = elmget(elem,elData(1),nodes(1))
             if (elData(EL_DEAD).eq.1) cycle
             call conget (elem,ival(1))
             icon = 0
             iseg = 0
             if (i.le.K_nMast) then
                 ifoun = abs(ival(3))
c --- fix for the slave
                 icon = abs(ival(2))
                 iseg = icnda(3,icon)
             else
                 ifoun = abs(ival(1))
c --- fix for the master
                 iseg = abs(ival(2))
                 icon = isegda(19,iseg)
             endif
c --- no element - no contact             
             if (ififo(1,ifoun).gt.0) cycle
             if (icon.gt.0) icnda(3,icon) = 0
             if (iseg.gt.0) isegda(19,iseg) = 0
             elData(EL_DEAD) = 1
             call elmput (elem,elData(1),nnod,nodes(1))
  230     continue


c --- sweep across database and make sure all unconnected links are deleted
          if (IsLinks) then
              itypold = 0      
              do 240 elem = 1,NumEl
                 nnod = elmget(elem,elData(1),nodes(1))
                 if (nnod.le.0) cycle
                 if (elData(EL_DEAD).eq.1) cycle
                 ityp = elData(EL_TYPE)
                 if (ityp .ne. itypold) then
                     itypold = ityp
                     ityp = etyget(ityp,elChar(1))
                     etyp = elChar(EL_TYPE)
                 endif
c --- only links are processed
                 if (etyp.eq.W_eLink) then
                     do 241 i = 1,nnod
                        nElemS  = ElemAttachedSet(1,nodes(i),Elptch(1))
                        do k = 1,nElemS
                           jelm = Elptch(k)
                           j = elmget(jelm,elData2(1),nodes2(1))
c --- probe the beams only
                           if (etyiqr(elData2(EL_TYPE),-JBEAM).eq.0) 
     &                         cycle
c --- is there any alive elements connected? go to the next node
                           if (elData2(EL_DEAD).eq.0) goto 241
                        enddo
c --- no alive element were detected, we need to delete the link
                        goto 242
  241                continue
c --- we came out clean, both ends are connected to some alive elements
                     goto 240 ! cycle
                 
  242                do k = 1,nElemS
                        jelm = Elptch(k)
                        j = elmget(jelm,elData2(1),nodes2(1))
                        if (elData2(EL_DEAD).eq.1) cycle
c --- process links only
                        if (etyiqr(elData2(EL_TYPE),-JETYP).ne.W_eLink) 
     &                      cycle
c --- is there any alive elements connected? go to the next node
                        nldele = nldele + 1
                        elData2(EL_DEAD) = 1
                        call elmput (jelm,elData2(1),j,nodes2(1))
                    enddo
                 endif
                       
  240         continue
          endif
  
          write (iott,2000) ndel
 2000     format('   -> Number of freely suspended beams eliminated',i6)
      endif
c      write (8,*) 'Enter Upd 15'     
      
      call MemFree(piElptch)
      call MemFree(prtBeamStat)
      call MemFree(pbcBeams)

      if (nldele.gt.0) then
          write (iott,2001) nldele
 2001     format('     > Database was reduced by',i9,' link(s)')
      endif

      if (ncdele.gt.0) then
          write (iott,2003) ncdele
 2003     format('     > Database was reduced by',i9,' contact pair(s)')
      endif      
      
  300 continue

c --- turn off update of SVRs
      if (.not. FirstRound) then
          do ifoun = 1,K_nFoun
c --- zero out the contact status and forces
             founda(34,ifoun) = ZERO !!!
             ififo(26,ifoun) = -2
          enddo
c --- line search does not require the matrix update, therefore we skip the connectivity update
          if (A_InLineSrch.eq.0) then
              kupsvr_save = A_KeySvrUpd
              A_KeySvrUpd = 0
              SVRfp_saveL = S_FilePos
c    _______________  this part is to be run in parallel  ______________
              parallelOff = 0
              if (numel.lt.InqSMP(2) .or. nprocs.eq.1) parallelOff = 1
              call InfoSMP (1,parallelOff)
              if (parallelOff.eq.0) then
                  call ResumeSMP
                  call ForkSMP4(UpdElemConn,NumEl,disp(1),SVRIndex(1,1),
     &                          ififo(1,1))
                  call SuspendSMP
              else
                  call UpdElemConn(NumEl,disp(1),SVRIndex(1,1),
     &                             ififo(1,1))
              endif
c    ___________________  end of parallel portion  _____________________           
              S_FilePos = SVRfp_saveL
c --- set update of SVRs back to what it was
              A_KeySvrUpd = kupsvr_save
          endif
      endif

      
 999  continue
    
c      write (8,*) 'Finnish Upd' 
      
#if defined(DEBUG)
      call SubEnd('UpdContConn')
#endif

      return
      end