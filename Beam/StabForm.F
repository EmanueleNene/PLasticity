      subroutine StabForm  (elId, nUsrDof, kelout, IncU, StbMat, 
     &                      elEnergy, cDamp, refEnergyEl, stbEnergyEl,
     &                      eStiff, fInt, fStb)
c --- primary function: form stabilization matrix and load vector
c
c  input arguments:
c     elId     (int,sc)             - element number
c     nUsrDof  (int,sc,in)          - matrix and LV size
c     kelout   (int,ar(10),in)      - matrix and load vector formed
c     IncU     (int,ar(*))          - incremental displacements
c     StbMat   (dp,ar(nr,nr))       - mass matrix with the uniform density
c
c  output arguments:
c     eStiff   (dp,ar(nr,nr),inout) - k matrix
c     fInt     (dp,ar(nr),out)      - n-r restoring f vector
c     elEnergy (dp,ar(10),out)      - element energies
c     cDamp    (dp,sc)              - elemental articifical damping
c     StbMat   (dp,ar(nr,nr))       - stabilization matrix
c     fStb     (dp,ar(nr),out)      - force due to stabilization
c     refEnergyEl (dp,ar(nr,nr))    - elemental energy used for calculating damping ratio
c     stbEnergyEl (dp,ar(nr,nr))    - energy dissipation due to stabilization
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
c
      double precision v_dot, v_sum, v_sumi

      integer          elId, nUsrDof, kelout(10)

      double precision IncU(*), StbMat(nUsrDof,*), elEnergy(10),
     &                 cDamp, refEnergyEl, stbEnergyEl,
     &                 eStiff(nUsrDof,*), fInt(*), fStb(*)



      integer          nUsrDof2, nUsrDof1

      double precision refEnergyEl0, stbEnergyEl0, refEnergyInc,
     &                 extraCoeff, deltTinv, Vel(nUsrDof), ratio, xw,
     &                 dperr(3), stbEnergyPred, cMode, yw, cDampAct


      integer          i, j

#if defined(DEBUG)
      call SubBeg('StabForm')
#endif

c --- initialization
      nUsrDof2 = nUsrDof*nUsrDof
      refEnergyEl0 = refEnergyEl
      stbEnergyEl0 = stbEnergyEl
      refEnergyEl = elEnergy(1)
      refEnergyInc = refEnergyEl - refEnergyEl0
      extraCoeff = (A_TimeEnd-A_TimeCur+A_TimeInc)/A_TimeInc
      deltTinv = 1.0d0/A_TimeInc
      call cv1_v2 (IncU(1), Vel(1), nUsrDof, deltTinv)

c --- handle constant stabilization factor for the time being
c --- assume the 1st substep is stable
      if (A_StabFirstStep.gt.0 .and. A_StabFirstStep.lt.3
     &                         .and. A_KeyMatUpd.eq.1) then
c          if (elId.eq.1) write (6,*) 'Branch 1', A_StabFirstStep
          if (A_StabFact.lt.0.0d0) then
              cDamp = abs(A_StabFact)
          else
              ratio = A_StabFact
              call Axv1_v2(StbMat(1,1),Vel(1),fStb(1),nUsrDof,nUsrDof)
              xw = v_dot (fStb(1),IncU(1),nUsrDof)
              if (xw.gt.0.0d0) then
                  stbEnergyPred = ratio*(refEnergyEl0
     &                          + extraCoeff*refEnergyInc)
                  if (stbEnergyPred.gt.stbEnergyEl0) then
                      cDamp = (stbEnergyPred-stbEnergyEl0)
     &                       /(extraCoeff*xw)
                  else
                      stbEnergyPred = ratio*refEnergyEl
                      if (stbEnergyPred.gt.stbEnergyEl0) then
                          cDamp = (stbEnergyPred-stbEnergyEl0)
     &                           /(extraCoeff*xw)
                      else
                         if (cDamp.lt.0.0d0) then
                             cDamp = 0.0d0
                         endif
                      endif
                  endif
              else
                  cDamp = 0.0d0
              endif
          endif
          stbEnergyEl = stbEnergyEl0
          call v_zero (fStb(1), nUsrDof)
          call v_zero (StbMat(1,1), nUsrDof2)
          cDampAct = cDamp

c --- the unstable first substep
      elseif (A_StabFirstStep.eq.3 .and. A_IterNum.eq.1
     &                             .and. A_KeyMatUpd.ne.1) then
c          if (elId.eq.1) write (6,*) 'Branch 2', A_StabFirstStep
          if (A_StabFact.lt.0.0d0) then
              cDamp = abs(A_StabFact)
          elseif (cDamp.le.0.0d0) then
              ratio = A_StabFact
              cMode = 0.9d0
              nUsrDof1 = nUsrDof - 1
              xw = cMode*cMode*StbMat(1,1) + cMode
     &            * (v_sum(StbMat(2,1), nUsrDof1)
     &            + v_sumi(StbMat(1,2), nUsrDof, nUsrDof1))
              yw = cMode*cMode*eStiff(1,1) + cMode
     &            * (v_sum(eStiff(2,1), nUsrDof1)
     &            + v_sumi(eStiff(1,2), nUsrDof, nUsrDof1))
              do i = 2, nUsrDof
                 do j = 2, nUsrDof
                    xw = xw + StbMat(i,j)
                    yw = yw + eStiff(i,j)
                 enddo
              enddo
             if (xw.le.0.0d0) then
                 dperr(1) = elId
                 dperr(2) = xw
                 call erhandler ('StabForm', 5000, 3, 'Element %I
     &                            got negative stabilization energy %I',
     &                            dperr(1),' ')
             endif
             yw = abs(yw)
             cDamp = ratio*A_TimeInc*yw/xw
          endif
          xw = cDamp/A_TimeInc
          call cv1_v1 (StbMat(1,1), nUsrDof*nUsrDof, xw)
          call Axv1_v2(StbMat(1,1), IncU(1), fStb(1), nUsrDof, nUsrDof)
          xw = v_dot(fStb(1), IncU(1), nUsrDof)
          stbEnergyEl = stbEnergyEl0
          if (kelout(1).eq.1) then
              call v1pv2_v1 (eStiff(1,1), StbMat(1,1), nUsrDof2)
          else
              call v_zero (StbMat(1,1), nUsrDof2)
          endif
          if (kelout(6).eq.1) then
              call v1pv2_v1 (fInt(1), fStb(1), nUsrDof)
          else
              call v_zero (fStb(1), nUsrDof)
          endif
          cDampAct = cDamp
      elseif (A_StabFirstStep.eq.3 .and. A_KeyMatUpd.eq.1) then
c          if (elId.eq.1) write (6,*) 'Branch 3', A_StabFirstStep
          call Axv1_v2(StbMat(1,1), Vel(1), fStb(1), nUsrDof, nUsrDof)
          xw = v_dot (fStb(1), IncU(1), nUsrDof)
          stbEnergyEl = stbEnergyEl0 + cDamp*xw
          call cv1_v1 (fStb(1),  nUsrDof, cDamp)
          call cv1_v1 (StbMat(1,1), nUsrDof*nUsrDof,cDamp)
          if (A_StabFact.gt.0.0d0) then
             ratio = A_StabFact
             stbEnergyPred = ratio*(refEnergyEl0
     &                     + extraCoeff*refEnergyInc)
             yw = stbEnergyEl0 + cDamp*extraCoeff*xw
             if (stbEnergyPred .gt. yw) then
             else
                if (xw.gt.0.0d0) then
                    if (stbEnergyPred.gt.stbEnergyEl0) then
                        cDamp = (stbEnergyPred-stbEnergyEl0)
     &                         /(extraCoeff*xw)
                    else
                        stbEnergyPred = ratio*refEnergyEl
                        if (stbEnergyPred.gt.stbEnergyEl0) then
                            cDamp = (stbEnergyPred-stbEnergyEl0)
     &                             /(extraCoeff*xw)
                        endif
                    endif
                else
                    cDamp = 0.0d0
                endif
             endif
          endif
          if (kelout(1).eq.1) then
              call v1pv2_v1 (eStiff(1,1), StbMat(1,1), nUsrDof2)
          else
              call v_zero (StbMat(1,1), nUsrDof2)
          endif
          if (kelout(6).eq.1) then
              call v1pv2_v1 (fInt(1), fStb(1), nUsrDof)
          else
              call v_zero (fStb(1), nUsrDof)
          endif
          cDampAct = cDamp

c --- substeps other than the first
      else
c          if (elId.eq.1) write (6,*) 'Branch 4', A_StabFirstStep
          if (A_KeyStab.eq.2) then
              cDampAct = cDamp*(1.0d0-A_CurLoadFact)
          else
              cDampAct = cDamp
          endif
          xw = cDampAct/A_TimeInc
          call cv1_v1 (StbMat(1,1), nUsrDof*nUsrDof, xw)
          call Axv1_v2(StbMat(1,1), IncU(1), fStb(1), nUsrDof, nUsrDof)
          xw = v_dot(fStb(1), IncU(1), nUsrDof)
          stbEnergyEl = stbEnergyEl0 + xw
          if (kelout(1).eq.1) then
              call v1pv2_v1 (eStiff(1,1), StbMat(1,1), nUsrDof2)
          else
              call v_zero (StbMat(1,1), nUsrDof2)
          endif
          if (kelout(6).eq.1) then
              call v1pv2_v1 (fInt(1), fStb(1), nUsrDof)
          else
              call v_zero (fStb(1), nUsrDof)
          endif
      endif
      elEnergy(6) = stbEnergyEl
      elEnergy(7) = refEnergyEl
      if (stbEnergyEl.lt.stbEnergyEl0) then
          dperr(1) = elId
          dperr(2) = stbEnergyEl - stbEnergyEl0
          call erhandler ('StabForm', 5000, 3, 'Element %I
     &                    got negative stabilization energy increment
     &                    %I ',dperr(1),' ')
      endif

#if defined(DEBUG)
      call SubEnd('StabForm')
#endif

      return
      end



