      subroutine EBStabFormA (elId, nUsrDof, kelout, IncU,
     &                        elEnergy, StabFact,
     &                        cDamp, refEnergyEl, stbEnergyEl,
     &                        eStiff, StbMat, fInt, add)
c --- primary function: form stabilization matrix and load vector
c
c  input arguments:
c     elId        (int,sc,in)           - element number
c     nUsrDof     (int,sc,in)           - matrix and load vector size
c     kelout      (int,ar(10),in)       - matrix and load vector formed
c     IncU        (int,ar(*),in)        - incremental displacements
c     StbMat      (dp,ar(nr,nr),in)     - mass matrix with the uniform density
c     add         (int,sc,in)           - add or not the stabilization
c
c  output arguments:
c     eStiff      (dp,ar(nr,nr),inout)  - stiffness matrix
c     fInt        (dp,ar(nr),out)       - n-r restoring f vector
c     elEnergy    (dp,ar(10),out)       - element energies
c     cDamp       (dp,sc,inout)         - elemental articifical damping
c     refEnergyEl (dp,ar(nr,nr),inout)  - elemental energy used for calculating damping ratio
c     stbEnergyEl (dp,ar(nr,nr),inout)  - energy dissipation due to stabilization
c
#include "impcom.inc"
#include "stepcm.inc"
#include "elecom.inc"
c
      double precision v_dot,v_sumi,v_sum

      integer          elId, nUsrDof, kelout(10), add

      double precision IncU(nUsrDof), StbMat(nUsrDof,nUsrDof), StabFact,
     &                 elEnergy(10), cDamp, refEnergyEl, stbEnergyEl,
     &                 eStiff(nUsrDof,nUsrDof), fInt(nUsrDof), fStb(12)

      integer          nUsrDof2, i, j, nUsrDof1

      double precision refEnergyEl0, stbEnergyEl0, refEnergyInc,
     &                 extraCoeff, deltTinv, Vel(12), ratio, xw, yw,
     &                 dperr(3), stbEnergyPred, cDampAct, cMode

c --- initialization
      nUsrDof2 = nUsrDof*nUsrDof
      refEnergyEl0 = refEnergyEl
      stbEnergyEl0 = stbEnergyEl
      refEnergyEl = elEnergy(1)
      refEnergyInc = refEnergyEl - refEnergyEl0
      deltTinv = 1.0d0/timinc
      extraCoeff = (timend-timval+timinc)*deltTinv

      call cv1_v2 (IncU(1), Vel(1), nUsrDof, deltTinv)
      
c --- assume stable solution
      if (add.eq.0 .and. kuphvr.eq.1) then
          if (StabFact.lt.0.0d0) then
              cDamp = abs(StabFact)
          else
              ratio = StabFact
c --- calculate the incremental energy rate
              call Axv1_v2(StbMat(1,1),Vel(1),fStb(1),nUsrDof,nUsrDof)
              xw = v_dot (fStb(1),IncU(1),nUsrDof)
c --- check if energy grows
              if (xw.gt.0.0d0) then
c --- energy predictor (the stabilization energy should be a fraction "ratio" of the total energy)
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
c --- assume unstable first step
      elseif (add.eq.1 .and. isubst.eq.1 .and. ieqitr.eq.1 
     &        .and. kuphvr.ne.1) then
          if (StabFact.lt.0.0d0) then
              cDamp = abs(StabFact)
          elseif (cDamp.le.0.0d0) then
              ratio = StabFact
              cMode = 0.9d0
              nUsrDof1 = nUsrDof - 1
              xw = cMode*cMode*StbMat(1,1) + cMode
     &            * (v_sum (StbMat(2,1), nUsrDof1)
     &            +  v_sumi(StbMat(1,2), nUsrDof, nUsrDof1))
              yw = cMode*cMode*eStiff(1,1) + cMode
     &            * (v_sum (eStiff(2,1), nUsrDof1)
     &            +  v_sumi(eStiff(1,2), nUsrDof, nUsrDof1))
              do i = 2, nUsrDof
                 do j = 2, nUsrDof
                    xw = xw + StbMat(i,j)
                    yw = yw + eStiff(i,j)
                 enddo
              enddo
             if (xw.le.0.0d0) then
                 dperr(1) = elId
                 dperr(2) = xw
                 call erhandler ('EBStabForm', 5000, 3, 'Element %I
     &                            got negative stabilization energy %I',
     &                            dperr(1),' ')
             endif
             yw = abs(yw)
             cDamp = ratio*timinc*yw/xw
          endif
          xw = cDamp*deltTinv
          call cv1_v1 (StbMat(1,1), nUsrDof*nUsrDof, xw)
          call Axv1_v2(StbMat(1,1), IncU(1), fStb(1), nUsrDof, nUsrDof)
          xw = v_dot(fStb(1), IncU(1), nUsrDof)
          stbEnergyEl = stbEnergyEl0
          if (kelout(1).eq.1) then
              call v1pv2_v1 (eStiff(1,1), StbMat(1,1), nUsrDof2)
          endif
          if (kelout(6).eq.1) then
              call v1pv2_v1 (fInt(1), fStb(1), nUsrDof)
          endif
          cDampAct = cDamp
c --- substeps other than the first
      elseif (add.eq.1) then
          cDampAct = cDamp
c          cDampAct = cDamp*lfactn**3          
          xw = cDampAct*deltTinv
          call cv1_v1 (StbMat(1,1), nUsrDof*nUsrDof, xw)
          call Axv1_v2(StbMat(1,1), IncU(1), fStb(1), nUsrDof, nUsrDof)
          xw = v_dot(fStb(1), IncU(1), nUsrDof)
          stbEnergyEl = stbEnergyEl0 + xw
          if (kelout(1).eq.1) then
              call v1pv2_v1 (eStiff(1,1), StbMat(1,1), nUsrDof2)
          endif
          if (kelout(6).eq.1) then
              call v1pv2_v1 (fInt(1), fStb(1), nUsrDof)
          endif
      endif
      if (stbEnergyEl.lt.stbEnergyEl0) then
          dperr(1) = elId
          dperr(2) = stbEnergyEl - stbEnergyEl0
          call erhandler ('EBStabForm', 5000, 3, 'Element %I
     &                    got negative stabilization energy increment
     &                    %I ',dperr(1),' ')
      endif

 999  continue

      return
      end