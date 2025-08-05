      subroutine BeamElement (elemId,elChar,elData,eomask,eConn,locsvrL,
     &                        kelreq,nUsrDof,xBeg,eleDisp,kelout,eStiff,
     &                        eMass,gStiff,fExt,fInt,elemVolume,
     &                        elemMass,elemCG,eEnergy,edindxL,lcerstL,
     &                        numElemNodes, LenMatRecMax, numWordsPt, 
     &                        nmatMisc,secDataI,RstNodeList,nCrnr,lwork,
     &                        work,twork,nResults,matStrain,matStress,
     &                        secDataR,lenSecDataR,numSectPts,
     &                        secResults,matRecord,pResults,postdata,
     &                        ple,peeq,pwrk,thStrSecPt,tSecPt,thStrSNDP,
     &                        hyStrSecNd,nDofNode,nGenStr,eU,incDisp,
     &                        secStiff,tGenTh,gradMat,tGenStress,
     &                        cre,creq,crwk,mRecIntPidx,nlkword,
     &                        miscdata)
c --- primary function: finite strain beam element
c
c  input arguments:
c     elemId (int,sc,in)             - element numbe
c     elChar(int,ar(IELCSZ),in)      - array of element type characteristics
c     elData (int,ar(EL_DIM),in)     - array of element data
c     eomask (int,sc,in)             - bit pattern of element output
c     eConn  (int,ar(4),in)          - array of element node numbers
c     locsvrL (intL,sc,in)           - location of this eles SVRs on SVR
c     kelreq (int,ar(10),in)         - matrix and load vector form requests
c     nUsrDof   (int,sc,in)          - matrix and lv size
c     xBeg  (dp,ar(6,4),in)          - nodal coords (orig) and rot angles
c     eleDisp(dp,ar(nUsrDof,5),in)   - element nodal solution values
c     numElemNodes(int,sc,in)        - number of element nodes
c
c  output arguments:
c     kelout (int,ar(10),out)        - keys indicating created mats and lvs
c     eStiff (dp,ar(nUsrDof,nUsrDof),inout) - k matrix (kelreq(1))
c     eMass  (dp,ar(nUsrDof,nUsrDof),inout) - m matrix (kelreq(2))
c     gStiff (dp,ar(nUsrDof,nUsrDof),inout) - s matrix (kelreq(4))
c     fExt   (dp,ar(nUsrDof),out)    - applied f vector (kelreq(5))
c     fInt   (dp,ar(nUsrDof),out)    - n-r restoring f vector (kelreq(6))
c     elemVolume  (dp,sc,out)        - element volume
c     elemMass (dp,sc,out)           - element mass
c     elemCG (dp,ar(3),out)          - centroid location
c     eEnergy(dp,ar(10),out)         - element energies
c     edindxL (intL,ar(25),out)      - element result data file indexes
c     lcerstL (intL,sc,inout)        - position on result file
c
c
c  keyopts:
c     keyopt(1) = force going thru MaterialMain even for elastic case
c     keyopt(2) = add warping upon torsion. this option is not yet debugged
c     keyopt(3) = include initial stresses - interface is needed
c
c  keyopts for debugging:
c     keyopt(4) = 1 ignore stress stiffness contributions to torsional mode from axial force
c               = 2 ignore all stress stiffness terms corresponding to axial force
c     keyopt(5) = 1 ignore stress stiffness from torque
c               = 2 ignore stress stiffness from moments
c               = 3 ignore ALL moment stress stiffness terms
c
c     keyopt(6) = 1 ignore stress stiffness terms that couple u and phi DOFs due to shear forces
c               = 2 ignore stress stiffness to all rot DOFs due to shear forces
c               = 3 ignore all transverse shear related stress stiffness terms
c     keyopt(7) = 1 keep the section data incore
c     keyopt(8) = 1 ignore volume-preserving assumption
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "RST_file.inc"
#include "sec_info.inc"
#if defined(MOIST)
#include "in_mem.inc"
#include "cont_info.inc"
#include "beam_hyg.inc"
#endif
c
      integer*8 locsvrL, edindxL(25), lcerstL
      integer
     &     elemId,elChar(IELCSZ),elData(EL_DIM),eomask,eConn(4),
     &     kelreq(10),nUsrDof,kelout(10),nreal,
     &     numElemNodes, LenMatRecMax, lwork,
     &     nDirect, nTens, mRecIntPidx(W_LenMatRecIndx)

      double precision
     &     xBeg(6,4),eleDisp(nUsrDof,8),elemCG(3),eEnergy(10),
     &     eStiff(nUsrDof,nUsrDof),eMass(nUsrDof,nUsrDof),
     &     gStiff(nUsrDof,nUsrDof),fExt(nUsrDof),fInt(nUsrDof),
     &     elemVolume,elemMass
      integer
     &     keyopt(12), matId, ireal,
     &     numBasicSaveVars,stateVarIndex(20)
c
c --- section properties summary array
c
      integer
     &     pArea,pCgy,pCgz,pIy,pIz,pIyy,pIyz,pIzz,pTorq,pIpp,
     &     pSww,pSy,pSz
      parameter (
     &     pArea =1, pIpp  =10,
     &     pCgy  =2,
     &     pCgz  =3,
     &     pIy   =4,
     &     pIz   =5, pSww  =14,
     &     pIyy  =6, pSy   =15,
     &     pIyz  =7, pSz   =16,
     &     pIzz  =8,
     &     pTorq =9)
      double precision sectSummary(20,1), mStrainInc(3),
     &                 eEnerSecPt(10), mStrain(3), mStress(3)

c --- section model data
c
      integer
     &     pNNODI,pSectionType,numIntPts,
     &     pNodes,pElems,pElemConn,lenSecDataR
      integer
     &     offY,offZ,nDofNode
      parameter (
     &     numIntPts    = 4,
     &     pNNODI       = 8,
     &     pSectionType = pNNODI+1,
     &     pNodes       = pSectionType+1,
     &     pElems       = pNodes+1,
     &     pElemConn    = pElems+1)

      integer  secDataI(*)
      double precision secDataR(lenSecDataR,*)

c --- section integration data
c
      integer pdA, pJ, pphiy, pphiz, pY, pZ, pM
      parameter (
     &     pdA  = 1, pJ = pdA+1 , pphiy = pJ+1 , pphiz=pphiy+1,
     &     pY = pphiz+1, pZ = pY+1 , pM = pZ+1)
c
      integer numSectPts, elasticSection, numSectNodes
      double precision secResults(B_NumResPts,numSectPts,1)
      double precision yOffset(1), zOffset(1), scaleSF1,
     &                 scaleSF2, aScale , yzScale(2), yzScale_t(2), GIP
c
      double precision ONE, ZERO, HALF, TWO, TWELVE, THREE
      parameter (ZERO = 0.0d0, HALF = 0.5d0, ONE = 1.0d0, TWO = 2.0d0,
     &           TWELVE = 12.0d0, THREE=3.0d0)

c --- local variables
c
      integer kelin(10),khygro,kpress,nod,i,j,jj,kspt,kk,
     &        numDomainIntPts, domInt, knod,ii,ll,nmisc,khygroIP,
     &        col,row, k1,k2,rlgetn,
     &        RstNodeList(*),inode,ibnode,NumRstNodes, NodeNum,
     &        nCrnr(*),ishape,nGenStr,maxGenStr, numWordsPt,
     &        nmatMisc, domIntP, nodeP, NumElems

      integer*8  nlkword(W_nMatCol,W_nMatRow)
      integer*8  keyplL, keycrL, keynlL
      parameter (maxGenStr = 7)

      double precision refmoist,moistbeg(9),moistcur(9),moistend(9),
     &     moistprev(9),prsBeg(14),pressure(14),prsEnd(14),dL_t0(2),
     &     GetMoistRef,work(lwork),xCur(3,4),
     &     eU(nDofNode*3),incDisp(nDofNode*3),
     &     basicSaveVars(182+20+704),histSaveVars(119+18),trv(3,3),
     &     tGenStrain_t(maxGenStr,2), shpe(3,2),Gauss2(2),
     &     basis(3,3,2),tGenStrain(maxGenStr,2),dxdg(3),tDomInt,tSecInt,
     &     dL,dLinv,dV,dV_t0(2),ncrd(3),secStiff(nGenStr,nGenStr),
     &     nBasis(3,3,3),nBasis_t(3,3,3),shpDer(3,2),eWarp(2),fWarp(2),
     &     y_pt,z_pt,dA_pt,phiy_pt,phiz_pt,matStiff(3,3),xIntp(3,3),
     &     modWt,tDomInt_t,tSecInt_t
      double precision
     &     n3sp(3,3),n2sp(3,3),tt(3,3),dxDx(3,3),n2xn3(3),prop(13),
     &     n3tp(3,3),n2tp(3,3),dpDp1(3,3),n2n3m(3,3),n2tm(3,3),
     &     n3tm(3,3),dpDp2(3,3),tn2(3,3),tn3(3,3),dxDp(3,3),
     &     YoungsMod,density,lambda,lambda2,tInt(4),alpha,
     &     ShearMxy,ShearMyz,ShearMxz,xPolate(4,4),twork(*),
     &     nResults(6,2,*),matStrain(3,2,*),bwarp(21),
     &     matStress(3,2,*),thStrain(2),nBasis_0(3,3,3),
     &     tGenTh(nGenStr),siglv(3),siglvi(3),coords_0(3),
     &     mvnd(4,4),cDensity,basis_t(3,3,2),basis_0(3,3,2),coordAx(3),
     &     gradMat(nGenStr,nDofNode*3),moffset(3,3,3),factor,
     &     tGenStress(nGenStr,2),siglvis(7),iniForce(21),tDer(3,2),
     &     hbasis(3,3,3),hbasis_t(3,3,3),densityDA,tInt_t(4),L,kof
      double precision matRecord(*),pResults(9,3,*),postdata(*),
     &     miscdata(*),tab(3)
c
      double precision Gauss3(3),Gaussw2(2),Gaussw3(3),n1(3),n2(3)
      double precision ple(3,2,*),peeq(3,2,*),rvr(6),
     &     pwrk(3,2,*),thStrSecPt(2,*),tSecPt(*),thStrSNDP(2,*),
     &     hyStrSecNd(3,*), cre(3,2,*),creq(3,2,*), crwk(3,2,*)
      
      
#if defined(MOIST)
      double precision eh,em(H_Row),xn,yn,zn,hight,width
#endif
      
#if defined(DEBUG)        
      double precision
     &     yzmxmn(4,2),semxmn(20)
#endif
c
      equivalence
     &     (basis_t(1,1,1),histSaveVars(1)),
     &     (tGenStrain_t(1,1),histSaveVars(19)),
     &     (nBasis_t(1,1,1),histSaveVars(33)),
     &     (dV_t0(1),histSaveVars(60)),
     &     (dL_t0(1),histSaveVars(62)),
     &     (yzScale_t(1),histSaveVars(64)),
     &     (nBasis_0(1,1,1),histSaveVars(66)),
     &     (hbasis_t(1,1,1),histSaveVars(93)),
     &     (basis_0(1,1,1),histSaveVars(120))
      
      equivalence
     &     (prsBeg(1),basicSaveVars(1)),
     &     (shpDer(1,1),basicSaveVars(15)),
     &     (histSaveVars(1),basicSaveVars(21))
c#if defined(MOIST)
      equivalence
     &     (moistbeg(1),basicSaveVars(158)),
     &     (moistprev(1),basicSaveVars(167)),
     &     (iniForce(1),basicSaveVars(176))      
c#else
c      equivalence
c     &     (iniForce(1),basicSaveVars(158))
c#endif
     
      integer tmp_svr
      logical WetMode, isWarp, iniFlag, compSSTIF

      integer
     &        iNumTerms


      double precision
     &     dEqPl(6), dEqCr(6)

      data xPolate/
     &      1.8660     ,-0.50000    , 0.13397    ,-0.50000    ,
     &     -0.50000    ,  1.8660    ,-0.50000    , 0.13397    ,
     &      0.13397    ,-0.50000    ,  1.8660    ,-0.50000    ,
     &     -0.50000    , 0.13397    ,-0.50000    ,  1.8660    /
      data mvnd/
     &     1.0, 0.0, 0.0, 0.0,
     &     0.0, 1.0, 0.0, 0.0,
     &     0.0, 0.0, 1.0, 0.0,
     &     0.0, 0.0, 0.0, 1.0/

      data Gauss2/-0.577350269189626D0,0.577350269189626D0/
      data Gauss3/-0.7745966692415D0,0.0d0,0.7745966692415D0/
      data Gaussw3/0.555555555556D0,0.888888888889D0,0.555555555556D0/
      data ncrd  /-1.0d0, 1.0d0, 0.0d0/
      data Gaussw2/1.0d0,1.0d0/

#if defined(DEBUG)
      call SubBeg('BeamElement')
#endif
      
      ishape = 0
      call vi_move (elChar(KYOP1),keyopt(1),12)
      call vi_move (kelreq(1),kelin(1),10)
      nTens = 3
      nDirect = 1
      WetMode = .FALSE.
      compSSTIF = .TRUE.
      if (WetMode) compSSTIF = .FALSE.
      
      call v_zero(tGenStrain(1,1),maxGenStr*2)
      thStrain(1) = ZERO
      thStrain(2) = ZERO
      eEnergy(1) = ZERO
      eEnergy(4) = ZERO
      eEnergy(5) = ZERO

c --- get element attributes
      matId = elData(EL_MAT)
      ireal = elData(EL_REAL)
      nreal = elChar(NMDRLC)
      call v_zero (rvr(1),nreal)
      i = rlgetn(ireal,1,nreal,rvr(1))
      
      if (i.lt.2) then
          call erhandler ('BeamElement',5000,3,
     &                    'Beam element %I has no properties defined',
     &                     dble(elemId),' ')
      endif
      
      numBasicSaveVars = elChar(NMSSVR)
      if (.not.(B_PtrRealN.ne.PTRFTNNULL .and. B_RealN(ireal).gt.0))then
          if (keyopt(7).eq.1) then
              numBasicSaveVars = numBasicSaveVars + 20
     &                         + B_NumResPts*numSectPts
          endif
      endif
          
c --- material nonlinear properties index

      keyplL = nlkword(1,1)
      keycrL = nlkword(1,2)
      keynlL = keyplL + keycrL

      call v_zero(siglvi(1),nTens)
      call v_zero(siglv(1),nTens)

c --- read real constants for the element
      isWarp = .FALSE.
      if (keyopt(2).eq.1) isWarp = .TRUE.

      iniFlag = .FALSE.
      if (keyopt(3).eq.1) then
          iniFlag = .TRUE.
      endif

c --- fetch section data
c
      sectSummary(1,1) = ZERO

      if (numElemNodes .eq. 2) then
          numDomainIntPts = 1
      elseif (numElemNodes .eq. 3) then
          numDomainIntPts = 2
      endif

      numSectNodes=secDataI(pNodes)
      ii = 6*numSectPts
      call v_zero(matStress(1,1,1),ii)
      call v_zero(matStrain(1,1,1),ii)

      if (A_KeyFirstPass.eq.1) then
         if (iniflag) call v_zero(iniForce(1),21)
         yzScale_t(1) = ONE
         yzScale_t(2) = ONE
      endif

c --- non-linear save variables
      call GetSVRIndx(elemId,locsvrL,stateVarIndex(1))
      tmp_svr = numBasicSaveVars
      call GetSVR (stateVarIndex(1),1,numBasicSaveVars,basicSaveVars(1))

      if (tmp_svr .lt. numBasicSaveVars) then
          call erhandler ('BeamElement',5000,3,
     &        'SVR record is mismatched in beam element',0.0d0,' ')
          goto 999
      endif
c --- only 1st column has data.
c --- both domIntP & domInt are identical
      domIntP = 1
      nodeP = 1
      call CalBeamSect (ireal,keyopt(1),secDataI(1),secDataR(1,1),
     &                  sectSummary(1,1),secResults(1,1,1),NumRstNodes,
     &                  basicSaveVars(elChar(NMSSVR)-2))

      yOffset(1) = sectSummary(pCgy,nodeP)
      zOffset(1) = sectSummary(pCgz,nodeP)

      if (LenMatRecMax.gt.0 .and. keyopt(1).ne.2) then
          call IncoreMatIndx(mRecIntPidx(1),W_LenMatRecIndx,numWordsPt)
      endif
c
      offY = 1
      offZ = secDataI(pNodes)+1
c
      if (kelin(1).eq.1) then
          call v_zero(eStiff(1,1), nUsrDof*nUsrDof)
      endif
      if (kelin(2).eq.1)then
          call v_zero(eMass(1,1) , nUsrDof*nUsrDof)
          kelout(2)=1
      endif
      if (kelin(5).eq.1) call v_zero (fExt(1), nUsrDof)
      if (kelin(6).eq.1) call v_zero (fInt(1), nUsrDof)
      if (kelin(4).eq.1 .or. kelin(1).eq.1)
     &    call v_zero (gStiff(1,1), nUsrDof*nUsrDof)

c --- get the element loads
c
      call v_zero (moistcur(1),9)

c#if defined(MOIST)
      refmoist = GetMoistRef(matId)
      call GetMoist (elemId,elChar(1),numElemNodes,eConn(1),
     &               refmoist,numElemNodes*3,moistbeg(1),
     &               moistcur(1),moistend(1),khygro,tab(1))
c#else
c      khygro = 0
c      refmoist = ZERO
c      call v_zero (moistbeg(1),9)
c      call v_zero (moistend(1),9)
c      call v_zero (tab(1),3)
c#endif

      if (A_KeyFirstPass.eq.1) then
          call v_move (moistbeg(1),moistprev(1),3*numElemNodes)
      else
#if defined(PRS)
          do i = 1,14
             if (prsBeg(1).ne.ZERO) goto 50
          enddo
          kpress = 0
          goto 60
#endif
      endif

 50   continue
#if defined(PRS)
      call prsget (elemId,elChar(1),5,14,
     &             prsBeg(1),pressure(1),prsEnd(1),kpress)
#else
      kpress = 0
#endif

 60   elasticSection  = 0
      if (keynlL.eq.0) then
          if  (moistcur(2).eq.moistcur(1).and.
     &         moistcur(3).eq.moistcur(1).and.
     &         moistcur(5).eq.moistcur(4).and.
     &         moistcur(6).eq.moistcur(4))then
               elasticSection = 1
          endif
          if (iniFlag) elasticSection = 0
          if (keyopt(1).eq.1) elasticSection = 0
      endif
      if (keyopt(1).ne.1) elasticSection = 1

      numBasicSaveVars  = tmp_svr
      tmp_svr = LenMatRecMax

      if (LenMatRecMax.gt.0 .and. keyopt(1).ne.2)
     &    call UpdSvrRec (0, elemId, stateVarIndex(1),LenMatRecMax)

      if (tmp_svr .ne. LenMatRecMax) then
          call erhandler ('BeamElement',5000,3,
     &        'SVR record is mismatched in beam element',0.0d0,' ')
          goto 999
      endif
      
      call v_move (eleDisp(1,A_RowDisp),eU(1),nUsrDof)

      if (WetMode) then
          call m_move (xBeg(1,1),xCur(1,1),6,3,3,numElemNodes)
          call v_zero(incDisp(1), nUsrDof)
      elseif (A_KeyFirstPass .eq.1) then
          call m_move (xBeg(1,1),xCur(1,1),6,3,3,numElemNodes)
          call v_move (eleDisp(1,A_RowIncS),incDisp(1),nUsrDof)
      else
          call UpdBeamXYZ (numElemNodes,eU(1),nUsrDof,
     &                    xBeg(1,1),6,xCur(1,1))
          call v_move (eleDisp(1,A_RowIncS),incDisp(1),nUsrDof)
      endif

c --- calculate nodal basis vectors
      call v_zero (nBasis(1,1,1),27)
      do knod = 1,numElemNodes
         if (numElemNodes .eq. 2) then
            shpe(1,2)=-HALF
            shpe(2,2)= HALF
         elseif (numElemNodes .eq. 3) then
            shpe(1,2)=ncrd(knod)-HALF
            shpe(2,2)=ncrd(knod)+HALF
            shpe(3,2)=-TWO*ncrd(knod)
         endif
         dL = ZERO
         do j = 1,3
            dxdg(j) = ZERO
            do nod=1,numElemNodes
               dxdg(j) = dxdg(j)+shpe(nod,2)*xCur(j,nod)
            enddo
            dL=dL+dxdg(j)*dxdg(j)
         enddo
         dL = sqrt(dL)
         if (dL .le. sqrt(TINY)) then
            ishape = 1
            goto 999
         endif
         dLinv = ONE/dL
         nBasis(1,1,knod)=dxdg(1)*dLinv
         nBasis(2,1,knod)=dxdg(2)*dLinv
         nBasis(3,1,knod)=dxdg(3)*dLinv
         if (A_KeyFirstPass.eq.1) then         
             if (eConn(numElemNodes+1).eq.0) then
                 call TransBeam (dxdg(1),trv(1,1))
                 nBasis(1,2,knod)=trv(2,1)
                 nBasis(2,2,knod)=trv(2,2)
                 nBasis(3,2,knod)=trv(2,3)
                 nBasis(1,3,knod)=trv(3,1)
                 nBasis(2,3,knod)=trv(3,2)
                 nBasis(3,3,knod)=trv(3,3)
                 call v_move(nBasis(1,1,knod),nBasis_t(1,1,knod),9)
                 call v_move(nBasis(1,1,knod),nBasis_0(1,1,knod),9)
             else
c --- beam orientation using 4-th node (under development)
                j=numElemNodes+1
                work(1)=xBeg(1,j)-xBeg(1,knod)
                work(2)=xBeg(2,j)-xBeg(2,knod)
                work(3)=xBeg(3,j)-xBeg(3,knod)
                work(4)=work(1)*work(1)+work(2)*work(2)+work(3)*work(3)
                if (work(4) .le. sqrt(TINY)) then
                    call erhandler ('BeamElement',5000, 3,
     &                   'Beam element %I has zero orientation vector',
     &                   dble(elemId),' ')
                    goto 999
                endif
                work(4)=ONE/sqrt(work(4))
                work(1)=work(4)*work(1)
                work(2)=work(4)*work(2)
                work(3)=work(4)*work(3)
                nBasis(1,2,knod)=
     &               -nBasis(2,1,knod)*work(3)+nBasis(3,1,knod)*work(2)
                nBasis(2,2,knod)=
     &               -nBasis(3,1,knod)*work(1)+nBasis(1,1,knod)*work(3)
                nBasis(3,2,knod)=
     &               -nBasis(1,1,knod)*work(2)+nBasis(2,1,knod)*work(1)
                work(4)=
     &               nBasis(1,2,knod)*nBasis(1,2,knod)+
     &               nBasis(2,2,knod)*nBasis(2,2,knod)+
     &               nBasis(3,2,knod)*nBasis(3,2,knod)
                if (work(4) .le. sqrt(TINY)) then
                    call erhandler ('BeamElement',5000, 3,
     &                   'Beam element %I has erroneous input.',
     &                   dble(elemId),' ')
                   goto 999
                endif
                work(4)=ONE/sqrt(work(4))
                nBasis(1,2,knod)=work(4)*nBasis(1,2,knod)
                nBasis(2,2,knod)=work(4)*nBasis(2,2,knod)
                nBasis(3,2,knod)=work(4)*nBasis(3,2,knod)
c
                nBasis(1,3,knod)=
     &               nBasis(2,1,knod)*nBasis(3,2,knod)-
     &               nBasis(3,1,knod)*nBasis(2,2,knod)
                nBasis(2,3,knod)=
     &               nBasis(3,1,knod)*nBasis(1,2,knod)-
     &               nBasis(1,1,knod)*nBasis(3,2,knod)
                nBasis(3,3,knod)=
     &               nBasis(1,1,knod)*nBasis(2,2,knod)-
     &               nBasis(2,1,knod)*nBasis(1,2,knod)
                call v_move(nBasis(1,1,knod),nBasis_t(1,1,knod),9)
                call v_move(nBasis(1,1,knod),nBasis_0(1,1,knod),9)
             endif
         endif

         if (WetMode) then
             call v_move(nBasis_t(1,2,knod),nBasis(1,2,knod),6)
         else
             call BasRot (nBasis_t(1,2,knod),nBasis(1,2,knod),
     &                    incDisp((knod-1)*nDofNode+4),2)
         endif

      enddo

      do j=1,numElemNodes
         work(1) = -(yOffset(nodeP)-sectSummary(pCgy,nodeP))
         work(2) = -(zOffset(nodeP)-sectSummary(pCgz,nodeP))
         moffset(1,1,j) = ZERO
         moffset(2,2,j) = ZERO
         moffset(3,3,j) = ZERO
         moffset(1,2,j) =  nBasis(3,2,j)*work(1)+nBasis(3,3,j)*work(2)
         moffset(1,3,j) = -nBasis(2,2,j)*work(1)-nBasis(2,3,j)*work(2)
         moffset(2,1,j) = -nBasis(3,2,j)*work(1)-nBasis(3,3,j)*work(2)
         moffset(2,3,j) =  nBasis(1,2,j)*work(1)+nBasis(1,3,j)*work(2)
         moffset(3,1,j) =  nBasis(2,2,j)*work(1)+nBasis(2,3,j)*work(2)
         moffset(3,2,j) = -nBasis(1,2,j)*work(1)-nBasis(1,3,j)*work(2)
      enddo

c --- calculate basis at integration points

      elemCG(1)=ZERO
      elemCG(2)=ZERO
      elemCG(3)=ZERO
      elemVolume = ZERO
      call v_zero (basis(1,1,1),18)
      yzScale(2) = ONE
      L = ZERO
      do domInt = 1,numDomainIntPts
         domIntP = 1
         if (numElemNodes .eq. 2) then
            shpe(1,1)= HALF
            shpe(2,1)= HALF
            shpe(1,2)=-HALF
            shpe(2,2)= HALF
         elseif (numElemNodes .eq. 3) then
            shpe(1,1)=HALF*(Gauss2(domInt)-ONE)*Gauss2(domInt)
            shpe(2,1)=HALF*(Gauss2(domInt)+ONE)*Gauss2(domInt)
            shpe(3,1)=ONE-Gauss2(domInt)*Gauss2(domInt)
            shpe(1,2)=Gauss2(domInt)-HALF
            shpe(2,2)=Gauss2(domInt)+HALF
            shpe(3,2)=-TWO*Gauss2(domInt)
         endif
         dL = ZERO
         do j=1,3
            dxdg(j) = ZERO
            do nod=1,numElemNodes
               dxdg(j) = dxdg(j)+shpe(nod,2)*xCur(j,nod)
            enddo
            dL=dL+dxdg(j)*dxdg(j)
         enddo
         dL    = sqrt(dL)
         dLinv = ONE/dL
         if (numDomainIntPts .eq. 1) then
             dV = TWO*dL
         else
             dV = dL
         endif
         L = L + dV
         basis(1,1,domInt)=dxdg(1)*dLinv
         basis(2,1,domInt)=dxdg(2)*dLinv
         basis(3,1,domInt)=dxdg(3)*dLinv
         if (A_KeyFirstPass.eq.1) then
             do j=1,numElemNodes
                shpDer(j,domInt)=shpe(j,2)*dLinv
             enddo
             dV_t0(domInt) = dV
             dL_t0(domInt) = dL
             if (eConn(numElemNodes+1).eq.0) then
                  call TransBeam (dxdg(1),trv(1,1))
                  basis(1,2,domInt)=trv(2,1)
                  basis(2,2,domInt)=trv(2,2)
                  basis(3,2,domInt)=trv(2,3)
                  basis(1,3,domInt)=trv(3,1)
                  basis(2,3,domInt)=trv(3,2)
                  basis(3,3,domInt)=trv(3,3)
                  call v_move (basis(1,1,domInt),basis_t(1,1,domInt),9)
                  call v_move (basis(1,1,domInt),basis_0(1,1,domInt),9)
             else
                  j = numElemNodes+1
                  if (numElemNodes .eq. 2) then
                      work(1)=xBeg(1,j)-HALF*(xBeg(1,1)+xBeg(1,2))
                      work(2)=xBeg(2,j)-HALF*(xBeg(2,1)+xBeg(2,2))
                      work(3)=xBeg(3,j)-HALF*(xBeg(3,1)+xBeg(3,2))
                  else
                      work(1)=xBeg(1,j)-xBeg(1,3)
                      work(2)=xBeg(2,j)-xBeg(2,3)
                      work(3)=xBeg(3,j)-xBeg(3,3)
                  endif
c --- beam orientation using Kth node (under development)
                  work(4)= work(1)*work(1)
     &                   + work(2)*work(2)
     &                   + work(3)*work(3)
                  work(4)=ONE/sqrt(work(4))
                  work(1)=work(4)*work(1)
                  work(2)=work(4)*work(2)
                  work(3)=work(4)*work(3)
                  basis(1,2,domInt)=
     &              -basis(2,1,domInt)*work(3)+basis(3,1,domInt)*work(2)
                  basis(2,2,domInt)=
     &              -basis(3,1,domInt)*work(1)+basis(1,1,domInt)*work(3)
                  basis(3,2,domInt)=
     &              -basis(1,1,domInt)*work(2)+basis(2,1,domInt)*work(1)

                  work(4)=
     &                basis(1,2,domInt)*basis(1,2,domInt)+
     &                basis(2,2,domInt)*basis(2,2,domInt)+
     &                basis(3,2,domInt)*basis(3,2,domInt)
                  work(4)=ONE/sqrt(work(4))
                      basis(1,2,domInt)=work(4)*basis(1,2,domInt)
                      basis(2,2,domInt)=work(4)*basis(2,2,domInt)
                      basis(3,2,domInt)=work(4)*basis(3,2,domInt)
c
                 basis(1,3,domInt)=
     &                basis(2,1,domInt)*basis(3,2,domInt)-
     &                basis(3,1,domInt)*basis(2,2,domInt)
                 basis(2,3,domInt)=
     &                basis(3,1,domInt)*basis(1,2,domInt)-
     &                basis(1,1,domInt)*basis(3,2,domInt)
                 basis(3,3,domInt)=
     &                basis(1,1,domInt)*basis(2,2,domInt)-
     &                basis(2,1,domInt)*basis(1,2,domInt)
            endif
            call v_move (basis(1,1,domInt),basis_t(1,1,domInt),9)
            call v_move (basis(1,1,domInt),basis_0(1,1,domInt),9)
         endif
c --- interpolated incremental rotation vector
         if (WetMode) then
            lambda = ONE
            call v_move (basis_t(1,2,domInt),basis(1,2,domInt),6)
            call v_move (basis_t(1,1,domInt),n2xn3(1),3)
         else
            work(1) = ZERO
            work(2) = ZERO
            work(3) = ZERO
            do nod = 1,numElemNodes
               j = (nod-1)*nDofNode
               work(1)=work(1)+shpe(nod,1)*incDisp(j+4)
               work(2)=work(2)+shpe(nod,1)*incDisp(j+5)
               work(3)=work(3)+shpe(nod,1)*incDisp(j+6)
            enddo
            call BasRot(basis_t(1,2,domInt),basis(1,2,domInt),work(1),2)
            lambda = ZERO
            do j=1,3
               work(j) = ZERO
               do nod=1,numElemNodes
                  work(j) = work(j)+shpDer(nod,domInt)*xCur(j,nod)
               enddo
               lambda=lambda+work(j)*work(j)
            enddo
            lambda = ONE/sqrt(lambda)
            n2xn3(1)=basis(2,2,domInt)*basis(3,3,domInt)-
     &           basis(3,2,domInt)*basis(2,3,domInt)
            n2xn3(2)=basis(3,2,domInt)*basis(1,3,domInt)-
     &           basis(1,2,domInt)*basis(3,3,domInt)
            n2xn3(3)=basis(1,2,domInt)*basis(2,3,domInt)-
     &           basis(2,2,domInt)*basis(1,3,domInt)
         endif
         
         do j=1,3
            xIntp(j,domInt)=ZERO
            do nod=1,numElemNodes
               xIntp(j,domInt)=xIntp(j,domInt)+shpe(nod,1)*xCur(j,nod)
            enddo
         enddo
         elemVolume=elemVolume+dV_t0(domInt)*sectSummary(pArea,domIntP)
         do j = 1,3
            elemCG(j)=elemCG(j)+
     &         dV_t0(domInt)*sectSummary(pArea,domIntP)*xIntp(j,domInt)
         enddo
         
         if (WetMode) then
            aScale = ONE
            yzScale(domInt) = ONE         
         else
            if (keyopt(8).eq.1) then
                aScale = ONE
                yzScale(domInt) = ONE
            else
                aScale  =dL_t0(domInt)/dL
                yzScale(domInt) =sqrt(aScale)
            endif
c --- dDphi/dS
            work(21) = ZERO
            work(22) = ZERO
            work(23) = ZERO
            do nod = 1,numElemNodes
               j = (nod-1)*nDofNode
               work(21)=work(21)+shpDer(nod,domInt)*incDisp(j+4)
               work(22)=work(22)+shpDer(nod,domInt)*incDisp(j+5)
               work(23)=work(23)+shpDer(nod,domInt)*incDisp(j+6)
            enddo
            work(1) = ZERO
            work(2) = ZERO
            work(3) = ZERO
            do nod = 1,numElemNodes
               j = (nod-1)*nDofNode
               work(1)=work(1)+shpe(nod,1)*incDisp(j+4)
               work(2)=work(2)+shpe(nod,1)*incDisp(j+5)
               work(3)=work(3)+shpe(nod,1)*incDisp(j+6)
            enddo
            work(17)=
     &           sqrt(work(1)*work(1)+work(2)*work(2)+work(3)*work(3))
            if (work(17).gt.1.0D-6) then
c --- p vector -- axis around which rotation is defined
               work(1)  = work(1)/work(17)
               work(2)  = work(2)/work(17)
               work(3)  = work(3)/work(17)
               work(18) = cos(work(17))
               work(19) = sin(work(17))
c
               work(18) = (ONE-work(18))/work(17)
               work(19) = work(19)/work(17)
c --- n.dDphi/dS
               work(20)=
     &              work(1)*work(21)+work(2)*work(22)+work(3)*work(23)
c
               work(24)=work(19)*work(21)+
     &              work(18)*(work(2)*work(23)-work(3)*work(22))
     &              +(ONE-work(19))*work(20)*work(1)
               work(25)=work(19)*work(22)+
     &              work(18)*(work(3)*work(21)-work(1)*work(23))
     &              +(ONE-work(19))*work(20)*work(2)
               work(26)=work(19)*work(23)+
     &              work(18)*(work(1)*work(22)-work(2)*work(21))
     &              +(ONE-work(19))*work(20)*work(3)
            else
               work(24) = work(21)
               work(25) = work(22)
               work(26) = work(23)
            endif
            if (A_KeyFirstPass.eq.1) then
                call v_zero(tGenStrain_t(1,domInt),nGenStr)
            endif

            tGenStrain(1,domInt) = log(dL/dL_t0(domInt))
            tGenStrain(2,domInt)=tGenStrain_t(2,domInt) +
     &           work(24)*basis(1,2,domInt)+
     &           work(25)*basis(2,2,domInt)+
     &           work(26)*basis(3,2,domInt)
            tGenStrain(3,domInt)=tGenStrain_t(3,domInt) +
     &           work(24)*basis(1,3,domInt)+
     &           work(25)*basis(2,3,domInt)+
     &           work(26)*basis(3,3,domInt)
            if (keyopt(2).eq.2) then
                tGenStrain(4,domInt) = ZERO
            else
                tGenStrain(4,domInt)=tGenStrain_t(4,domInt) +
     &               work(24)*n2xn3(1)+
     &               work(25)*n2xn3(2)+
     &               work(26)*n2xn3(3)
            endif
            tGenStrain(5,domInt)=
     &           basis(1,1,domInt)*basis(1,3,domInt)+
     &           basis(2,1,domInt)*basis(2,3,domInt)+
     &           basis(3,1,domInt)*basis(3,3,domInt)
            tGenStrain(6,domInt)=
     &           basis(1,1,domInt)*basis(1,2,domInt)+
     &           basis(2,1,domInt)*basis(2,2,domInt)+
     &           basis(3,1,domInt)*basis(3,2,domInt)
            tGenStrain(1,domInt)=tGenStrain(1,domInt)+
     &           HALF*aScale*(lambda**2)*
     &           (sectSummary(pIpp,domIntP)/sectSummary(pArea,domIntP))*
     &           (tGenStrain(4,domInt)**2-tGenStrain_t(4,domInt)**2)
            work(1)=
     &           -yzScale(domInt)*
     &           (sectSummary(pCgy,domIntP)-sectSummary(pSy,domIntP))
            work(2)=
     &           yzScale(domInt)*
     &           (sectSummary(pCgz,domIntP)-sectSummary(pSz,domIntP))
            tGenStrain(5,domInt)=tGenStrain(5,domInt)+
     &           work(1)*tGenStrain(4,domInt)
            tGenStrain(6,domInt)=tGenStrain(6,domInt)+
     &           work(2)*tGenStrain(4,domInt)
            work(1) =
     &           (yOffset(domIntP)-
     &           sectSummary(pCgy,domIntP))*yzScale(domInt)
            work(2) =
     &           (zOffset(domIntP)-
     &           sectSummary(pCgz,domIntP))*yzScale(domInt)
            tGenStrain(1,domInt)=tGenStrain(1,domInt)+
     &           work(1)*tGenStrain(3,domInt)-
     &           work(2)*tGenStrain(2,domInt)
            if (isWarp) then
                tGenStrain(7,domInt) = ZERO
                do knod=1,numElemNodes
                   j=(knod-1)*nDofNode
                   tGenStrain(7,domInt)=tGenStrain(7,domInt)+
     &                        shpDer(knod,domInt)*eU(j+nDofNode)
                enddo 
            endif
         endif
            
         call v_zero(gradMat(1,1),nGenStr*nUsrDof)
         do knod=1,numElemNodes
            j=(knod-1)*nDofNode
            do jj=1,3
               gradMat(1,j+jj  )= lambda*
     &              basis(jj,1,domInt)*shpDer(knod,domInt)
               gradMat(2,j+jj+3)=basis(jj,2,domInt)*shpDer(knod,domInt)
               gradMat(3,j+jj+3)=basis(jj,3,domInt)*shpDer(knod,domInt)
               gradMat(4,j+jj+3)=n2xn3(jj)*shpDer(knod,domInt)
               gradMat(5,j+jj)  = lambda*
     &              basis(jj,3,domInt)*shpDer(knod,domInt)
               gradMat(5,j+jj+3)=basis(jj,2,domInt)*shpe(knod,1)
               gradMat(6,j+jj)  = lambda*
     &              basis(jj,2,domInt)*shpDer(knod,domInt)
               gradMat(6,j+jj+3)=-basis(jj,3,domInt)*shpe(knod,1)
            enddo
            if (isWarp) then
                gradMat(7,j+nDofNode)=shpDer(knod,domInt)
            endif
         enddo

c --- transverse shear will cause a torsional moment unless it is
c     passing thru shear center.
c
         work(1) = -yzScale(domInt)*
     &        (sectSummary(pCgy,domIntP)-sectSummary(pSy,domIntP))
         work(2) =  yzScale(domInt)*
     &        (sectSummary(pCgz,domIntP)-sectSummary(pSz,domIntP))
         if (WetMode) then
             work(3) = ZERO
         else
             work(3) = aScale*(lambda**2)*
     &             sectSummary(pIpp,domIntP)/sectSummary(pArea,domIntP)*
     &             tGenStrain(4,domInt)
         endif
         do jj = 1,nUsrDof
            gradMat(1,jj) = gradMat(1,jj)+work(3)*gradMat(4,jj)
            gradMat(5,jj) = gradMat(5,jj)+work(1)*gradMat(4,jj)
            gradMat(6,jj) = gradMat(6,jj)+work(2)*gradMat(4,jj)
         enddo
c
         do j=1,numElemNodes
            jj=(j-1)*nDofNode
            do ll = 1,6
               do kk = 1,3
                  gradMat(ll,jj+kk+3)=gradMat(ll,jj+kk+3)+
     &                 gradMat(ll,jj+1)*moffset(1,kk,j)*yzScale(domInt)+
     &                 gradMat(ll,jj+2)*moffset(2,kk,j)*yzScale(domInt)+
     &                 gradMat(ll,jj+3)*moffset(3,kk,j)*yzScale(domInt)
               enddo
            enddo
         enddo
         if (isWarp) then
             do j=1,nUsrDof
                bwarp(j)=gradMat(4,j)
             enddo
             do j=1,numElemNodes
                bwarp((j-1)*nDofNode+nDofNode)=-shpe(j,1)
             enddo
             eWarp(domInt)=ZERO
             do j=1,nUsrDof
                eWarp(domInt)=eWarp(domInt)+bwarp(j)*eU(j)
             enddo
         endif

         if (WetMode) then
             call Axv1_v2(gradMat(1,1),eU(1),tGenStrain(1,domInt),
     &                    nGenStr,nUsrDof)
         endif

         call v_zero(secStiff(1,1),nGenStr*nGenStr)
         call v_zero(tGenStress(1,domInt),nGenStr)
         if (numElemNodes .eq. 2) then
             scaleSF1 = 0.25d0
             scaleSF2 = 0.25d0
          else
             scaleSF1 = 0.25d-4
             scaleSF2 = 0.25d-4
          endif
          do j=17,19,2
             if (sectSummary(j,domIntP).le.TINY)
     &           sectSummary(j,domIntP)=1.0d0
          enddo
          scaleSF1 = sectSummary(19,domIntP)/
     &           (ONE+
     &           (scaleSF1*dL_t0(domInt)*dL_t0(domInt)
     &           *sectSummary(pArea,domIntP)
     &           /(TWELVE * sectSummary(pIyy,domIntP))))
          scaleSF2 = sectSummary(17,domIntP)/
     &           (ONE+
     &           (scaleSF2*dL_t0(domInt)*dL_t0(domInt)*
     &           sectSummary(pArea,domIntP)
     &           /(TWELVE * sectSummary(pIzz,domIntP))))
          
      
         if (numElemNodes.eq.2) then
             tInt(1)=HALF
             tInt(2)=HALF
         elseif(numElemNodes.eq.3)then
             tInt(1)=HALF*(ONE-Gauss2(domInt))
             tInt(2)=HALF*(ONE+Gauss2(domInt))
         endif
         tDomInt   = tInt(1)*moistcur(1) + tInt(2)*moistcur(4)
         tDomInt_t = tInt(1)*moistprev(1) + tInt(2)*moistprev(4)

#if defined(MOIST)
         eh = ZERO
         if (K_nMast.gt.0) then
             ii = (elemId-1)*H_Row*H_Col
             do jj = 1,H_Row
                em(jj) = M_HygStrains(ii+jj)
             enddo
         endif

         if (numDomainIntPts .eq. 1) then
             xn = ZERO
         else
             xn = Gauss2(domInt)
         endif      
      
       
!         eh = ZERO  !!!! to disable hygro strains
c#else
c         tDomInt   = ZERO
c         tDomInt_t = ZERO
#endif defined(MOIST)
         
         call v_zero(tGenTh(1),nGenStr)
         if (elasticSection.eq.1) then
            call GetElasticProps (elemId,matId,tDomInt,prop(1))
            secStiff(1,1)= sectSummary(pArea,domIntP)*aScale
            secStiff(1,2)= sectSummary(pIz,domIntP)*aScale
     &           *yzScale(domInt)
            secStiff(1,3)=-sectSummary(pIy,domIntP)*aScale*
     &           yzScale(domInt)
            secStiff(2,1)= sectSummary(pIz,domIntP)*aScale*
     &           yzScale(domInt)
            secStiff(2,2)= sectSummary(pIyy,domIntP)*aScale*aScale
            secStiff(2,3)=-sectSummary(pIyz,domIntP)*aScale*aScale
            secStiff(3,1)=-sectSummary(pIy,domIntP)*aScale*
     &           yzScale(domInt)
            secStiff(3,2)=-sectSummary(pIyz,domIntP)*aScale*aScale
            secStiff(3,3)= sectSummary(pIzz,domIntP)*aScale*aScale
            secStiff(4,4)= sectSummary(pTorq,domIntP)*aScale*aScale
            secStiff(5,5)= sectSummary(pArea,domIntP)*aScale
            secStiff(6,6)= sectSummary(pArea,domIntP)*aScale
            YoungsMod = prop(1)
            ShearMxy  = prop(5)
            ShearMyz  = prop(6)
            ShearMxz  = prop(4)
            alpha     = prop(10)
            density   = prop(13)
            do ii=1,3
               do jj=1,3
                  secStiff(ii,jj) = YoungsMod*secStiff(ii,jj)
               enddo
            enddo
            secStiff(4,4) = ShearMxy * secStiff(4,4)
            secStiff(5,5) = ShearMyz * secStiff(5,5)
            secStiff(6,6) = ShearMxz * secStiff(6,6)
            if (isWarp) secStiff(7,7)=
     &          YoungsMod*sectSummary(pSww,domIntP)
            thStrain(domInt) = ZERO
#if defined(MOIST)
            if (khygro.eq.1) then
                thStrain(domInt) = alpha*(tDomInt-refmoist)
                if (K_nMast.gt.0) then
                    call InterpHygStrn(xn,ZERO,ZERO,em(1),eh)
                    thStrain(domInt) = thStrain(domInt) + eh
                endif
                
                if (A_KeyFirstIter.eq.1) then
                    tGenTh(1)=thStrain(domInt)*secStiff(1,1)
                else
                    tGenStrain(1,domInt)=
     &                         tGenStrain(1,domInt)-thStrain(domInt)
                endif
            endif
#endif defined(MOIST)            
            
            call Axv1_v2(secStiff(1,1), tGenStrain(1,domInt),
     &                   tGenStress(1,domInt),nGenStr,nGenStr)
            do ii=1,4
               eEnergy(1)=eEnergy(1)+HALF*dV*
     &                   (tGenStress(ii,domInt)*tGenStrain(ii,domInt))
            enddo
            do kspt=1,numSectPts
               tSecPt((domInt-1)*numSectPts+kspt)=tDomInt
            enddo
         else
c            tInt(3)   = tInt(1)*moistcur(2)+tInt(2)*moistcur(5)
c            tInt(4)   = tInt(1)*moistcur(3)+tInt(2)*moistcur(6)
c            tInt_t(3) = tInt(1)*moistprev(2)+tInt(2)*moistprev(5)
c            tInt_t(4) = tInt(1)*moistprev(3)+tInt(2)*moistprev(6)
            tInt(3)   = tInt(1)*moistcur(1)+tInt(2)*moistcur(4)
            tInt(4)   = tInt(3)
            tInt_t(3) = tInt(1)*moistprev(1)+tInt(2)*moistprev(4)
            tInt_t(4) = tInt_t(3)
             
            call v_zero(siglvis(1),7)
            do ii = 1,3
               coordAx(ii) = 0.0d0
               do jj = 1,numElemNodes
                  coordAx(ii) = coordAx(ii) + shpe(jj,1)*xBeg(ii,jj)
               enddo
            enddo
            densityDA=0.0d0
            khygro = 0
            do k1=1,secDataI(pElems)
               do k2=1,4
                  kspt=(k1-1)*4+k2
                  if (secResults(pM,kspt,domIntP).gt.ZERO) then
                      
                  if (A_IterNum.gt.1) then
                      ii = 0
                  endif
                  if (A_SubStepNum.eq.2) then
                      i = 0
                  endif
                  do ii = 1,3
                     coords_0(ii) = coordAx(ii)
     &                + secResults(pY,kspt,domIntP)*basis_0(ii,2,domInt)
     &                + secResults(pZ,kspt,domIntP)*basis_0(ii,3,domInt)
                  enddo
                  y_pt=secResults(pY,kspt,domIntP)*yzScale(domInt)
                  z_pt=secResults(pZ,kspt,domIntP)*yzScale(domInt)
                  
                  refmoist = GetMoistRef(matId)

                  if (tab(domInt).ne.ZERO) then
                      call GetTabMoist(tab(domInt),coords_0(1),tSecInt)
                  else
                     tSecInt =
     &                    tDomInt +
     &                    (tInt(3)-tDomInt)*y_pt+
     &                    (tInt(4)-tDomInt)*z_pt
c                     tSecInt_t =
c     &                    tDomInt_t +
c     &                    (tInt_t(3)-tDomInt_t)*y_pt+
c     &                    (tInt_t(4)-tDomInt_t)*z_pt
                  endif
                  
#if defined(MOIST)
                  width = rvr(2)
                  if (nint(rvr(1)).le.2) then
c --- rectangular cross section - output height and not the radius
                      hight = rvr(3)      
                  else
                      hight = width
                  endif

c --- the original coordinates without scaling
                  yn = secResults(pY,kspt,domIntP)/width*TWO
                  zn = secResults(pZ,kspt,domIntP)/hight*TWO
                  if (K_nMast.gt.0) then                  
                      call InterpHygStrn(xn,yn,zn,em(1),eh)
                  endif
                  
#endif defined(MOIST)            
                  
                  tSecPt((domInt-1)*numSectPts+kspt)=tSecInt
                  dA_pt=secResults(pdA,kspt,domIntP)*aScale

c --- change from stress function to warping function
                  twork(1)=y_pt-sectSummary(pSy,domIntP)+
     &                 sectSummary(pCgy,domIntP)
                  twork(2)=z_pt-sectSummary(pSz,domIntP)+
     &                 sectSummary(pCgz,domIntP)
                  phiy_pt=(secResults(pphiy,kspt,domIntP)-twork(2))*
     &                 yzScale(domInt)
                  phiz_pt=(secResults(pphiz,kspt,domIntP)+twork(1))*
     &                 yzScale(domInt)
                  matStrain(1,domInt,kspt)=
     &                 tGenStrain(1,domInt)+
     &                 z_pt*tGenStrain(2,domInt)-
     &                 y_pt*tGenStrain(3,domInt)
                  if (isWarp) then
                      matStrain(1,domInt,kspt)=matStrain(1,domInt,kspt)+
     &                   tGenStrain(7,domInt)*secResults(2,kspt,domIntP)
                  endif
c --- change from stress function to warping function
                  matStrain(2,domInt,kspt)=tGenStrain(4,domInt)*phiy_pt
                  matStrain(3,domInt,kspt)=tGenStrain(4,domInt)*phiz_pt
                  mStrain(1) = matStrain(1,domInt,kspt)
                  mStrain(2) = matStrain(2,domInt,kspt)
                  mStrain(3) = matStrain(3,domInt,kspt)
                  ii=(domInt-1)*numSectPts+kspt

c --- update beam incremental strain
                  call RetMatRec0 (W_iStrain, ii, nTens, work(1))
                  call v1mv2_v3(mStrain(1),work(1),mStrainInc(1),nTens)

c --- k1 is the beam cell index
c --- k2 is the section integration point within the cell
                  call MaterialMain (elemId, ii, matId, ireal, 
     &                               nTens, nlkword(1,1), numWordsPt,
     &                               mRecIntPidx(1), domInt, kspt, ii,
     &                               tSecInt, mStrainInc(1),
     &                               lwork, work(1), LenMatRecMax,
     &                               refmoist, khygroIP,
     &                               coords_0(1),basis(1,1,domInt),
     &                               iniFlag, matRecord(1),
     &                               mStrain(1),mStress(1),
     &                               matStiff(1,1),prop(1),
     &                               siglv(1), siglvi(1), eEnerSecPt(1)
#if defined(MOIST)
     &                               ,eh)
#else
     &                               )
#endif



                  if (khygroIP .eq. 1) khygro = 1
                  matStress(1,domInt,kspt) = mStress(1)
                  matStress(2,domInt,kspt) = mStress(2)
                  matStress(3,domInt,kspt) = mStress(3)
                  YoungsMod = prop(1)
                  ShearMxy = prop(2)
                  ShearMxz = prop(3)
                  density  = prop(5)
                  densityDA = densityDA+density*dA_pt
                  eEnergy(1) = eEnergy(1)
     &                       + eEnerSecPt(1)*dV*dA_pt
                  if (keyplL.gt.0) then
                      eEnergy(1) = eEnergy(1)
     &                     + eEnerSecPt(4)*dV*dA_pt
                      eEnergy(4) = eEnergy(4)
     &                     + eEnerSecPt(4)*dV*dA_pt
                  endif
                  if (keycrL.gt.0) then
                      eEnergy(1) = eEnergy(1)
     &                     + eEnerSecPt(5)*dV*dA_pt
                      eEnergy(5) = eEnergy(5)
     &                     + eEnerSecPt(5)*dV*dA_pt
                  endif
                  tGenTh(1) = tGenTh(1)+siglv(1)*dA_pt
                  tGenTh(2) = tGenTh(2)+siglv(1)*z_pt*dA_pt
                  tGenTh(3) = tGenTh(3)-siglv(1)*y_pt*dA_pt
                  
                  tGenStress(1,domInt)=tGenStress(1,domInt)+
     &                 matStress(1,domInt,kspt)*dA_pt
                  tGenStress(2,domInt)=tGenStress(2,domInt)+
     &                 matStress(1,domInt,kspt)*z_pt*dA_pt
                  tGenStress(3,domInt)=tGenStress(3,domInt)-
     &                 matStress(1,domInt,kspt)*y_pt*dA_pt
                  tGenStress(4,domInt)=tGenStress(4,domInt)+
     &                 dA_pt*(-matStress(2,domInt,kspt)*z_pt +
     &                 matStress(3,domInt,kspt)*y_pt)
                  if (iniFlag) then
                      siglvis(1)=siglvis(1)+siglvi(1)*dA_pt
                      siglvis(2)=siglvis(2)+siglvi(1)*z_pt*dA_pt
                      siglvis(3)=siglvis(3)-siglvi(1)*y_pt*dA_pt
                      siglvis(4)=siglvis(4)+dA_pt*
     &                         (-siglvi(2)*z_pt +siglvi(3)*y_pt)
                  endif
                  secStiff(1,1)=
     &                 secStiff(1,1)+matStiff(1,1)*dA_pt
                  secStiff(1,2)=
     &                 secStiff(1,2)+matStiff(1,1)*z_pt*dA_pt
                  secStiff(1,3)=
     &                 secStiff(1,3)-matStiff(1,1)*y_pt*dA_pt
                  secStiff(2,1)=
     &                 secStiff(2,1)+matStiff(1,1)*z_pt*dA_pt
                  secStiff(2,2)=
     &                 secStiff(2,2)+matStiff(1,1)*z_pt**2*dA_pt
                  secStiff(2,3)=
     &                 secStiff(2,3)-matStiff(1,1)*y_pt*z_pt*dA_pt
                  secStiff(3,1)=
     &                 secStiff(3,1)-matStiff(1,1)*y_pt*dA_pt
                  secStiff(3,2)=
     &                 secStiff(3,2)-matStiff(1,1)*z_pt*y_pt*dA_pt
                  secStiff(3,3)=
     &                 secStiff(3,3)+matStiff(1,1)*y_pt**2*dA_pt
                  secStiff(4,4)=secStiff(4,4)+
     &                 HALF*(matStiff(2,2)+matStiff(3,3))*dA_pt*
     &                 (y_pt**2+z_pt**2-z_pt*
     &                 secResults(pphiy,kspt,domIntP)+
     &                 y_pt*secResults(pphiz,kspt,domIntP))
                  secStiff(5,5)=secStiff(5,5)+matStiff(2,2)*dA_pt
                  secStiff(6,6)=secStiff(6,6)+matStiff(3,3)*dA_pt
                  if (isWarp) then
                      tGenStress(7,domInt)=tGenStress(7,domInt)+
     &                     matStress(1,domInt,kspt)*dA_pt*
     &                     secResults(2,kspt,domIntP)
                      secStiff(7,7)=secStiff(7,7)+matStiff(1,1)*dA_pt*
     &                    (secResults(2,kspt,domIntP)**2)
                  endif

                  endif
               enddo
            enddo
         endif

         GIP = secStiff(5,5)
         
         if (i.gt.0 .and. (abs(rvr(5))+abs(rvr(6)).gt.ZERO) ) then
             
             if (rvr(5).gt.ZERO) then
                 scaleSF1 = scaleSF1*rvr(5)
                 secStiff(5,5)=secStiff(5,5)*scaleSF1                 
             elseif (rvr(5).lt.ZERO) then
                 secStiff(5,5) = abs(rvr(5))
             endif

             if (rvr(6).gt.ZERO) then
                 scaleSF2 = scaleSF2*rvr(6)
                 secStiff(6,6)=secStiff(6,6)*scaleSF2                 
             elseif (rvr(6).lt.ZERO) then
                 secStiff(6,6) = abs(rvr(6))
             endif
             
         else
             secStiff(5,5)=secStiff(5,5)*scaleSF1
             secStiff(6,6)=secStiff(6,6)*scaleSF2
         endif
         
         tGenStress(5,domInt)=secStiff(5,5)*tGenStrain(5,domInt)
         tGenStress(6,domInt)=secStiff(6,6)*tGenStrain(6,domInt)
         eEnergy(1)=eEnergy(1)+HALF*dV*(
     &              tGenStress(5,domInt)*tGenStrain(5,domInt)+
     &              tGenStress(6,domInt)*tGenStrain(6,domInt))
         if (isWarp) then
             eEnergy(1)=eEnergy(1)+HALF*dV*
     &            (tGenStress(7,domInt)*tGenStrain(7,domInt))
             fWarp(domInt) = eWarp(domInt) * GIP*
     &             sectSummary(pIpp,domIntP)/sectSummary(pArea,domIntP)
         endif

c --- stiffness matrix
c
         if (kelin(1).eq.1) then
            call cAtxBxA_C (gradMat(1,1),secStiff(1,1),eStiff(1,1),
     &           nGenStr,nGenStr,nUsrDof,nGenStr,nUsrDof,work(1),dV)
            kelout(1) = 1
            if (isWarp) then
               work(1) = GIP*dV*sectSummary(pIpp,domIntP)/
     &                          sectSummary(pArea,domIntP)
               do kk=1,nUsrDof
                  do jj=kk,nUsrDof
                     eStiff(jj,kk)=eStiff(jj,kk)+
     &                    bwarp(kk)*bwarp(jj)*work(1)
                  enddo
               enddo
            endif
         endif
c --- initial stress matrix
c
         if (compSSTIF) then
            lambda2 = lambda**2
            work(1) = tGenStress(1,domInt)*lambda2*dV
            work(7) = tGenStress(1,domInt)*lambda2*dV*
     &           (sectSummary(pIpp,domIntP)/
     &            sectSummary(pArea,domIntP))
            if (keyopt(4) .eq. 1) then
                work(7) = ZERO
            elseif (keyopt(4) .eq. 2) then
                work(1) = ZERO
                work(7) = ZERO
            endif
            work(10)= dV*tGenStress(4,domInt)
            work(11)= dV*tGenStress(2,domInt)*lambda*HALF
            work(12)= dV*tGenStress(3,domInt)*lambda*HALF
            if (keyopt(5).eq.1) then
                work(10) = ZERO
            elseif (keyopt(5).eq.2) then
                work(11) = ZERO
                work(12) = ZERO
            elseif (keyopt(5).eq.3) then
                work(10) = ZERO
                work(11) = ZERO
                work(12) = ZERO
            endif
            work(5) = tGenStress(5,domInt)*lambda2*dV
            work(6) = tGenStress(6,domInt)*lambda2*dV
            work(8) = HALF*dV*tGenStress(5,domInt)
            work(9) = HALF*dV*tGenStress(6,domInt)
            work(13)= dV*tGenStress(5,domInt)*lambda
            work(14)= dV*tGenStress(6,domInt)*lambda
            if (keyopt(6) .eq. 1) then
               work(13)= ZERO
               work(14)= ZERO
            elseif (keyopt(6) .eq. 2) then
               work(8) = ZERO
               work(9) = ZERO
               work(13)= ZERO
               work(14)= ZERO
            elseif (keyopt(6) .eq. 3) then
               work(5) = ZERO
               work(6) = ZERO
               work(8) = ZERO
               work(9) = ZERO
               work(13)= ZERO
               work(14)= ZERO
            endif
            do i=1,3
               do j=1,3
                  n3sp(i,j)=basis(i,3,domInt)*basis(j,1,domInt)+
     &                 basis(i,1,domInt)*basis(j,3,domInt)
                  n2sp(i,j)=basis(i,2,domInt)*basis(j,1,domInt)+
     &                 basis(i,1,domInt)*basis(j,2,domInt)
                  tt(i,j)=n2xn3(i)*n2xn3(j)
                  n3tp(i,j)=basis(i,3,domInt)*n2xn3(j)+
     &                 n2xn3(i)*basis(j,3,domInt)
                  n2tp(i,j)=basis(i,2,domInt)*n2xn3(j)+
     &                 n2xn3(i)*basis(j,2,domInt)
                  n2n3m(i,j)=basis(i,2,domInt)*basis(j,3,domInt)-
     &                 basis(i,3,domInt)*basis(j,2,domInt)
                  n3tm(i,j)=basis(i,3,domInt)*n2xn3(j)-
     &                 n2xn3(i)*basis(j,3,domInt)
                  n2tm(i,j)=basis(i,2,domInt)*n2xn3(j)-
     &                 n2xn3(i)*basis(j,2,domInt)
                  tn3(i,j)=basis(i,3,domInt)*n2xn3(j)
                  tn2(i,j)=basis(i,2,domInt)*n2xn3(j)
               enddo
            enddo
            call v_zero(dxDx(1,1),9)
            do i=1,3
               dxDx(i,i) = work(1)
               do j=1,3
                  dxDx(i,j)=dxDx(i,j)-
     &                 work(5)*n3sp(i,j)-work(6)*n2sp(i,j)
                  dpDp1(i,j)=work(8)*n3tp(i,j)+work(9)*n2tp(i,j)
                  dpDp2(i,j)=work(10)*n2n3m(i,j)+
     &                 work(11)*n3tm(i,j)-work(12)*n2tm(i,j)
                  dxDp(i,j)=-work(13)*tn2(i,j)+work(14)*tn3(i,j)
               enddo
            enddo
            do j=1,numElemNodes
               jj=(j-1)*nDofNode
               work(15)=shpe(j,1)
               work(17)=shpDer(j,domInt)
               do i=j,numElemNodes
                  ii=(i-1)*nDofNode
                  work(16)=shpe(i,1)
                  work(18)=shpDer(i,domInt)
                  do kk=1,3
                     col=jj+kk
                     if (i.ne.j) then
                        do ll=1,3
                           row=ii+ll
                           gStiff(row,col)=gStiff(row,col)+
     &                          work(17)*work(18)*dxDx(ll,kk)
                           gStiff(row+3,col+3)=gStiff(row+3,col+3)+
     &                          work(15)*work(16)*dpDp1(ll,kk)+
     &                          (work(18)*work(15)-
     &                          work(17)*work(16))*dpDp2(ll,kk)+
     &                          work(7)*work(17)*work(18)*tt(ll,kk)
                        enddo
                     else
                        do ll=kk,3
                           row=ii+ll
                           gStiff(row,col)=gStiff(row,col)+
     &                          work(17)*work(18)*dxDx(ll,kk)
                           gStiff(row+3,col+3)=gStiff(row+3,col+3)+
     &                          work(15)*work(16)*dpDp1(ll,kk)+
     &                          (work(18)*work(15)-
     &                          work(17)*work(16))*dpDp2(ll,kk)+
     &                          work(7)*work(17)*work(18)*tt(ll,kk)
                        enddo
                     endif
                     do ll=1,3
                        row=ii+ll
                        gStiff(row+3,col)=gStiff(row+3,col)+
     &                       work(16)*work(17)*dxDp(kk,ll)
                        if (i.gt.j) then
                           gStiff(row,col+3)=gStiff(row,col+3)+
     &                          work(15)*work(18)*dxDp(ll,kk)
                        endif
                     enddo
                  enddo
               enddo
            enddo
         endif

c --- newton-raphson restoring force
         if (kelin(6).eq.1) then
            call v_move(tGenStress(1,domInt), work(1), nGenStr)
            call cv1_v1(work(1),nGenStr,dV)
            call Atxv1_v2(gradMat(1,1),work(1),fInt(1),nGenStr,-nUsrDof)
            if (isWarp) then
               do ii=1,nUsrDof
                  fInt(ii)=fInt(ii)+bwarp(ii)*fWarp(domInt)*dV
               enddo
            endif
            kelout(6) = 1
         endif
         if (iniFlag.and.A_KeyFirstPass.eq.1) then
            call v_move (siglvis(1), work(1), nGenStr)
            call cv1_v1(work(1),nGenStr,dV)
            call Atxv1_v2 (gradMat(1,1),work(1),
     &           iniForce(1),nGenStr,-nUsrDof)
            if (isWarp) then
               do ii=1,nUsrDof
                  iniForce(ii)=iniForce(ii)+bwarp(ii)*fWarp(domInt)*dV
               enddo
            endif
         endif

c --- hygro load vector
         if (kelin(5).eq.1) then
            if (khygro.eq.1 .and. A_KeyOut.eq.0) then
               call v_move (tGenTh(1), work(1),nGenStr)
               call cv1_v1(work(1),nGenStr,dV)
               call Atxv1_v2 (gradMat(1,1),work(1),fExt(1),nGenStr,
     &                        -nUsrDof)
               kelout(5)=1
            endif
         endif

      enddo


      if (A_KeyOut .ne. 1 .and. iniFlag) then
         if (A_LoadStepNum .eq. 1) then
            factor =  1.d0 - A_CurLoadFact
            call v1pcv2_v1 (fExt(1),iniForce(1),nUsrDof,factor)
            kelout(5) = 1
         endif
      endif

      do j=1,3
         elemCG(j)=elemCG(j)/elemVolume
      enddo

c --- calculate basis directions at a higher order integration rule

      do domInt = 1,numDomainIntPts+1
         domIntP = 1
         if (numElemNodes.eq.2) then
            shpe(1,1)=HALF*(ONE-Gauss2(domInt))
            shpe(2,1)=HALF*(ONE+Gauss2(domInt))
            shpe(1,2)=-HALF
            shpe(2,2)= HALF
         elseif(numElemNodes.eq.3) then
            shpe(1,1)=HALF*(Gauss3(domInt)-ONE)*Gauss3(domInt)
            shpe(2,1)=HALF*(Gauss3(domInt)+ONE)*Gauss3(domInt)
            shpe(3,1)=ONE-Gauss3(domInt)*Gauss3(domInt)
            shpe(1,2)=Gauss3(domInt)-HALF
            shpe(2,2)=Gauss3(domInt)+HALF
            shpe(3,2)=-TWO*Gauss3(domInt)
         endif
         dL = ZERO
         do j=1,3
            dxdg(j) = ZERO
            do nod=1,numElemNodes
               dxdg(j) = dxdg(j)+shpe(nod,2)*xCur(j,nod)
            enddo
            dL=dL+dxdg(j)*dxdg(j)
         enddo
         dL    = sqrt(dL)
         dLinv = ONE/dL
         hbasis(1,1,domInt)=dxdg(1)*dLinv
         hbasis(2,1,domInt)=dxdg(2)*dLinv
         hbasis(3,1,domInt)=dxdg(3)*dLinv
         if (eConn(numElemNodes+1).eq.0) then
            if(A_KeyFirstPass.eq.1) then
               call TransBeam (dxdg(1),trv(1,1))
               hbasis(1,2,domInt)=trv(2,1)
               hbasis(2,2,domInt)=trv(2,2)
               hbasis(3,2,domInt)=trv(2,3)
               hbasis(1,3,domInt)=trv(3,1)
               hbasis(2,3,domInt)=trv(3,2)
               hbasis(3,3,domInt)=trv(3,3)
               call v_move(hbasis(1,1,domInt),hbasis_t(1,1,domInt),9)
            endif
         elseif (A_KeyFirstPass.eq.1) then
            j = numElemNodes+1
            if(numElemNodes .eq. 2) then
               work(1)=xBeg(1,j)-HALF*(xBeg(1,1)+xBeg(1,2))
               work(2)=xBeg(2,j)-HALF*(xBeg(2,1)+xBeg(2,2))
               work(3)=xBeg(3,j)-HALF*(xBeg(3,1)+xBeg(3,2))
            else
               work(1)=xBeg(1,j)-xBeg(1,3)
               work(2)=xBeg(2,j)-xBeg(2,3)
               work(3)=xBeg(3,j)-xBeg(3,3)
            endif
c --- beam orientation using Kth node
            work(4)=work(1)*work(1)+work(2)*work(2)+work(3)*work(3)
            work(4)=ONE/sqrt(work(4))
            work(1)=work(4)*work(1)
            work(2)=work(4)*work(2)
            work(3)=work(4)*work(3)
            hbasis(1,2,domInt)=
     &           -hbasis(2,1,domInt)*work(3)+hbasis(3,1,domInt)*work(2)
            hbasis(2,2,domInt)=
     &           -hbasis(3,1,domInt)*work(1)+hbasis(1,1,domInt)*work(3)
            hbasis(3,2,domInt)=
     &           -hbasis(1,1,domInt)*work(2)+hbasis(2,1,domInt)*work(1)
c
            work(4)=
     &           hbasis(1,2,domInt)*hbasis(1,2,domInt)+
     &           hbasis(2,2,domInt)*hbasis(2,2,domInt)+
     &           hbasis(3,2,domInt)*hbasis(3,2,domInt)
            work(4)=ONE/sqrt(work(4))
            hbasis(1,2,domInt)=work(4)*hbasis(1,2,domInt)
            hbasis(2,2,domInt)=work(4)*hbasis(2,2,domInt)
            hbasis(3,2,domInt)=work(4)*hbasis(3,2,domInt)
c
            hbasis(1,3,domInt)=
     &           hbasis(2,1,domInt)*hbasis(3,2,domInt)-
     &           hbasis(3,1,domInt)*hbasis(2,2,domInt)
            hbasis(2,3,domInt)=
     &           hbasis(3,1,domInt)*hbasis(1,2,domInt)-
     &           hbasis(1,1,domInt)*hbasis(3,2,domInt)
            hbasis(3,3,domInt)=
     &           hbasis(1,1,domInt)*hbasis(2,2,domInt)-
     &           hbasis(2,1,domInt)*hbasis(1,2,domInt)

            work(4)=
     &           hbasis(1,3,domInt)*hbasis(1,3,domInt)+
     &           hbasis(2,3,domInt)*hbasis(2,3,domInt)+
     &           hbasis(3,3,domInt)*hbasis(3,3,domInt)
            work(4)=ONE/sqrt(work(4))
            hbasis(1,3,domInt)=work(4)*hbasis(1,3,domInt)
            hbasis(2,3,domInt)=work(4)*hbasis(2,3,domInt)
            hbasis(3,3,domInt)=work(4)*hbasis(3,3,domInt)
            call v_move(hbasis(1,1,domInt),hbasis_t(1,1,domInt),9)
         endif

         if (WetMode) then
             call v_move(hbasis_t(1,1,domInt),hbasis(1,1,domInt),9)
         else   
c --- interpolated incremental rotation vector
             work(1) = ZERO
             work(2) = ZERO
             work(3) = ZERO
             do nod = 1,numElemNodes
                j = (nod-1)*nDofNode
                work(1)=work(1)+shpe(nod,1)*incDisp(j+4)
                work(2)=work(2)+shpe(nod,1)*incDisp(j+5)
                work(3)=work(3)+shpe(nod,1)*incDisp(j+6)
             enddo
             call BasRot(hbasis_t(1,2,domInt),hbasis(1,2,domInt),
     &                   work(1),2)
          endif
      enddo
c --- account for element added mass per unit length
      if (kelin(2).eq.1) then
          domIntP = 1
          cDensity = density
      else
          cDensity=0.0d0
      endif
      elemMass=elemVolume*cDensity
      if (cDensity.gt.ZERO.and.kelin(2).eq.1) then
          call BeamMass(numDomainIntPts, numElemNodes, nDofNode,
     &                  Gauss2(1), Gauss3(1), xBeg(1,1), Gaussw2(1),
     &                  Gaussw3(1), hbasis(1,1,1), sectSummary(1,1),
     &                  pIzz, pIyy, pIyz, pArea, nGenStr, cDensity,
     &                  density, moffset(1,1,1), A_MassLump, 0,
     &                  nUsrDof, secStiff(1,1), gradMat(1,1),
     &                  work(1), eMass(1,1))
         kelout(2)=1
      endif
      if (compSSTIF)then
         call m_symLtoU(gStiff(1,1),nUsrDof,nUsrDof)
      endif
      if (kpress.eq.1) then
         if (kelin(5).eq.1) then
            if (pressure(13).ne.ZERO) then
               do ii=1,3
                  fExt(ii) = fExt(ii)+nBasis(ii,1,1)*pressure(13)
               enddo
            endif
            if (pressure(14).ne.ZERO) then
               do ii=1,3
                  fExt(ii+nDofNode)=fExt(ii+nDofNode)
     &                 -nBasis(ii,1,2)*pressure(14)
               enddo
            endif
         endif
         do ii = 1,3
            jj = (ii-1)*4
            if (pressure(jj+1).ne.ZERO .or.
     &           pressure(jj+2).ne.ZERO)then
               do domInt=1,numDomainIntPts+1
                  domIntP = 1
                  if (numElemNodes.eq.2) then
                     shpe(1,1)=HALF*(ONE-Gauss2(domInt))
                     shpe(2,1)=HALF*(ONE+Gauss2(domInt))
                     shpe(1,2)=-HALF
                     shpe(2,2)= HALF
                  elseif(numElemNodes.eq.3) then
                     shpe(1,1)=HALF*(Gauss3(domInt)-ONE)*Gauss3(domInt)
                     shpe(2,1)=HALF*(Gauss3(domInt)+ONE)*Gauss3(domInt)
                     shpe(3,1)=ONE-Gauss3(domInt)*Gauss3(domInt)
                     shpe(1,2)=Gauss3(domInt)-HALF
                     shpe(2,2)=Gauss3(domInt)+HALF
                     shpe(3,2)=-TWO*Gauss3(domInt)
                  endif
                  dL = ZERO
                  do j=1,3
                     dxdg(j) = ZERO
                     do nod=1,numElemNodes
                        dxdg(j) = dxdg(j)+shpe(nod,2)*xCur(j,nod)
                     enddo
                     dL=dL+dxdg(j)*dxdg(j)
                  enddo
                  dL=sqrt(dL)
                  if (numElemNodes.eq.2) then
                      dV = dL*Gaussw2(domInt)
                  elseif (numElemNodes.eq.3) then
                      dV = dL*Gaussw3(domInt)
                  endif
                  dL = ONE/dL
                  dxdg(1)=dxdg(1)*dL
                  dxdg(2)=dxdg(2)*dL
                  dxdg(3)=dxdg(3)*dL
                  n1(1) = hbasis(1,2,domInt)
                  n1(2) = hbasis(2,2,domInt)
                  n1(3) = hbasis(3,2,domInt)
                  n2(1) = hbasis(1,3,domInt)
                  n2(2) = hbasis(2,3,domInt)
                  n2(3) = hbasis(3,3,domInt)
                  do nod = 1,numElemNodes
                     tDer(nod,1)=shpe(nod,1)
                     tDer(nod,2)=shpe(nod,2)*dL
                  enddo
                  if (numElemNodes.eq.2) then
                     work(1)=HALF*(ONE-Gauss2(domInt))
                     work(2)=HALF*(ONE+Gauss2(domInt))
                  elseif(numElemNodes.eq.3)then
                     work(1)=HALF*(ONE-Gauss3(domInt))
                     work(2)=HALF*(ONE+Gauss3(domInt))
                  endif
                  work(3) = (pressure(jj+1)*work(1) +
     &                 pressure(jj+2)*work(2))*dV
                  if (ii.eq.1) then
                     work(4)=-n2(1)*work(3)
                     work(5)=-n2(2)*work(3)
                     work(6)=-n2(3)*work(3)
                     work(7)=-n1(1)*work(3)
                     work(8)=-n1(2)*work(3)
                     work(9)=-n1(3)*work(3)
                  elseif(ii.eq.2) then
                     work(4)=-n1(1)*work(3)
                     work(5)=-n1(2)*work(3)
                     work(6)=-n1(3)*work(3)
                     work(7)= n2(1)*work(3)
                     work(8)= n2(2)*work(3)
                     work(9)= n2(3)*work(3)
                  else
                     work(4)=dxdg(1)*work(3)
                     work(5)=dxdg(2)*work(3)
                     work(6)=dxdg(3)*work(3)
                  endif
                  if (kelin(5).eq.1) then
                     do nod=1,numElemNodes
                        kk=(nod-1)*nDofNode
                        fExt(kk+1) = fExt(kk+1) + shpe(nod,1)*work(4)
                        fExt(kk+2) = fExt(kk+2) + shpe(nod,1)*work(5)
                        fExt(kk+3) = fExt(kk+3) + shpe(nod,1)*work(6)
                     enddo
                  endif
                  if (ii.le.2 .and. WetMode) then
                      call PresLoadStiff(nUsrDof,numElemNodes,
     &                                   nDofNode,3,0,work(7),tDer(1,1),
     &                                   work(11),gStiff(1,1))
                  endif
               enddo
            endif
         enddo
         kelout(5) = 1
      endif

c --- finish off matrices
c
      if (kelin(1).eq.1)  then
          call m_symLtoU (eStiff(1,1),nUsrDof,nUsrDof)
      endif
      if (compSSTIF) then
          call v1pv2_v1(eStiff(1,1),gStiff(1,1),nUsrDof*nUsrDof)
          call v_zero(gStiff(1,1),nUsrDof*nUsrDof)
      endif

c --- update past boundary conditions
c --- set bc at beginning of next load step to end of this step
      if (A_KeyBCsUpd.eq.1)then
          if (khygro.eq.1) 
     &        call v_move (moistcur(1),moistbeg(1),3*numElemNodes)
          if (kpress.eq.1) 
     &        call v_move (pressure(1),prsBeg(1),14)
      endif
c --- set fields at beginning of next substep to end of this substep
      if (A_KeyMoistUpd.eq.1 .and. khygro.eq.1)
     &    call v_move (moistcur(1),moistprev(1),3*numElemNodes)

      if (A_KeyMatUpd.eq.1) then
          call v_move(basis(1,1,1), basis_t(1,1,1), 18)
          call v_move(hbasis(1,1,1), hbasis_t(1,1,1),27)
          call v_move(tGenStrain(1,1),tGenStrain_t(1,1),2*maxGenStr)
          call v_move(nBasis(1,1,1),nBasis_t(1,1,1),27)
          call v_move(yzScale(1),yzScale_t(1),2)
      endif

      if (A_KeyStab.gt.0) then
          call StabDriver (elemId, elChar(1), numDomainIntPts,
     &                     numElemNodes, nDofNode,
     &                     Gauss2(1), Gauss3(1), xBeg(1,1),
     &                     Gaussw2(1), Gaussw3(1), hbasis(1,1,1),
     &                     sectSummary(1,1), pIzz, pIyy, pIyz, pArea,
     &                     nGenStr, moffset(1,1,1), A_MassLump,
     &                     nUsrDof, work(1), eleDisp(1,1),kelout(1),
     &                     eStiff(1,1),fInt(1),eEnergy(1),
     &                     basicSaveVars(1), numBasicSaveVars,
     &                     secStiff(1,1),gradMat(1,1))
      endif

c --- handle SVR files
      call StoreSVR (stateVarIndex(1),1,numBasicSaveVars,
     &               basicSaveVars(1))
      
      if (LenMatRecMax.gt.0 .and. keyopt(1).ne.2)
     &    call UpdSvrRec (1, elemId, stateVarIndex(1), LenMatRecMax)

      call PutSVRIndx(elemId,locsvrL,stateVarIndex(1))

c --- output results

      if (A_KeyOut.eq.1) then
         if (isWarp) then
            do ii=1,numDomainIntPts
               tGenStress(4,ii)=tGenStress(4,ii)+fWarp(ii)
            enddo
         endif
         call v_zero(thStrSecPt(1,1),2*numSectPts)
         if (khygro .gt. 0) then
            if (elasticSection.eq.1) then
               do jj=1,numDomainIntPts
                  do kk=1,numSectPts
                     k1=  (jj-1)*numSectPts+kk
                     thStrSecPt(jj,kk)=thStrain(jj)
                  enddo
               enddo
            else
               do jj=1,numDomainIntPts
                  do kk=1,numSectPts
                     k1=  (jj-1) * numSectPts + kk
                     iNumTerms = 1
                     call RetMatRec (W_iStrainHyg, k1,
     &                               iNumTerms, thStrSecPt(jj,kk))
                  enddo
               enddo
            endif
         endif
         if (elasticSection.eq.0) then
             if (mRecIntPidx(W_iStrainPl).gt.0) then
                 call v_zero(ple (1,1,1),6*numSectPts)
                 call v_zero(peeq(1,1,1),6*numSectPts)
                 call v_zero(pwrk(1,1,1),6*numSectPts)
                 do jj=1,numDomainIntPts
                    do kk=1,numSectPts
                       k1 = (jj-1) * numSectPts + kk
                       iNumTerms = 3
                       call RetMatRec(W_iStrainPl, k1,
     &                                iNumTerms, ple(1,jj,kk) )
                       iNumTerms = 2
                       call RetMatRec(W_iEqvPl, k1,
     &                                iNumTerms , dEqPl )
                       call v_move (dEqPl(1),peeq(1,jj,kk),1)
                       call v_move (dEqPl(2),pwrk(1,jj,kk),1)
                       call v_move (dEqPl(1),peeq(2,jj,kk),1)
                       call v_move (dEqPl(2),pwrk(2,jj,kk),1)
                       call v_move (dEqPl(1),peeq(3,jj,kk),1)
                       call v_move (dEqPl(2),pwrk(3,jj,kk),1)
                    enddo
                 enddo
             endif
c --- retrive creep strain
             if (mRecIntPidx(W_iStrainCr).gt.0) then
                 do jj=1,numDomainIntPts
                    do kk=1,numSectPts
                       k1 = (jj-1) * numSectPts + kk
                       iNumTerms = 3
                       call RetMatRec(W_iStrainCr, k1,
     &                                iNumTerms, cre(1,jj,kk))
                       iNumTerms = 2
                       call RetMatRec (W_iEqvCr,k1,iNumTerms,dEqCr)
                       call v_move (dEqCr(1),creq(1,jj,kk),1)
                       call v_move (dEqCr(2),crwk(1,jj,kk),1)
                       call v_move (dEqCr(1),creq(2,jj,kk),1)
                       call v_move (dEqCr(2),crwk(2,jj,kk),1)
                       call v_move (dEqCr(1),creq(3,jj,kk),1)
                       call v_move (dEqCr(2),crwk(3,jj,kk),1)
                    enddo
                 enddo
             endif
         elseif (elasticSection.eq.1) then
            do domInt = 1,numDomainIntPts
               domIntP = 1
               if (numElemNodes.eq.2) then
                   tInt(1)=HALF
                   tInt(2)=HALF
               elseif(numElemNodes.eq.3)then
                   tInt(1)=HALF*(ONE-Gauss2(domInt))
                   tInt(2)=HALF*(ONE+Gauss2(domInt))
               endif
               tDomInt=tInt(1)*moistcur(1)+tInt(2)*moistcur(4)
               call GetElasticProps (elemId,matId,tDomInt,prop(1))
               YoungsMod = prop(1)
               ShearMxy  = prop(5)
               ShearMyz  = prop(6)
               ShearMxz  = prop(4)
               alpha     = prop(10)
               density   = prop(13)
               do kspt=1,numSectPts
                  if(secResults(pM,kspt,domIntP).gt.ZERO) then
                     y_pt=secResults(pY,kspt,domIntP)*yzScale(domInt)
                     z_pt=secResults(pZ,kspt,domIntP)*yzScale(domInt)

                     twork(1)=y_pt-sectSummary(pSy,domIntP)+
     &                    sectSummary(pCgy,domIntP)
                     twork(2)=z_pt-sectSummary(pSz,domIntP)+
     &                    sectSummary(pCgz,domIntP)
                     phiy_pt=(secResults(pphiy,kspt,domIntP)-twork(2))*
     &                    yzScale(domInt)
                     phiz_pt=(secResults(pphiz,kspt,domIntP)+twork(1))*
     &                    yzScale(domInt)
                     matStrain(1,domInt,kspt)=
     &                    tGenStrain(1,domInt)+
     &                    z_pt*tGenStrain(2,domInt)-
     &                    y_pt*tGenStrain(3,domInt)
                     if (isWarp) then
                        matStrain(1,domInt,kspt)=
     &                       matStrain(1,domInt,kspt)+
     &                       tGenStrain(7,domInt)*
     &                       secResults(2,kspt,domIntP)
                     endif

c --- change from stress function to warping function (if needed)

                     matStrain(2,domInt,kspt)=
     &                    tGenStrain(4,domInt)*phiy_pt
                     matStrain(3,domInt,kspt)=
     &                    tGenStrain(4,domInt)*phiz_pt
                     matStress(1,domInt,kspt) =
     &                    YoungsMod * matStrain(1,domInt,kspt)
                     matStress(2,domInt,kspt) =
     &                    ShearMyz  * matStrain(2,domInt,kspt)
                     matStress(3,domInt,kspt) =
     &                    ShearMxz  * matStrain(3,domInt,kspt)

c --- add hygro strain only for output of total strain
                     if (khygro.eq.1) then
                         if (A_KeyFirstIter.ne.1) then
                              matStrain(1,domInt,kspt)=
     &                                  matStrain(1,domInt,kspt) +
     &                                  thStrain(domInt)
                        endif
                     endif
                  endif
               enddo
            enddo
         endif
c --- calculate shear stresses for output
         do domInt = 1,numDomainIntPts
               domIntP = 1
               do kspt=1,numSectPts
                  if(secResults(pM,kspt,domIntP).gt.ZERO) then
                     modWt = secResults(pM,kspt,domIntP)/
     &                       secResults(pdA,kspt,domIntP)
                     matStress(2,domInt,kspt)=
     &                 matStress(2,domInt,kspt)+
     &                 modWt*(
     &                 tGenStress(5,domInt)*secResults(9,kspt,domIntP)+
     &                 tGenStress(6,domInt)*secResults(11,kspt,domIntP))
                     matStress(3,domInt,kspt)=
     &                 matStress(3,domInt,kspt)+
     &                 modWt*(
     &                 tGenStress(5,domInt)*secResults(8,kspt,domIntP)+
     &                 tGenStress(6,domInt)*secResults(10,kspt,domIntP))
                  endif
               enddo
         enddo
         if (khygro.eq.1.and.elasticSection.eq.1) then
            if (A_KeyFirstIter.ne.1) then
                tGenStrain(1,1)=
     &              tGenStrain(1,1)+thStrain(1)
                tGenStrain(1,2)=
     &              tGenStrain(1,2)+thStrain(2)
            endif
         endif
         if (btest(eomask,W_OutElIndx)) then
            if (btest(eomask,W_OutMisc)) then
#if defined(DEBUG)
               call v_move(secDataR(1,1),twork(1),secDataI(pNodes)*2)
               do ii=1,secDataI(pNodes)
                  twork(offY+ii-1)=twork(offY+ii-1)
     &                   -sectSummary(pCgy,1)
                  twork(offZ+ii-1)=twork(offZ+ii-1)
     &                   -sectSummary(pCgz,1)
               enddo
               yzmxmn(1,1)= -1.0d15
               yzmxmn(2,1)=  1.0d15
               yzmxmn(3,1)= -1.0d15
               yzmxmn(4,1)=  1.0d15
               yzmxmn(1,2)=yzmxmn(1,1)
               yzmxmn(2,2)=yzmxmn(2,1)
               yzmxmn(3,2)=yzmxmn(3,1)
               yzmxmn(4,2)=yzmxmn(4,1)
               do ii=1,secDataI(pNodes)
                  yzmxmn(1,1)=max(yzmxmn(1,1),twork(offY+ii-1))
                  yzmxmn(2,1)=min(yzmxmn(2,1),twork(offY+ii-1))
                  yzmxmn(3,1)=max(yzmxmn(3,1),twork(offZ+ii-1))
                  yzmxmn(4,1)=min(yzmxmn(4,1),twork(offZ+ii-1))
               enddo
               yzmxmn(1,2)=yzmxmn(1,1)
               yzmxmn(2,2)=yzmxmn(2,1)
               yzmxmn(3,2)=yzmxmn(3,1)
               yzmxmn(4,2)=yzmxmn(4,1)
#endif
c --- miscellaneous data extrapolated to element nodes
c
               if (numElemNodes.eq.2) then
                   work(101) = ONE
                   work(102) = ONE
                   work(103) = ONE
                   work(104) = ONE
                   kof = HALF
               else
                   if (keyplL.eq.0) then
                       work(101) = HALF*(ONE+sqrt(THREE))
                       work(102) = HALF*(ONE-sqrt(THREE))
                   else
                       work(101) = 1.0d0
                       work(102) = 0.0d0
                   endif
                   work(103) = work(102)
                   work(104) = work(101)
                   work(105) = HALF
                   work(106) = HALF
                   kof = 1/THREE
               endif
               call v_zero(work(1),numElemNodes*(2*nGenStr+1))
               do nod = 1,numElemNodes
                  ii=(nod-1)*13
                  kk=(nod-1)*2+100
                  do domInt=1,numDomainIntPts
c --- resultant data
                     do jj=1,6
                        work(ii+jj)=work(ii+jj)+
     &                       work(kk+domInt)*tGenStress(jj,domInt)
                        work(ii+jj+6)=work(ii+jj+6)+
     &                       work(kk+domInt)*tGenStrain(jj,domInt)
                     enddo
                     work(ii+13)=work(ii+13)+
     &                    work(kk+domInt)*yzScale(domInt)
#if defined(DEBUG)
c --- bending data
                     if(nGenStr.eq.7) then
                        jj=26+(nod-1)*2
                        work(jj+1)=work(jj+1)+
     &                       work(kk+domInt)*tGenStress(7,domInt)
                        work(jj+2)=work(jj+2)+
     &                       work(kk+domInt)*tGenStrain(7,domInt)
                     endif
#endif
                  enddo
               enddo
c --- cross section data
               call v_zero (miscdata(28),3)
               do nod = 1,numElemNodes
                  domIntP = 1
                  ii = (nod-1)*13
                  aScale = work(ii+13)**2
                  work(ii+13) = sectSummary(pArea,domIntP)*aScale
#if !defined(DEBUG)               
c --- stretching energy (averaged over nodes, thus may not be 100% accurate)
                  miscdata(28) = miscdata(28) + 
     &              kof*HALF*L/(YoungsMod*work(ii+13))*
     &              work(ii+1)**2
c --- bending energy
                  miscdata(29) = miscdata(29) + 
     &              kof*HALF*L/YoungsMod*(
     &              work(ii+2)**2/(sectSummary(pIyy,domIntP)*aScale**2)+
     &              work(ii+3)**2/(sectSummary(pIzz,domIntP)*aScale**2))
c --- torsion energy
                  miscdata(30) = miscdata(30) + 
     &              kof*HALF*L/(sectSummary(pTorq,domIntP)*aScale**2)*
     &              work(ii+4)**2/ShearMxy
#endif               
               enddo
               
c --- the missilanious variable is the total number of
c     section intg. forces which are stored
c     as the first set of numbers on the result file
c     under the record EDEMS
c
               if (numElemNodes.eq.2) nmisc = 30
               if (numElemNodes.eq.3) nmisc = 45
               nmisc = nmisc + 30
#if defined(DEBUG)               
               call v_move(work(1),miscdata(1),30)
               semxmn( 1)=work( 1)/sectSummary(pArea,1)
               semxmn( 2)=-work( 3)*yzmxmn(1,1)/sectSummary(pIzz,1)
               semxmn( 3)=-work( 3)*yzmxmn(2,1)/sectSummary(pIzz,1)
               semxmn( 4)=work( 2)*yzmxmn(3,1)/sectSummary(pIyy,1)
               semxmn( 5)=work( 2)*yzmxmn(4,1)/sectSummary(pIyy,1)
               semxmn( 6)=work(14)/sectSummary(pArea,1)
               semxmn( 7)=-work(16)*yzmxmn(1,2)/sectSummary(pIzz,1)
               semxmn( 8)=-work(16)*yzmxmn(2,2)/sectSummary(pIzz,1)
               semxmn( 9)=work(15)*yzmxmn(3,2)/sectSummary(pIyy,1)
               semxmn(10)=work(15)*yzmxmn(4,2)/sectSummary(pIyy,1)

               semxmn(11)=work( 7)
               semxmn(12)=-work( 9)*yzmxmn(1,1)
               semxmn(13)=-work( 9)*yzmxmn(2,1)
               semxmn(14)=work( 8)*yzmxmn(3,1)
               semxmn(15)=work( 8)*yzmxmn(4,1)
               semxmn(16)=work(20)
               semxmn(17)=-work(22)*yzmxmn(1,2)
               semxmn(18)=-work(22)*yzmxmn(2,2)
               semxmn(19)=work(21)*yzmxmn(3,2)
               semxmn(20)=work(21)*yzmxmn(4,2)
               do jj=1,20
                  miscdata(30+jj)=semxmn(jj)
               enddo
               do jj=1,9
                  miscdata(50+jj)=moistcur(jj)
               enddo
               call v_zero(miscdata(60),nmisc-59)
#else
c --- move data only for two nodes (2*13 = 26 variables)
               call v_move(work(1),miscdata(1),26)
c --- store energy (not stored during debug phase)
               miscdata(27) = eEnergy(1)
c --- shear energy (hard to calculate otherwise)
               miscdata(31) = eEnergy(1) - eEnergy(4) - eEnergy(5)
     &                      - miscdata(28) - miscdata(29) - miscdata(30)
c --- plastic energy
               miscdata(32) = eEnergy(4)
c --- creep energy
               miscdata(33) = eEnergy(5)
c --- element length
               miscdata(34) = L
               
               call v_zero(miscdata(35),nmisc-34)
#endif
            endif
         endif
            ii = 6 * numSectNodes * 2
            call v_zero(nResults(1,1,1),ii)
            call v_zero(thStrSNDP(1,1),numSectNodes*2)
            call v_zero(work(1),secDataI(pNodes))
            do ii=1,secDataI(pElems)
               if (secDataR(2*secDataI(pNodes)+ii,1).gt.ZERO) then
                   j=pElemConn+(ii-1)*9
                   do i=1,4
                      nod = secDataI(j+i-1)
                      work(nod)=work(nod)+ONE
                   enddo
               endif
            enddo
            do ii=1,secDataI(pNodes)
               nCrnr(ii) = 0
               if (work(ii).gt.ZERO) then
                   work(ii) = ONE/work(ii)
               endif
            enddo
            do domInt=1,numDomainIntPts
               do ii=1,secDataI(pElems)
                  j=pElemConn+(ii-1)*9
                  jj=(ii-1)*4
                  do i=1,4
                     nod = secDataI(j+i-1)
                     nCrnr(nod) = 1
                     do ll=1,3
                        if (keyplL.eq.0) then
                           nResults(ll,domInt,nod)=
     &                          nResults(ll,domInt,nod)+
     &                          matStress(ll,domInt,jj+1)*xPolate(1,i)+
     &                          matStress(ll,domInt,jj+2)*xPolate(2,i)+
     &                          matStress(ll,domInt,jj+3)*xPolate(3,i)+
     &                          matStress(ll,domInt,jj+4)*xPolate(4,i)
                           nResults(ll+3,domInt,nod)=
     &                          nResults(ll+3,domInt,nod)+
     &                          matStrain(ll,domInt,jj+1)*xPolate(1,i)+
     &                          matStrain(ll,domInt,jj+2)*xPolate(2,i)+
     &                          matStrain(ll,domInt,jj+3)*xPolate(3,i)+
     &                          matStrain(ll,domInt,jj+4)*xPolate(4,i)
                        else
                           nResults(ll,domInt,nod)=
     &                          nResults(ll,domInt,nod)+
     &                          matStress(ll,domInt,jj+1)*mvnd(1,i)+
     &                          matStress(ll,domInt,jj+2)*mvnd(2,i)+
     &                          matStress(ll,domInt,jj+3)*mvnd(3,i)+
     &                          matStress(ll,domInt,jj+4)*mvnd(4,i)
                           nResults(ll+3,domInt,nod)=
     &                          nResults(ll+3,domInt,nod)+
     &                          matStrain(ll,domInt,jj+1)*mvnd(1,i)+
     &                          matStrain(ll,domInt,jj+2)*mvnd(2,i)+
     &                          matStrain(ll,domInt,jj+3)*mvnd(3,i)+
     &                          matStrain(ll,domInt,jj+4)*mvnd(4,i)
                        endif
                     enddo
                     if (khygro.gt.0) then
                        thStrSNDP(domInt,nod)=
     &                       thStrSNDP(domInt,nod)+
     &                       thStrSecPt(domInt,jj+1)*xPolate(1,i)+
     &                       thStrSecPt(domInt,jj+2)*xPolate(2,i)+
     &                       thStrSecPt(domInt,jj+3)*xPolate(3,i)+
     &                       thStrSecPt(domInt,jj+4)*xPolate(4,i)
                     endif
                  enddo
               enddo
            enddo
            do domInt=1,numDomainIntPts
               do ii=1,secDataI(pNodes)
                  if (nCrnr(ii).eq.1) then
                     do ll=1,6
                        nResults(ll,domInt,ii)=
     &                       nResults(ll,domInt,ii)*work(ii)
                     enddo
                     if (khygro.gt.0) then
                         thStrSNDP(domInt,ii)=
     &                        thStrSNDP(domInt,ii)*work(ii)
                     endif
                  endif
               enddo
            enddo

            if (btest(eomask,W_OutElIndx)) then
               call v_move(secDataR(1,1),twork(1),secDataI(pNodes)*2)
               do ii=1,secDataI(pNodes)
                  twork(offY+ii-1)=twork(offY+ii-1)-yOffset(1)
                  twork(offZ+ii-1)=twork(offZ+ii-1)-zOffset(1)
               enddo
               jj = secDataI(pNodes)*2+secDataI(pElems)
               do ii=1,secDataI(pNodes)
                  work(ii)=secDataR(jj+ii,1)
               enddo
               call ResuToNodes(
     &              numElemNodes,numDomainIntPts,secDataI(pNodes),
     &              nDofNode,twork(1),offY,offZ,nResults(1,1,1),eU(1),
     &              nBasis(1,1,1),nBasis_0(1,1,1),
     &              yzScale(1),pResults(1,1,1),nCrnr(1),
     &              khygro, thStrSNDP(1,1),hyStrSecNd(1,1),
     &              keyplL,work(1),lenSecDataR)

c --- gather all the corner node stresses and strains to be
c     placed for post processing.

               NumRstNodes = 0
               NumElems =  numSectPts / numIntPts ! change that
               do ibnode = 1, NumElems
                 do inode = 1, 4
                     NodeNum = secDataI(7+4+(ibnode-1)*9+inode) ! change that
                     do i = 1, NumRstNodes
                        if (NodeNum .eq. RstNodeList(i)) goto 888
                     enddo
                     NumRstNodes = NumRstNodes + 1
                     RstNodeList(NumRstNodes) = NodeNum
 888             enddo
               enddo

               i = 1
               do ibnode = 1, numElemNodes
                  do inode = 1, NumRstNodes
c --- gather the secton node stresses
                     postdata(i) = pResults(4,ibnode,RstNodeList(inode))
                     postdata(i+1) =
     &                      pResults(5,ibnode,RstNodeList(inode))
                     postdata(i+2) =
     &                      pResults(6,ibnode,RstNodeList(inode))

c --- gather the secton node strains
                     postdata(numElemNodes*NumRstNodes*3+i) =
     &                      pResults(7,ibnode,RstNodeList(inode))
                     postdata(numElemNodes*NumRstNodes*3+i+1) =
     &                      pResults(8,ibnode,RstNodeList(inode))
                     postdata(numElemNodes*NumRstNodes*3+i+2) =
     &                      pResults(9,ibnode,RstNodeList(inode))

c --- gather the secton node deformed locations
                     postdata(numElemNodes*NumRstNodes*3*2+i) =
     &                     pResults(1,ibnode,RstNodeList(inode))
                     postdata(numElemNodes*NumRstNodes*3*2+i+1) =
     &                     pResults(2,ibnode,RstNodeList(inode))
                     postdata(numElemNodes*NumRstNodes*3*2+i+2) =
     &                     pResults(3,ibnode,RstNodeList(inode))

                     i = i + 3
                  enddo
                  i = 3*NumRstNodes*ibnode + 1
               enddo
c --- write out stresses
               if (btest(eomask,W_OutStress)) then
                   call WriteToRST(EDENS,lcerstL,edindxL(1),
     &                          3*NumRstNodes*numElemNodes,postdata(1))
               endif
c --- write out strains
               if (btest(eomask,W_OutElaStra)) then
c --- we are storing total strains here i.e., elastic+plastic+hygro
                   call WriteToRST(EDEEL,lcerstL,edindxL(1),
     &                           3*NumRstNodes*numElemNodes,
     &                           postdata(3*numElemNodes*NumRstNodes+1))
               endif

c --- here we are storing the deformed location of the
c     section nodes to be used later for plotting deformed
c     shapes in the postprocessor
c
               if (btest(eomask,W_OutMisc)) then
                   do ii = 1, 3*NumRstNodes*numElemNodes
                      miscdata(nmisc+ii) =
     &                    postdata(numElemNodes*NumRstNodes*3*2+ii)
                   enddo
                   nmisc = nmisc+3*NumRstNodes*numElemNodes
                   call WriteToRST (EDEMS,lcerstL,edindxL(1),nmisc,
     &                              miscdata(1))
               endif
               if (btest(eomask,W_OutPlStra) .and. keyplL.gt.0) then
                   kk = 0
                   ll = numSectPts*3*numDomainIntPts
                   do ii = 1, numDomainIntPts
                      do jj= 1, numSectPts
                         postdata(kk+1)=ple(1,ii,jj)
                         postdata(kk+2)=ple(2,ii,jj)
                         postdata(kk+3)=ple(3,ii,jj)
                         postdata(ll+1)=peeq(1,ii,jj)
                         postdata(ll+2)=peeq(2,ii,jj)
                         postdata(ll+3)=peeq(3,ii,jj)
                         postdata(ll+4)=pwrk(1,ii,jj)
                         postdata(ll+5)=pwrk(2,ii,jj)
                         postdata(ll+6)=pwrk(3,ii,jj)
                         kk = kk+3
                         ll = ll+6
                      enddo
                   enddo
                   call WriteToRST(EDEPL,lcerstL,edindxL(1),
     &                  numSectPts*3*numDomainIntPts,postdata(1))
                   call WriteToRST(EDENL,lcerstL,edindxL(1),
     &                  numSectPts*6*numDomainIntPts,
     &                  postdata(1+numSectPts*3*numDomainIntPts))
               endif

               if (btest(eomask,W_OutCrStra) .and. keycrL.gt.0) then
                   kk = 0
c                   ll = numSectPts*3*numDomainIntPts

                   do ii = 1, numDomainIntPts
                      do jj = 1, numSectPts
                         postdata(kk+1)=cre(1,ii,jj)
                         postdata(kk+2)=cre(2,ii,jj)
                         postdata(kk+3)=cre(3,ii,jj)
c                         postdata(ll+1)=creq(1,ii,jj)
c                         postdata(ll+2)=creq(2,ii,jj)
c                         postdata(ll+3)=creq(3,ii,jj)
c                         postdata(ll+4)=crwk(1,ii,jj)
c                         postdata(ll+5)=crwk(2,ii,jj)
c                         postdata(ll+6)=crwk(3,ii,jj)
                         kk = kk+3
c                         ll = ll+6

                      enddo
                  enddo
                  call WriteToRST(EDECR,lcerstL,edindxL(1),
     &                            numSectPts*3*numDomainIntPts,
     &                            postdata(1))
               endif

               if (btest(eomask,W_OutGyStra)) then
                   i=1
                   do ibnode = 1, numElemNodes
                      do inode = 1, NumRstNodes
                         postdata(i)=
     &                        hyStrSecNd(ibnode,RstNodeList(inode))
                         i = i + 1
                      enddo
                      i = NumRstNodes*ibnode + 1
                   enddo
                   call WriteToRST(EDETH,lcerstL,edindxL(1),
     &                  NumRstNodes*numElemNodes,postdata(1))
                endif
            endif
      endif


 999  continue
      call SetCutBackFact(ishape)
#if defined(DEBUG)
      call SubEnd('BeamElement')
#endif

      return
      end
