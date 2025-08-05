      subroutine LinkElement (elemId,elChar,elData,eomask,eConn,locsvrL,
     &                        kelreq,nUsrdof,xBeg,eleDisp,kelout,
     &                        eStiff,eMass,gStiff,fExt,fInt,elemVolume,
     &                        elemMass,elemCG,eEnergy,edindxL,lcerstL)
c --- primary function: finite strain link element
c     
c  input arguments:
c     elemId  (int,sc,in)             - element number
c     elChar  (int,ar(IELCSZ),in)     - array of element type characteristics
c     elData  (int,ar(EL_DIM),in)     - array of element data
c     eomask  (int,sc,in)             - bit pattern of element output
c     eConn   (int,ar(2),in)          - array of element node numbers
c     locsvrL (int*8,sc,in)           - location of this eles svrs on esav
c     kelreq  (int,ar(10),in)         - matrix and load vector form requests
c     nUsrdof (int,sc,in)             - matrix and LV size
c     xBeg    (dp,ar(6,2),in)         - nodal coords and rot angles
c     eleDisp (dp,ar(nUsrdof,5),in)   - element nodal solution values
c     
c  output arguments:
c     kelout (int,ar(10),out)        - keys indicating created mats and lvs
c     eStiff (dp,ar(nUsrdof,nUsrdof),inout) - k matrix (kelreq(1))
c     eMass  (dp,ar(nUsrdof,nUsrdof),inout) - m matrix (kelreq(2))
c     gStiff (dp,ar(nUsrdof,nUsrdof),inout) - s matrix (kelreq(4))
c     fExt   (dp,ar(nUsrdof),out)    - applied f vector (kelreq(5))
c     fInt   (dp,ar(nUsrdof),out)    - n-r restoring f vector (kelreq(6))
c     elemVolume  (dp,sc,out)        - element volume
c     elemMass (dp,sc,out)           - element mass
c     elemCG (dp,ar(3),out)          - centroid location
c     eEnergy(dp,ar(10),out)         - element energies
c     edindxL(int*8,ar(25),out)      - element result data file indexes
c     lcerstL(int*8,sc,inout)        - position on result file
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
c     
      integer*8  locsvrL, lcerstL, edindxL(25)
      integer
     &     elemId,elChar(*),elData(EL_DIM),eomask,eConn(2),
     &     kelreq(10),nUsrdof,kelout(10),rlgetn,v_last

      double precision
     &     xBeg(6,2),eleDisp(nUsrdof,8),elemCG(5),eEnergy(10),
     &     eStiff(nUsrdof,nUsrdof),eMass(nUsrdof,nUsrdof),
     &     gStiff(nUsrdof,nUsrdof),fExt(nUsrdof),fInt(nUsrdof),
     &     elemVolume,elemMass,aScale,CalDist2P
     
      integer
     &     keyopt(12), matId, ireal, nrvr
c     
      double precision 
     &     mStrain(1), mStress(1)
c     
      double precision ONE, ZERO, HALF, TWO
      parameter (ZERO = 0.0d0, HALF = 0.5d0, ONE = 1.0d0, TWO = 2.0d0)

      integer lwork
      parameter (lwork = 3)
      integer kelin(10),nod,i,ii,j,jj,ishape,lZero
      double precision 
     &     Lengv,work(lwork),postdata(10),avemoist,FailForc,
     &     xCur(3,2),eU(6),incDisp(6),rvr(W_NumContParam),
     &     shpe(2,2),basis(3),tGenStrain,dxdg(3),Leng,
     &     Linv,Linv0,secStiff,Area,Areav,matStiff(1,1),
     &     prop(2),Ex,density,lambda,gradMat(6),Force,
     &     Stiff,Disp,tr(3,3)
      
      equivalence
     &   (Areav, rvr(1)),
     &   (FailForc, rvr(2))
          
      double precision
     &    shpDer(2,2),Length,Lengthv,ltol

      parameter (ltol = 1.0d-9)

c#if defined(MOIST)
      integer khygro,numBasicSaveVars,stateVarIndex(20)
      double precision basicSaveVars(2),refmoist,GetMoistRef,
     &                 moistcur(2),moistbeg(2),moistend(2)
      equivalence
     &   (moistbeg(1), basicSaveVars(1))
c#endif

      integer lp(2)
      data lp /1,13/

      equivalence
     &   (Ex, prop(1)),
     &   (density, prop(2))

#if defined(DEBUG)
      call SubBeg('LinkElement')
#endif
      
c --- basic options   
c     keyopt(2) = 0 section is scaled wrt stretch 
c     keyopt(2) = 1 section is assumed rigid
c
c     rvr(1) = Areav
c
      ishape = 0
      tGenStrain = ZERO
      eEnergy(1) = ZERO
      matId = elData(EL_MAT)
      ireal = elData(EL_REAL)
      nrvr = elChar(NMDRLC)

c --- material properties
c#if defined(MOIST)
      call v_zero(moistcur(1),2)
      refmoist = GetMoistRef(matId)
      call GetMoist (elemId,elChar(1),2,eConn(1),
     &               refmoist,2,moistbeg(1),
     &               moistcur(1),moistend(1),khygro,work(1))
      avemoist = (moistcur(1) + moistcur(2))*HALF
c#else
c      avemoist = ZERO
c#endif      

      call vi_move (elChar(KYOP1),keyopt(1),12)
     
c --- read real constants for the element
      call v_zero (rvr(1),nrvr)
      i = rlgetn(ireal,1,nrvr,rvr(1)) 

c --- read basic, non-linear save variables
c#if defined(MOIST)
      numBasicSaveVars = elChar(NMSSVR)

      call GetSVRIndx(elemId,locsvrL,stateVarIndex(1))

      i = numBasicSaveVars
      call GetSVR(stateVarIndex(1),1,i,basicSaveVars(1))
      
      if (i .ne. numBasicSaveVars) then
          call erhandler ('link',5000,3,
     &    ' SVR record is mismatched in link element.',0.0d0,' ')
          goto 999
      endif
c#endif
      call vi_move (kelreq(1),kelin(1),10)
           
c --- get the element loads
c     
      call v_move (eleDisp(1,A_RowDisp), eU(1), nUsrdof)
      if (A_KeyFirstPass.eq.1) then
          call m_move (xBeg(1,1),xCur(1,1),6,3,3,2)
          call v_move (eleDisp(1,A_RowIncS),incDisp(1),nUsrdof)
      else
          call UpdBeamXYZ (2,eU(1),nUsrdof,xBeg(1,1),6,xCur(1,1))
          call v_move (eleDisp(1,A_RowIncS),incDisp(1),nUsrdof)
      endif

      if (keyopt(1).eq.1) then
c --- spring element (allows for zero length)

c --- compute current length
          Length = CalDist2P (xCur(1,1),xCur(1,2))

c --- compute the original length
          if (A_KeyFirstPass.eq.1) then
              Lengthv = Length
          else
              Lengthv = CalDist2P (xBeg(1,1),xBeg(1,2))
          endif

c --- transformation matrix
          
          if (Length.lt.ltol) then
              if (Lengthv.lt.ltol) then
                  lZero = 2
              else
                  lZero = 1
                  call LineGeomTrans (xBeg(1,1),3,tr(1,1))
              endif
          else
              lZero = 0
              call LineGeomTrans (xCur(1,1),3,tr(1,1))
          endif
          
          Disp = Length - Lengthv

          Stiff = rvr(1)
          Force = Stiff*Disp
          if (keyopt(3).ne.0 .and. v_last(rvr(2),3).ne.0) then
              call UserForce (Disp,rvr(2),Stiff,Force,keyopt(3))
          endif
          
          eEnergy(1) = HALF*Force*Disp
          elemMass = ZERO

c --- stiffness matrix
          if (kelin(1).eq.1) then
              call v_zero (eStiff(1,1),nUsrdof*nUsrdof)
              if (Stiff.ne.0.0d0) then
                  eStiff(1,1) = Stiff
                  if (lZero.lt.2) then
                      call BtxAxB_A (eStiff(1,1),nUsrdof,tr(1,1),3,3,3)
                  else
                      do i = 2,3
                         eStiff(i,i) = eStiff(1,1)
                      enddo
                  endif
                  do i = 1,3
                     ii = i + 3
                     do j = 1,3
                        jj = j + 3                   
                        eStiff(ii,j)  = -eStiff(i,j)
                        eStiff(i,jj)  = -eStiff(i,j)
                        eStiff(ii,jj) =  eStiff(i,j)
                      enddo
                  enddo

              endif
              kelout(1) = 1
          endif

c --- stress stiffness matrix
          if (kelin(4).eq.1) then
              call v_zero (gStiff(1,1),nUsrdof*nUsrdof)
              if (Force.ne.ZERO .and. lZero.eq.0) then
                  gStiff(2,2) = Force/Length
                  gStiff(3,3) = gStiff(2,2)
                  call BtxAxB_A(gStiff(1,1),nUsrdof,tr(1,1),3,3,3)
                  do i = 1,3
                     ii = i + 3
                     do j = 1,3
                        jj = j + 3
                        gStiff(ii,j)  = -gStiff(i,j)
                        gStiff(i,jj)  = -gStiff(i,j)
                        gStiff(ii,jj) =  gStiff(i,j)
                     enddo
                  enddo
              endif
              kelout(4) = 1
          endif

c --- restoring force
          if (kelin(6).eq.1) then
              call v_zero (fInt(1),nUsrdof)
              if (Force.ne.ZERO  .and. lZero.lt.2) then
                  work(1) = -Force
                  work(2) =  0.0d0
                  work(3) =  0.0d0
                  call AtxB_C(tr(1,1),work(1),fInt(1),3,3,nUsrdof,3,1,3)
                  do i = 1,3
                     ii = i + 3
                     fInt(ii) = -fInt(i)
                  enddo
                  kelout(6) = 1
              endif
          endif
          
       else
c --- link element

c --- retrive material data
          call GetMatProps(elemId,matId,lp(1),avemoist,prop(1),2)

          shpe(1,1) = HALF
          shpe(2,1) = HALF
          shpe(1,2) = -HALF
          shpe(2,2) = HALF
c --- current length
          Leng = ZERO    
          do j = 1,3
             dxdg(j) = ZERO
             do nod = 1,2
                dxdg(j) = dxdg(j) + shpe(nod,2)*xCur(j,nod)
             enddo
             Leng = Leng + dxdg(j)*dxdg(j)
          enddo
          Leng = sqrt(Leng)

c --- original length
          if (A_KeyFirstPass.eq.1) then
              Lengv = Leng
          else
              Lengv = ZERO
              do j = 1,3
                 work(1) = ZERO
                 do nod = 1,2
                    work(1) = work(1) + shpe(nod,2)*xBeg(j,nod)
                 enddo
                 Lengv = Lengv + work(1)*work(1)
              enddo
              Lengv = sqrt(Lengv)
          endif

          if (Leng.le.TINY) then
              ishape = 1
              call erhandler('LinkElement',5000,3,
     &                       'Link element %I has zero length.',
     &                        dble(elemId),' ')
              goto 999
          endif

          Linv = ONE/Leng
          Length = TWO*Leng
          basis(1) = dxdg(1)*Linv
          basis(2) = dxdg(2)*Linv
          basis(3) = dxdg(3)*Linv
          Lengthv = TWO*Lengv
          Linv0 = ONE/Lengv
          do nod =1,2
             shpDer(nod,1) = shpe(nod,2)*Linv0
          enddo

c --- define lambda
          lambda = ZERO
          do j = 1,3
             work(j) = ZERO
             do nod = 1,2
                work(j) = work(j)+shpDer(nod,1)*xCur(j,nod)
             enddo
             lambda = lambda+work(j)*work(j)
          enddo
          lambda = ONE/sqrt(lambda)
          elemVolume =  Lengthv*Areav
          aScale = Lengv/Leng
          if (keyopt(2).eq.1) then
              aScale = ONE
          endif
          tGenStrain = log(Leng/Lengv)

          call v_zero (gradMat(1),nUsrdof)

          do nod = 1,2
             j = (nod-1)*3
             do jj = 1,3
                gradMat(j+jj) = lambda*basis(jj)*shpDer(nod,1)
             enddo
          enddo

          Area = Areav*aScale
          mStrain(1) = tGenStrain

c --- do material computations, only elastic
          matStiff(1,1) = Ex
          mStress(1) = matStiff(1,1) * mStrain(1)
          Force = mStress(1)*Area
          
c --- check the failure criterion only if the convergence is achieved
          if (FailForc.lt.ZERO) then
              if (Force.gt.-FailForc) then
                  Force = -FailForc
c --- recalculate actual elastic stress and strain for elastic energy calculations
                  mStress(1) = Force / Area
                  mStrain(1) = mStress(1) / matStiff(1,1)
c --- behave as ideal plastisity, so zero tangent stiffness
                  matStiff(1,1) = TINY
              endif
          elseif (FailForc.gt.ZERO) then
              if (A_KeyConv.eq.1 .and.  A_KeySvrUpd.eq.1) then
                  if (Force.ge.FailForc) then
                      elData(EL_DEAD) = 1
                      call elmput (elemId,elData(1),2,eConn(1))
                  endif
              endif    
          endif
          eEnergy(1) = HALF*mStress(1)*mStrain(1)*Length*Area
          secStiff = matStiff(1,1)*Area      

c --- compute mass matrix
          elemMass = elemVolume*density
          if (kelin(2) .eq. 1) then
              call v_zero (eMass(1,1),nUsrdof*nUsrdof)
              if (elemMass .gt. ZERO) then
                  work(1) = elemMass*HALF
                  do i = 1,6
                     eMass(i,i) = work(1)
                  enddo
              endif
              kelout(2)=1
          endif
          
c --- stiffness matrix
c     
          if (kelin(1).eq.1) then
              call v_zero (eStiff(1,1), nUsrdof*nUsrdof)
              call cAtxA_C (gradMat(1),eStiff(1,1),nUsrdof,nUsrdof,
     &                      Length*secStiff)
              call m_symLtoU (eStiff(1,1),nUsrdof,nUsrdof)
              kelout(1) = 1
          endif
     
c --- stress-stiffness matrix
c     
          call v_zero (gStiff(1,1), nUsrdof*nUsrdof)
          if (Force.ne.ZERO) then
              work(1) = Force*lambda**2*Length/Lengthv**2
              do i = 1,3
                 gStiff(i  ,i  ) =  work(1)
                 gStiff(i+3,i+3) =  work(1)
                 gStiff(i  ,i+3) = -work(1)
                 gStiff(i+3,i  ) = -work(1)
             enddo
c --- newton-raphson restoring force
c
             if (kelin(6).eq.1) then
                 call cv1_v2 (gradMat(1),fInt(1),nUsrdof,Length*Force)
                 kelout(6) = 1
             endif
          elseif (kelin(6).eq.1) then
             call v_zero (fInt(1), nUsrdof)
          endif
      endif

         
c --- external load vector (nothing here yet)
c
      if (kelin(5).eq.1) then
          call v_zero (fExt(1),nUsrdof)
      endif

c --- compute current centroid
       do j=1,3
          elemCG(j) = (xCur(j,1) + xCur(j,2))*HALF
       enddo
      
c --- produce output
      if (A_KeyOut.eq.1) then
          if (btest(eomask,W_OutElIndx)) then
              if (btest(eomask,W_OutMisc)) then 
                  i = elChar(NMSMIS)
                  postdata(1) = Force
                  if (keyopt(1).eq.1) then
                      postdata(2) = Disp
                  else
                      postdata(2) = Area
                  endif
                  postdata(3) = eEnergy(1)
                  call WriteToRST (EDEMS,lcerstL,edindxL(1),i,
     &                             postdata(1))
             endif
          endif
      endif             

c --- finish off matrices
c     
      if (kelout(1).eq.1 .and. kelout(4).eq.1) then
          call v1pv2_v1 (eStiff(1,1),gStiff(1,1),nUsrdof*nUsrdof)
          call v_zero (gStiff(1,1),nUsrdof*nUsrdof)
      endif


c#if defined(MOIST)
c --- handle SVR files
      call StoreSVR(stateVarIndex(1),1,numBasicSaveVars,
     &              basicSaveVars(1))
      call PutSVRIndx(elemId,locsvrL,stateVarIndex(1))
c#endif
      
 999  continue
      call SetCutBackFact(ishape)


#if defined(DEBUG)
      call SubEnd('LinkElement')
#endif

      return
      end


