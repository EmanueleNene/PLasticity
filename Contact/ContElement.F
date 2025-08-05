      subroutine ContElement (elem,elChar,elData,eomask,nodes,
     &                        kelrqq,nr,nvmax,npi,ival,nSeg,nFoun,
     &                        nPerSeg,nPerFounI,nPerFounDp,nPerMast,
     &                        nMast,isegda,ififo,founda,icnda,xyz,u,
     &                        coor0,nodbac,qa,qb,qab,kelout,
     &                        Kmtrx,Fext,Fnr,
     &                        elvol,center,eEnergy,ContStat,
     &                        edindxL,lcerstL)
c --- primary function: compute element matrices, load vectors, and results for beam-beam contact
c
c  input arguments:
c       elem   (int,sc,in)         - element number
c       elChar (int,ar(IELCSZ),in) - array of element type characteristics
c       elData (int,ar(EL_DIM),in) - array of element data
c       eomask (int,sc,in)         - bit pattern of element output
c       nodes  (int,ar(nnod),in)   - array of element node numbers
c       kelrqq (int,ar(10),in)     - matrix and load vector form requests
c       nr     (int,sc,in)         - matrix and lv size
c       nvmax  (int,sc,in)         - number of state variables on each i.p.
c       npi    (int,sc,in)         - number of intergration points
c       ival   (int,ar(10),in)     - element options
c       keyopt (int,ar(12),in)     - element key options
c       nSeg   (int,sc,in)         - total number of segments
c       nFoun  (int,sc,in)         - number of foundations
c       nPerSeg  (int,sc,in)       - number of items per segment database entry
c       nPerFounI (int,sc,in)      - number of integer items per foundation
c       isegda (int,ar(nPerSeg,nSeg),in)
c                                  - segment information
c          isegda(1,iseg)        = code of segment type (icode)
c                                = 1   : straight line
c                                = 2   : parabolic arc
c                                = 7   : quadrilateral
c          isegda(2,iseg)        = fondation number           
c          isegda(j,iseg)        = 4 nodes (j=3-10) on seg. iseg
c          isegda(j,iseg)        = 4 adjacent segments (j=11-14)
c                                = 2 adjacent segments (j=11,12)
c                                     for seg. iseg
c                                = 0  (icode > 10)
c                                  --> none adjacent segment
c       coor0  (dp,ar(*),in)         - initial nodal coordinate coordinate
c       nodbac (int,ar(D_BacLen),out)- internal to exter. list
c
c  output arguments:
c       qb     (dp,ar(nvmax,npi),out)-stresses and state variables at updated
c                                     configuration for each i.p.
c       qb   (1) = presb  - normal contact pressure
c       qb   (2) = taurb  - tangent contact stress 1
c       qb   (3) = tausb  - tangent contact stress 2
c       qb   (5) = ield   -  1: if state is sticking contact (elastic)
c                         -  2: if state is sliding contact (elastoplastic)
c                         -  3: open but near contact
c                         -  4: open and not near contact
c       qb   (6) = ddisp  - energy dissipated per unit area due to friction
c       qb   (8) = iseg   - attached segment number calculated by CalPeni
c       qb   (9) = db     - total penetration distance calculated by CalPeni
c       qb  (10) = nocot  - calculated by CalPeni for initial situation
c       qb  (11) = dajct  - jocoba value used for hygro
c       qb  (12) = vtfrb  - 1st relative tangential velocity
c       qb  (13) = vtfsb  - 2nd relative tangential velocity
c       qb  (14) = gapi   - initial gap
c       qb  (15) = iconv  - convergent key
c       qb  (16) =        - total energy density
c       qb  (17) =        - sliding energy density
c       qb  (18) = rigp   - penetration coefficient
c       qb  (19) = toln   - penetration tolerance
c       qb  (20) = dba    - interpenetration variation
c       qb  (24) = rigt   - tangent penetration coefficient
c       qb  (25) =        -
c       qb  (26) =        - adapted penalty coefficient
c       qb  (27) =        - normal presure at which slip started
c       qb  (28) =        - 0: contact never happen
c                         - 1: contact has happened
c       qb  (29) =        - accumulated sliding in 1st tangent direction
c       qb  (30) =        - accumulated sliding in 2nd tangent direction
c       qb  (31) =        - absolute accumulated sliding in 1st tangent direction
c       qb  (32) =        - absolute accumulated sliding in 2nd tangent direction
c       qb  (33) =        - surface (fluid) pressure
c       qb  (34) =        - moisture
c       qb  (35-43) =     - components of the normal vector
c       qb(nvmax-13)=     - energy released during debonding in mode II
c       qb(nvmax-12)=     - energy released during debonding in mode I
c       qb(nvmax-11)=     - visc. regularized damage parameter for debonding
c       qb(nvmax-10)=     - failure mode for debonding
c                         - 0: no debonding
c                         - 1: mode I (normal mode)
c                         - 2: mode II (tangential mode)
c                         - 3: mixed mode (mode I + mode II)
c       qb(nvmax-9) =     - damage parameter for debonding
c       qb(nvmax-8) =     - debonding fraction, time
c       qb(nvmax-7) =     - initial interference
c       qb(nvmax-6) =     - old status
c       qb(nvmax-5) =     - old gap 
c       qb(nvmax-4) =     - old normal pressure 
c       qb(nvmax-3) =     - old 1st tangent velocity 
c       qb(nvmax-2) =     - old 2nd tangent velocity 
c       qb(nvmax-1) =     - old 1st friction stress 
c       qb(nvmax) =       - old 2nd friction stress 
c       qab     (dp,ar(*),inout) - geometrical variablesc                                
c          qab(1) = rigp     - relative penalty
c          qab(2) = cradi    - contact radius
c          qab(3) = pinb     - pinball
c          qab(4) = length   - contact length
c          qab(5) = ebeam    - Young's modulu
c          qab(7) = elvol    - contact element volume
c          qab(9) =          - attached element shape factor
c          qab(12) =         - igtol
c          qab(13) =         - avegage area for contact stiffness adjustements
c
c       kelout (int,ar(10),out)     - keys indicating created mats and LVs
c       Kmtrx  (dp,ar(nr,nr),inout) - k matrix (kelrqq(1))
c       Fnr    (dp,ar(nr),out)      - n-r restoring force vector
c       elvol  (dp,sc,out)          - element volume
c       center (dp,ar(3),out)       - centroid location
c       eEnergy(dp,ar(5),out)       - element energies
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
#include "solu_info.inc"
#include "conv_stat.inc"
#if defined(MOIST)
#include "in_mem.inc"
#include "beam_hyg.inc"
#endif

c
      integer*8 MakeLongInt, keyL, edindxL(25), lcerstL

      integer elmgct,rlget,ContAdjStat,CheckHyg

      double precision  CalDist2P,GetMoistRef

      integer
     &  elem,elChar(IELCSZ),elData(EL_DIM),nodes(*),nPerSeg,nSeg,
     &  nPerFounI,nFoun,kelrqq(10),nr,kelout(10),i2key,elmgcn,
     &  nvmax,npi,keyopt(12),ival(3),eomask,ContStat,nPerMast,nMast,
     &  isegda(nPerSeg,nSeg),methi,
     &  nodbac(*),ififo(nPerFounI,nFoun),nPerFounDp,kerr,intfer,
     &  icnda(nPerMast,nMast),icon,ind

      double precision
     &  xyz(6,*),u(nr,5),Kmtrx(nr*nr),Fext(nr),Fnr(nr),eEnergy(5),
     &  center(3),qa(nvmax,npi),qb(nvmax,npi),elvol,coor0(*),refmoist,
     &  founda(nPerFounDp,nFoun),qab(40),hygro(16,9),mudat(9),pinb2

      logical lOrthoFric,lTrans,ifDebo,ifFrac,ifConnected

      integer kelin(10),nofo(72),kk(12),nnode,lp(2)
            
      parameter (nnode = 7)
      data lp /16,17/

      integer
     &  mat,ireal,nrvr,nssvr,i,khygro,ifoun,ibeam,mbeam,sbeam,irigf,
     &  nmax,lintv,nnodefu,ipi,ield,lny,lm1,ickcel,isegn,sbeamn,
     &  iseg,nocon,lmx,nno,noeud,k,j,lmy,lnx,imod,
     &  ivolu,kpred,nnodr,isega,nnodmst,cnnode(2),
     &  lmz,lnz,ik,jk,lma,ipass,nrtv,nnodtv,
     &  ksmsld,mramp,ipene,ndin,
     &  kpress,elcnup0,inode,kdsp,kfarcn,iseg0,ickcel1,
     &  keyopt_3,keyopt_9,numpt,option,ipair,krot,
     &  lintvr,lmr,keyopt_1,aupd

      double precision
     &  ZERO,HALF,ONE,avemoist,length,gapi,vect(20),
     &  mu,cradi,rigp,rigt,ebeam,pinb,gapcri,vecnor(2),
     &  da,dajct,x,y,dajc,tsb,tfb,db,presb,
     &  weight,fx,fy,findj,vtfra,vtfrb,vtfsa,vtfsb,
     &  xmax,ymax,xmin,ymin,fpin,zmax,zmin,z,
     &  taurb,tausb,fz,dum1,dum2,dum3,findi,dajcw,elength,
     &  mux,muy,smallmu,TWO,bestU,ta(3),tb(3),ac(3),bc(3)

      parameter (smallmu=0.001d0)

      double precision
     &  rvr(W_NumContParam),moistbeg(nnode),moistend(nnode),
     &  moistcur(nnode),xyt(6,7),xytp(6,7),vxyt(3,7),fifoa(7),
     &  fifob(7),epsil(7),cinem(43),c(4,4),dvet(11),
     &  fi(7),ww(npi),dperr(4),fcell(6),r(3,3),u0(W_NrMax),
     &  prs(4),prsbeg(4),prsend(4),c1(3,3),cc(3,3),
     &  slides,lfact0,coffs,tols,dir(3),tol,
     &  qab1,norm, trgshft,dpres,props(48),MH,mx,my,mz,wid(2),
     &  vax(7),vay(7),vaz(7),xa(7),ya(7),za(7),
     &  rxfb,ryfb,rzfb,rxsb,rysb,rzsb,rigpr,rigpr0,ONED

#if defined(MOIST)      
      double precision
     & bp(11,2),em(H_Row,H_Col+1),cp(2),ang
#endif

      parameter (ZERO=0.0d0,ONE=1.0d0,HALF=0.5d0,TWO=2.0d0,ONED=1.0d0)

#if defined(DEBUG)
      call SubBeg('ContElement')
#endif
c --- define initial data
c
c     nnode =  number of nodes
c     ifoun =  foundation number
c     isold =  attached element number
c     intyp =  type of integration
c     mat   =  material number
      icon = abs(ival(2))
      iseg = 0
      ibeam = icnda(2,icon)
      sbeam = 0
      ifoun = abs(ival(3))
      ield = 4
c      if (ififo(12,ifoun).eq.28177) then
c          write (6,*) '----- > Enter ELEM', elem
c      endif

c --- define initial data
      mat = elData(EL_MAT)
      ireal = elData(EL_REAL)
      nrvr = elChar(NMDRLC)
      nssvr = elChar(NMSSVR)
c --- initialize
      call v_zero(cinem(1),43)
      tsb = ZERO
      tfb = ZERO
c --- get keyopts
      call vi_move (elChar(KYOP1),keyopt(1),12)
c --- get the element real constant data
      ipene = nvmax - 7
      call v_zero (rvr(1),W_NumContParam)
      j = rlget(ireal,rvr(1))
      keyopt_3 = keyopt(3)
      keyopt_9 = keyopt(9)
      if (keyopt_3.gt.2) keyopt_3 = keyopt_3 - 3
c --- adjust contact stiffness due to bisection for the first sub-step
      if (A_SubStepNum.eq.1.and.A_LoadStepNum.eq.1.and.
     &    rvr(3).ge.-TINY) then
          if (keyopt_3.eq.1.or.(keyopt_3.eq.2.and.A_IterNum.eq.1)) then
              if (rvr(3).le.TINY) rvr(3) = ONE
              rvr(3) = rvr(3)*(0.2d0)**A_BisecNum
          endif
      endif

c --- contact pair ID
      ipair = ififo(16,ifoun)
      if (ififo(15,ifoun).lt.1.and.C_ConnUpd.eq.1.and.keyopt_3.lt.2)
     &    goto 999


c --- first pass key
      if (A_KeyMatUpd.eq.1) then
c --- average stress
          qab(16) = founda(12,ifoun)
c --- plastic occurs update contact stiffness
          qab(17) = founda(29,ifoun)
      endif
      if (A_IterNum.eq.1.or.keyopt_3.eq.2) then
          if (A_KeyFirstPass.eq.1) then
              qab(16) = ZERO
              qab(17) = ONE
              qab(18) = ONE
          endif
          founda(29,ifoun) = qab(17)
          if (keyopt_3.eq.1) then
              founda(29,ifoun) = qab(17)
              founda(12,ifoun) = qab(16)
          endif
          if (A_KeyFirstPass.ne.1.and.founda(17,ifoun).lt.TINY)
     &       founda(17,ifoun) = qab(2)
      endif

      call v_zero (u0(1),nr)
      if (A_RowIncS.ne.0) then
          call v_move (u(1,A_RowIncS),u0(1),nr)
      endif

c --- element node number
      nmax = elChar(NMNDST)
      nnodmst = 3
      nnodr = 3
      if (nodes(3).eq.0) nnodr = 2

c --- DOF compoments
      kdsp = 1
      methi = keyopt(8)
c --- integration type
c --- slave type: irigf=1 rigid slave, irigf=2 deformable slave
      irigf = 0

c --- interface option
      intfer = 0
      if ((methi.eq.2.or.methi.eq.4).and.A_LoadStepNum.eq.1) intfer = 1
      mramp = 0
      if (ififo(15,ifoun).ge.1) then
          irigf = 2
          if (intfer.eq.1) mramp = 1
          nmax = elChar(NMNDST)
      endif

c --- error key
      kerr = 2
      kpred = keyopt(7)
      if (A_KeyAutoStep.eq.1.and.A_TimeInc.gt.2.d0*A_MinTimeStep.and.
     &    kpred.ge.1) kerr = 3
      if (A_SubStepNum.eq.1.and.intfer.eq.1) kerr = 2
      if (A_KeyOut.eq.1.or.A_KeyMatUpd.eq.1) kerr = 2

c --- contact configuration update key
      imod = 1
      if (A_LoadStepNum.eq.1.and.mramp.eq.1) imod = 1
      if (A_IterNum.ne.1) imod = 0
      if (A_KeyOut.eq.1.or.A_KeyMatUpd.eq.1) imod = 0
      ksmsld = 0
      if (A_KeyAutoStep.eq.1.and.kpred.ne.0) ksmsld = 1
      if (C_Disp*C_Forc.eq.1) ksmsld = 0
      if (A_IterNum.eq.1) ksmsld = 0
      if (imod.eq.1) ksmsld = 0

c --- slave edge tolerence to prevent slipping off
      if (rvr(25).lt.TINY) then
          tols = ONE  ! in % of the segment length
      else
          tols = rvr(25)
      endif
      
      tol = 0.01d0*tols
      
      call v_zero (hygro(1,1),144)
      if (A_IterNum.eq.1) then
          do ipi = 1,npi
             qb(44,ipi) = ZERO
          enddo
      endif
      call v_zero (mudat(1),npi)

c --- rotational DOFs
      lintv = 3
      if (iand(elChar(KDOFS),ROTZ).le.0) then
          krot = 0
          lintvr = 0
      else
          krot = 1
          lintvr = 3
          lintv = lintv + lintvr
      endif

      nnodtv = lintv*nnode
      nrtv = nr*lintv

c --- update coordinates

      call UpdContXYZ (nnode,u(1,A_RowDisp),u0(1),nnodtv,xyz(1,1),6,
     &                 xyt(1,1),xytp(1,1),krot)

c --- calculate element length and centroid location

c --- element center and length based on corner nodes
      length = CalDist2P(xytp(1,1),xytp(1,2))
      center(1) = (xytp(1,1)+xytp(1,2))*0.5d0
      center(2) = (xytp(2,1)+xytp(2,2))*0.5d0
      center(3) = (xytp(3,1)+xytp(3,2))*0.5d0

      khygro = 0
      refmoist = ZERO
      avemoist = ZERO

c --- reset stiffness formation key if needed
      call vi_move (kelrqq(1),kelin(1),10)

c --- initialize matrices
      if (kelin(1).eq.1) then
          call v_zero (Kmtrx(1),nr*nr)
          kelout(1) = 1
      endif
      if (kelin(6).eq.1) then
          call v_zero (Fnr(1),nr)
          kelout(6) = 1
      endif

#if defined(MOIST)
      if (kelin(5).eq.1) then
          call v_zero (Fext(1),nr)
          kelout(5) = 1
      endif
      call v_zero(wid(1),2)
c --- get moisture
      call v_zero (moistcur(1),nnode)
      refmoist = GetMoistRef(mat)
      
      call GetMoist (elem,elChar(1),nnodmst,nodes(1),refmoist,
     &               nnodmst,moistbeg(1),moistcur(1),moistend(1),khygro,
     &               vect(1))
      
c --- average nodal moisture
      avemoist = ZERO
      ipi = 0
      do 01 i = 1,nnode
         if (nodes(i).eq.0) goto 01
         if (i.le.nnodmst) then
             avemoist = avemoist + moistcur(i)
             ipi = ipi + 1
         endif
  01  continue
      avemoist = avemoist/ipi
      do ipi = 1,npi
         qb(34,ipi) = avemoist
      enddo
#endif


      elvol = ZERO
      aupd = 0
c --- check if it is right after contact update
      if  (C_ConnUpd.eq.0 .and. A_KeyOut.ne.1 .and.
     &     A_KeyFirstIter.eq.0 .and. A_KeySvrUpd.eq.1) aupd = 1
      
      if (icnda(3,icon).eq.0 .and. aupd.eq.1) then
c --- the element status checked above comes from the connectivity update
c          write (8,*) elem
          goto 995
      endif

         
c --- pressures
      kpress = 0

#if defined(PRS)
c --- pressures field
      do ipi = 1,npi
         prsbeg(ipi) = qa(33,ipi)
      enddo
      call prsget (elem,elChar(1),1,4,prsbeg(1),prs(1),prsend(1),
     &             kpress)
#endif

      if (kpress.eq.0) then
          call v_zero (prs(1),npi)
      else
          if (prs(1).lt.TINY) then
              kpress = 0
              call v_zero(prs(1),npi)
          endif
      endif

      kfarcn = kpress

c --- update contact dates for each sub-step
      if (A_IterNum.eq.1.or.keyopt_3.eq.2) then
         qab(6) = ZERO

c --- first pass
         if (A_KeyFirstPass.eq.1) then
            qab(2) = founda(17,ifoun)
            qab(3) = founda(20,ifoun)
            qab(4) = founda(16,ifoun)
            qab(5) = founda(18,ifoun)
            qab(11) = founda(18,ifoun)*1.d-3/founda(17,ifoun)
            qab(12) = founda(22,ifoun)
            qab(13) = founda(21,ifoun)/founda(19,ifoun)
            qab(1) = qab(5)*qab(2)*10.d0

            if (keyopt(3).gt.2) then
c --- using local stiffness
               ivolu = elmgct (elem,dvet(1))
               qab(5) = dvet(5)
               if (ififo(15,ifoun).eq.3) qab(5) = 0.5d0*qab(5)
            endif
         elseif (founda(17,ifoun).gt.TINY.and.A_KeyFirstPass.eq.0.and.
     &           keyopt_3.ge.1) then
            qab1 = 100.d0*founda(12,ifoun)*founda(17,ifoun)
            qab1 = max (qab1,0.1d0*qab(1))
            qab1 = min (qab1,10.d0*qab(1))
            qab(1) = qab1
         endif

         if (keyopt(3).gt.2) then
            qab(4) = length
            cradi = min(abs(rvr(1)),abs(rvr(2)))
            if (cradi.lt.TINY) cradi = max(abs(rvr(1)),abs(rvr(2)))
            qab(2) = cradi

c --- automatically calculate penetration coefficient
            rigp = qab(5)*qab(2)*10.d0
            if (ibeam.ne.elem.and.ibeam.gt.0.and.keyopt(3).ge.4) then
                ivolu = elmgct (ibeam,dvet(1))
c --- using small contact stiffness
                dpres = 1.d-3
                if (keyopt(1).lt.2) then
                    if (keyopt(6).eq.1) dpres = 1.d-6
                    if (keyopt(6).ge.2) dpres = 1.d-9
                endif
                dvet(6) = max(dvet(6),dpres*qab(5))
                qab1 = 100.d0*dvet(6)*qab(2)*qab(18)
                qab1 = max (qab1,0.1d0*rigp)
                qab1 = min (qab1,10.d0*rigp)
                rigp = qab1
            endif
            qab(1) = rigp

c --- calculate contact pinball range inside which contact may occur
            if (rvr(6).gt.TINY) then
                fpin = rvr(6)
            else
                fpin = 3.d0
                if (irigf.eq.2) fpin = 2.0d0
                if (intfer.eq.1) fpin = 4.d0
                if (ififo(15,ifoun).eq.3) fpin = fpin/2.d0
                if (keyopt(1).eq.5.or.keyopt(1).eq.4) fpin = fpin/4.d0
            endif
            pinb = fpin*qab(2)
            pinb = max(pinb,1.1d0*abs(rvr(10)))
            pinb = pinb + abs(rvr(1)+rvr(2))
            if (rvr(6).lt.ZERO) pinb = -rvr(6)
            if (A_KeyFirstPass.ne.1) pinb = max(pinb,qab(3))
            qab(3) = pinb
         endif
      endif

c --- get stress of underlying element
      if (ibeam.ne.elem.and.ibeam.gt.0.and.keyopt_3.ge.1) then
         if (qab(18).ne.0.1d0) then
             qab(18) = ONE
         else
             founda(28,ifoun) = 0.1d0
             qab(17) = 0.1d0
         endif
         ivolu = elmgct (ibeam,dvet(1))

c --- using small contact stiffness
         dpres = 1.d-3
         if (keyopt(1).lt.2) then
            if (keyopt(6).eq.1) dpres = 1.d-6
            if (keyopt(6).ge.2) dpres = 1.d-9
         endif
         dvet(6) = max(dvet(6),dpres*qab(5))
         if (keyopt_3.eq.1)
     &       founda(11,ifoun) = max (founda(11,ifoun),dvet(6))
         if (keyopt(3).eq.4) qab(16) = dvet(6)
      endif
      if (ififo(15,ifoun).lt.1.and.C_ConnUpd.eq.1) goto 999

c --- small sliding case, it may be unstable
      if (ksmsld.eq.1.and.kfarcn.ne.1) then
          do 11 ipi = 1,npi
             if (nint(qa(5,ipi)).ne.4) goto 12
             qb(5,ipi) = 4.d0
             qb(9,ipi) = founda(20,ifoun)
  11      continue
          elvol = qab(7)
          if (irigf.eq.2.and.nodes(nnode+1).ne.0.and.A_KeyOut.ne.1)then
              call vi_zero (nodes(nnode+1),nmax-nnode)
              ContStat = 1
              call lockSMP(W_LockAnsys)
              call elmpnd (elem,nmax,nodes(1))
              call unlockSMP(W_LockAnsys)
          endif
          call v_zero (ww(1),npi)
          goto 995
      endif
  12  continue

c --- re-issue all parameters back

      rigp = qab(1)
      pinb = qab(3)
      cradi = qab(2)
      length = qab(4)
      ebeam = qab(5)

      if (krot.eq.1) then
          if (rvr(3).ge.TINY) then
              rigpr0 = rigp*rvr(3)
          elseif (rvr(3).lt.-TINY) then
              rigpr0 = -rvr(3)
          else
              rigpr0 = rigp
          endif
          rigpr0 = rigpr0*length
      endif

      if (A_KeyMatUpd.eq.1.and.keyopt(3).gt.2) then
          ivolu = elmgct (elem,dvet(1))
          dvet(6) = cradi
          call elmpct (elem,dvet(1))
      endif

c --- igtol
      gapcri = rvr(5)
      if (abs(gapcri).le.TINY.and.keyopt(5).eq.0) gapcri = qab(12)
      if (abs(gapcri).le.TINY.and.keyopt(5).eq.4.and.keyopt(1).eq.6)
     &    gapcri = qab(12)
      coffs = rvr(10)
      if (keyopt(5).gt.0.and.keyopt(5).lt.4) coffs = rvr(10)+qab(12)

      call vi_zero (kk(1),12)
      call vi_zero (nofo(1),npi*12)
      dvet(1) = ZERO

c --- check for debonding in case of bonded contact;
c     store cohesive zone material data for debonding in
c     cinem(36-43)
c
      ifConnected = .true.
c --- check if they are connected
      do ipi = 1,npi
         ifConnected = (ifConnected.or.abs(qa(nvmax-9,ipi)).lt.ONED)
      enddo

      if (keyopt(1).ge.2 .and. keyopt(1).le.6 .and. ifConnected) then
          call v_zero(props(1),6)
          keyL = MakeLongInt(mat,NLMP_USERMAT)
          call GetTabParam(keyL,0.0d0,numpt,props(1),option)
          if (rvr(27).ne.ZERO) then
              call cv1_v1 (props(1),min(numpt,4),abs(rvr(27))) 
          endif
          if (option.eq.1) then
              call v_move(props(1),cinem(36),min(numpt,6))
          elseif (option.eq.2) then
              call v_move(props(1),cinem(36),min(numpt,6))
              if (abs(cinem(36)).gt.TINY) then
                  cinem(37) = cinem(37)*TWO/cinem(36)
              else
                  cinem(37) = ZERO
              endif
              if (abs(cinem(38)).gt.TINY) then
                  cinem(39) = cinem(39)*TWO/cinem(38)
              else
                  cinem(39) = ZERO
              endif
          elseif (option.eq.3) then
              call v_move(props(1),cinem(36),min(numpt,6))              
              if (rvr(3).le.-TINY) then
                  rigp = -rvr(3)
                  if (abs(cinem(36)).gt.TINY) then
                      dum1 = HALF*rigp/cinem(36)
                      cinem(42) = cinem(36)
                      cinem(36) = rigp*dum1-cinem(36)*dum1**2
                      cinem(37) = max(cinem(37),dum1)
                  else
                      cinem(37) = ZERO
                  endif
              endif
              if (rvr(12).le.-TINY) then
                  rigt = -rvr(12)
                  if (abs(cinem(38)).gt.TINY) then
                      dum1 = HALF*rigt/cinem(38)
                      cinem(43) = cinem(38)
                      cinem(38) = rigt*dum1-cinem(38)*dum1**2
                      cinem(39) = max(cinem(39),dum1)
                  else
                      cinem(39) = ZERO
                  endif
              endif
          elseif (option.eq.4) then
              call v_move(props(1),cinem(36),min(numpt,6))              
              if (rvr(3).le.-TINY) then
                  rigp = -rvr(3)
                  if (cinem(36).gt.TINY .and. cinem(37).gt.TINY) then
                      cinem(42)=(cinem(36)-rigp*cinem(37))/cinem(37)**2
                      cinem(42) = min(ZERO,cinem(42))
                  else
                      cinem(37) = ZERO
                  endif
              endif
              if (rvr(12).le.-TINY) then
                  rigt = -rvr(12)
                  if (cinem(38).gt.TINY .and. cinem(39).gt.TINY) then
                      cinem(43)=(cinem(38)-rigt*cinem(39))/cinem(39)**2
                      cinem(43) = min(ZERO,cinem(43))
                  else
                      cinem(39) = ZERO
                  endif
              endif
          elseif (option.gt.0) then
              dperr(1) = elem
              dperr(2) = ireal
              call erhandler ('ContElement',5000,3,'Invalid
     &                         model option for debonding
     &                         %I (real ID %I)',dperr(1),' ')
          endif
      endif
      
c --- no friction for bonded contact      
      if ((keyopt(1).eq.3 .or. keyopt(1).ge.5) .and. ifConnected) then
          lOrthoFric = .false.
          lTrans = .false.
          mu = 0.0d0
          goto 100
      endif
      
c --- get friction coefficient
      call GetMatProps (elem,mat,14,avemoist,mu,1)
      if (mu.lt.0.0d0) mu = 0.0d0
c --- get friction coefficient
      lOrthoFric = .false.
      lTrans = .false.
      if (mu.eq.0.0d0) then
         option = -1
         call GetMatProps (elem,mat,lp(1),avemoist,props(1),2)
         if (props(1).ge.0.0d0 .and. props(2).gt.0.0d0) then
             option = 2
         elseif (props(1) .gt. 0.0d0) then
             option = 1
         endif
         if (option.eq.1) then
c --- isotropic friction
            mu = props(1)
         elseif (option.eq.2) then
c --- orthotropic friction
            mux = props(1)
            muy = props(2)
            if (mux.lt.ZERO) mux = ZERO
            if (muy.lt.ZERO) muy = ZERO
            if (mux.eq.ZERO .and. muy.eq.ZERO) then
                mu = ZERO
            else
               if (mux.eq.ZERO) then
                   mux = smallmu*muy
               elseif (muy.eq.ZERO) then
                   muy = smallmu*mux
               endif
               mu = sqrt(HALF*(mux*mux+muy*muy))
               lOrthoFric = .true.
            endif
         endif
      endif
100   continue
      cinem(19) = ZERO
      cinem(20) = ZERO
      cinem(21) = ZERO
      cinem(22) = ZERO
      cinem(23) = ZERO
      cinem(25) = ZERO
      cinem(32) = ONE
      ww(1) = ONE
      fi(1) = ONE

c --- set element depth options
      if (ififo(16,ifoun).eq.0.and.ififo(17,ifoun).eq.0.and.
     &    ((C_ConnUpd.eq.0.and.
     &     ififo(15,ifoun).ge.1.and.A_KeyMatUpd.ne.1).or.
     &     ksmsld.eq.1)) then
           ickcel = 1
           goto 03
      endif

c --- check whether the current element is in the foundation cell
      ickcel = 0
      if (ifoun .le. 0) then
          ickcel = 1
      else
          if (A_KeyFirstPass.ne.1) then
              if (keyopt(6).eq.3) then
                  ickcel = 3
                  goto 03
              endif
              if (A_KeyMatUpd.ne.1) then
                  do ipi=1,npi
                     if (qa(5,ipi).lt.4.0d0) then
                         goto 13
                     endif
                  enddo
                  ickcel = 1
                  goto 03
              endif
 13           continue
              do ipi=1,npi
                 if (qb(5,ipi).lt.4.0d0) then
                     ickcel = 2
                     goto 03 
                 endif
              enddo
          endif
          xmax = xytp(1,1)
          ymax = xytp(2,1)
          zmax = xytp(3,1)
          xmin = xytp(1,1)
          ymin = xytp(2,1)
          zmin = xytp(3,1)
          do 02 i = 2,nnode
             if (nodes(i).eq.0) goto 02
             xmax = max (xytp(1,i),xmax)
             xmin = min (xytp(1,i),xmin)
             ymax = max (xytp(2,i),ymax)
             ymin = min (xytp(2,i),ymin)
             zmax = max (xytp(3,i),zmax)
             zmin = min (xytp(3,i),zmin)
  02      continue
          pinb2 = pinb
          if (pinb2.le.length) pinb2 = length
          fcell(1) = founda(1,ifoun) - pinb2
          fcell(2) = founda(2,ifoun) + pinb2
          fcell(3) = founda(3,ifoun) - pinb2
          fcell(4) = founda(4,ifoun) + pinb2
          fcell(5) = founda(5,ifoun) - pinb2
          fcell(6) = founda(6,ifoun) + pinb2
          if (xmax.lt.fcell(1)) ickcel = 1
          if (xmin.gt.fcell(2)) ickcel = 1
          if (ymax.lt.fcell(3)) ickcel = 1
          if (ymin.gt.fcell(4)) ickcel = 1
          if (zmax.lt.fcell(5)) ickcel = 1
          if (zmin.gt.fcell(6)) ickcel = 1
      endif

  03  continue

c --- out of the cell region, calculate element volume and return
      if (ickcel.eq.1 .and. kfarcn.ne.1) then
          do 05 ipi = 1,npi
             if (A_KeyFirstPass.ne.1) then
                 elvol = qab(7)
                 qb(11,ipi) = qa(11,ipi)
                 goto 42
             endif
c --- coordinates and velocities in initial configuration 
             dajc = ONE
             elvol = elvol + dajc
             qb(11,ipi) = dajc
  42         continue
             qb(1,ipi) = ZERO
             qb(2,ipi) = ZERO
             qb(3,ipi) = ZERO
             qb(21,ipi) = ZERO
             if (C_ConnUpd.eq.0) then
                 ield = 4
                 if (A_IterNum.eq.1) then
                     if (qa(5,ipi).eq.ZERO) qa(5,ipi) = 4.0d0
                     qb(44,ipi) = abs(abs(ield)-abs(qa(5,ipi)))
                 elseif (A_IterNum.gt.1) then
                     if (qb(5,ipi).eq.ZERO) qb(5,ipi) = 4.0d0
                     qb(44,ipi) = qb(44,ipi) +
     &                            abs(abs(ield)-abs(qb(5,ipi)))
                 endif
             endif
             qb(5,ipi) = 4.0d0
             qb(6,ipi) = ZERO
             qb(8,ipi) = ZERO
             qb(9,ipi) = pinb
             qa(18,ipi) = -ONE
             qb(29,ipi) = qa(29,ipi)
             qb(30,ipi) = qa(30,ipi)
             qb(31,ipi) = qa(31,ipi)
             qb(32,ipi) = qa(32,ipi)
             if (abs(qa(5,ipi)).lt.HALF) qa(5,ipi) = 4.0d0
             if (abs(nint(qa(5,ipi))).le.2 .and.
     &          (C_ConnUpd.eq.0.or.kerr.eq.3)) then
                 if (C_nConF.eq.0) then
                     call lockSMP (W_LockConvVar)
                     C_nConF = C_nConF + 1
                     call unlockSMP (W_LockConvVar)
                     if (A_Debug.eq.1) then
                         dperr(1) = elem
                         dperr(2) = ireal
                         dperr(3) = isegda(15,nint(qa(8,ipi)))
                         call erhandler ('ContElement',5000,kerr,
     &                   'Contact element %I (ID %I) status changed
     &                    suddenly from contact (with slave element %I)
     &                    to no-contact.',dperr(1),' ')
                     endif
                 endif
             endif
 05       continue
          if (irigf.eq.2.and.nodes(nnode+1).ne.0.and.A_KeyOut.ne.1) then
              call vi_zero (nodes(nnode+1),nmax-nnode)
              ContStat = 1
              call lockSMP(W_LockAnsys)              
              call elmpnd (elem,nmax,nodes(1))
              call unlockSMP(W_LockAnsys)
          endif 
          goto 995
      endif
c
      do 06 i = 1,nnode
         if (nodes(i).eq.0) goto 06
         vxyt(1,i) = xytp(1,i) - xyt(1,i)
         vxyt(2,i) = xytp(2,i) - xyt(2,i)
         vxyt(3,i) = xytp(3,i) - xyt(3,i)
 06   continue

      eEnergy(1) = ZERO
      eEnergy(4) = ZERO
      elvol     = ZERO
      cinem(8)  = ZERO
      cinem(9)  = ONE
      cinem(10) = mu
      cinem(11) = qab(1)
      cinem(12) = qab(2)
      cinem(13) = qab(4)
      cinem(14) = elem
      cinem(15) = pinb
      if (lOrthoFric) then
          cinem(34) = mux/mu
          cinem(35) = muy/mu
      endif

c --- initial contact detection
c
      if (A_KeyFirstPass.eq.1) then
         if (gapcri.gt.ZERO) then
            gapcri = min(gapcri*qab(4),gapcri*qab(2))
            if (gapcri.ge.0.9d0*qab(2)) then
               call erhandler ('ContElement',5000,
     &          1,'Specified initial contact closure is too large for
     &          element %i. Default value is used.',dble(elem),' ')
               gapcri = min(0.05d0*qab(4),0.05d0*qab(2))
            endif
         else
            gapcri = - gapcri
         endif
      endif

c --- loop on integration points

      nnodefu = 0
      if (ickcel.eq.1) imod = 0

      dajc = ONE
      dajct = dajc
      
      do 900 ipi=1,npi

         ickcel1 = ickcel

         isega  = nint(qa(8,ipi))
         iseg  = nint(qb(8,ipi))

         if (ipi.eq.1) then
             iseg0 = iseg
         else
             if (isega.eq.0) isega = iseg0
         endif
         dajct = qa(11,ipi)
         cinem(6) = ZERO
         call vi_move (nodes(1),cnnode(1),2)

c --- small sliding case
         if (ksmsld.eq.1.and.kfarcn.ne.1) then
            if (nint(qa(5,ipi)).eq.4) then
               elvol = elvol + ww(ipi)*dajct
               qb(21,ipi) = ZERO
               qb(5,ipi)  = 4.d0
               qb(11,ipi) = dajct
               goto 900
            endif
         endif

         if (ickcel1.eq.3) then
             if (qb(10,ipi).ne.2.or.iseg.eq.0) ickcel1 = 1
         endif

c --- initial penetration

c --- set initial contact closure
         if (A_KeyFirstPass.eq.1) then
             qb(14,ipi) = ZERO
             qb(ipene,ipi) = ZERO
         endif

         if (dajct.eq.ZERO) imod = 1
         if (imod.eq.1.and.ickcel1.ne.1) then
c
            do inode = 1,nnode
               if (nodes(inode).gt.0) then
                   xa(inode) = xyt(1,inode)
                   ya(inode) = xyt(2,inode)
                   za(inode) = xyt(3,inode)
                   vax(inode) = vxyt(1,inode)
                   vay(inode) = vxyt(2,inode)
                   vaz(inode) = vxyt(3,inode)
               endif
            enddo
            
            if (A_KeyFirstPass.eq.0) then
                trgshft = qb(ipene,ipi) + coffs
            else
                trgshft = ZERO
            endif

            call v_move (qa(35,ipi),r(1,1),9)
            call CalPeni(nodes(1),ickcel1,x,y,z,xa(1),ya(1),za(1),
     &                   vax(1),vay(1),vaz(1),r(1,1),
     &                   pinb,nSeg,ifoun,nPerSeg,nPerFounI,nPerFounDp,0,
     &                   ififo(1,ifoun),coor0(1),isegda(1,1),
     &                   founda(1,ifoun),nocon,isega,iseg,epsil(1),
     &                   fi(1),fifoa(1),tols,trgshft,cnnode(1),krot,
     &                   qab(1),ta(1),tb(1),ac(1),bc(1),elem)
            call v_move (r(1,1),qa(35,ipi),9)
c
            do inode = 1,2
               if (cnnode(inode).lt.0) then
                   cnnode(inode) = abs(cnnode(inode))
                   epsil(3) = -TINY
               endif
            enddo
            call v_move (r(1,1),qb(35,ipi),9)
            da = epsil(3)

c --- initial contact closure
            if (A_KeyFirstPass.eq.1.and.
     &          nocon.eq.0.and.abs(da).lt.pinb) then
              
                if ((da.gt.ZERO.and.da.le.gapcri).or.
     &             ((da.lt.ZERO.or.keyopt(1).eq.4.or.keyopt(1).eq.5)
     &             .and.(methi.eq.3.or.methi.eq.4.or.
     &              methi.eq.0))) then

                    qb(ipene,ipi) = da
                    da = ZERO

c         if (ififo(12,ifoun).eq.28177) then
c             write (6,*) '-- > qb(ipene,ipi) ini', elem, qb(ipene,ipi)
c         endif

                endif
            endif

            if ((methi.eq.3.or.methi.eq.4).and.A_KeyFirstPass.eq.1) then
               trgshft = qb(ipene,ipi) + coffs

               call v_move (qa(35,ipi),r(1,1),9)
               call CalPeni (nodes(1),ickcel1,x,y,z,xa(1),ya(1),za(1),
     &                       vax(1),vay(1),vaz(1),r(1,1),
     &                       pinb,nSeg,ifoun,nPerSeg,nPerFounI,
     &                       nPerFounDp,0,ififo(1,ifoun),coor0(1),
     &                       isegda(1,1),founda(1,ifoun),nocon,isega,
     &                       iseg,epsil(1),fi(1),fifoa(1),tols,trgshft,
     &                       cnnode(1),krot,qab(1),ta(1),tb(1),
     &                       ac(1),bc(1),elem)
               do inode=1,2
                  if (cnnode(inode).lt.0) then
                      cnnode(inode) = abs(cnnode(inode))
                      epsil(3) = -TINY
                  endif
                  da    = epsil(3)
               enddo
            endif

            do inode=1,2
               cnnode(inode) = abs(cnnode(inode))
            enddo
            vtfra = epsil(1)
            vtfsa = epsil(2)

            if (nocon.eq.2) then
                iseg = 0
            elseif (nocon.eq.0) then
                iseg0 = iseg
            endif

            qa(8,ipi) = iseg
            qa(11,ipi) = dajc
            qa(10,ipi)= nocon + 2
            cinem(3) = 0

            if (nint(qa(5,ipi)).eq.0) qa(5,ipi) = 4.0d0
            if (nint(qa(5,ipi)).eq.4.and.abs(da).lt.pinb) then
                if (A_KeyFirstPass.eq.1.or.da.ge.ZERO) then
                    cinem(3) = 0
                    if (nocon.eq.0.and.ififo(15,ifoun).eq.3.
     &                  and.iseg.gt.0) then
                        ndin = 2
                        if (ibeam.eq.abs(isegda(16,iseg))) cinem(3) = 2
                        if (cinem(3).lt.0.5d0)
     &                      cinem(3) = ContAdjStat(nodes(1),ndin,
     &                                             isegda(1,1),iseg,
     &                                             nodbac(1),nPerSeg,
     &                                             nSeg)
                    endif
                    if (cinem(3).eq.0) then
                        qa(5,ipi) = 3.0d0
                        if (nocon.eq.0.and.da.le.ZERO) qa(5,ipi) = 2.d0
                    endif
                endif
            endif

            qa(9,ipi) = da

c --- ramp contact option only for first load step
            if (nocon.eq.0.and.da.lt.ZERO.and.intfer.eq.1) then
                if (A_KeyFirstPass.eq.1.or.mramp.eq.1) qb(14,ipi) = da
            endif

            if ((A_KeyFirstPass.eq.1).and.nocon.eq.0.and.
     &         -da.ge.0.8d0*qab(2).and.methi.le.1.and.cinem(3).eq.0.and.
     &         (C_nPenIni.eq.0.or.C_nPenIni.eq.2).and.
     &         (C_ConnUpd.eq.0.or.kerr.eq.3)) then
               dperr(1) = -da
               dperr(2) = elem
               dperr(3) = isegda(15,iseg)
               call lockSMP (W_LockConvVar)
               if (C_nPenIni.eq.0) C_nPenIni = 1
               if (C_nPenIni.eq.2) C_nPenIni = 3
               call unlockSMP (W_LockConvVar)
               call erhandler ('ContElement',5000,
     &                2,'Initial penetration %g may be too
     &                large for contact element %i (with slave
     &                element %i).',dperr(1),' ')
            endif

            if ((A_KeyFirstPass.eq.1).and.nocon.eq.0.and.
     &         -da.ge.0.6d0*qab(3).and.intfer.eq.1.and.
     &         (C_nPenIni.eq.0.or.C_nPenIni.eq.1).and.
     &         (C_ConnUpd.eq.0.or.kerr.eq.3)) then
               dperr(1) = -da
               dperr(2) = elem
               dperr(3) = isegda(15,iseg)
               call lockSMP (W_LockConvVar)
               if (C_nPenIni.eq.0) C_nPenIni = 2
               if (C_nPenIni.eq.1) C_nPenIni = 3
               call unlockSMP (W_LockConvVar)
               call erhandler ('ContElement',5000,
     &                2,'Initial penetration %g may be too
     &                large for contact element %i (with slave
     &                element %I).',dperr(1),' ')
            endif

            qa(12,ipi) = vtfra
            qa(13,ipi) = vtfsa
            cinem(4) = ZERO
            cinem(5) = nocon
            cinem(7) = 1
            cinem(10) = mu
            if (lOrthoFric) then
                cinem(34) = mux/mu
                cinem(35) = muy/mu
            endif
         endif

         do inode = 1,nnode
            if (nodes(inode).gt.0) then
                xa(inode) = xytp(1,inode)
                ya(inode) = xytp(2,inode)
                za(inode) = xytp(3,inode)
                vax(inode) = vxyt(1,inode)
                vay(inode) = vxyt(2,inode)
                vaz(inode) = vxyt(3,inode)
            endif
         enddo

c --- set ramping contact interfernce
         gapi = ZERO
         if (intfer.eq.1.and.qb(14,ipi).lt.ZERO) then
             gapi = (ONE-A_CurLoadFact)*qb(14,ipi)
             if (mramp.eq.1.and.A_TimeEnd.gt.A_TimeBeg) then
                 lfact0 = (A_TimeCur - A_TimeInc - A_TimeBeg)/
     &                    (A_TimeEnd - A_TimeBeg)
                 gapi = (ONE-A_CurLoadFact)/(ONE-lfact0)*qb(14,ipi)
             endif
         endif

         trgshft = gapi + qb(ipene,ipi) + coffs
         cinem(17) = 0
         if (isega.gt.0) cinem(17) = isegda(15,isega)
         call v_move (qa(35,ipi),r(1,1),9)

c         if (ififo(12,ifoun).eq.28177) then
c             write (6,*) '----- > Enter CalPen', elem, qb(ipene,ipi)
c         endif
      
         call CalPeni(nodes(1),ickcel1,x,y,z,xa(1),ya(1),za(1),
     &                vax(1),vay(1),vaz(1),r(1,1),
     &                pinb,nSeg,ifoun,nPerSeg,nPerFounI,nPerFounDp,
     &                1,ififo(1,ifoun),coor0(1),isegda(1,1),
     &                founda(1,ifoun),nocon,isega,iseg,epsil(1),
     &                fi(1),fifob(1),tols,trgshft,cnnode(1),krot,
     &                qab(1),ta(1),tb(1),ac(1),bc(1),elem)


         qb(ipene,ipi) = trgshft - gapi - coffs

c         if (ififo(12,ifoun).eq.28177) then
c             write (6,*) 'trgshft, gapi, coffs', trgshft, gapi, coffs
c             write (6,*) 'qb(ipene,ipi)',qb(ipene,ipi)
c             write (6,*) 'qab(19:21)', qab(19:21)
c             write (6,*) '<-----  Exit ELEM', elem, nocon, epsil(3)
c         endif

         
         i2key = 0
         if (nocon.eq.0) then
             i2key = elmgcn(elem,vecnor(1))
         else
             call v_zero (fi(1),nnode)
             if (nodes(3).eq.0) then
                 fi(1) = 0.5d0
                 fi(2) = 0.5d0
             else
                 fi(3) = 1.0d0
             endif
         endif
         
         bestU = qab(19)
         cinem(6) = ZERO
                  
         if (i2key.gt.0) then
             if (abs(bestU+ONE).le.tol) then
                 if (nint(vecnor(1)).ne.elem) then
                     cinem(6) = ONE 
                 endif
             endif
             if (abs(bestU-ONE).le.tol) then
                 if (nint(vecnor(2)).ne.elem) then
                     cinem(6) = ONE
                 endif
             endif
         endif
         if (A_KeyFirstPass.eq.1) then
             call v_move (r(1,1),qa(35,ipi),9)
             goto 105
         endif

     
         if (nocon.ne.0) goto 105

         if (qa(10,ipi).eq.2) then
             dir(1) = qa(37,ipi)
             dir(2) = qa(40,ipi)
             dir(3) = qa(43,ipi)
         elseif (qb(10,ipi).eq.2) then
             dir(1) = qb(37,ipi)
             dir(2) = qb(40,ipi)
             dir(3) = qb(43,ipi)
         else
             goto 105
         endif

         norm = dir(1)*r(3,1)+dir(2)*r(3,2)+dir(3)*r(3,3)

c         if (norm.le.TINY .and. A_Debug.eq.1) then
         if (norm.le.TINY. and. A_KeySvrUpd.eq.1 
     &       .and. qa(nvmax-9,ipi).gt.TINY) then

             dperr(1) = elem
             dperr(3) = ireal
             dperr(2) = isegda(15,iseg)
             ififo(22,ifoun) = ififo(22,ifoun) + 4
             call PutErr(1)
             call erhandler ('ContElement',5000,kerr,'Contact element %I
     &                        with slave element %I (real %I) changes 
     &                        normal sign. The underlying beams cross.',
     &                        dperr(1),' ')

             if (A_IterNum.le.2 .or. A_BisecNum.gt.1) then
c --- mark as killed
                 elData(EL_DEAD) = 1
                 call vi_zero (nodes(nnode+1),nmax-nnode)
                 call lockSMP(W_LockAnsys)
                 call elmput (elem,elData(1),nnode,nodes(1))
                 call unlockSMP(W_LockAnsys)
c --- in this case we simply kill the element, so the update flag can stay intack
c                 ContStat = 1
                 goto 999
             endif

         endif

 105     call v_move (r(1,1),qb(35,ipi),9)

         do inode = 1,2
            if (cnnode(inode).lt.0) then
                cnnode(inode) = abs(cnnode(inode))
                epsil(3) = -TINY
            endif
         enddo
         vtfrb = epsil(1)
         vtfsb = epsil(2)
         db    = epsil(3)

         if (abs(nint(qa(5,ipi))).eq.1 .or. abs(nint(qa(5,ipi))).eq.2)
     &       ififo(25,ifoun) = 0
c --- bisection for first time near field to contact
         if (kerr.eq.3.and.db.le.ZERO.and.kpred.ge.2.and.
     &       ififo(25,ifoun).lt.0) then
             if (C_nConF.eq.0) then
                 call lockSMP (W_LockConvVar)
                 C_nConF = C_nConF + 1
                 call unlockSMP (W_LockConvVar)
                 if (A_Debug.eq.1) then
                     dperr(1) = elem
                     dperr(2) = ireal
                     call erhandler ('ContElement',5000,kerr,'Contact
     &                element %I (ID %I) status changes suddenly.',
     &                dperr(1),' ')
                 endif
             endif
         endif

         if (A_KeyFirstPass.eq.1) qa(ipene,ipi) = db

         tsb = ZERO
         rxsb = ZERO
         rysb = ZERO
         rzsb = ZERO
         if (krot.eq.1) then
             do inode = 1,nnode
                if (nodes(inode).gt.0) then
                    rxsb = rxsb + fi(inode)*xytp(4,inode)
                    rysb = rysb + fi(inode)*xytp(5,inode)
                    rzsb = rzsb + fi(inode)*xytp(6,inode)
                endif
             enddo
         endif

         tfb   = epsil(7)

         if (krot.eq.1) then
             rxfb = epsil(4)
             ryfb = epsil(5)
             rzfb = epsil(6)
         endif

         cinem(1) = tsb
         cinem(2) = tfb
         cinem(5) = nocon
         cinem(7) = 0
         cinem(11) = qab(1)
         if (nocon.eq.0.and.iseg.gt.0) then
             cinem(16) = isegda(15,iseg)
         else
            if (nint(qa(8,ipi)).eq.0) then
                cinem(16) = ZERO
            else
                cinem(16) = isegda(15,nint(qa(8,ipi)))
            endif
         endif

         if (nocon.eq.2) then
             iseg = 0
         elseif (nocon.eq.0) then
             iseg0 = iseg
         endif


         qb(8,ipi) = iseg
         qb(9,ipi) = db
         qb(11,ipi) = dajc
         qb(10,ipi) = nocon + 2
         
c         if (qb(10,ipi).ne.qa(10,ipi)) then
c             write (6,*) 'Strange 1', elem, qa(10,ipi), qb(10,ipi)
c         endif
         
c         if (qb(5,ipi).ne.qa(5,ipi)) then
c             write (6,*) 'Strange 2', elem, qa(5,ipi), qb(5,ipi)
c         endif
         
         
         qb(12,ipi) = vtfrb
         qb(13,ipi) = vtfsb

         ipass = 0
         cinem(25) = ZERO
         if (ififo(1,ifoun).eq.0) cinem(25) = ONE
 1569    continue
         if (keyopt(3).lt.3) goto 65
         if (A_IterNum.eq.1.and.nocon.eq.0.and.iseg.ne.0.and.
     &       rvr(3).ge.ZERO.and.irigf.eq.2) then
             ebeam = ZERO
             elength = ZERO
             if (isegda(16,iseg).ne.0) then
c --- find slave stiffness if deformable-deformable contact
                 ivolu = elmgct (isegda(15,iseg),dvet(1))
                 if (ivolu.ne.0) ebeam = dvet(4)
                 if (ififo(15,ifoun).eq.3) ebeam = 0.5d0*ebeam
             endif
c --- select the softer one for penalty coefficient
             if (ebeam.gt.0) then
                 if (elength.le.ZERO) elength = qab(2)
                 if (qab(1)*elength.gt.ebeam) cinem(11) = ebeam*elength
                 qab(1) = cinem(11)
             endif
         endif
  65     continue

         cinem(3) = 0
         if (nocon.eq.0.and.ififo(15,ifoun).eq.3.and.
     &       qa(5,ipi).eq.3.0d0 .and. -qb(9,ipi) .gt. 0.1d0*pinb) then
             ndin = 2
             if (ibeam.eq.abs(isegda(16,iseg))) cinem(3) = 2
             if (cinem(3).lt.0.5d0)
     &           cinem(3) = ContAdjStat(nodes(1),ndin,isegda(1,1),
     &                                  iseg,nodbac(1),nPerSeg,nSeg)
         endif

         weight = ww(ipi)
         if (ipass.eq.1) weight = HALF*weight

c --- compute element volume and energy

         dajcw = weight*dajc

         cinem(10) = mu
         cinem(24) = ONE
         if (lOrthoFric) then
             cinem(34) = mux/mu
             cinem(35) = muy/mu
         endif

         qb(33,ipi) = prs(ipi)

         if (lOrthoFric) then
c --- transform tangential stress and sliding for
c --- orthotropic friction
             lTrans = .true.
             call OrthoToIsoFric(nvmax,cinem(1),qa(1,ipi),qb(1,ipi),
     &                           c(1,1),lTrans,keyopt(1))
         endif

c         if (ififo(12,ifoun).eq.28177) then
c             write (6,*) 'FrictionLaw', db, qa(9,ipi), qb(9,ipi)
c         endif

         call FrictionLaw (nvmax,keyopt(1),cinem(1),rvr(3),qa(1,ipi),
     &                     A_TimeInc,qb(1,ipi),c(1,1),intfer,kerr,
     &                     hygro(1,ipi),ififo(1,ifoun),
     &                     founda(1,ifoun),qab(1),ONED,keyopt_1)
         
c         if (kelin(1).eq.1)  write (6,*) 'Matrix c',c(1,1),c(2,2),c(3,3)
         
         if (qb(5,ipi).eq.5.0d0) then
c --- mark as killed
             elData(EL_DEAD) = 1
             call vi_zero (nodes(nnode+1),nmax-nnode)
             call lockSMP(W_LockAnsys)
             call elmput(elem,elData(1),nnode,nodes(1))
             call unlockSMP(W_LockAnsys)
c --- in this case we simply kill the element, so the update flag can stay intack
c             ContStat = 1
             goto 999
         endif

         if (lOrthoFric) then
c --- transform tangential stress and sliding back
c
             lTrans = .false.
             call OrthoToIsoFric(nvmax,cinem(1),qa(1,ipi),qb(1,ipi),
     &                           c(1,1),lTrans,keyopt(1))
         endif
         ield  = abs(nint(qb(5,ipi)))
         mudat(ipi) = cinem(10)

#if defined(MOIST)
         if (ield.le.2) then
            if (A_KeyFirstPass.eq.1) then

                dir(1) = r(3,1)
                dir(2) = r(3,2)
                dir(3) = r(3,3)
                
                call v_move(dir(1),qab(23),3)
                call v_move(ta(1),qab(26),3)
                call v_move(tb(1),qab(29),3)
                call v_move(ac(1),qab(32),3)
                call v_move(bc(1),qab(35),3)
                
                cp(1) = qab(19)
                cp(2) = qab(20)                
                qab(38) = cp(1)
                qab(39) = cp(2)
                ang = epsil(7)
                qab(40) = ang
                
            else

                call v_move(qab(23),dir(1),3)
                call v_move(qab(26),ta(1),3)
                call v_move(qab(29),tb(1),3)
                call v_move(qab(32),ac(1),3)
                call v_move(qab(35),bc(1),3)
                cp(1) = qab(38)
                cp(2) = qab(39)
                ang = qab(40)
                
                
            endif         
         endif
         
#endif defined(MOIST)
c --- skip if only updating connectivity
         if (C_ConnUpd.eq.1 .or. A_KeyOut.eq.1) then
             if (ield.ge.3.and.kfarcn.ne.1) goto 900
#if defined(MOIST)
             if ( (keyopt(1).eq.3 .or. keyopt(1).ge.5)
c --- only during first iteration
     &           .and. (A_IterNum.eq.1 .or. 
c --- or for the output
     &           (A_KeyOut.eq.1 .and. btest(eomask,W_OutMisc)) )) then
                 
                 mbeam = ibeam
                 sbeam = abs(isegda(16,iseg))
                 if (mbeam.eq.0.or.sbeam.eq.0) goto 500
c --- limit the angle
                 if (abs(ang).le.0.05d0 .or. 
     &               abs(ang-W_Pi).le.0.05d0) goto 500

c             if ( (keyopt(1).eq.3 .or. keyopt(1).ge.5) ) then
                 call GetBeamData(mbeam,cp(1),ac(1),bp(1,1))
                 call GetBeamData(sbeam,cp(2),bc(1),bp(1,2))
                 
                 if (abs(bp(9,1))+abs(bp(9,2)).le.TINY) goto 500

c                 if (mbeam.eq.522 .or. sbeam.eq.522) then
c                     write (8,*) 'elem,mbeam,sbeam',elem,mbeam,sbeam
c                     write (8,*) 'bp(9,1),bp(9,2)',bp(9,1),bp(9,2)
c                     write (8,*) 'xyz',x,y,z
c                 endif
                
c     & (L,beamprop(1)), (H,beamprop(2)), (W,beamprop(3)),
c     & (Ex,beamprop(4)), (Ez,beamprop(5)), (Gx,beamprop(6)),
c     & (alpx,beamprop(7)), (alpz,beamprop(8)), 
c     & (dh,beamprop(9))
                 
                 if (CheckHyg(bp(1,1),ang) .eq. 0) then
                     if (H_key(1) .eq. 1) then
c --- 1D                        
                         call BndHyg1D (bp(1,1),epsil(7),
     &                                  cp(1),em(1,1),MH)
                         
                     else
c --- 3D
                        if (dir(3).lt.ZERO) then
c --- reverse the sides if the normal is not right                            
                            call cv1_v1 (dir(1),3,-ONE)
                            call cv1_v1 (ta(1),3,-ONE)
                            call cv1_v1 (tb(1),3,-ONE)
                            call v_move(bp(1,1),vect(1),11)
                            call v_move(bp(1,2),bp(1,1),11)
                            call v_move(vect(1),bp(1,2),11)
                            call v_move(ta(1),vect(1),3)
                            call v_move(tb(1),ta(1),3)
                            call v_move(vect(1),tb(1),3)
                            call v_move(ac(1),vect(1),3)
                            call v_move(bc(1),ac(1),3)
                            call v_move(vect(1),bc(1),3)
                            
                            vect(1) = cp(1)
                            cp(1) = -cp(2)
                            cp(2) = -vect(1)
                            i = mbeam
                            mbeam = sbeam
                            sbeam = i
                         endif

c                         write(8,*) 'elem',elem,mbeam,sbeam
c                         write(8,*) 'bp',bp
c                         write(8,*) 'ang',ang
c                         write(8,*) 'ac',ac
c                         write(8,*) 'bc',bc
c                         write(8,*) 'dir',dir
c                         write(8,*) 'ta',tb
c                         write(8,*) 'tb',tb
c                         write(8,*) 'cp',cp
c                         write(8,*) 'elem',elem

c --- OBS
c                         dir(1) = 0.0d0
c                         dir(2) = 0.0d0
c                         dir(3) = 1.0d0
                         
                         call BndHygStrn (bp(1,1),ang,ac(1),bc(1),
     &                                    dir(1),ta(1),tb(1),cp(1),
     &                                    em(1,1),MH,wid(1))
c                         write(8,*) 'em',em
                     endif
                        
                 else
                     goto 500
                 endif

                 
                 if (A_KeyOut.eq.1) goto 501
                 if (tol.gt.ZERO .and. abs(ififo(21,ifoun)).gt.1) then
                     dum1 = abs(cp(1)) - ONE
                     if (abs(dum1).le.tol) then
c --- this is the case in which the contact points is on the border between two contact elements
                         dum2 = HALF*(tol-dum1)/tol
                         MH = MH * dum2
                         call cv1_v1(em(1,1),H_Col,dum2)
                     endif
                 endif
                 ind = (M_Beams(mbeam)-1)*H_Col*H_Row

                 do i=1,H_Row
!$omp atomic update                     
                    M_HygStrains(ind+i) = M_HygStrains(ind+i) + 
     &                                    em(i,1)
c                    if (abs(em(i,1)).gt.5.0d-2) then
c                        write (6,*) 'elemID 1', elem, em(i,1)
c                        write (6,*) 'bp', bp
c                        write (6,*) 'ang', ang
c                        write (6,*) 'stop'
c                        stop
c                    endif
                 enddo
                 sbeamn = 0
                 if (tol.gt.ZERO  .and. ififo(1,ifoun).gt.1) then
                     dum1 = abs(cp(2)) - ONE
                     if (abs(dum1).le.tol) then
c --- this is the case in which the contact points is on the border between two target elements
                         if (cp(2).gt.ZERO) then
c --- target element to the right
                             isegn = isegda(12,iseg)
                         else
c --- target element to the left
                             isegn = isegda(11,iseg)
                         endif
c --- get the element connected to it
                         if (isegn.gt.0) sbeamn = abs(isegda(16,isegn))

                         if (sbeamn.gt.0) then
c --- split between two beams
                             dum2 = HALF*(tol-dum1)/tol
                             dum3 = ONE - dum2
c                             write (6,*) 'case 2',dum2,dum3
                             call cv1_v2(em(1,2),em(1,3),H_Row,dum3)
                             call cv1_v1(em(1,2),H_Row,dum2)
                         endif
                     endif
                 endif
                 
                 ind = (M_Beams(sbeam)-1)*H_Col*H_Row
                 do i=1,H_Row
!$omp atomic update                     
                    M_HygStrains(ind+i) = M_HygStrains(ind+i) + 
     &                                    em(i,2)
c                    if (abs(em(i,2)).gt.5.0d-2) then
c                        write (6,*) 'elemID 2', elem, em(i,2)
c                        write (6,*) 'bp', bp
c                        write (6,*) 'ang', ang
c                        write (6,*) 'stop'
c                        stop
c                    endif
                 enddo
                 
                 if (sbeamn.gt.0) then
                     ind = (M_Beams(sbeamn)-1)*H_Col*H_Row
                     do i=1,H_Row
!$omp atomic update                         
                        M_HygStrains(ind+i) = M_HygStrains(ind+i) + 
     &                                        em(i,3)
c                        if (abs(em(i,3)).gt.5.0d-2) then
c                            write (6,*) 'elemID 3', elem, em(i,3)
c                            write (6,*) 'bp', bp
c                            write (6,*) 'ang', ang
c                            write (6,*) 'stop'
c                            stop
c                        endif

                     enddo
                 endif
                 if (keyopt_9.eq.2) then
                     ivolu = elmgct (elem,dvet(1))
                     dvet(7) = MH
                     call elmpct (elem,dvet(1))
                 endif

             endif
            
c             if (A_IterNum.eq.1) write(8,*) elem, em
c             if (A_IterNum.eq.2) write(9,*) elem, em
#endif defined(MOIST)
             goto 500
         endif
  501    continue
         
         qb(33,ipi) = ZERO
         if (kpress.eq.1) then
             if (ield.ge.3) then
                 qb(1,ipi) = qb(1,ipi) + prs(ipi)
                 qb(33,ipi) = prs(ipi)
             endif
         endif

         presb  = qb(1,ipi)
         taurb  = qb(2,ipi)
         tausb  = qb(3,ipi)

         elvol = elvol + dajcw
         if (ipass.eq.1) elvol = elvol + dajcw
         eEnergy(1) = eEnergy(1) + qb(16,ipi)*dajcw
         eEnergy(4) = eEnergy(4) + qb(17,ipi)*dajcw

c         if (ififo(12,ifoun).eq.28177) then
c             write (6,*) 'cinem(6),ield', cinem(6),ield
c         endif

         if (cinem(6).eq.ONE) goto 900
         if (ield.ge.3.and.kfarcn.ne.1) goto 900

         
c --- compute contact stresses in global axes
         fx = r(1,1)*taurb+r(2,1)*tausb+r(3,1)*presb
         fy = r(1,2)*taurb+r(2,2)*tausb+r(3,2)*presb
         fz = r(1,3)*taurb+r(2,3)*tausb+r(3,3)*presb
#if defined(MOIST)   
         if (keyopt_9.eq.2) then
              ivolu = elmgct (elem,dvet(1))
              MH = dvet(7)
              mx = r(3,1)*MH
              my = r(3,2)*MH
              mz = r(3,3)*MH
         endif
#endif defined(MOIST)   
c         if (ififo(12,ifoun).eq.28177) then
c             write (6,*) 'taurb, tausb, presb', taurb, tausb, presb
c         endif


c --- stiffness matrix in global axes

         if (kelin(1).ne.1) goto 500

         call v_zero (c1(1,1),9)
         if (ield.ge.3) goto 500

c --- product c * r
         call v_zero (cc(1,1),9)
         do 80 i = 1, 3
            do 75 j = 1, 3
               do 70 k = 1, 3
                  cc(i,j) = cc(i,j) + c(i,k) * r(k,j)
  70           continue
  75        continue
  80     continue

c                 t
c --- product r  d  r
         do 95 i = 1, 3
            do 90 j = 1, 3
               do 85 k = 1, 3
                  c1(i,j) = c1(i,j) + r(k,i) * cc(k,j)
  85           continue
  90        continue
  95     continue

 500     continue

         nno = 0

         if (irigf.ne.0 .and. (nocon.eq.0.or.kfarcn.eq.1)) then

            nno = 7
c --- in the case the slave nodes are reached
            if (nocon.ne.0.and.kfarcn.eq.1) then
                if (ififo(2,ifoun).gt.0) then
                    iseg = ififo(10,ifoun)
                    fifob(1) = ONE
                else
                    iseg = 0
                endif
            endif

            do 570 i = 1,nno
               noeud = isegda(2+i,iseg)
               if (noeud.le.0) goto 570
               k = 0
               if (nnodefu.eq.0) goto 565
               do 560 j = 1,nnodefu
                  if (noeud.eq.nofo(j)) then
                      k = j
                      goto 567
                  endif
  560          continue

  565          nnodefu = nnodefu + 1
               k = nnodefu
               nofo(k) = noeud
  567          kk(i) = k
  570       continue

            if (nnodefu+nnode.gt.W_NdMax.and.irigf.eq.2) then
               dperr(1) = nnodefu
               dperr(2) = elem
               dperr(3) = W_NdMax - nnode
               call erhandler('ContElement',5000,3,
     &           'We have found %i nodes of slave element are in
     &           contact with the element %i, which exceeds allowable
     &           number %i. ',dperr(1),' ')
            endif
         endif

         if (C_ConnUpd.eq.1) goto 900

         do 700 j = 1,nnode+nno
#if !defined(MOIST)
            if (kelin(6).ne.1) goto 575
#endif
            if (j.le.nnode.and.nodes(j).eq.0) goto 700
            if (j.gt.nnode) then
                if (kk(j-nnode).eq.0) goto 700
            endif
            if (j.le.nnode) then

c --- internal forces for beam 1
               findj = fi (j)
               dum1 = findj*dajcw
               
               if (kelin(6).ne.1) goto 571

c --- stress stiffness matrix

c --- internal forces for beam 1
               lmx = 1 + (j-1)*lintv
               lmy = lmx + kdsp
               lmz = lmy + kdsp

               Fnr(lmx) = Fnr(lmx) + fx*dum1
               Fnr(lmy) = Fnr(lmy) + fy*dum1
               Fnr(lmz) = Fnr(lmz) + fz*dum1

               if (krot.eq.1.and.ield.eq.1.and.keyopt_1.gt.0) then
                   rigpr = -rigpr0*dum1
                   Fnr(lmx+3) = Fnr(lmx+3) + rigpr*(rxfb-rxsb)
                   Fnr(lmx+4) = Fnr(lmx+4) + rigpr*(ryfb-rysb)
                   Fnr(lmx+5) = Fnr(lmx+5) + rigpr*(rzfb-rzsb)
               endif

 571           continue

#if defined(MOIST)
               if (kelin(5).ne.1) goto 671
c --- external forces for beam 1
               if (keyopt_9.eq.2 .and. krot.eq.1 .and. 
     &                            ield.eq.1.and.keyopt_1.gt.0) then
                   lmx = 1 + (j-1)*lintv
                   Fext(lmx+3) = Fext(lmx+3) + mx*dum1
                   Fext(lmx+4) = Fext(lmx+4) + my*dum1
                   Fext(lmx+5) = Fext(lmx+5) + mz*dum1
               endif
               
 671           continue               
#endif defined(MOIST)
            elseif (nocon.eq.0.or.kfarcn.eq.1) then

c --- internal forces for beam 2

               findj =  fifob(j-nnode)
               dum2  = findj*dajcw
               if (kelin(6).ne.1) goto 572
               if (irigf.eq.2) then
                   lnx =  nnodtv + (kk(j-nnode)-1)*lintv + 1
               else
                   lnx = nnodtv + 1
               endif
               lny = lnx + kdsp
               lnz= lny + kdsp
               Fnr(lnx)= Fnr(lnx) - fx*dum2
               Fnr(lny)= Fnr(lny) - fy*dum2
               Fnr(lnz)= Fnr(lnz) - fz*dum2

               if (krot.eq.1.and.ield.eq.1.and.keyopt_1.gt.0) then
                   rigpr = -rigpr0*dum2
                   Fnr(lnx+3) = Fnr(lnx+3) - rigpr*(rxfb-rxsb)
                   Fnr(lnx+4) = Fnr(lnx+4) - rigpr*(ryfb-rysb)
                   Fnr(lnx+5) = Fnr(lnx+5) - rigpr*(rzfb-rzsb)
               endif

 572           continue
#if defined(MOIST)
               if (kelin(5).ne.1) goto 672
c --- external forces for beam 2
               if (keyopt_9.eq.2 .and. krot.eq.1 .and. 
     &                            ield.eq.1.and.keyopt_1.gt.0) then
                   if (irigf.eq.2) then
                       lnx =  nnodtv + (kk(j-nnode)-1)*lintv + 1
                   else
                       lnx = nnodtv + 1
                   endif
                   Fext(lnx+3) = Fext(lnx+3) - mx*dum2
                   Fext(lnx+4) = Fext(lnx+4) - my*dum2
                   Fext(lnx+5) = Fext(lnx+5) - mz*dum2
               endif
               
 672           continue               
#endif defined(MOIST)

                  
            endif

 575        continue

            if (kelin(1).ne.1) goto 700
c --- stiffness matrix

            do 600 i = 1,nnode+nno
               
               if (i.le.nnode.and.nodes(i).eq.0) goto 600
               if (i.gt.nnode) then
                   if (kk(i-nnode).eq.0) goto 600
               endif
               if (j.le.nnode.and.i.le.nnode) then

                  findi = fi(i)
                  findj = fi(j)
                  dum1 = findi*findj*dajcw
                  dum2 = findi*weight

c --- mechanical...

                  lm1 = (j-1)*nrtv+lintv*(i-1)
                  if (ield.ge.3) goto 580
                  do 675 ik=1,3
                     do 673 jk=1,3
                        c(ik,jk) = c1(ik,jk)*dum1
  673                continue
  675             continue

                  do 679 jk = 1,3
                     lma = lm1 + nr*(jk-1)
                     do 678 ik = 1,3
                        lma = lma + kdsp
                        Kmtrx(lma) = c(ik,jk) + Kmtrx(lma)
  678                continue
  679             continue

c --- rotation
                  if (krot.eq.1.and.ield.eq.1.and.keyopt_1.gt.0) then
                      rigpr = rigpr0*dum1
                      do jk = 3,5
                         lmr = lm1 + 1 + jk*(1+nr)
                         Kmtrx(lmr) = Kmtrx(lmr) + rigpr
                      enddo
                  endif

  580             continue
               endif
               if (nocon.ne.0.and.kfarcn.ne.1) goto 600

               if (j.gt.nnode.and.i.le.nnode) then

                  findi = fi (i)
                  findj = fifob(j-nnode)
                  dum1 = findi*findj*dajcw

                  if (irigf.eq.2) then
                      lm1 = nnodtv*nr + (kk(j-nnode)-1)*nrtv
     &                                + lintv*(i-1)
                  else
                      lm1 = nnodtv*nr + lintv*(i-1)
                  endif
                  if (ield .ge. 3) goto 583

                  do 685 ik=1,3
                     do 683 jk=1,3
                        c(ik,jk) = c1(ik,jk)*dum1
  683                continue
  685             continue

                  do 689 jk = 1,3
                     lma = lm1 + nr*(jk-1)
                     do 688 ik = 1,3
                        lma = lma + kdsp
                        Kmtrx(lma) = - c(ik,jk) + Kmtrx(lma)
  688                continue
  689             continue

c --- rotation
                  if (krot.eq.1.and.ield.eq.1.and.keyopt_1.gt.0) then
                      rigpr = -rigpr0*dum1
                      do jk = 3,5
                         lmr = lm1 + 1 + jk*(1+nr)
                         Kmtrx(lmr) = Kmtrx(lmr) + rigpr
                      enddo
                  endif

 583              continue
               endif

               if (j.le.nnode.and.i.gt.nnode) then

                  findi = fifob(i-nnode)
                  findj = fi (j)
                  dum1 = findi*findj*dajcw
                  dum2 = findi*weight
                  if (irigf.eq.2) then
                     lm1 = (j-1)*nrtv+nnodtv +
     &                     (kk(i-nnode)-1)*lintv
                  else
                     lm1 = (j-1)*nrtv + nnodtv
                  endif
                  if (ield .ge. 3) goto 584

                  do 693 ik=1,3
                     do 691 jk=1,3
                        c(ik,jk) = c1(ik,jk)*dum1
  691                continue
  693             continue


                  do 695 jk = 1,3
                     lma = lm1 + nr*(jk-1)
                     do 694 ik = 1,3
                        lma = lma + kdsp
                        Kmtrx(lma) = - c(ik,jk) + Kmtrx(lma)
  694                continue
  695             continue

c --- rotation
                  if (krot.eq.1.and.ield.eq.1.and.keyopt_1.gt.0) then
                      rigpr = -rigpr0*dum1
                      do jk = 3,5
                         lmr = lm1 + 1 + jk*(1+nr)
                         Kmtrx(lmr) = Kmtrx(lmr) + rigpr
                      enddo
                  endif
                  
 584              continue
               endif
               if (j.gt.nnode.and.i.gt.nnode) then

                  findi = fifob(i-nnode)
                  findj = fifob(j-nnode)
                  dum1 = findi*findj*dajcw
                  dum2 = findj*dajcw
                  if (i.ne.(nnode+1)) dum2 = ZERO
                  if (irigf.eq.2) then
                     lm1 = nnodtv*nr + (kk(j-nnode)-1)*nrtv
     &                   + nnodtv    + (kk(i-nnode)-1)*lintv
                  else
                     lm1 = nnodtv*nr + nnodtv
                  endif
                  if (ield .ge. 3) goto 680

                  do 697 ik=1,3
                     do 696 jk=1,3
                        c(ik,jk) = c1(ik,jk)*dum1
  696                continue
  697             continue

                  do 699 jk = 1,3
                     lma = lm1 + nr*(jk-1)
                     do 698 ik = 1,3
                        lma = lma + kdsp
                        Kmtrx(lma) = c(ik,jk) + Kmtrx(lma)
  698                continue
  699             continue

c --- rotation
                  if (krot.eq.1.and.ield.eq.1.and.keyopt_1.gt.0) then
                      rigpr = rigpr0*dum1
                      do jk = 3,5
                         lmr = lm1 + 1 + jk*(1+nr)
                         Kmtrx(lmr) = Kmtrx(lmr) + rigpr
                      enddo
                  endif

  680             continue
c
               endif
  600       continue
  700    continue

c --- do not set contact states if contact damping is activated
         if (qb(5,ipi).eq.-2.d0) then
             qb(5,ipi) = 3.d0
         endif

  900 continue
      
      elcnup0 = 0

      if (nnodefu.gt.0.and.irigf.eq.2) then
         do 910 j = 1,nnodefu
            nofo(j) = nodbac((nofo(j)-1)/D_NumDof +1)
            if (nofo(j).ne.nodes(nnode+j)) then
                nodes(nnode+j) = nofo(j)
                elcnup0 = 1
            endif
 910     continue
         if (nnodefu+nnode.lt.nmax) then
             if (nodes(nnode+nnodefu+1).ne.0) then
                 elcnup0 = 1
                 call vi_zero(nodes(nnode+nnodefu+1),nmax-nnode-nnodefu)
             endif
         endif

         if (elcnup0.eq.1.and.A_KeyOut.ne.1) then
c             if (A_SubStepNum.gt.1) then
c                 write (6,*) 'update 3', elem
c             endif
             ContStat = 1
             call lockSMP(W_LockAnsys)
             call elmpnd (elem,nmax,nodes(1))
             call unlockSMP(W_LockAnsys)             
         endif
      endif

      if (kelin(1).eq.1) then
!dir$ ivdep
          do 994 i = 1,nr-1
              do 993 j = i+1,nr
                Kmtrx(nr*(j-1)+i) = HALF*(Kmtrx(nr*(j-1)+i) +
     &                                    Kmtrx(nr*(i-1)+j))
                Kmtrx(nr*(i-1)+j) = Kmtrx(nr*(j-1)+i)
  993        continue
  994     continue
      endif

  995 qab(7) = elvol
c --- chattering measurement (only when primary computations and not line search)
      if (C_ConnUpd.ne.1) then
         ifDebo = .False.
         ifFrac = .False.
         if (founda(33,ifoun).lt.qab(3)) founda(33,ifoun) = qab(3)
         do ipi =1,npi
            ififo(22,ifoun) = ififo(22,ifoun) + nint(qb(44,ipi))
            if (qb(nvmax-9,ipi).ge.ONED.and.qa(nvmax-9,ipi).lt.ONED)then
                ifDebo = .True.
            elseif (qb(nvmax-9,ipi).gt.TINY) then
                ifFrac = .True.
            endif
            if (abs(nint(qb(5,ipi))).le.2) then
                if (founda(32,ifoun).lt.qb(nvmax-4,ipi))
     &              founda(32,ifoun) = qb(nvmax-4,ipi)
                if (qb(nvmax-4,ipi).gt.TINY.and.
     &              founda(36,ifoun).gt.qb(nvmax-4,ipi))
     &              founda(36,ifoun) = qb(nvmax-4,ipi)
                if (founda(37,ifoun).lt.qb(nvmax-3,ipi))
     &              founda(37,ifoun) = qb(nvmax-3,ipi)
                if (qb(nvmax-3,ipi).gt.TINY.and.
     &              founda(28,ifoun).gt.qb(nvmax-3,ipi))
     &              founda(28,ifoun) = qb(nvmax-3,ipi)
                if (abs(founda(34,ifoun)).lt.abs(qb(1,ipi)))
     &              founda(34,ifoun) = qb(1,ipi)
                slides = sqrt(qb(2,ipi)*qb(2,ipi)+qb(3,ipi)*qb(3,ipi))
                if (abs(founda(35,ifoun)).lt.abs(slides)) then
                    founda(35,ifoun) = slides
                endif
                slides = sqrt(qb(29,ipi)*qb(29,ipi) + 
     &                        qb(30,ipi)*qb(30,ipi))
                if (founda(31,ifoun).lt.slides) then
                    founda(31,ifoun) = slides
                endif
                if (abs(nint(qb(5,ipi))).eq.1) then
                    if (ififo(23,ifoun).eq.0) founda(30,ifoun) = ZERO
                    slides = sqrt(qb(12,ipi)*qb(12,ipi) +
     &                            qb(13,ipi)*qb(13,ipi))
                    if (founda(30,ifoun).lt.slides) then
                        founda(30,ifoun) = slides
                        ififo(23,ifoun) = elem
                    endif
                endif
            endif
         enddo
         if (ifDebo .or. ifFrac) then
             call lockSMP (W_LockConvVar)
             if (ifFrac) then
                 C_nFraBo = C_nFraBo + 1
             elseif (ifDebo) then
                 C_nDeBo = C_nDeBo + 1
c                 write (8,*) 'ELEM', elem, qb(1:3,1)
c                 write (8,*) 'qb(nvmax)', qb(nvmax-13:nvmax,1)
c                 write (8,*) 'qb(nvmax)', qb(nvmax-13:nvmax,1)
c                 write (8,*) 'Fnr(nr)', Fnr(1:nr)
             endif             
             call unlockSMP (W_LockConvVar)
         endif    
      endif

c --- save final contact status for removal of loose elements
      ield = 4
      icnda(3,icon) = 0
      do ipi = 1,npi
         ield = min(ield,nint(qb(5,ipi)))
         icnda(3,icon) = max(icnda(3,icon),nint(qb(8,ipi)))
      enddo

      if (ield.le.2) then
c --- mark contact
          if (ififo(26,ifoun).lt.0) then
              call lockSMP (W_LockConvVar)
              ififo(26,ifoun) = 0
              call unlockSMP (W_LockConvVar)
          endif
          iseg = icnda(3,icon)
          isegda(19,iseg) = icon
      else
          if (ield.eq.3) then
c --- mark near contact only if the pair is not in contact
              if (ififo(26,ifoun).eq.-2) then
                  call lockSMP (W_LockConvVar)
                  ififo(26,ifoun) = max(ififo(26,ifoun),-1)
                  call unlockSMP (W_LockConvVar)
              endif
          endif
          iseg = icnda(3,icon)
c --- update the pair connection only if solution converged
          if (A_KeyConv.eq.1 .and. A_KeySvrUpd.eq.1) then
              icnda(3,icon) = 0
c --- zero out segment info, if it has the current contact element in
              if (iseg.gt.0) then
                  if (isegda(19,iseg).eq.icon) isegda(19,iseg) = 0
              endif
          endif
      endif

c --- perform time predictions
      if(A_KeyAutoStep.eq.1 .and. A_KeyMatUpd.eq.1 .and. kpred.ge.2)then
         call ContPredi (elem,kpred-1,nvmax,npi,qa(1,1),qb(1,1),mu,
     &                   ipair,cinem(13))
      endif

c --- update history variables
      if (A_KeyMatUpd .eq. 1) then
          do 998 ipi = 1,npi
             qa(nvmax-6,ipi) = qa(5,ipi)
             qa(nvmax-5,ipi) = qa(9,ipi)
             qa(nvmax-4,ipi) = qa(1,ipi)
             qa(nvmax-3,ipi) = qa(12,ipi)
             qa(nvmax-2,ipi) = qa(13,ipi)
             qa(nvmax-1,ipi) = qa(2,ipi)
             qa(nvmax,ipi) = qa(3,ipi)
             do 997 i = 1,nvmax-8
                qa(i,ipi) = qb(i,ipi)
  997        continue
  998     continue
      endif
#if defined(MOIST)
c --- update past moisture
      if (A_KeyBCsUpd.eq.1) then
          do ipi = 1,npi
             qa(34,ipi) = moistend(ipi)
          enddo
       else
          do ipi = 1,npi
             qa(34,ipi) = moistbeg(ipi)
          enddo
       endif
#endif       
#if defined(PRS)
c --- update past pressure
      if (kpress.ne.0) then
          if (A_KeyBCsUpd.eq.1) then
              do ipi = 1,npi
                 qa(33,ipi) = qb(33,ipi)
              enddo
              if (nnodmst.eq.3) qb(33,4) = qb(33,3)
          else
              do ipi = 1,npi
                 qa(33,ipi) = prsbeg(ipi)
              enddo
          endif
      endif
#endif
c               if (ield.le.2) then
c                   write (8,*) ififo(12,ifoun), ield, x,y,z
c               endif
               
c --- produce output
      if (A_KeyOut.eq.1 .and. ield.le.3 .and. cinem(6).eq.ZERO) then

         if (btest(eomask,W_OutElIndx)) then
             if (btest(eomask,W_OutMisc)) then
c --- output contact status through misc data
               call v_zero (vect(1),elChar(NMSMIS))
               vect(1) = dble(ield)
c --- these is damage output
               vect(2) = qb(nvmax-9,1)
#if !defined(DEBUG)
               if (ield.le.2) then
                   vect(3) = x
                   vect(4) = y
                   vect(5) = z
               endif    
c --- save energy
               vect(6) = eEnergy(1)
               vect(7) = ibeam
               if (iseg.gt.0) vect(8) = abs(isegda(16,iseg))

c --- save slide distances for Hanna
               vect(9) = sqrt(qb(29,npi)**2 + qb(30,npi)**2)
c --- contact forces for Svetlana
               vect(10) = qb(1,1)
               vect(11) = sqrt(qb(2,1)*qb(2,1)+qb(3,1)*qb(3,1))
c --- widths for Kurosh
#if defined(MOIST)
               vect(12) = wid(1)
               vect(13) = wid(2)
               if (ield.le.2) then
                   vect(14) = qab(19)
                   vect(15) = qab(20)              
               endif
#endif               
c            if (qb(nvmax-9,ipi).ge.ONED .or. qa(nvmax-9,ipi).ge.ONED) 
c     &          write (6,*) 'Elem sepa', elem,ield,vect(10:11)
               
#else
c --- debug output
               vect(3) = qb(nvmax-8,1)
               vect(4) = qb(nvmax-11,1)
               do ipi = 2,npi
                  vect(3) = max(vect(2),qb(nvmax-8,npi))
                  vect(4) = max(vect(3),qb(nvmax-11,npi))
               enddo

c --- max. chattering level
               vect(5) = dble(ififo(22,ifoun))

c --- max. penetration / min. gap
               if (ififo(16,ifoun).eq.0) then
                   if (founda(33,ifoun).gt.TINY) then
                       vect(6) = -founda(33,ifoun)
                   else
                       vect(6) = ZERO
                   endif
               else
                   vect(6) = -founda(9,ifoun)
               endif

c --- max. normal stiffness
               if (founda(32,ifoun).gt.TINY) then
                   vect(7) = founda(32,ifoun)
               else
                   vect(7) = ZERO
               endif


c --- min. normal stiffness
               if (founda(36,ifoun).gt.TINY .and. 
     &             founda(36,ifoun).lt.1.d25) then
                   vect(8) = founda(36,ifoun)
               else
                   vect(8) = ZERO
               endif

c --- max. resulting pinball
               if (founda(33,ifoun).gt.TINY) then
                   vect(9) = founda(33,ifoun)
               else
                   vect(9) = ZERO
               endif
c --- max. elastic slip distance
               if (founda(30,ifoun).gt.TINY) then
                   vect(10) = founda(30,ifoun)
               else
                   vect(10) = ZERO
               endif
c --- max. tangential stiffness
               if (founda(37,ifoun).gt.TINY) then
                   vect(11) = founda(37,ifoun)
               else
                   vect(11)= ZERO
               endif

c --- min. tangential stiffness
               if (founda(28,ifoun).gt.TINY .and. 
     &             founda(28,ifoun).lt.1.d25) then
                   vect(12) = founda(28,ifoun)
               else
                   vect(12) = ZERO
               endif
c --- max. sliding distance
               if (founda(31,ifoun).gt.TINY) then
                   vect(13) = founda(31,ifoun)
               else
                   vect(13) = ZERO
               endif
c --- max. contact pressure
               if (abs(founda(34,ifoun)).gt.TINY) then
                   vect(14) = founda(34,ifoun)
               else
                   vect(14) = ZERO
               endif
c --- max. friction stress
               if (founda(35,ifoun).gt.TINY) then
                   vect(15) = founda(35,ifoun)
               else
                   vect(15) = ZERO
               endif
c --- average. contact depth
               if (founda(17,ifoun).gt.TINY) then
                   vect(16) = founda(17,ifoun)
               else
                   vect(16) = ZERO
               endif
#endif
               i = elChar(NMSMIS)
               call WriteToRST (EDEMS,lcerstL,edindxL(1),i,vect(1))
             endif
         endif
      endif

  999 continue

      if (A_KeyFirstPass.eq.1 .and. keyopt(1).eq.6 .and. ield.gt.2) then
c --- mark as killed
          elData(EL_DEAD) = 1
          call vi_zero (nodes(nnode+1),nmax-nnode)
          call lockSMP(W_LockAnsys)
          call elmput(elem,elData(1),nnode,nodes(1))
          call unlockSMP(W_LockAnsys)
          ContStat = 1
      endif

c      write (8,*) 'elem', elem, epsil(3), nocon, nnode, nno
c      write (8,*) 'nodes', nodes(1:nnode+nno)

c      if (kelin(6).eq.1) then
c          write (6,*) 'Force'
c          write (6,*) Fnr(1:nr)
c      endif

c      if (kelin(1).eq.1) then
c          write (6,*) 'Matrix'
c          write (6,*) Kmtrx(1:nr*nr)
c      endif
c      if (ContStat.eq.1) write (6,*) 'Changed', elem

#if defined(DEBUG)
      call SubEnd('ContElement')
#endif

      return
      end
