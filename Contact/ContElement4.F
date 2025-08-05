      subroutine ContElement4 (elem,elChar,elData,eomask,nodes,
     &                         kelrqq,nr,nvmax,npi,ival,nSeg,nFoun,
     &                         nPerSeg,nPerFounI,nPerFounDp,nPerMast,
     &                         nMast,isegda,ififo,founda,icnda,xyz,u,
     &                         dispn,coor0,nodbac,qa,qb,qab,kelout,
     &                         Kmtrx,Fnr,elvol,center,eEnergy,ContStat,
     &                         edindxL,lcerstL)
c --- primary function: compute element matrices, load vectors, and results for beam-surface element
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
c
c                                  --> nONE adjacent segment
c       coor0  (dp,ar(*),in)       - initial nodal coordinate coordinate
c       dispn  (dp,ar(*),in)       - current nodal displacement
c       nodbac (int,ar(lenbac),out)- internal to external list
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
c       qb   (7) =        - information concerning mechanical dissipation
c                           into heat flow
c       qb   (8) = iseg   - attached segment number calculated by CalPeni4
c       qb   (9) = db     - total penetration distance calculated by CalPeni4
c       qb  (10) = nocot  - calculated by CalPeni4 for initial situation
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
c          qab(2) = thickl   - contact thickness
c          qab(3) = pinb     - pinball
c          qab(4) = length   - contact length
c          qab(5) = Ebeam    - Young's modulu for attached beam element
c          qab(7) = elvol    - contact element volume
c          qab(9) =          - attached element shape factor
c          qab(13) =         - avegage area for contact stiffness adjustements
c
c
c       kelout (int,ar(10),out)     - keys indicating created mats and lvs
c       Kmtrx  (dp,ar(nr,nr),inout) - k matrix (kelrqq(1))
c       Fnr  (dp,ar(nr),out)        - n-r restoring force vector
c       elvol  (dp,sc,out)          - element volume
c       center (dp,ar(3),out)       - centroid location
c       eEnergy (dp,ar(5),out)      - element energies
c       edindxL(LONG,ar(25),out)    - element result data file indexes
c       lcerstL(LONG,sc,inout)      - position on result file
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
#include "solu_info.inc"
#include "conv_stat.inc"
#include "in_mem.inc"
#if defined(MOIST2)
#include "beam_hyg.inc"
#endif
c
      integer  elmgct,elmget,ContAdjStat,rlget,CheckHyg

      integer*8  lcerstL,edindxL(25),MakeLongInt,keyL
      integer
     &  elem,elChar(IELCSZ),elData(EL_DIM),eomask,nodes(*),nPerSeg,nSeg,
     &  nPerFounI,nFoun,kelrqq(10),nr,kelout(10),nPerMast,nMast,
     &  nvmax,ContStat,npi,keyopt(12),ival(3),
     &  isegda(nPerSeg,nSeg),methi,nodbac(*),ififo(nPerFounI,nFoun),
     &  nod,nPerFounDp,kerr,i2key,elmgcn,icnda(nPerMast,nMast),icon,
     &  elcnup0,nmax,hi

      double precision
     &  xyz(6,*),u(nr,5),Kmtrx(nr*nr),Fnr(nr),eEnergy(5),center(3),
     &  qa(nvmax,npi),qb(nvmax,npi),elvol,coor0(*),dispn(*),
     &  founda(nPerFounDp,nFoun),qab(21),mudat(9),pinb2,
     &  dx,curvature,dfb(10),
     &  bp(9,2),em(2,3)
      
      logical lOrthoFric,lTrans,ifDebo,ifConnected

      integer kelin(10),nofo(72),kk(8),elDatab(EL_DIM)

      integer
     &  mat,ireal,irealb,nrvr,nssvr,i,khygro,
     &  nnode,ifoun,ibeam,sbeam,irigf,intfer,
     &  inopil,iropil,kddlx,kddly,kddlz,lintv,nnodefu,ipi,ield,
     &  iseg,nocon,lmx,nno,icode,noeud,k,j,n,lmy,lnx,lny,imod,
     &  lna,lm1,inod,iddl,ickcel,lp(2),
     &  ieqnr,ivolu,kpred,isega,nnodmst,cnnode,
     &  lmz,lnz,ik,jk,lma,indic,ipass,
     &  lmr,lmf,lmo,lmh,nrtv,nnodtv,
     &  nodesb(W_NdMax),ksmsld,ipene,ndin,
     &  kpress,inode,ipair,kdsp,kfarcn,iseg0,ickcel1,
     &  keyopt_3,krot,lintvr,numpt,option,keyopt_1,aupd,
     &  isegn,ind,sbeamn
           
      double precision
     &  ZERO,HALF,ONE,TWO,
     &  avemoist,length,gapi,xpilo,ypilo,zpilo,
     &  mu,thickl,rigp,Ebeam,pinb,gapcri,
     &  da,dajct,x,y,z,vx0,vy0,vz0,refmoist,
     &  dajc,db,presb,weight,fx,fy,
     &  findj,vtfra,vtfrb,vtfsa,vtfsb,vect(20),
     &  dxj,dyj,dzi,dzj,ddum,dxi,dyi,xmax,ymax,
     &  xmin,ymin,zmax,zmin,fpin,GetFibRad,
     &  taurb,tausb,fz,dum1,dum2,dum3,findi,dajcw,
     &  mux,muy,smallmu,GetMoistRef,crad,rad(2),radi

      parameter (smallmu=0.001d0)

      double precision
     &  rvr(W_NumContParam),rvrb(W_NumContParam),
     &  moistbeg(8),moistend(8),moistcur(8),
     &  xyt(6,7),xytp(6,7),vxyt(3,8),fifoa(7),fifob(7),df(6),
     &  epsil(8),cinem(43),c(4,4),dvet(11),fi(3,3),u0(W_NrMax),
     &  dperr(4),fcell(6),r(3,3),dispx,dispy,dispz,prs(9),prsbeg(9),
     &  c1(3,3),c2(3,3),cr(3,3),ct(3,3),prsend(9),
     &  cc(3,3),slides,dxk,dyk,dzk,cnof,
     &  tols,hygro(16,9),tol,MH,wid(2),
     &  qab1,trgshft,rxsb,rysb,rzsb,rigpr,rxfb,ryfb,rzfb,
     &  dpres,props(6),natcor(5),vecnor(2),rigpr0,ONED

      data lp /16,17/
      
      parameter (ZERO=0.0d0,ONE=1.0d0,HALF=0.5d0,TWO=2.0d0,ONED=1.0d0)

#if defined(DEBUG)
      call SubBeg('ContElement4')
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
      
c --- shift to account for beam radius
      crad = ZERO     

c --- define initial data
      mat = elData(EL_MAT)
      ireal = elData(EL_REAL)
      nrvr = elChar(NMDRLC)
      nssvr = elChar(NMSSVR)
c --- initialize
      call v_zero(cinem(1),43)
      
c --- get keyopts
      call vi_move (elChar(KYOP1),keyopt(1),12)

      nnode = 3      
      if (ififo(15,ifoun).ge.1) nnode = 7
c --- get the element real constant data
      ipene = nvmax - 7
      call v_zero (rvr(1),W_NumContParam)
      j = rlget(ireal,rvr(1))
      keyopt_3 = keyopt(3)
      if (keyopt_3.gt.2) keyopt_3 = keyopt_3 - 3
      hi = 1
      if (keyopt(12).eq.1) hi = 0
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
     &        founda(17,ifoun) = qab(2) 
      endif
      founda(29,ifoun) = ONE 

      call v_zero (u0(1),nr)
      if (A_RowIncS.ne.0) then
          call v_move (u(1,A_RowIncS),u0(1),nr)
      endif

c --- element node number
      nnodmst = nnode
      nmax = elChar(NMNDST)
      
c --- DOF compoments
      kdsp = 1
      methi = keyopt(8)

c --- slave type: irigf=1 rigid slave
      irigf = 0
      if (ififo(15,ifoun).eq.-1.or.ififo(15,ifoun).le.-3) irigf = 1
      if (irigf.eq.1.and.nodes(nnode+1).le.0) irigf = 0
      if (ififo(15,ifoun).ge.1) irigf = 2

c --- pilot node info
      if (ifoun .le. 0) then
         inopil = 1 - D_NumDof
         iropil = 0
      else
         inopil = ififo(2,ifoun)
         iropil = ififo(3,ifoun)
      endif

      intfer = 0
      if ((methi.eq.2.or.methi.eq.4).and.A_LoadStepNum.eq.1) intfer = 1

c --- error key
      kerr = 2
      kpred = keyopt(7)

      if (A_KeyAutoStep.eq.1.and.A_TimeInc.gt.TWO*A_MinTimeStep.and.
     &    kpred.ge.1) kerr = 3
      if (A_KeyOut.eq.1.or.A_KeyMatUpd.eq.1) kerr = 2

      if (A_KeyFirstPass.eq.1) qab(10) = 0

      if (ibeam.gt.0.and.ibeam.ne.elem) then
c --- get attached beam element radius (only for surface contact)
          if (keyopt(11).eq.1 .and. keyopt(10).eq.1) then
              crad = GetFibRad(ibeam)
          endif
      endif

c --- contact configuration update key
      imod = 1
      if (A_IterNum.ne.1) imod = 0
      if (A_KeyOut.eq.1.or.A_KeyMatUpd.eq.1) imod = 0

      ksmsld = 0
      if (A_KeyAutoStep.eq.1.and.kpred.ne.0) ksmsld = 1
      if (C_Disp*C_Forc.eq.1) ksmsld = 0
      if (A_IterNum.eq.1) ksmsld = 0
      if (imod.eq.1) ksmsld = 0  

c --- slave edge tolerence to prevent slipping off
      if (rvr(25).lt.TINY) then
          tols = TWO
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
     
c --- element center and length based on corner nodes

      center(1) = (xytp(1,1)+xytp(1,2))*0.5d0
      center(2) = (xytp(2,1)+xytp(2,2))*0.5d0
      center(3) = (xytp(3,1)+xytp(3,2))*0.5d0

c --- get moistures
      khygro = 0
      refmoist = ZERO
      avemoist = ZERO
#if defined(MOIST)
      call v_zero (moistcur(1),7)
      moistbeg(1) = qa(34,1)
      moistbeg(2) = qa(34,2)
      moistbeg(3) = qa(34,1)
      moistbeg(4) = qa(34,1)
      moistbeg(5) = qa(34,1)
      moistbeg(6) = qa(34,1)
      moistbeg(7) = qa(34,1)
      moistbeg(8) = qa(34,1)
      refmoist = GetMoistRef(mat)         
      call GetMoist (elem,elChar(1),npi,nodes(1),refmoist,
     &               1,moistbeg(1),moistcur(1),moistend(1),khygro,
     &               vect(1))
      do ipi =1,npi
         qb(34,ipi) = moistcur(ipi)
      enddo
c --- average nodal moisture
      avemoist = ZERO
      do 01 i = 1,nnode
         if (nodes(i).eq.0) goto 01
         if (i.le.nnodmst) avemoist = avemoist + moistcur(i)
  01  continue
      avemoist = avemoist/nnodmst
#endif

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
      elvol = ZERO
c      write (8,*) elem,icnda(3,icon),C_ConnUpd,A_KeySvrUpd
      aupd = 0
c --- check if it is right after contact update
      if  (C_ConnUpd.eq.0 .and. A_KeyOut.ne.1 .and.
     &     A_KeyFirstIter.eq.0 .and. A_KeySvrUpd.eq.1) aupd = 1
      if (icnda(3,icon).eq.0 .and. aupd.eq.1) then
c --- the element status checked above comes from the connectivity update
c          write (8,*) elem
          goto 995
      endif
      
      kpress = 0
#if defined(PRS)
c --- pressures field
      do ipi = 1,npi
         prsbeg(ipi) = qa(33,ipi)
      enddo
      call prsget (elem,elChar(1),1,4,prsbeg(1),prs(1),prsend(1),kpress)
#endif
      if (kpress.eq.0) then
          call v_zero (prs(1),npi)
      else
          ddum = ZERO
          do ipi = 1,npi
             if (ipi.le.4) then
                 ddum = ddum + prs(ipi)
             else
                 prs(ipi) = ddum/4.d0
             endif
          enddo         
          if (ddum.lt.TINY) then
              kpress = 0
              call v_zero(prs(1),npi)
          endif          
      endif

      kfarcn = kpress

c --- update contact for each sub-step
      if (A_IterNum.eq.1.or.keyopt_3.eq.2) then
          qab(6) = ZERO

c --- first pass
         if (A_KeyFirstPass.eq.1) then
             qab(2) = founda(17,ifoun)
             qab(3) = founda(20,ifoun)
             qab(4) = founda(16,ifoun)
             qab(5) = founda(18,ifoun)
             qab(11) = founda(18,ifoun)*1.0d-3/founda(17,ifoun)
             qab(12) = founda(22,ifoun)
             qab(13) = founda(21,ifoun)/founda(19,ifoun)
             
             if (0.eq.0) then
                qab(1) = qab(5)*qab(2)*10.d0
             else
c               qab(1) = qab(5)*qab(2)/qab(13)
                qab(1) = qab(5)/qab(2)
             endif

             if (keyopt(3).gt.2) then
c --- using local stiffness
                 ivolu = elmgct (elem,dvet(1))
                 qab(5) = dvet(5)
                 if (ififo(15,ifoun).eq.3) qab(5) = 0.5d0*qab(5)
             endif
          elseif (founda(17,ifoun).gt.TINY.and.A_KeyFirstPass.eq.0.and.
     &           keyopt_3.ge.1) then

             qab1 = 100.d0*founda(12,ifoun)/founda(17,ifoun)
             if (0.eq.0) then
                 qab1 = qab1*qab(13)
             endif
             qab1 = max (qab1,0.1d0*qab(1))
             qab1 = min (qab1,10.d0*qab(1))
             qab(1) = qab1
         endif

         if (keyopt(3).gt.2) then

            thickl = min(abs(rvr(1)),abs(rvr(2)))
            if (thickl.lt.TINY)thickl = max(abs(rvr(1)),abs(rvr(2))) 
            qab(2) = thickl
c --- automatically calculate penetration coefficient
            if (0.eq.0) then
                rigp = qab(5)*qab(2)*10.d0
            else
                rigp = qab(5)*qab(2)/qab(13)
            endif            
            if (ibeam.ne.elem.and.ibeam.gt.0.and.keyopt(3).ge.4) then
                ivolu = elmgct (ibeam,dvet(1))
c --- using small contact stiffness
                dpres = 1.d-3
                if (keyopt(1).lt.2) then
                    if (keyopt(6).eq.1) dpres = 1.d-6
                    if (keyopt(6).ge.2) dpres = 1.d-9
                endif
                dvet(6) = max(dvet(6),dpres*qab(5))
                qab1 = 100.d0*founda(12,ifoun)/qab(2)*qab(18)
                if (0.eq.0) then
                    qab1 = qab1*qab(13)
                endif
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
              if (keyopt(6).eq.1) dpres = 1.0d-6
              if (keyopt(6).ge.2) dpres = 1.0d-9
          endif
          dvet(6) = max(dvet(6),dpres*qab(5))
          if (keyopt_3.eq.1) 
     &        founda(11,ifoun) = max(founda(11,ifoun),dvet(6))
          if (keyopt(3).eq.4) qab(16) = dvet(6)
      endif
      if (ififo(15,ifoun).lt.1.and.C_ConnUpd.eq.1) goto 999

c --- small sliding case
      if (ksmsld.eq.1.and.kfarcn.ne.1) then
          do 11 ipi = 1,npi
             if (nint(qa(5,ipi)).ne.4) goto 12
             qb(5,ipi) = 4.d0
             qb(9,ipi) = founda(20,ifoun)
  11      continue
          elvol = qab(7)
          if (irigf.eq.2.and.nodes(nnode+1).ne.0.and.A_KeyOut.ne.1) then
              call vi_zero (nodes(nnode+1),nmax-nnode)
              ContStat = 1
              call lockSMP(W_LockAnsys)
              call elmpnd (elem,nmax,nodes(1))
              call unlockSMP(W_LockAnsys)
          endif
          goto 995 
      endif
  12  continue

c --- re-issue all parameters back
      rigp = qab(1)
      pinb = qab(3)
      thickl = qab(2)
      length = qab(4) 
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
      
      Ebeam = qab(5)

      if (A_KeyMatUpd.eq.1.and.keyopt(3).gt.2) then      
          ivolu = elmgct (elem,dvet(1))
          dvet(6) = thickl
          call elmpct (elem,dvet(1))
      endif

      gapcri = rvr(5)
      if (abs(gapcri).le.TINY.and.keyopt(5).eq.0) gapcri = qab(12)
      if (abs(gapcri).le.TINY.and.keyopt(5).eq.4.and.keyopt(1).eq.6)
     &    gapcri = qab(12)
      cnof = rvr(10)
      if (keyopt(5).gt.0.and.keyopt(5).lt.4) cnof = rvr(10)+qab(12)

      call vi_zero (kk(1),8)
      call vi_zero (nofo(1),npi*8)
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
      
c --- check for debonding in case of bonded contact;
c     store cohesive zONE material data for debonding in
c     cinem(36-41)
c
      if (keyopt(1).ge.2 .and. keyopt(1).le.6 .and. ifConnected) then
         option = -1
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
         elseif (option.gt.0) then
            dperr(1) = elem
            dperr(2) = ireal
            call erhandler ('ContElement4',5000,3,'Invalid
     &                       model option for debonding
     &                       %I (real ID %I)',dperr(1),' ')
         endif
      endif

      if ((keyopt(1).eq.3 .or. keyopt(1).ge.5) .and. ifConnected) then
          lOrthoFric = .false.
          lTrans = .false.
          mu = 0.0d0
          goto 100
      endif
      
c --- get friction coefficient--isotropic friction
      call GetMatProps (elem,mat,14,avemoist,mu,1)
      if (mu .lt. 0.0d0) mu = 0.0d0
c --- get friction coefficient--orthotropic friction
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
      cinem(32) = ONE

c --- initiation of shape functions , their derivates ,
c     and integration weights

      call v_zero (fi(1,1),9)

      fi(1,1) = ONE
      fi(2,2) = ONE
      fi(3,3) = ONE

      if (ififo(16,ifoun).eq.0.and.ififo(17,ifoun).eq.0.and.
     &    (C_ConnUpd.eq.0.and.
     &     ififo(15,ifoun).ge.1.and.A_KeyMatUpd.ne.1)) then
          ickcel = 1
          goto 03
      endif

c --- check whether the current element is in the foundation cell
      ickcel = 0
      if (ifoun .le. 0) then
          call v_zero (fcell(1),6)
          ickcel = 1
      else
          xmax = xytp(1,1)
          ymax = xytp(2,1)
          zmax = xytp(3,1) 
          xmin = xytp(1,1)
          ymin = xytp(2,1) 
          zmin = xytp(3,1) 
          do 02 i = 2,3
             if (nodes(i).eq.0) goto 02
             xmax = max (xytp(1,i),xmax)
             xmin = min (xytp(1,i),xmin)
             ymax = max (xytp(2,i),ymax)
             ymin = min (xytp(2,i),ymin)
             zmax = max (xytp(3,i),zmax)
             zmin = min (xytp(3,i),zmin)
  02      continue
          pinb2 = pinb
          if (pinb2.le.TWO*length) pinb2 = TWO*length
          fcell(1) = founda(1,ifoun) - pinb2
          fcell(2) = founda(2,ifoun) + pinb2
          fcell(3) = founda(3,ifoun) - pinb2
          fcell(4) = founda(4,ifoun) + pinb2
          fcell(5) = founda(5,ifoun) - pinb2
          fcell(6) = founda(6,ifoun) + pinb2
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
          if (xmax .lt. fcell(1)) ickcel = 1
          if (xmin .gt. fcell(2)) ickcel = 1
          if (ymax .lt. fcell(3)) ickcel = 1
          if (ymin .gt. fcell(4)) ickcel = 1
          if (zmax .lt. fcell(5)) ickcel = 1
          if (zmin .gt. fcell(6)) ickcel = 1 
      endif

  03  continue

c --- out of the cell region, calculate element volume and return
      if (ickcel.eq.1 .and. kfarcn.ne.1) then
         do 05 ipi = 1,npi
            if (A_KeyFirstPass.ne.1) then
                elvol = qab(7) 
                qb(11,ipi) = qa(11,ipi)
                goto 04
            endif
c --- coordinates and velocities in the initial configuration 
            dajc = ONE
            elvol = length
            qb(11,ipi) = dajc
  04        continue 
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
     &                         abs(abs(ield)-abs(qb(5,ipi)))
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
                        call erhandler ('ContElement',5000,kerr,'Contact
     &                  element %I (ID %I) status changed suddenly from
     &                  contact (with slave element %I) to no-contact.',
     &                  dperr(1),' ')
                    endif
                endif
            endif
 05      continue
         if (irigf.eq.2.and.nodes(nnode+1).ne.0.and.A_KeyOut.ne.1) then
             call vi_zero (nodes(nnode+1),nmax-nnode)
             ContStat = 1
             call lockSMP(W_LockAnsys)
             call elmpnd (elem,nmax,nodes(1))
             call unlockSMP(W_LockAnsys)
         endif 
         goto 995
      endif

      do 06 i = 1,3
         if (nodes(i).eq.0) goto 06
         vxyt(1,i) = xytp(1,i) - xyt(1,i)
         vxyt(2,i) = xytp(2,i) - xyt(2,i) 
         vxyt(3,i) = xytp(3,i) - xyt(3,i)
 06   continue 

      eEnergy(1) = ZERO
      eEnergy(4) = ZERO
      elvol = ZERO
      cinem(8) = ZERO
      cinem(9) = ONE
      cinem(10)= mu
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

      if (inopil .gt. 0) then
          kddlx = inopil
          kddly = kddlx + kdsp
          kddlz = kddly + kdsp
          if (iropil.le.0) goto 07
          dispx = dispn(kddlx)
          dispy = dispn(kddly)
          dispz = dispn(kddlz)
          xpilo = dispx + coor0(kddlx)
          ypilo = dispy + coor0(kddly)
          zpilo = dispz + coor0(kddlz)
 07       continue
      endif

c --- loop on integration points

      nnodefu = 0
      if (ickcel.eq.1) imod = 0

      i2key = elmgcn(elem,vecnor(1))
            
c --- get contact normal
      do 900 ipi=1,npi

         ickcel1 = ickcel
         inode = ipi
         if (nodes(inode).eq.0) goto 900

         isega  = nint(qa(8,ipi))
         iseg  = nint(qb(8,ipi))
         if (ipi.eq.1) then
             iseg0 = iseg
         else
             if (isega.eq.0) isega = iseg0
         endif
         dajct = qa(11,ipi)
         cinem(6) = ZERO
         cnnode = nodes(inode)
         if (ipi.eq.3) then
             cinem(6) = ZERO
         elseif (i2key.gt.0) then
             if (nint(vecnor(ipi)).ne.elem) cinem(6) = ONE
         endif

c --- small sliding case
         if (ksmsld.eq.1.and.kfarcn.ne.1) then
            if (nint(qa(5,ipi)).eq.4) then
                elvol = length
                qb(21,ipi) = ZERO
                qb(5,ipi) = 4.0d0
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
            x = xyt(1,ipi)
            y = xyt(2,ipi)
            z = xyt(3,ipi)

            vx0 = vxyt(1,ipi)
            vy0 = vxyt(2,ipi)
            vz0 = vxyt(3,ipi)

            dajc = ONE
            dajct = dajc
            trgshft = qb(ipene,ipi) + cnof + crad
            if (A_KeyFirstPass.eq.1) then
                if (methi.eq.3.or.methi.eq.4) trgshft = crad
                if (methi.eq.0) trgshft = cnof + crad
            endif
            call v_move (qa(35,ipi),r(1,1),9)
            call CalPeni4 (ickcel1,x,y,z,vx0,vy0,vz0,dx,r(1,1),
     &                     pinb,nSeg,ifoun,nPerSeg,nPerFounI,nPerFounDp,
     &                     0,ififo(1,ifoun),coor0(1),
     &                     isegda(1,1),fcell(1),founda(1,ifoun),
     &                     nocon,isega,iseg,epsil(1),fifoa(1),tols,
     &                     trgshft,elem,cnnode,krot,natcor(1),
     &                     dfb(1),ibeam,rad(1),hi)
            call v_move (r(1,1),qa(35,ipi),9)
c
            if (cnnode.lt.0) then
                cnnode = abs(cnnode)
                epsil(3) = -TINY
            endif
            call v_move (r(1,1),qb(35,ipi),9)
            da = epsil(3)

c --- initial contact closure
            if (A_KeyFirstPass.eq.1.and.nocon.eq.0.and.abs(da).lt.pinb) 
     &         then 
               if ((da.gt.ZERO.and.da.le.gapcri).or.
     &            ((da.lt.ZERO.or.keyopt(1).eq.4.or.keyopt(1).eq.5)
     &            .and.(methi.eq.3.or.methi.eq.4.or.methi.eq.0))) then
                  qb(ipene,ipi) = da           
                  da = ZERO
               endif
            endif
            if ((methi.eq.3.or.methi.eq.4).and.A_KeyFirstPass.eq.1) then
               trgshft = qb(ipene,ipi) + cnof 
               call v_move (qa(35,ipi),r(1,1),9)
               call CalPeni4 (ickcel1,x,y,z,vx0,vy0,vz0,dx,r(1,1),
     &                        pinb,nSeg,ifoun,nPerSeg,nPerFounI,
     &                        nPerFounDp,0,ififo(1,ifoun),coor0(1),
     &                        isegda(1,1),fcell(1),
     &                        founda(1,ifoun),nocon,isega,iseg,epsil(1),
     &                        fifoa(1),tols,trgshft,elem,cnnode,krot,
     &                        natcor(1),dfb(1),ibeam,
     &                        rad(1),hi)

               if (cnnode.lt.0) then
                  cnnode = abs(cnnode)
                  epsil(3) = -TINY
               endif
               da = epsil(3)
            endif
c
            cnnode = abs(cnnode)
            vtfra = epsil(1)
            vtfsa = epsil(2)
            if (nocon.eq.2) then
                iseg = 0
            elseif (nocon.eq.0) then
                iseg0 = iseg
            endif

            qa(8,ipi) = iseg
            qa(11,ipi) = dajc
            qa(10,ipi) = nocon + 2
            cinem(3) = 0

            if (nint(qa(5,ipi)).eq.0) qa(5,ipi) = 4.0d0
            if (nint(qa(5,ipi)).eq.4.and.abs(da).lt.pinb) then
               if (A_KeyFirstPass.eq.1.or.da.ge.ZERO) then
                  cinem(3) = 0
                  if (nocon.eq.0.and.ififo(15,ifoun).eq.3.
     &               and.qab(10).eq.0.and.iseg.gt.0) then
                     ndin = 2
                     if (ibeam.eq.abs(isegda(16,iseg))) cinem(3) = 2
                     if (cinem(3).lt.0.5d0)
     &                   cinem(3) = ContAdjStat(nodes(1),ndin,
     &                                          isegda(1,1),iseg,
     &                                          nodbac(1),nPerSeg,nSeg)
                  endif
                  if (cinem(3).eq.0) then
                      qa(5,ipi) = 3.0d0
                      if (nocon.eq.0.and.da.le.ZERO) qa(5,ipi) = 2.d0
                  endif
               endif
            endif

            qa(9,ipi) = da

            if (nocon.eq.0.and.da.lt.ZERO.and.intfer.eq.1) then
                if (A_KeyFirstPass.eq.1) qb(14,ipi) = da
            endif

            if (A_KeyFirstPass.eq.1.and.nocon.eq.0.and.
     &         -da.ge.0.8d0*qab(2).and.methi.le.1.and.cinem(3).eq.0.and.
     &         (C_nPenIni.eq.0.or.C_nPenIni.eq.2).and.
     &          (C_ConnUpd.eq.0.or.kerr.eq.3))
     &         then
               dperr(1) = -da
               dperr(2) = elem
               dperr(3) = isegda(15,iseg) 
               call lockSMP (W_LockConvVar)
               if (C_nPenIni.eq.0) C_nPenIni = 1
               if (C_nPenIni.eq.2) C_nPenIni = 3
               call unlockSMP (W_LockConvVar)
               call erhandler ('ContElement4',5000,
     &                2,'Initial penetration %g may be too
     &                large for contact element %i (with slave
     &                element %i).',dperr(1),' ')
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

         x = xytp(1,ipi)
         y = xytp(2,ipi)
         z = xytp(3,ipi)
         if (krot.eq.1) then
             rxsb = xytp(4,ipi)
             rysb = xytp(5,ipi)
             rzsb = xytp(6,ipi)
         endif

         vx0 = vxyt(1,ipi)
         vy0 = vxyt(2,ipi)
         vz0 = vxyt(3,ipi)

         dajc = ONE
         dajct  = dajc

c --- set ramping contact interfernce
         gapi = ZERO
         if (intfer.eq.1.and.qb(14,ipi).lt.ZERO) then
             gapi = (ONE-A_CurLoadFact)*qb(14,ipi)
         endif
         
         trgshft = gapi + qb(ipene,ipi) + cnof + crad
         cinem(17) = 0
         if (isega.gt.0) cinem(17) = isegda(15,isega)
         call v_move (qa(35,ipi),r(1,1),9)

c         if (elem.eq.65998) then
c             call lockSMP (W_LockAnsys)
c             write (6,*) '------->', ipi
c             call unlockSMP (W_LockAnsys)
c         endif

         call CalPeni4 (ickcel1,x,y,z,vx0,vy0,vz0,dx,r(1,1),
     &                  pinb,nSeg,ifoun,nPerSeg,nPerFounI,nPerFounDp,1,
     &                  ififo(1,ifoun),coor0(1),isegda(1,1),
     &                  fcell(1),founda(1,ifoun),nocon,isega,iseg,
     &                  epsil(1),fifob(1),tols,trgshft,elem,cnnode,krot,
     &                  natcor(1),dfb(1),ibeam,rad(1),
     &                  hi)
 
         if (A_KeyFirstPass.eq.1) call v_move (r(1,1),qa(35,ipi),9)     
         call v_move (r(1,1),qb(35,ipi),9)

c         if (elem.eq.65998) then
c             call lockSMP (W_LockAnsys)
c             write (6,*) 'do', epsil(3),iseg,int(qa(5,ipi))
c             call unlockSMP (W_LockAnsys)
c         endif
         
         if (cnnode.lt.0) then
             cnnode = abs(cnnode)
             epsil(3) = -TINY
         endif
         vtfrb = epsil(1)
         vtfsb = epsil(2)
         db = epsil(3)

         if (abs(nint(qa(5,ipi))).eq.1 .or. abs(nint(qa(5,ipi))).eq.2)
     &       ififo(25,ifoun) = 0
c --- bisection for first time near field to contact
         if (kerr.eq.3.and.db.le.ZERO.and.kpred.ge.2.and.
     &      ififo(25,ifoun).lt.0) then
            if (C_nConF.eq.0) then
               call lockSMP (W_LockConvVar)
               C_nConF = C_nConF + 1
               call unlockSMP (W_LockConvVar)
               if (A_Debug.eq.1) then
                   dperr(1) = elem
                   dperr(2) = ireal
                   call erhandler ('ContElement',5000,kerr,'Contact
     &              element %I (ID %I) status changes suddenly.',
     &              dperr(1),' ')
               endif
            endif
         endif

         if (A_KeyFirstPass.eq.1) qa(ipene,ipi) = db
         if (krot.eq.1) then
             rxfb = epsil(6)
             ryfb = epsil(7)
             rzfb = epsil(8)
         endif

         cinem(1) = ZERO
         cinem(2) = ZERO
         cinem(26) = ZERO
         cinem(27) = ZERO
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

         cinem(30) = ONE
         cinem(31) = ONE
         cinem(32) = ONE
         qb(8,ipi) = iseg
         qb(9,ipi) = db
         qb(11,ipi) = dajc
         qb(10,ipi) = nocon + 2       
         qb(12,ipi) = vtfrb 
         qb(13,ipi) = vtfsb

         ipass = 0

         if (nocon.eq.0.and.iseg.gt.0) then
             if (isegda(17,iseg).gt.0.and.irigf.eq.2) then
c --- symmetric contact pair
                 if (ififo(15,ifoun).eq.3) ipass = 1
             endif
         endif
 65    continue


         cinem(18) = qab(10)
         cinem(3) = 0
         if (nocon.eq.0.and.ififo(15,ifoun).eq.3.and.
     &       qa(5,ipi).eq.3.0d0 .and. -qb(9,ipi) .gt. 0.1d0*pinb) then
             ndin = 2
             if (ibeam.eq.abs(isegda(16,iseg))) cinem(3) = 2
             if (cinem(3).lt.0.5d0)
     &           cinem(3) = ContAdjStat(nodes(1),ndin,isegda(1,1),
     &                                  iseg,nodbac(1),nPerSeg,nSeg)
         endif
         weight = ONE
         if (ipass.eq.1) weight = HALF*weight

c --- compute element volume and energy

         dajc = ONE
         dajcw = weight*dajc

         cinem(10) = mu
         if (lOrthoFric) then
             cinem(34) = mux/mu
             cinem(35) = muy/mu
         endif
         if (0.eq.0) then
             cinem(24) = ONE
         else
             cinem(24) = qab(13)
         endif
         
         qb(33,ipi) = prs(ipi)

c
         if (lOrthoFric) then
c --- transform tangential stress and sliding for
c --- orthotropic friction
c
             lTrans = .true.
             call OrthoToIsoFric(nvmax,cinem(1),qa(1,ipi),qb(1,ipi),
     &                           c(1,1),lTrans,keyopt(1))
         endif


         call FrictionLaw (nvmax,keyopt(1),cinem(1),rvr(3),qa(1,ipi),
     &                     A_TimeInc,qb(1,ipi),c(1,1),intfer,kerr,
     &                     hygro(1,ipi),ififo(1,ifoun),founda(1,ifoun),
     &                     qab(1),ONED,keyopt_1)

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
     &           c(1,1),lTrans,keyopt(1))
         endif
         ield  = abs(nint(qb(5,ipi)))

c         if (elem.eq.65998) then
c             call lockSMP (W_LockAnsys)
c             write (6,*) 'po', qa(9,ipi), iseg, ield
c             write (6,*) '< ------'
c             call unlockSMP (W_LockAnsys)
c         endif

         
         mudat(ipi) = cinem(10)

        
c --- skip if only updating connectivity
         if (C_ConnUpd.eq.1) then
             if (ield.ge.3.and.kfarcn.ne.1) goto 900
#if defined(MOIST2)
             if ( (keyopt(1).eq.3 .or. keyopt(1).ge.5)
     &           .and. A_IterNum.eq.1) then
                 sbeam = abs(isegda(16,iseg))
                 if (sbeam.eq.0) goto 500
                 goto 500
c ---> epsil(7) not ready !!!
                 call GetBeamData(ibeam,natcor(1),bp(1,1))
                 call GetBeamData(sbeam,natcor(2),bp(1,2))
                 if (abs(bp(9,1))+abs(bp(9,2)).le.TINY) goto 500
c     & (L,beamprop(1)), (H,beamprop(2)), (W,beamprop(3)),
c     & (Ex,beamprop(4)), (Ez,beamprop(5)), (Gx,beamprop(6)),
c     & (alpx,beamprop(7)), (alpz,beamprop(8)), 
c     & (dh,beamprop(9))
                 if (CheckHyg(bp(1,1),epsil(7)) .eq. 0) then
                     if (H_key(1) .eq. 1) then                                              1D element
                        call BndHyg1D (bp(1,1),epsil(7),
     &                                 qab(19),em(1,1),MH)
                     else
                        call BndHygStrn (bp(1,1),epsil(7),
     &                                   qab(19),em(1,1),MH,wid(1))
                     endif
                 endif

                 if (tol.gt.ZERO) then
                     dum1 = abs(natcor(1)) - ONE
                     if (abs(dum1).le.tol) then
                         dum2 = HALF*(tol-dum1)/tol
                         call cv1_v1(em(1,1),4,dum2)
                     endif
                 endif
                 ind = (M_Beams(ibeam)-1)*H_Col+1
!$omp atomic update
                 M_HygStrains(ind) = M_HygStrains(ind) + em(1,1)
!$omp atomic update
                 M_HygStrains(ind+1) = M_HygStrains(ind+1) + em(2,1)
                 sbeamn = 0
                 if (tol.gt.ZERO) then
                     dum1 = abs(natcor(2)) - ONE
                     if (abs(dum1).le.tol) then
                         if (natcor(2).gt.ZERO) then
                             isegn = isegda(12,iseg)
                         else
                             isegn = isegda(11,iseg)
                         endif
                         if (isegn.gt.0) sbeamn = abs(isegda(16,isegn))
                         if (sbeamn.gt.0) then
                             dum2 = HALF*(tol-dum1)/tol
                             dum3 = ONE - dum2
                             call cv1_v2(em(1,2),em(1,3),2,dum3)
                             dum3 = em(1,3)
                             em(1,3) = em(2,3)
                             em(2,3) = dum3
                             call cv1_v1(em(1,2),2,dum2)
                        endif
                     endif
                 endif
                 ind = (M_Beams(sbeam)-1)*H_Col+1
!$omp atomic update
                 M_HygStrains(ind) = M_HygStrains(ind) + em(1,2)
!$omp atomic update
                 M_HygStrains(ind+1) = M_HygStrains(ind+1) + em(2,2)
                 
                 if (sbeamn.gt.0) then
                     ind = (M_Beams(sbeamn)-1)*H_Col+1
!$omp atomic update
                     M_HygStrains(ind) = M_HygStrains(ind) + em(1,3)
!$omp atomic update
                     M_HygStrains(ind+1) = M_HygStrains(ind+1) + em(2,3)
                 endif
             endif
#endif defined(MOIST)
             goto 500
         endif
         
         
         qb(33,ipi) = ZERO
         if (kpress.eq.1) then
             if (ield.ge.3.or.qb(1,ipi).le.ZERO) then
                 qb(1,ipi) = qb(1,ipi) + prs(ipi)
                 qb(33,ipi) = prs(ipi)
             endif
         endif

         presb = qb(1,ipi)                    
         taurb = qb(2,ipi)
         tausb = qb(3,ipi)

         if (cinem(6).eq.1) goto 900

         elvol = length
         eEnergy(1) = eEnergy(1) + qb(16,ipi)*dajcw
         eEnergy(4) = eEnergy(4) + qb(17,ipi)*dajcw
         if (iropil.gt.0) then
             dxk = x - xpilo - founda(13,ifoun)
             dyk = y - ypilo - founda(14,ifoun)
             dzk = z - zpilo - founda(15,ifoun)
         endif

         if (ield.ge.3.and.kfarcn.ne.1) goto 900

         if (rad(1).eq.ZERO .or. rad(2).eq.ZERO) then
             call v_move(rvr(1),rad(1),2)
         endif
         if (rad(1).ge.ZERO) then
             radi = rad(1)+abs(rad(2))
         else
             radi = -abs((abs(rad(1))-abs(rad(2))))
         endif
         
c --- compute contact stresses in global axes
         curvature = ONE
         if (abs(dfb(9)).le.TINY) then
             fx = r(1,1)*taurb+r(2,1)*tausb+r(3,1)*presb
             fy = r(1,2)*taurb+r(2,2)*tausb+r(3,2)*presb
             fz = r(1,3)*taurb+r(2,3)*tausb+r(3,3)*presb
         else
             curvature = ONE/(ONE-dfb(9)*radi)
             fx = r(1,1)*taurb*curvature+r(2,1)*tausb+r(3,1)*presb
             fy = r(1,2)*taurb*curvature+r(2,2)*tausb+r(3,2)*presb
             fz = r(1,3)*taurb*curvature+r(2,3)*tausb+r(3,3)*presb
         endif

c --- stiffness matrix in global axes

         if (kelin(1).ne.1) goto 500

         call v_zero (c1(1,1),9)
         if (ield.ge.3) goto 500

c --- product c * r
         call v_zero (cc(1,1),9)
         do 80 i = 1, 3
            do 75 j = 1, 3
               do 70 k = 1, 3
                  if (k.eq.1.and.abs(dfb(9)).gt.TINY) then
                      cc(i,j) = cc(i,j) + c(i,k)*r(k,j)*curvature
                  else
                      cc(i,j) = cc(i,j) + c(i,k) * r(k,j)
                  endif
  70           continue
  75        continue
  80     continue

c                 t
c --- product r  d  r
         do 95 i = 1, 3
            do 90 j = 1, 3
               do 85 k = 1, 3
                  if (k.eq.1.and.abs(dfb(9)).gt.TINY) then
                      c1(i,j) = c1(i,j) + r(k,i)*cc(k,j)*curvature
                  else
                      c1(i,j) = c1(i,j) + r(k,i) * cc(k,j)
                  endif
  85           continue
  90        continue
  95     continue

 500     continue

         if (cinem(6).eq.1) goto 900

         nno = 0
         if (irigf.ne.0 .and. (nocon.eq.0.or.kfarcn.eq.1)) then
            nno = 1
c --- in the case the foundation nodes are reached
            if (nocon.ne.0.and.kfarcn.eq.1) then
                if (ififo(2,ifoun).gt.0) then
                    iseg = ififo(10,ifoun)
                    fifoa(1) = ONE
                    fifob(1) = ONE
                else
                    iseg = 0
                    nno = 0
                    icode = 0
                endif
            endif
            if (iseg.gt.0) then
                icode = isegda(1,iseg)
            else
                nno = 0
                icode = 0
            endif               

            if (icode.eq.1) then
                nno = 3
                if (nint(dx).eq.0) then
                    nno = 2
                elseif (nint(dx).eq.-1) then
                    if (isegda(8,iseg).gt.0) nno = 4 
                elseif (nint(dx).eq.1) then
                    if (isegda(9,iseg).gt.0) nno = 4 
                endif
            elseif (icode.eq.2) then
                nno = 5
                if (nint(dx).eq.0) then
                    nno = 3
                elseif (nint(dx).eq.-1) then
                    if (isegda(8,iseg).le.0) nno = 4 
                elseif (nint(dx).eq.1) then
                    if (isegda(9,iseg).le.0) nno = 4 
                endif
            elseif (icode.eq.7) then
                nno = 4
            endif

            do 570 i = 1,nno
               noeud = isegda(2+i,iseg)

               if (icode.eq.1) then
                   if (i.eq.2) then
                       noeud = isegda(5,iseg)
                       fifoa(2) = fifoa(3) 
                       fifob(2) = fifob(3) 
                       dfb(2) = dfb(3)
                   endif
                   if (i.eq.3) then
                       if (nint(dx).eq.-1) then
                           noeud = isegda(6,iseg)
                           fifoa(3) = fifoa(4) 
                           fifob(3) = fifob(4) 
                           dfb(3) = dfb(4)
                       elseif (nint(dx).eq.1) then
                           noeud = isegda(7,iseg)
                           fifoa(3) = fifoa(5) 
                           fifob(3) = fifob(5) 
                           dfb(3) = dfb(5)
                       else
                           noeud = 0
                           fifob(3) = 0.d0
                           dfb(3) = 0.d0
                       endif
                   endif
                   if (i.eq.4) then
                       if (nint(dx).eq.-1) then
                           noeud = isegda(8,iseg)
                           fifoa(4) = fifoa(6) 
                           fifob(4) = fifob(6) 
                           dfb(4) = dfb(6)
                       elseif (nint(dx).eq.1) then
                           noeud = isegda(9,iseg)
                           fifoa(4) = fifoa(7) 
                           fifob(4) = fifob(7) 
                           dfb(4) = dfb(7)
                       else
                           noeud = 0
                           fifob(4)= 0.d0
                           dfb(4) = 0.d0
                       endif
                   endif
               elseif (icode.eq.2) then
                   if (i.eq.4) then
                       if (nint(dx).eq.-1) then
                           noeud = isegda(6,iseg)
                       elseif (nint(dx).eq.1) then
                           noeud = isegda(7,iseg)
                           fifoa(4) = fifoa(5)
                           fifob(4) = fifob(5) 
                           dfb(4) = dfb(5)
                       else
                           noeud = 0
                           fifob(4) = 0.d0
                           dfb(4) = 0.d0
                       endif
                    endif
                    if (i.eq.5) then
                        if (nint(dx).eq.-1) then
                            noeud = isegda(8,iseg)
                            fifoa(5) = fifoa(6) 
                            fifob(5) = fifob(6) 
                            dfb(5) = dfb(6)
                        elseif (nint(dx).eq.1) then
                            noeud = isegda(9,iseg)
                            fifoa(5) = fifoa(7) 
                            fifob(5) = fifob(7) 
                            dfb(5) = dfb(7)
                        else
                            noeud = 0
                            fifob(5) = 0.d0
                            dfb(5) = 0.d0
                        endif
                    endif
               endif
               
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
     &            'We have found %i nodes of slave element are in
     &            contact with the element %i, which exceeds allowable
     &            number %i. ',dperr(1),' ')
            endif
         endif

         if (C_ConnUpd.eq.1) goto 900

         do 700 j = 1,nnode + nno

            call v_zero (c2(1,1),9)
            if (j.le.nnode) then
                if (j.gt.3) goto 700
                if (nodes(j).eq.0) goto 700
                if (j.ne.ipi) goto 700
            endif
            if (kelin(6).ne.1) goto 575

            if (j.le.nnode) then

c --- internal forces for beam1

               findj = fi(j,ipi)
               dum1 = findj*dajcw
               if (kelin(6).ne.1) goto 571

c --- internal forces for beam1

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
              
            else

               if (nocon.eq.0.or.kfarcn.eq.1) then

c --- internal forces for surface2
c               
                  findj =  fifob(j-nnode)
                  dum2  = findj*dajcw
                  if (kelin(6).ne.1) goto 572
                  if (irigf.eq.2) then
                      lnx =  nnodtv + (kk(j-nnode)-1)*lintv + 1
                  else
                      lnx = nnodtv + 1
                  endif
                  lny = lnx + kdsp             
                  lnz = lny + kdsp
                  Fnr(lnx)= Fnr(lnx) - fx*dum2
                  Fnr(lny)= Fnr(lny) - fy*dum2
                  Fnr(lnz)= Fnr(lnz) - fz*dum2
                  if (abs(dfb(9)).gt.TINY) then
                      Fnr(lnx)= Fnr(lnx) +
     &                          curvature*radi*r(3,1)*dfb(j-nnode)
                      Fnr(lny)= Fnr(lny) +
     &                          curvature*radi*r(3,2)*dfb(j-nnode)
                      Fnr(lnz)= Fnr(lnz) +
     &                          curvature*radi*r(3,3)*dfb(j-nnode)
                  endif
                  if (iropil.gt.0) then
                      nod = nofo(kk(j-nnode))
                      kddlx = nod
                      kddly = kddlx + kdsp
                      kddlz = kddly + kdsp
                      dispx = dispn(kddlx)
                      dispy = dispn(kddly)
                      dispz = dispn(kddlz)
                      dxj = dispx + coor0(kddlx) - xpilo
                      dyj = dispy + coor0(kddly) - ypilo
                      dzj = dispz + coor0(kddlz) - zpilo
                      if (nno.eq.1) then
                          dxj = dxk
                          dyj = dyk
                          dzj = dzk
                      endif

                      lna = lnx + lintv

                      Fnr(lna) = Fnr(lna) - (fz*dyj-fy*dzj)*dum2
                      Fnr(lna+1) = Fnr(lna+1) - (fx*dzj - fz*dxj)*dum2
                      Fnr(lna+2) = Fnr(lna+2) - (fy*dxj - fx*dyj)*dum2
                   endif

                   if (krot.eq.1.and.ield.eq.1.and.keyopt_1.gt.0) then
                       rigpr = -rigpr0*dum2
                       Fnr(lnx+3) = Fnr(lnx+3) - rigpr*(rxfb-rxsb)
                       Fnr(lnx+4) = Fnr(lnx+4) - rigpr*(ryfb-rysb)
                       Fnr(lnx+5) = Fnr(lnx+5) - rigpr*(rzfb-rzsb)
                   endif

 572              continue
               endif

            endif

 575        continue

            if (kelin(1).ne.1) goto 700

c --- stiffness matrix

            do 600 i = 1,nnode+nno

               if (i.le.nnode) then
                   if (i.gt.3) goto 600
                   if (nodes(i).eq.0) goto 600
                   if (i.ne.ipi) goto 600
               endif
               if (j.le.nnode.and.i.le.nnode) then

                  findi = fi(i,ipi)
                  findj = fi(j,ipi)
                  dum1 = findi*findj*dajcw
                  dum2 = findi*weight

c --- mechanical...

                  lm1 = (j-1)*nrtv + lintv*(i-1)
                  if (ield.ge.3) goto 580
                  do 675 ik=1,3
                     do 673 jk=1,3
                        c(ik,jk) = c1(ik,jk)*dum1+c2(ik,jk)*dum2
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

                  findi = fi(i,ipi)
                  findj = fifob(j-nnode)
                  dum1 = findi*findj*dajcw

                  if (irigf.eq.2) then
                      lm1 = nnodtv*nr + (kk(j-nnode)-1)*nrtv 
     &                    + lintv*(i-1)
                  else   
                      lm1 = nnodtv*nr + lintv*(i-1)
                  endif

                  if (ield .ge. 3) goto 583

                  do 685 ik=1,3
                     do 683 jk=1,3
                        c(ik,jk) = c1(ik,jk)*dum1
  683                continue
  685             continue
                  indic = 0
  687             continue
                  do 689 jk = 1,3
                     lma = lm1 + nr*(jk-1)
                     do 688 ik = 1,3
                        lma = lma + kdsp
                        Kmtrx(lma) = - c(ik,jk) + Kmtrx(lma)
  688                continue
  689             continue

                  if (indic.ne.0) goto 583

                  if (iropil.gt.0) then
                      nod = nofo(kk(j-nnode))
                      kddlx = nod
                      kddly = kddlx + kdsp
                      kddlz = kddly + kdsp
                      dispx = dispn(kddlx)
                      dispy = dispn(kddly)
                      dispz = dispn(kddlz)
                      dxj = dispx + coor0(kddlx) - xpilo
                      dyj = dispy + coor0(kddly) - ypilo
                      dzj = dispz + coor0(kddlz) - zpilo
                      if (nno.eq.1) then
                          dxj = dxk
                          dyj = dyk
                          dzj = dzk
                      endif
                      lm1 = (nnode+1)*nrtv + lintv*(i-1)
                      c(1,1) = -(c1(1,2)*dzj-c1(1,3)*dyj)*dum1
                      c(2,1) = -(c1(2,2)*dzj-c1(2,3)*dyj)*dum1
                      c(3,1) = -(c1(3,2)*dzj-c1(3,3)*dyj)*dum1
                      c(1,2) = -(-c1(1,1)*dzj+c1(1,3)*dxj)*dum1
                      c(2,2) = -(-c1(2,1)*dzj+c1(2,3)*dxj)*dum1
                      c(3,2) = -(-c1(3,1)*dzj+c1(3,3)*dxj)*dum1
                      c(1,3) = -(c1(1,1)*dyj-c1(1,2)*dxj)*dum1
                      c(2,3) = -(c1(2,1)*dyj-c1(2,2)*dxj)*dum1
                      c(3,3) = -(c1(3,1)*dyj-c1(3,2)*dxj)*dum1
                      indic = 1
                      goto 687
                  endif

c --- rotation stiffness
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
                  findj = fi(j,ipi)
                  dum1 = findi*findj*dajcw
                  dum2 = findi*weight
                  if (irigf.eq.2) then
                      lm1 = (j-1)*nrtv+nnodtv +
     &                      (kk(i-nnode)-1)*lintv
                  else
                      lm1 = (j-1)*nrtv + nnodtv
                  endif
                  if (ield .ge. 3) goto 584

                  do 693 ik=1,3
                     do 691 jk=1,3
                        c(ik,jk) = c1(ik,jk)*dum1+c2(ik,jk)*dum2
  691                continue
  693             continue


                  do 695 jk = 1,3
                     lma = lm1 + nr*(jk-1)
                     do 694 ik = 1,3
                        lma = lma + kdsp
                        Kmtrx(lma) = - c(ik,jk) + Kmtrx(lma)
  694                continue
  695             continue

                  if (iropil.gt.0) then
                      nod = nofo(kk(i-nnode))
                      kddlx = nod
                      kddly = kddlx + kdsp
                      kddlz = kddly + kdsp
                      dispx = dispn(kddlx)
                      dispy = dispn(kddly)
                      dispz = dispn(kddlz)
                      dxi = dispx + coor0(kddlx) - xpilo
                      dyi = dispy + coor0(kddly) - ypilo
                      dzi = dispz + coor0(kddlz) - zpilo
                      if (nno.eq.1) then
                          dxi = dxk
                          dyi = dyk
                          dzi = dzk
                      endif
                      lmr = (nnode+1)*lintv + nrtv*(j-1)
                      cr(1,1) = -c(2,1)*dzi+c(3,1)*dyi
                      cr(1,2) = -c(2,2)*dzi+c(3,2)*dyi
                      cr(1,3) = -c(2,3)*dzi+c(3,3)*dyi
                      cr(2,1) =  c(1,1)*dzi-c(3,1)*dxi 
                      cr(2,2) =  c(1,2)*dzi-c(3,2)*dxi
                      cr(2,3) =  c(1,3)*dzi-c(3,3)*dxi
                      cr(3,1) = -c(1,1)*dyi+c(2,1)*dxi
                      cr(3,2) = -c(1,2)*dyi+c(2,2)*dxi
                      cr(3,3) = -c(1,3)*dyi+c(2,3)*dxi
                      do 2695 jk = 1,3
                         lma = lmr + nr*(jk-1)
                         do 1694 ik = 1,3
                            lma = lma + kdsp
                            Kmtrx(lma) = - cr(ik,jk) + Kmtrx(lma)
 1694                    continue
2695                 continue
                  endif                      

c --- rotation stiffness
                  if (krot.eq.1.and.ield.eq.1.and.keyopt_1.gt.0) then
                      rigpr = -rigpr0*dum1
                      do jk = 3,5
                         lmr = lm1 + 1 + jk*(1+nr)
                         Kmtrx(lmr) = Kmtrx(lmr) + rigpr
                      enddo
                  endif
c                 
 584              continue
               endif

               if (j.gt.nnode.and.i.gt.nnode) then     

                  findi = fifob(i-nnode)
                  findj = fifob(j-nnode)
                  dum1 = findi*findj*dajcw
                  if (irigf.eq.2) then
                      lm1 = nnodtv*nr + (kk(j-nnode)-1)*nrtv
     &                    + nnodtv    + (kk(i-nnode)-1)*lintv
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

c --- rotation stiffness

                  if (iropil.gt.0) then
                      nod = nofo(kk(i-nnode))
                      kddlx = nod
                      kddly = kddlx + kdsp
                      kddlz = kddly + kdsp
                      dispx = dispn(kddlx)
                      dispy = dispn(kddly)
                      dispz = dispn(kddlz)
                      dxi = dispx + coor0(kddlx) - xpilo
                      dyi = dispy + coor0(kddly) - ypilo
                      dzi = dispz + coor0(kddlz) - zpilo
 
                      nod = nofo(kk(j-nnode))
                      kddlx = nod
                      kddly = kddlx + kdsp
                      kddlz = kddly + kdsp
                      dispx = dispn(kddlx)
                      dispy = dispn(kddly)
                      dispz = dispn(kddlz)
                      dxj = dispx + coor0(kddlx) - xpilo
                      dyj = dispy + coor0(kddly) - ypilo
                      dzj = dispz + coor0(kddlz) - zpilo
                      if (nno.eq.1) then
                          dxi = dxk
                          dyi = dyk
                          dzi = dzk
                          dxj = dxk
                          dyj = dyk
                          dzj = dzk
                      endif
                      lmh = nnodtv*nr + lintv*(nnode+1)
                      lmf = (nnode+1)*nrtv + nnodtv
                      lmo = (nnode+1)*nrtv + lintv*(nnode+1)

                      indic = 1
                      lm1 = lmh
                      cr(1,1) = -c(2,1)*dzi+c(3,1)*dyi
                      cr(1,2) = -c(2,2)*dzi+c(3,2)*dyi
                      cr(1,3) = -c(2,3)*dzi+c(3,3)*dyi
                      cr(2,1) =  c(1,1)*dzi-c(3,1)*dxi 
                      cr(2,2) =  c(1,2)*dzi-c(3,2)*dxi
                      cr(2,3) =  c(1,3)*dzi-c(3,3)*dxi
                      cr(3,1) = -c(1,1)*dyi+c(2,1)*dxi
                      cr(3,2) = -c(1,2)*dyi+c(2,2)*dxi
                      cr(3,3) = -c(1,3)*dyi+c(2,3)*dxi
                      goto 1699
 1600                 indic = 2
                      lm1 = lmf
                      cr(1,1) = -(c(1,2)*dzj-c(1,3)*dyj)
                      cr(2,1) = -(c(2,2)*dzj-c(2,3)*dyj)
                      cr(3,1) = -(c(3,2)*dzj-c(3,3)*dyj)
                      cr(1,2) = -(-c(1,1)*dzj+c(1,3)*dxj)
                      cr(2,2) = -(-c(2,1)*dzj+c(2,3)*dxj)
                      cr(3,2) = -(-c(3,1)*dzj+c(3,3)*dxj)
                      cr(1,3) = -(c(1,1)*dyj-c(1,2)*dxj)
                      cr(2,3) = -(c(2,1)*dyj-c(2,2)*dxj)
                      cr(3,3) = -(c(3,1)*dyj-c(3,2)*dxj)

 1699                 do 2700 jk = 1,3
                         lma = lm1 + nr*(jk-1)
                         do 1700 ik = 1,3
                            lma = lma + kdsp
                            Kmtrx(lma) = cr(ik,jk) + Kmtrx(lma)
 1700                    continue
 2700                 continue
                      if (indic.eq.1) goto 1600

                      ct(1,1) = -cr(2,1)*dzi+cr(3,1)*dyi
                      ct(1,2) = -cr(2,2)*dzi+cr(3,2)*dyi
                      ct(1,3) = -cr(2,3)*dzi+cr(3,3)*dyi
                      ct(2,1) =  cr(1,1)*dzi-cr(3,1)*dxi 
                      ct(2,2) =  cr(1,2)*dzi-cr(3,2)*dxi
                      ct(2,3) =  cr(1,3)*dzi-cr(3,3)*dxi
                      ct(3,1) = -cr(1,1)*dyi+cr(2,1)*dxi
                      ct(3,2) = -cr(1,2)*dyi+cr(2,2)*dxi
                      ct(3,3) = -cr(1,3)*dyi+cr(2,3)*dxi

                      do 4700 jk = 1,3
                         lma = lmo + nr*(jk-1)
                         do 3700 ik = 1,3
                            lma = lma + kdsp
                            Kmtrx(lma) = ct(ik,jk) + Kmtrx(lma)
 3700                    continue
 4700                 continue
                  endif                      

c --- rotation stiffness
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

c --- do not set contact states
         if (qb(5,ipi).eq.-2.d0) then
             qb(5,ipi) = 3.d0
         endif

c --- reaction applied to the foundation nodes 
         if (irigf.ne.0 .or. iseg.le.0) goto 900
         
         icode = isegda(1,iseg)
         nno = 1
         if (icode.eq.7) nno = 4
         do 800 i=1,nno

            inod  = isegda(2+i,iseg)
            if (inod.le.0) goto 800
            if (inopil.eq.inod) goto 800

            df(1) = fx * fifob(i)*dajcw   
            df(2) = fy * fifob(i)*dajcw
            df(3) = fz * fifob(i)*dajcw

            call lockSMP (W_LockLoadVect)
            do iddl = 1,lintv            
               if (inopil.le.0) then
                   ieqnr = inod + iddl - 1
               else
                   ieqnr = inopil + iddl - 1
               endif
               if (ieqnr.le.0) cycle
c --- force applied on the foundation
               M_ForcNR(ieqnr) = M_ForcNR(ieqnr) - df(iddl)
               M_ForcTot(ieqnr) = M_ForcTot(ieqnr) + df(iddl)
            enddo
            call unlockSMP (W_LockLoadVect)

  800    continue
c
  900 continue

      kelout(1) = 1
      kelout(6) = 1
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
             ContStat = 1
             call lockSMP(W_LockAnsys)
             call elmpnd (elem,nmax,nodes(1))
             call unlockSMP(W_LockAnsys)
         endif
      endif


      if (kelin(1).eq.1) then
         do 994 i = 1,nr-1
            do 993 j = i+1,nr 
               Kmtrx(nr*(j-1)+i) = HALF*(Kmtrx(nr*(j-1)+i)
     &                           + Kmtrx(nr*(i-1)+j))
               Kmtrx(nr*(i-1)+j) = Kmtrx(nr*(j-1)+i)
  993       continue
  994    continue
      endif

 995  qab(7) = elvol

c --- charttering measurement (only when primary computations and not line search)
      if (C_ConnUpd.ne.1) then
         ifDebo = .False.
         if (founda(33,ifoun).lt.qab(3)) founda(33,ifoun) = qab(3)
         do ipi =1,npi
            ififo(22,ifoun) = ififo(22,ifoun) + nint(qb(44,ipi))
            if (qb(nvmax-9,ipi).ge.ONE) ifDebo = .True.
            if (abs(nint(qb(5,ipi))).le.2) then
               if (founda(32,ifoun).lt.qb(nvmax-4,ipi))
     &            founda(32,ifoun) = qb(nvmax-4,ipi)
               if (qb(nvmax-4,ipi).gt.TINY.and.
     &            founda(36,ifoun).gt.qb(nvmax-4,ipi))
     &            founda(36,ifoun) = qb(nvmax-4,ipi)
               if (founda(37,ifoun).lt.qb(nvmax-3,ipi))
     &            founda(37,ifoun) = qb(nvmax-3,ipi)
               if (qb(nvmax-3,ipi).gt.TINY.and.
     &             founda(28,ifoun).gt.qb(nvmax-3,ipi))
     &            founda(28,ifoun) = qb(nvmax-3,ipi)
               if (abs(founda(34,ifoun)).lt.abs(qb(1,ipi)))
     &            founda(34,ifoun) = qb(1,ipi)
               slides =sqrt(qb(2,ipi)*qb(2,ipi)+qb(3,ipi)*qb(3,ipi))
               if (abs(founda(35,ifoun)).lt.abs(slides)) then
                  founda(35,ifoun) = slides
               endif
               slides =sqrt(qb(29,ipi)*qb(29,ipi)+qb(30,ipi)*qb(30,ipi))
               if (founda(31,ifoun).lt.slides) then
                   founda(31,ifoun) = slides
               endif
               if (abs(nint(qb(5,ipi))).eq.1) then
                  if (ififo(23,ifoun).eq.0) founda(30,ifoun) = ZERO
                  slides =
     &                 sqrt(qb(12,ipi)*qb(12,ipi)+qb(13,ipi)*qb(13,ipi))
                  if (founda(30,ifoun).lt.slides) then
                      founda(30,ifoun) = slides
                      ififo(23,ifoun) = elem
                  endif
               endif
            endif
         enddo
         if (ifDebo) then
             call lockSMP (W_LockConvVar)
             C_nDeBo = C_nDeBo + 1
             call unlockSMP (W_LockConvVar)
         endif    
      endif

c --- save final contact status for removal of loose elements
      ield = 4
      icnda(3,icon) = 0
      do ipi = 1,npi
         ield = min(ield,abs(nint(qb(5,ipi))))
         icnda(3,icon) = max(icnda(3,icon),nint(qb(8,ipi)))
      enddo
      
      if (ield .le. 2) then
c --- mark contact
          if (ififo(26,ifoun).lt.0) then
              call lockSMP (W_LockConvVar)
              ififo(26,ifoun) = 0
              call unlockSMP (W_LockConvVar)
          endif
c          isegda(19,iseg) = icon
      else
          if (ield .eq. 3) then
c --- mark near contact only if the pair is not in contact
              if (ififo(26,ifoun).lt.-1) then
                  call lockSMP (W_LockConvVar)
                  ififo(26,ifoun) = max(ififo(26,ifoun),-1)
                  call unlockSMP (W_LockConvVar)
              endif
          endif
c          iseg = icnda(3,icon)
c --- update the pair connection only if solution converged
          if (A_KeyConv.eq.1 .and. A_KeySvrUpd.eq.1) then
              icnda(3,icon) = 0
c --- zero out segment info, if it has the current contact element in
c              if (iseg.gt.0) then
c                  if (isegda(19,iseg).eq.icon) isegda(19,iseg) = 0
c              endif
          endif
      endif
      
c --- perform time predictions
      if(A_KeyAutoStep.eq.1 .and. A_KeyMatUpd.eq.1 .and. kpred.ge.2)then
         call ContPredi (elem,kpred-1,nvmax,npi,qa(1,1),qb(1,1),mu,
     &                   ipair,cinem(13))
      endif

c --- update history variables
      if (A_KeyMatUpd .eq. 1) then
          qab(10) = qab(8)
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
          if (nnodmst.eq.3) moistend(4) = moistend(3)
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
  999 continue

c --- produce output
      if (A_KeyOut.eq.1) then
          if (btest(eomask,W_OutElIndx)) then
              if (btest(eomask,W_OutMisc)) then
c --- output contact status through misc data
                  call v_zero (vect(1),elChar(NMSMIS))
                  vect(1) = dble(ield)
                  if (nodes(3).eq.0) then
                      n = 2
                  else
                      n = 3
                  endif
                  do ipi = 1,n
c --- these is damage output
                     vect(2) = max(vect(2),qb(nvmax-9,ipi))
c --- contact area
                     if (qb(5,ipi).le.2) then
                         vect(7) = vect(7) + ONE
                     endif
                  enddo
                  if (vect(7).gt.ZERO) then
                      do j = 1,n
                         lmx = 1 + (j-1)*lintv
                         lmy = lmx + kdsp
                         lmz = lmy + kdsp
                         vect(3) = vect(3) + Fnr(lmx)
                         vect(4) = vect(4) + Fnr(lmy)
                         vect(5) = vect(5) + Fnr(lmz)
                      enddo
c --- save energy
                      vect(6) = eEnergy(1)
c --- save contact length
                      vect(7) = vect(7)/dble(n)*length
                  endif
                  i = elChar(NMSMIS)
                  call WriteToRST (EDEMS,lcerstL,edindxL(1),i,vect(1))
              endif
          endif
      endif
      
#if defined(DEBUG)
      call SubEnd('ContElement4')
#endif
c      if (A_SubStepNum.gt.1) stop
      return
      end
