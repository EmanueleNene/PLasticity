      subroutine EBeam3D (elem,ielc,elmdat,eomask,nodes,locsvrL,kelreq,
     &                    nr,xyz,u,kelout,zs,zass,gstif,
     &                    zsc,zscnr,elvol,elmass,center,elener,
     &                    edindxL,lcerstL)
c --- primary function:
c        1. compute element matrices, load vectors, and results
c --- secondary functions:
c        2. maintain element solution data
c
c  input arguments:
c       elem   (int,sc,in)          - element label (number)
c       ielc   (int,ar(IELCSZ),in)  - array of element type characteristics
c       elmdat (int,ar(EL_DIM),in)  - array of element data
c       eomask (int,sc,in)          - bit pattern for element output
c       nodes  (int,ar(nnod),in)    - array of element node numbers
c       locsvrL (LONGINT,sc,in)     - location of this eles svrs on file esav
c       kelreq (int,ar(10),in)      - matrix and load vector form requests
c       nr     (int,sc,in)          - matrix and lv size
c       xyz    (dp,ar(6,nnod),in)   - nodal coords (orig) and rotation angles
c       u      (dp,ar(nr,5),in)     - element nodal solution values
c
c  output arguments:
c       kelout (int,ar(10),out)     - keys indicating created mats and lvs
c       zs     (dp,ar(nr,nr),inout) - k matrix (kelreq(1))
c       zass   (dp,ar(nr,nr),inout) - m matrix (kelreq(2))
c       gstif  (dp,ar(nr,nr),inout) - s matrix (kelreq(4))
c       zsc    (dp,ar(nr),out)      - applied f vector (kelreq(5))
c       zscnr  (dp,ar(nr),out)      - n-r restoring f vector (kelreq(6))
c       elvol  (dp,sc,out)          - element volume
c       elmass (dp,sc,out)          - element mass
c       center (dp,ar(3),out)       - centroid location
c       elener (dp,ar(5),out)       - element energies
c       edindxL(LONG,ar(25),out)    - element result data file indexes
c       lcerstL(LONG,sc,inout)      - position on result file
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
c
      double precision
     & v_dot,CalDist2P

      integer*8 locsvrL, lcerstL, edindxL(25)

      integer nnod
      parameter (nnod=2)
      
      integer
     & elem,ielc(IELCSZ),elmdat(EL_DIM),eomask,nodes(nnod),kelreq(10),
     & nr,kelout(10),nlgeom,nrkey

      integer
     & mat,ireal,nrvr,kelin(10),keyopt(12),i,num,kcompt,sstif_c,add,
     & SectionType,kstabnl,nnrsvr,svindx(20),istar,jstar,n,j,kerr

      double precision
     & xyz(6,3),u(nr,5),zs(nr,nr),zass(nr,nr),gstif(nr,nr),
     & zsc(nr),zscnr(nr),elvol,elmass,center(3),elener(5)

      double precision
     & rvr(7),rOuter,rInner,StabFact,StbDpVar(3),forx(6,2),
     & area,aiz,aiy,thick,thickz,thicky,epinit,aip,shearz,sheary

      double precision
     & ex,nuxy,gxy,dens

      double precision
     & prop(13),con,vect(33),alenv,uth(12),utho(12),
     & trv(3,3),trr(3,3),tr(3,3),udef(12),force,ulv(12),forl(12),
     & forg(12),ustr(12),uthoo(12),stress(1),
     & xyzh(3,2),qlc(3),gam(12,3),alen,keba(12,12),feba(12),fe(12),
     & hthb(3,6),pte(12,12),vk1(12),wk1(12,12),deba(12),fba(12,3),ioa,
     & eleRot(3),nrsvr(4),StbMat(12,12),udef0(12),cDamp

      equivalence
     & (ex,prop(1)), (nuxy,prop(7)),
     & (gxy,prop(4)),(dens,prop(13))

c --- controls formation of consistent stress-stiffness
      equivalence
     & (kstabnl,keyopt(3))
      
c --- define initial data
      mat = elmdat(EL_MAT)
      ireal = elmdat(EL_REAL)
      nrvr = ielc(NMDRLC)
      nnrsvr = ielc(NMSSVR)
      num = nr*nr
      elener(1) = 0.0d0
      kerr = 0
      nlgeom = 1
      nrkey = 1
      call vi_move (ielc(KYOP1),keyopt(1),12)
      area = 0.0d0

      if (nrkey.eq.1) then
          sstif_c = keyopt(2)
      else
          sstif_c = 0
      endif

c --- consistent tangent is undesirable in transient
      if (A_Trans.gt.0) sstif_c = 0

c --- get the element real constant data
      call rvrget (elem,ireal,ielc(1),nrvr,rvr(1))
     
      if (kstabnl.gt.0 .and. nrkey.eq.1) then
          StabFact = rvr(6)
c --- get the SVR index vector
          call GetSVRIndx (elem,locsvrL,svindx(1))
c --- retrieve the nonlinear SVRs
          call GetSVR (svindx(1),1,nnrsvr,nrsvr(1))
      endif

c --- section type 
      SectionType = nint(abs(rvr(1)))
c --- initial strain
      epinit = rvr(5)

      if (SectionType.eq.1) then
c --- rectangular
          thickz = abs(rvr(2))
          thicky = abs(rvr(3))
          area = thickz*thicky
          aiz = thicky*thickz**3/12.0d0
          aiy = thickz*thicky**3/12.0d0
          shearz = 1.2d0 ! 6/5
          sheary = shearz
      elseif (SectionType.eq.2)then
c --- thin wall hollow rectangular
          thickz = abs(rvr(2))
          thicky = abs(rvr(3))
          thick = abs(rvr(4))
          if (thick.eq.0.0d0) then
              thick = (0.5d0-1.0d-3)*min(thickz,thicky)
          endif
          area = 2.0d0*thick*(thickz + thicky - 2.0d0*thick)
          aiz = thick*thickz**3/6.d0 + 0.5d0*thick*thicky*thickz**2
          aiy = thick*thicky**3/6.d0 + 0.5d0*thick*thickz*thicky**2
          shearz = 2.4d0 ! 12/5
          sheary = shearz
      elseif (SectionType.eq.3) then
c --- circular
          rOuter = abs(rvr(2))
          area = W_Pi*rOuter**2
          aiz = W_Pi*rOuter**4*0.25d0
          aiy = aiz
          shearz = 10.0d0/9.0d0
          sheary = shearz
      elseif (SectionType.eq.4) then
c --- thin-wall pipe
          rOuter = abs(rvr(2))
          rInner = abs(rvr(3))
          if (rInner.eq.0.0d0) rInner = rOuter*1.0d-3
          area = W_Pi*(rOuter**2-rInner**2)
          aiz = W_Pi*(rOuter**4-rInner**4)*0.25d0
          aiy = aiz
          shearz = 2.0d0
          sheary = shearz
      endif 

      if (area.le.TINY) then
          call erhandler ('EBeam3D',5000,3,
     &                    'Beam element %I has zero area!',
     &                     dble(elem), ' ')
          goto 999
      endif
c --- polar moment of inertia
      aip = aiz + aiy      
      
c --- compute materials for every pass
      call GetElasticProps (elem,mat,0.0d0,prop(1))

c --- geometry computations

c --- geometry computations
c --- compute original length
      alenv = CalDist2P (xyz(1,1),xyz(1,2))

c --- compute original transformation
      call LineGeomTrans (xyz(1,1),6,trv(1,1))
      call v_move (trv(1,1),tr(1,1),9)
c --- compute updated transformation for large rotations
      if (nlgeom.eq.1) then
          if (A_KeyFirstPass.eq.0) then
c --- compute transformation matrix from average rotation
              call CalAveRot (u(4,A_RowDisp),u(10,A_RowDisp),eleRot(1))
              call CalRotTr (eleRot(1),trr(1,1))
              call AxB_C (trv(1,1),trr(1,1),tr(1,1), 3,3,3, 3,3,3)
          endif
      endif

      if (nlgeom.eq.1.and.sstif_c.eq.1.and.A_KeyFirstPass.ne.1) then
          call UpdBeamXYZ (2,u(1,A_RowDisp),nr,xyz(1,1),6,xyzh(1,1))
          call CalTransMat (xyzh(1,1),qlc(1),tr(1,1),alen)
          call AtxB_C (trv(1,1),tr(1,1),trr(1,1), 3,3,3, 3,3,3)
          call CalGamma (gam(1,1),qlc(1),tr(1,1),alen)
      endif

c --- set up deformational and non-linear displacements
      if (nlgeom.eq.0 .or. A_KeyFirstPass.eq.1) then
          call v_move (u(1,A_RowDisp),udef(1),12)
          call AxB_C (tr(1,1),udef(1),deba(1), 3,3,3, 3,4,3)
      else
          call CalUdef (xyz(1,1),6,2,trr(1,1),u(1,A_RowDisp),
     &                  udef(1),kerr)
c          if (i.eq.1) write (8,*) 'error', elem
          call AxB_C (tr(1,1),udef(1),deba(1), 3,3,3, 3,4,3)
          if (nlgeom.eq.1.and.sstif_c.eq.1) then
              call CalHtb (hthb(1,1),deba(1),2)
              call CalProjMat (u(1,A_RowDisp),xyz(1,1),tr(1,1),gam(1,1),
     &                         pte(1,1),2,6,12)
          endif
      endif

c --- non-linear displacement vector relative to node i in local c.s.
c --- assumes a simply supported beam
      call v_zero (vect(1),12)
      vect(7) = -epinit*alenv

      call AtxB_C (tr(1,1),vect(1),utho(1), 3,3,3, 3,4,3)

c --- compute element mass and centroid
      elvol = alenv*(1.0d0 - epinit)*area
      elmass = dens*elvol
      center(1) = (xyz(1,1) + xyz(1,2))*0.5d0
      center(2) = (xyz(2,1) + xyz(2,2))*0.5d0
      center(3) = (xyz(3,1) + xyz(3,2))*0.5d0

      call vi_move (kelreq(1),kelin(1),10)
      kelin(1) = max(kelin(1),kelin(5),kelin(6))
      if (A_KeyOut.eq.1) kelin(1) = 1

c --- compute forces
      kcompt = 0
      if (A_KeyOut.eq.1 .or. nrkey.eq.1) then
          if (nlgeom.eq.1.and.sstif_c.eq.1) then
              call v_zero (zs(1,1),144)
              call BeamEStif (ex,gxy,area,aiz,aiy,aip,
     &                        shearz,sheary,alenv,nr,zs(1,1))
              call v_move (zs(1,1),keba(1,1),144)
              call MatRotTran (zs(1,1),tr(1,1),nr,3)
              kcompt = 1
          endif
          call v_zero (forg(1),12)
c --- set up stress-causing displacements
          if (nrkey.eq.0) then
              call v_move (udef(1),ustr(1),12)
          else
c --- remove previous non-linear effects
              call v_zero (vect(1),12)
              vect(7) = -epinit*alenv
              call AtxB_C (tr(1,1),vect(1),uthoo(1), 3,3,3, 3,4,3)
              call v1mv2_v3 (udef(1),uthoo(1),ustr(1),12)
          endif
c --- strain energy calculation
          call v_move (udef(1),udef0(1),12)
          call v1mv2_v1 (udef(1),utho(1),12)
      endif

c --- stiffness matrix
      if (kelin(1).eq.1) then
          if (kcompt.eq.0) then
              call v_zero (zs(1,1),144)
              call BeamEStif (ex,gxy,area,aiz,aiy,aip,
     &                        shearz,sheary,alenv,nr,zs(1,1))
              call v_move (zs(1,1),keba(1,1),num)
          else
c ---- local stiffness matrix has been formed
              call v_move (keba(1,1),zs(1,1),num)
          endif

c --- compute f-e overbar and no overbar which are used for convergence
c     checkings and calculations of stress stiffness
          if (A_KeyFirstPass.eq.1) call v_zero (fe(1),12)
c --- ref: nour-omid&rankin 'finite rotation analysis ...'
c          CMANE, vol,93, pp. 353-384;1991; eqns. (76-78)
          if (nlgeom.eq.1.and.sstif_c.eq.1.and.A_KeyFirstPass.ne.1) then
              call AxB_C (tr(1,1),udef(1),vk1(1),3,3,3,3,4,3)
              call Axv1_v2 (keba(1,1),vk1(1),feba(1),12,12)
              call v_move (feba(1),fe(1),12)
              call BlockDiagMult (hthb(1,1),feba(1),vk1(1),2)
              call Axv1_v2 (pte(1,1),vk1(1),fe(1),12,12)

c --- material stiffness matrix for finite rotations
              call MatHTran (zs(1,1),hthb(1,1),nr,3)
              call CalCmt (wk1(1,1),feba(1),deba(1),2,12)
              call v1pv2_v1i (zs(1,1),nr,wk1(1,1),nr,nr)
              call AxB_C (pte(1,1),zs(1,1),wk1(1,1),12,12,12,12,12,12)
              call AxBt_C (wk1(1,1),pte(1,1),zs(1,1),12,12,12,12,12,12)

c --- stress stiffness matrix for finite rotations
              call CalFb (fba(1,1),fe(1),2,12)
              call AxBt_C(gam(1,1),fba(1,1),wk1(1,1),12,12,12,12,12,3)
              call v_move (wk1(1,1),gstif(1,1),144)
              call AxBt_C (pte(1,1),wk1(1,1),gstif(1,1),12,12,12,
     &                     12,12,-12)
              call v1pv2_v1 (zs(1,1),gstif(1,1),num)
          endif

c --- calculate k x u in local c.s.
          if (kelin(4).eq.1 .or. A_KeyOut.eq.1)
     &        call Axv1_v2 (zs(1,1),deba(1),forl(1),12,12)

          call MatRotTran (zs(1,1),tr(1,1),nr,3)

c --- calculate axial force
          call Axv1_v2 (zs(1,1),udef(1),forg(1), 12,12)

          kelout(1) = 1
      endif

c --- mass matrix
      if (kelin(2).eq.1) then
          call v_zero (zass(1,1),144)
          con = dens
          if (con.ne.0.0d0) then
              call BeamEMass (con,ex,gxy,area,aiz,aiy,shearz,sheary,
     &                        alenv,epinit,0,nr,zass(1,1),A_MassLump)
              if (A_MassLump .eq. 0) 
     &            call MatRotTran (zass(1,1),tr(1,1),nr,3)
          endif
          kelout(2) = 1
      endif

c --- stress stiffness matrix
      if (kelin(4).eq.1.and.sstif_c.eq.0) then
          call v_zero (gstif(1,1),144)
          if (A_KeyFirstPass.eq.1) then
c --- first step force
              force = area*ex*epinit
          else
              force = 0.5d0*(-forl(1) + forl(7))
          endif

c --- stress stiffness matrix
          if (force.ne.0.0d0) then
              ioa = 0.0d0 ! (aiz + aiy)/area
              call BeamSStif (force,alenv,nr,ioa,gstif(1,1))
              call MatRotTran (gstif(1,1),tr(1,1),nr,3)
          endif
          kelout(4) = 1
      endif

c --- applied load vector
      if (kelin(5).eq.1) then
          call v_zero (zsc(1),12)
          if (epinit .ne. 0.0d0) then
c --- set up initial strain displacement in local c.s.
              call v_zero (vect(1),12)
              vect(7) = -epinit*alenv
              call AtxB_C (tr(1,1),vect(1),uth(1), 3,3,3, 3,4,3)
              if (nrkey.eq.0) then
                  call v_move (uth(1),ulv(1),12)
              else
                  call v1mv2_v3 (uth(1),utho(1),ulv(1),12)
              endif
              if (nlgeom.eq.1.and.sstif_c.eq.1) then
                  call AxB_C (tr(1,1),ulv(1),vect(1), 3,3,3, 3,4,3)
                  call Axv1_v2 (keba(1,1),vect(1),vk1(1), 12,12)
                  call AtxB_C (tr(1,1),vk1(1),zsc(1), 3,3,3, 3,4,3)
              else
                  call Axv1_v2 (zs(1,1),ulv(1),zsc(1), 12,12)
              endif
              kelout(5) = 1
          endif
      endif

c --- restoring force for newton-raphson
      if (kelin(6).eq.1) then
          if (nlgeom.eq.1.and.sstif_c.eq.1) then
              if (kelin(1).eq.1) then
                  call AtxB_C (tr(1,1),fe(1),zscnr(1), 3,3,3, 3,4,3)
              else
                  call AxB_C (tr(1,1),udef(1),vect(1), 3,3,3, 3,4,3)
                  call Axv1_v2 (keba(1,1),vect(1),vk1(1), 12,12)
                  call AtxB_C (tr(1,1),vk1(1),zscnr(1), 3,3,3, 3,4,3)
              endif
          else
              call Axv1_v2 (zs(1,1),udef(1),zscnr(1), 12,12)
          endif
          kelout(6) = 1
      endif

      if (kelout(1).eq.1) then
c --- compute elastic energy
          elener(1) = 0.5d0*v_dot(udef(1),forg(1),12)
      endif

      if (kstabnl.gt.0 .and. nrkey.eq.1) then
c --- apply non-linear stabilization but skip if it is the first substep
          add = 0
          call v_zero (StbDpVar(1),3)
          if ((A_KeyMatUpd.eq.0 .and. A_SubStepNum.eq.1 
     &         .and. kstabnl.lt.3) .or. (abs(StabFact).le.TINY)) 
     &        goto 100
     
          con = 1.0d0
          call v_zero (StbMat(1,1),num)
          call BeamEMass (con,ex,gxy,area,aiz,aiy,shearz,sheary,
     &                    alenv,epinit,0,nr,StbMat(1,1),A_MassLump)
          if (A_MassLump .eq. 0) 
     &        call MatRotTran (StbMat(1,1),tr(1,1),nr,3)
c --- neglect any translation and rotation coupling
          if (A_MassLump.ne.1) then
              do n = 1,2
                 istar = (n-1)*6
                 jstar = 0
                 do i = 1,2
                    do j = 1,3
                       StbMat(istar+1, jstar+j+3) = 0.0d0
                       StbMat(istar+2, jstar+j+3) = 0.0d0
                       StbMat(istar+3, jstar+j+3) = 0.0d0
                       StbMat(istar+4, jstar+j) = 0.0d0
                       StbMat(istar+5, jstar+j) = 0.0d0
                       StbMat(istar+6, jstar+j) = 0.0d0
                    enddo
                    jstar = jstar + 6
                 enddo
              enddo
          endif          

          call v_move (nrsvr(1),StbDpVar(1),3)
          add = nint(nrsvr(4))
          cDamp = 0.0d0          
          if (add.eq.0) then
              if (kstabnl.eq.3) then
c --- stabilization already during the first substep
                  add = 1
              elseif (kstabnl.eq.2 .and. A_SubStepNum.gt.1) then
c --- stabilization after the first substep
                  add = 1
              elseif (A_BisecNum.gt.0 .and. 
     &                A_TimeInc.eq.A_MinTimeStep) then
c --- stabilization after the convergence failure
                  add = 1
              endif
          else
              if (kstabnl.eq.3) then
                  if (A_BisecNum.gt.1 .and. 
     &                A_TimeInc.eq.A_MinTimeStep .and. 
     &                A_TimeIno.eq.A_MinTimeStep) then
                      cDamp = StbDpVar(1)
                      StbDpVar(1) = 
     &                             max(1.1d0,abs(rvr(7)))*StbDpVar(1)
     &                             *dble(A_BisecNum-1.0d0)
                  endif
     
              elseif (kstabnl.eq.1) then
                  if (A_BisecNum.eq.0) then
                      add = 0
                  else
                      cDamp = StbDpVar(1)
                      StbDpVar(1) = 
     &                             max(1.1d0,abs(rvr(7)))*StbDpVar(1)
     &                             *dble(A_BisecNum-1.0d0)
                  endif
              endif
          endif

          call EBStabForm (elem,nr,kelout(1),u(1,A_RowIncS),elener(1),
     &                     StabFact,StbDpVar(1),StbDpVar(2),StbDpVar(3),
     &                     zs(1,1),StbMat(1,1),zscnr(1),add)

c --- check if there is a saved factor
          if (cDamp.ne.0.0d0) StbDpVar(1) = cDamp
          
100       continue
c --- save stabilization variables
          if (A_KeyMatUpd.eq.1) then
              call v_move (StbDpVar(1),nrsvr(1),3)
              nrsvr(4) = dble(add)
c --- save non-linear SVRs
              call StoreSVR (svindx(1),1,nnrsvr,nrsvr(1))
          elseif (A_SubStepNum.eq.1 .and. A_IterNum.eq.1 
     &            .and. kstabnl.eq.3) then
              nrsvr(1) = StbDpVar(1)
              nrsvr(4) = dble(add)
c --- handle SVR files
              call StoreSVR (svindx(1),1,nnrsvr,nrsvr(1))
          endif
c --- write out the SVR index vector
          call PutSVRIndx (elem,locsvrL,svindx(1))  
      endif

      if (A_KeyOut.eq.1) then
          call cv1_v2 (forl(1),forx(1,1),6,-1.0d0)
          call v_move (forl(7),forx(1,2),6)
          stress(1) = 0.5d0*(forx(1,1)+forx(1,2))/area
          vect(1) = 0.5d0*stress(1)**2*alenv*area/ex
c --- bending energy along element y axis (account for linear distribution of moment)
          vect(2) = (forx(6,1)**2+forx(6,1)*forx(6,2)+forx(6,2)**2)
c --- add bending energy along element z axis
          vect(2) = vect(2) +
     &              (forx(5,1)**2+forx(5,1)*forx(5,2)+forx(5,2)**2)
c --- multiply by the coefficients
          vect(2) = 0.5d0*vect(2)*alenv/3.0d0/ex/aiy
c --- torsion energy
          vect(3) = 0.5d0*((forx(4,1)+forx(4,2))*0.5d0)**2*alenv/gxy/aip
c --- shear energy
          vect(4) = max(TINY,elener(1) - vect(1) - vect(2) - vect(3))
          vect(5) = stress(1)
c --- save to RST
           if (btest(eomask,W_OutMisc)) then
               call WriteToRST(EDEMS,lcerstL,edindxL(1),
     &                         ielc(NMSMIS),vect(1))
          endif
      endif

 999  continue

      return
      end
      
