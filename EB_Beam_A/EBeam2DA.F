      subroutine EBeam2DA (elem,ielc,elmdat,eomask,nodes,locsvrL,kelreq,
     &                     kelfil,nr,xyz,u,kelout,zs,zass,gstif,
     &                     zsc,zscnr,elvol,elmass,center,elener,
     &                     edindxL,lcerstL)
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
c       kelfil (int,ar(10),in)      - keys indicating incoming mats and lvs
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
#include "ansysdef.inc"
#include "echprm.inc"
#include "soptcm.inc"
#include "stepcm.inc"
#include "elucom.inc"
#include "elecom.inc"
c
      double precision
     & v_dot,CalDist2P

      LONGINT  locsvrL, lcerstL, edindxL(25)

      integer
     & elem,ielc(IELCSZ),elmdat(EL_DIM),eomask,nodes(2),kelreq(10),
     & kelfil(10),nr,kelout(10)

      integer
     & mat,ireal,nrvr,nnrsvr,kelin(10),svindx(20),keyopt(12),
     & kstabnl,add,istar,jstar,n,i,j

      double precision
     & xyz(6,2),u(nr,5),zs(nr,nr),zass(nr,nr),gstif(nr,nr),
     & zsc(nr),zscnr(nr),elvol,elmass,center(3),elener(5)

      double precision
     & rvr(7),area,aiz,shearz,epinit,nrsvr(4),StabFact

      double precision
     & ex,nuxy,gxy,dens

      double precision
     & prop(13),vect(28),alenv,con,averot,uth(6),utho(6),
     & trv(3,3),trr(3,3),tr(3,3),udef(6),force,ulv(6),
     & forl(6),forg(6),uthoo(6),ustr(6),deba(6),udef0(6),
     & StbMat(6,6),StbDpVar(3)

      equivalence
     & (area,rvr(1)), (aiz,rvr(2)),
     & (shearz,rvr(3)), (epinit,rvr(4))

      equivalence
     & (ex,prop(1)), (nuxy,prop(7)),
     & (gxy,prop(4)),(dens,prop(13))

c --- controls formation of consistent stress-stiffness
      equivalence
     & (kstabnl,keyopt(3))

c --- define initial data
      mat = elmdat(EL_MAT)
      ireal = elmdat(EL_REAL)
      nrvr = ielc(NMTRLC)
      nnrsvr = ielc(NMNSVR)
      elener(1) = 0.0d0
      call vi_move (ielc(KYOP1),keyopt(1),12)

c --- get the element real constant data
      call rvrget (elem,ireal,ielc(1),nrvr,rvr(1))

      if (rvr(1).le.TINY .or. rvr(2).le.TINY) then
          call erhandler ('EBeam2DA',5000,3,
     &                    'Beam element %I has zero area!',
     &                     dble(elem), ' ')
          goto 999
      endif

      if (kstabnl.gt.0 .and. nrkey.eq.1) then
          StabFact = rvr(6)
c --- get the SVR index vector
          call svgidx (locsvrL,svindx(1))
c --- retrieve the nonlinear SVRs
          call svrget (svindx(1),4,nnrsvr,nrsvr(1))
      endif

c --- compute materials for every pass
      call GetElasticProps (elem,mat,0.0d0,prop(1))

c --- geometry computations

c --- compute original length
      alenv = CalDist2P (xyz(1,1),xyz(1,2))

c --- compute original transformation
      call v_zero (trv(1,1),9)
      trv(1,1) = (xyz(1,2) - xyz(1,1))/alenv
      trv(1,2) = (xyz(2,2) - xyz(2,1))/alenv
      trv(2,1) = -trv(1,2)
      trv(2,2) =  trv(1,1)
      trv(3,3) =  1.0d0
      call v_move (trv(1,1),tr(1,1),9)
c --- compute updated transformation for large rotations
      if (nlgeom.eq.1) then
          call v_zero (trr(1,1),9)
          call CalAveRot2D (u(3,putot),u(6,putot),averot)
          trr(1,1) =  cos(averot)
          trr(1,2) =  sin(averot)
          trr(2,2) =  trr(1,1)
          trr(2,1) = -trr(1,2)
          trr(3,3) =  1.0d0
          call AxB_C (trv(1,1),trr(1,1),tr(1,1), 3,3,3, 3,3,3)
      endif

c --- set up deformational displacements
      if (nlgeom.eq.0) then
          call v_move (u(1,putot),udef(1),6)
          call AxB_C (tr(1,1),udef(1),deba(1), 3,3,3, 3,2,3)
      else
          call CalUdef2D (xyz(1,1),6,averot,u(1,putot),udef(1))
          call AxB_C (tr(1,1),udef(1),deba(1), 3,3,3, 3,2,3)
      endif

c --- initial displacement vector relative to node i in global c.s.
      call v_zero (utho(1),6)
      con = - epinit*alenv
      call cv1_v2i (tr(1,1),3,utho(4),1, 2,con)

c --- compute element mass and centroid
      elvol = alenv*(1.0d0 - epinit)*area
      elmass = dens*elvol
      center(1) = (xyz(1,1) + xyz(1,2))*0.5d0
      center(2) = (xyz(2,1) + xyz(2,2))*0.5d0

c --- set form keys
      call vi_move (kelreq(1),kelin(1),10)
      if (ielc(NOEMAT).eq.1) kelin(1) = max(kelin(1),kelin(5),kelin(6))

c --- compute forces
      if (nrkey.eq.1 .or. sstif.ne.0) then
          call v_zero (forg(1),6)
c --- set up stress-causing displacements
          if (nrkey.eq.0) then
              call v_move (udef(1),ustr(1),6)
          else
c --- remove previous thermal effects
              call v_zero (uthoo(1),6)
              con = - epinit*alenv
              call cv1_v2i (tr(1,1),3,uthoo(4),1, 2,con)
              call v1mv2_v3 (udef(1),uthoo(1),ustr(1),6)
          endif
c --- calculate k x u
          if (kelfil(1).eq.1) call Axv1_v2 (zs(1,1),ustr(1),forg(1),6,6)
c --- subtract off applied loads
          if (kelfil(5).eq.1) call v1mv2_v1 (forg(1),zsc(1),6)
c --- rotate to element c.s.
          if (kelin(4).eq.1) 
     &        call AxB_C (tr(1,1),forg(1),forl(1), 3,3,3, 3,2,3)
c --- strain energy calculation
          call v_move (udef(1),udef0(1),6)
          call v1mv2_v1 (udef(1),utho(1),6)
          if (kelfil(1).eq.1) then
              call Axv1_v2 (zs(1,1),udef(1),forg(1), 6,6)
              elener(1) = 0.5d0*v_dot (udef(1),forg(1),6)
          endif
      endif

c --- stiffness matrix
      if (kelin(1).eq.1) then
          call v_zero (zs(1,1),36)
          call BeamEStif2D (ex,gxy,area,aiz,shearz,alenv,nr,zs(1,1))
c --- calculate k x u in local c.s.
          if (kelin(4).eq.1 .and. ielc(NOEMAT).eq.1)
     &        call Axv1_v2 (zs(1,1),deba(1),forl(1),6,6)
          call MatRotTran (zs(1,1),tr(1,1),nr,3)
          kelout(1) = 1
      endif

c --- mass matrix
      if (kelin(2).eq.1) then
          call v_zero (zass(1,1),36)
          con = dens
          if (con.ne.0.0d0) then
              call BeamEMass2D (con,ex,gxy,area,aiz,shearz,alenv,epinit,
     &                          0,nr,zass(1,1),lumpm)
              if (lumpm .eq. 0) 
     &            call MatRotTran (zass(1,1),tr(1,1),nr,3)
          endif
          kelout(2) = 1
      endif

c --- stress stiffness matrix
      if (kelin(4).eq.1) then
          call v_zero (gstif(1,1),36)
          if (kfstps.eq.1) then
c --- first step force
              force = area*ex*epinit
          else
              force = 0.5d0*(-forl(1) + forl(4))
          endif

          if (force.ne.0.0d0) then
              call BeamSStif2D (force,alenv,nr,gstif(1,1))
              call MatRotTran (gstif(1,1),tr(1,1),nr,3)
          endif
          kelout(4) = 1
      endif

c --- applied load vector
      if (kelin(5).eq.1) then
          call v_zero (zsc(1),6)
          if (epinit .ne. 0.0d0) then
c --- set up initial strain displacement in global c.s.
              call v_zero (uth(1),6)
              con = - epinit*alenv
              call cv1_v2i (tr(1,1),3,uth(4),1,2,con)
              if (nrkey.eq.0) then
                  call v_move (uth(1),ulv(1),6)
              else
                  call v1mv2_v3 (uth(1),utho(1),ulv(1),6)
              endif
              call Axv1_v2 (zs(1,1),ulv(1),zsc(1), 6,6)
              kelout(5) = 1
          endif
      endif

c --- restoring force for newton-raphson
      if (kelin(6).eq.1) then
          call Axv1_v2 (zs(1,1),udef(1),zscnr(1), 6,6)
          kelout(6) = 1
      endif

      if (kelout(1).eq.1 .and. ielc(NOEMAT).eq.1) then
c --- calculate elastic energy
          call Axv1_v2 (zs(1,1),udef(1),forg(1), 6,6)
          elener(1) = 0.5d0*v_dot(udef(1),forg(1),6)
      endif

      if (kstabnl.gt.0 .and. nrkey.eq.1) then
c --- apply non-linear stabilization but skip if it is the first substep
          add = 0
          call v_zero (StbDpVar(1),3)
          if ((kuphvr.eq.0 .and. isubst.eq.1 .and. kstabnl.lt.3) 
     &         .or. (abs(StabFact).le.TINY)) goto 100
          con = 1.0d0
          call v_zero (StbMat(1,1),36)
          call BeamEMass2D (con,ex,gxy,area,aiz,shearz,alenv,epinit,
     &                      0,nr,StbMat(1,1),lumpm)
          if (lumpm .eq. 0) 
     &        call MatRotTran (StbMat(1,1),tr(1,1),nr,3)
c --- neglect any translation and rotation coupling
          if (lumpm.ne.1) then
              do n = 1,2
                 istar = (n-1)*3
                 jstar = 0
                 do i = 1,2
                    do j = 1,2
                       StbMat(istar+j, jstar+3) = 0.0d0;
                       StbMat(istar+3, jstar+j) = 0.0d0;
                    enddo
                    jstar = jstar + 3
                 enddo
              enddo
          endif          

          call v_move (nrsvr(1),StbDpVar(1),3)
          add = nint(nrsvr(4))

          if (add.eq.0) then
              if (kstabnl.eq.3) then
c --- stabilization already during the first substep
                  add = 1
              elseif (kstabnl.eq.2 .and. isubst.gt.1) then
c --- stabilization after the first substep
                  add = 1
              elseif (ibisec.gt.0 .and. timinc.eq.deltmn) then
c --- stabilization after the convergence failure
                  add = 1
              endif
          else
              if (kstabnl.eq.3) then
                  if (timinc.eq.deltmn .and. 
     &               (ibisec.gt.0 .or. timino.eq.deltmn)) then
                      StbDpVar(1) = max(1.1d0,abs(rvr(7)))*StbDpVar(1)
     &                              *dble(ibisec)
                  endif
              endif
          endif

          call EBStabFormA(elem,nr,kelout(1),u(1,puincn),elener(1),
     &                     StabFact,StbDpVar(1),StbDpVar(2),StbDpVar(3),
     &                     zs(1,1),StbMat(1,1),zscnr(1),add)

 100      continue
c --- save stabilization variables
          if (kuphvr.eq.1) then
              call v_move (StbDpVar(1),nrsvr(1),3)
              nrsvr(4) = dble(add)
          endif
c --- save non-linear SVRs
          call svrput (svindx(1),4,nnrsvr,nrsvr(1))

c --- write out the svr index vector
          call svpidx (locsvrL,svindx(1))
      endif

 999  return
      end