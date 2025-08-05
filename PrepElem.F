      subroutine PrepElem (elem,kelout,nr,elData,elChar,
     &                     ls,Kmtrx,Fnr,uelm,uLen,
     &                     force,forcnr,forceStb,kReform)
c --- primary function: prepare an element matrix and load vector for solution
c
c  input arguments:
c     elem     (int,sc,in)         - current element number
c     nr       (int,sc,in)         - number of rows in element matrix
c     elData   (int,ar(EL_DIM),in) - element data vector
c     elChar   (int,ar(*),in)      - element type data vector
c     ls       (int,ar(nr),in)     - dof vector for this element
c     Kmtrx    (dp,ar(nr,nr,4),in) - element matrices
c     Fnr      (dp,ar(nr,2),in)    - element load vector(s)
c     uLen     (int,sc,in)         - length of disp vector (number of DOFs)
c     kReform  (int,sc,in)         - stiffness matrix reform flag
c
c  output arguments:
c     force    (dp,ar(uLen),out)   - accumulated global applied force vector
c     forcnr   (dp,ar(uLen),out)   - accumulated global n-r restoring force
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "dof_info.inc"
#include "in_mem.inc"
c
      integer  elem,kelout(10),nr,elData(EL_DIM),elChar(*),
     &         ls(nr),uLen,lw(W_NrMax),
     &         nr2,i,j,KeyStabLoc,kReform
      double precision  Kmtrx(nr,nr,4),Fnr(nr,3),uelm(nr,5),
     &                  force(uLen),forcnr(uLen),forceStb(uLen),
     &                  Fstb(nr),kdstif,Fnrc(2*W_NrMax)


#if defined(DEBUG)
      call SubBeg('PrepElem')
#endif

c --- adjust the matrices and load vector as appropriate
      nr2 = nr*nr
c --- zero load vectors if they do not exist - they may may modified later
      if (kelout(5).eq.0) then
          call v_zero (Fnr(1,1),nr)
          kelout(5) = 1
      endif
      if (kelout(6).eq.0) then
          call v_zero (Fnr(1,2),nr)
          kelout(6) = 1
      endif
c --- process dead elements (just killed)
      if (elData(EL_DEAD).eq.1) then
          kdstif=1.d-6
          if (kelout(1).ne.0) then
              if (elChar(JBEAM).gt.0) then
                  call v_zero (Kmtrx(1,1,1),nr2)
                  call cv1_v1 (Kmtrx(1,1,1),nr2,kdstif)
              else
                  call v_zero (Kmtrx(1,1,1),nr2)
              endif
          endif
          if (kelout(2).ne.0) call v_zero (Kmtrx(1,1,2),nr2)
          if (kelout(3).ne.0) call v_zero (Kmtrx(1,1,3),nr2)
          if (kelout(4).ne.0) call v_zero (Kmtrx(1,1,4),nr2)
          if (kelout(5).ne.0) call v_zero (Fnr(1,1),nr)
          if (kelout(6).ne.0) call v_zero (Fnr(1,2),nr)
      endif

c --- get stabilization force
      KeyStabLoc = 0
      if (A_KeyStab.gt.0.and.elChar(NMSTVR).gt.0) KeyStabLoc = 1
      if (KeyStabLoc.gt.0)
     &    call GetStabData (elem,elChar(1),nr,2,nr,Fstb(1))

c --- adjust matrices and load vector for transient
      if (A_Trans.gt.0) then
         call DynMatAdj (nr,kelout(1),uelm(1,A_RowIncS),
     &                   uelm(1,A_RowVelo),uelm(1,A_RowAcce),
     &                   Kmtrx(1,1,2),Kmtrx(1,1,1),
     &                   Fnr(1,1),Fnr(1,2),kReform)
      endif

      if (A_Gravity.gt.0 .and. kelout(2).ne.0) then
          call AddGrav (nr,ls(1),Kmtrx(1,1,2),kelout(5),Fnr(1,1))
      endif

c --- further adjustments for reformed global matrix
c --- copy element dof vector
      call vi_move (ls(1),lw(1),nr)
c --- copy element load vectors: used for summing forces
      if (kelout(5).eq.1) then
          call v_move(Fnr(1,1),Fnrc(1),nr)
      else
          call v_zero(Fnrc(1),nr)
      endif
      if (kelout(6).eq.1) then
          call v_move(Fnr(1,2),Fnrc(nr+1),nr)
      else
          call v_zero(Fnrc(nr+1),nr)
      endif

c --- adjust for the same dofs
      if (D_NumCp .gt. 0) call SetCoup (nr,lw(1))
c --- adjust for specified displacements
      if (kelout(1).eq.1) !.and. A_IterNum.eq.1) 
     &    call SetDisp (nr,lw(1),kelout(1),
     &                  Kmtrx(1,1,1),Fnrc(1),M_DofBits(1))

c --- add element load vector to total if required
      call lockSMP(W_LockLoadVect)
!dir$ ivdep
      do i = 1,nr
         j = lw(i)
         if (j.gt.0) then
             force(j) = force(j) + Fnrc(i)
             forcnr(j) = forcnr(j) + Fnrc(i+nr)
             if (KeyStabLoc.gt.0) forceStb(j) = forceStb(j) + Fstb(i)
         endif
      enddo
      call unlockSMP(W_LockLoadVect)
c --- sum the applied loads and the n-r restoring force for the solver (REDUNDANT NOW)
c      if (kelout(6).eq.1) then
c          call v1mv2_v1 (Fnr(1,1),Fnr(1,2),nr)
c          kelout(5) = 1
c      endif

      kelout(1) = 1
      kelout(3) = 0
      kelout(4) = 0

#if defined(DEBUG)
      call SubEnd('PrepElem')
#endif

      return
      end
