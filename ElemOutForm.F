      subroutine ElemOutForm (elem,elChar,elData,eomask,nodes,nr,
     &                        ls,rforce,kelfil,xyzang,uelm,Kmtrx,Fnr,
     &                        elvol,elener,lcSVRL,EleResultIdxL,lcerstL)
c --- primary function: set up the call to the element routines needed for output
c
c  input arguments:
c       elem    (int,sc,in)             - element number
c       elChar  (int,ar(*),in)          - array of element type characteristics
c       elData  (int,ar(EL_DIM),in)     - array of element data
c       eomask  (int,sc,in)             - bit pattern for element output
c       nodes   (int,ar(*),in)          - array of element node numbers
c       nr      (int,sc,in)             - matrix and lv size
c       kelfil  (int,ar(10),in)         - keys indicating incoming mats and LV's
c       xyzang  (dp,ar(6,ndim),in)      - nodal coords (orig) and rotation angles
c       uelm    (dp,ar(nr,5),in)        - element nodal solution values
c       Kmtrx   (dp,ar(nr,nr,4),inout)  - matrices
c       Fnr     (dp,ar(nr,2),inout)     - load vectors
c       lcSVRL  (int*8,sc,inout)        - SVR file pointer
c       lcerstL (int*8,sc,inout)        - result file pointer
c
c  output arguments:
c       EleResultIdxL(intL,ar(25),out)  - index to element result file data
c       elvol  (dp,sc,out)              - element volume
c       elener (dp,ar(10),out)          - element energies
c
#include "impcom.inc"
#include "in_mem.inc"
#include "constants.inc"
#include "solu_info.inc"
c
      integer*8  lcSVRL, lcerstL, EleResultIdxL(25)
      integer elem,elChar(*),elData(EL_DIM),nodes(*),nr,i,j,
     & kelfil(10),kelout(10),eomask,kelreqLoc(10),ls(*),idummy
      double precision  xyzang(6,*),uelm(nr,*),Kmtrx(nr,nr,4),Fnr(nr,2),
     & elvol,elener(10),elmass,center(3),rforce(*),forc(2*W_NrMax)


#if defined(DEBUG)
      call SubBeg('ElemOutForm')
#endif

c --- initialize request keys
      call vi_zero (kelreqLoc(1),10)

      if (A_Trans.gt.0) then
c --- transient always needs stiffness and mass
          kelreqLoc(1) = 1
          kelreqLoc(2) = 1
          kelreqLoc(3) = 1
      endif

      if (A_Gravity.eq.1) then
c --- gravity needs mass
          kelreqLoc(2) = 1
      endif 
            
c --- request appropriate matrices and LVs to compute forces, etc.
      kelreqLoc(5) = 1
      kelreqLoc(6) = 1

c --- initialize the output form keys
      call vi_zero (kelout(1),10)

c --- call the element routines
      call ElemDrivers (elem,elChar(1),elData(1),eomask,nodes(1),lcSVRL,
     &                  kelreqLoc(1),nr,xyzang(1,1),uelm(1,1),kelout(1),
     &                  Kmtrx(1,1,1),Kmtrx(1,1,2),
     &                  Kmtrx(1,1,3),Kmtrx(1,1,4),
     &                  Fnr(1,1),elmass,center(1),elvol,elener(1),
     &                  EleResultIdxL(1),lcerstL,idummy)
      
      if (elData(EL_DEAD) .eq. 0) then

c --- add to the nodal force sum
          call v_zero (forc(1),2*nr)

c --- generate inertial load vector
          if (A_Gravity.gt.0 .and. kelout(2).eq.1) then
              call AddGrav (nr,ls(1),Kmtrx(1,1,2),kelout(5),Fnr(1,1))
          endif

c --- element calculated K x u
          if (kelout(5).eq.1) call v1pv2_v1 (forc(1),Fnr(1,1),nr)
          if (kelout(6).eq.1) call v1mv2_v1 (forc(1),Fnr(1,2),nr)

c --- check if conditions are appropriate to calculate inertia force
          if (kelout(2).eq.1 .and. A_Trans.gt.0) then
c --- calculate forces
              call Axv1_v2 (Kmtrx(1,1,2),uelm(1,A_RowAcce),forc(1+nr),
     &                      nr,nr)
          endif
      
          call lockSMP (W_LockLoadVect)

!dir$ ivdep
          do i = 1,nr
             j = ls(i)
             if (j.gt.0) rforce(j) = rforce(j) - forc(i) - forc(i+nr)
          enddo

          call unlockSMP (W_LockLoadVect)

      endif

      call vi_move (kelout(1),kelfil(1),10)

#if defined(DEBUG)
      call SubEnd('ElemOutForm')
#endif

      return
      end

