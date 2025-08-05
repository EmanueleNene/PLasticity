      function user04()
c --- primary function: calculate the volume of all fibers and NFCs
c
#include "impcom.inc"
#include "constants.inc"
c
      PTRFTN  MemAlloc
      
      integer user04, ndgxyz, ndinqr, elmiqr, elmget, etyiqr, rlget
      double precision BeamLength
            
      integer i, j, node, elem, ndmax, elmax, nnoda, sect,
     &        elData(EL_DIM), nodes(W_NdMax), nnod, etyp, 
     &        nSubs, nLeng

      double precision rvr(W_NumContParam), Sub(5,3), th,
     &                 Ro, Ri, Volu_B, Volu_L, L, Area, 
     &                 XYZ(3,3), ZERO
      
      parameter (ZERO = 0.0d0)
      
      pointer (PtrNdXYZ,NdXYZ)
      double precision NdXYZ(3,*)
      character*(32) string, cName1, cName2
      data string /'    '/       
      
      
      user04 = -999
      elmax = elmiqr(0,DB_MAXDEFINED)
      ndmax = ndinqr(0,DB_MAXDEFINED)

      if (ndmax.le.1 .or. elmax.le.1) then
          call erhandler ('user04',5000,3,'No elements or nodes',
     &                     ZERO,' ') 
          goto 999
      endif

      PtrNdXYZ = MemAlloc(ndmax*3,W_LenDouble,'PtrNdXYZ        ')
      call v_zero (NdXYZ(1,1),ndmax*3)
      do node = 1,ndmax
         i = ndgxyz (node,NdXYZ(1,node))
      enddo

      Volu_B = ZERO
      Volu_L = ZERO
c --- fill in the fiber number for each node
      do elem = 1,elmax
         nnod = elmget(elem,elData(1),nodes(1))
         if (elData(EL_DEAD).eq.1) cycle
         etyp = etyiqr(elData(EL_TYPE),-JETYP)
c --- only beams and links go through
         if (etyp.ne.W_eBeam3nd .and. etyp.ne.W_eBeam2nd
     &       .and. etyp.ne.W_eLink .and. etyp.ne.W_eBeamEB) cycle
         
c -- count number of active nodes
         nnoda = 0
         do i = 1,nnod
            if (nodes(i).le.0) cycle
            nnoda = nnoda + 1
            call v_move (NdXYZ(1,nodes(i)),XYZ(1,i),3)
         enddo   
         
         
c --- skip unselected components
         if (nnoda.le.1) cycle

         L = BeamLength(nnoda,XYZ(1,1),0,3)

         call v_zero (rvr(1),W_NumContParam)
         i = rlget (elData(EL_REAL),rvr(1))
         
         if (etyp.eq.W_eBeam3nd .or. etyp.eq.W_eBeam2nd 
     &       .or. etyp.eq.W_eBeamEB) then
             sect = nint(rvr(1))
             Ro = rvr(2)
             Ri = rvr(3)
             if (sect.eq.4 .and. Ri.gt.ZERO) then
                 Area = W_Pi*(Ro**2-Ri**2)
             elseif (sect.eq.3) then
                 Area = W_Pi*Ro**2
             elseif (sect.eq.2) then
                 th = rvr(4)
                 Area = Ro*Ri - (Ro-2.0d0*th)*(Ri-2.0d0*th)
             elseif (sect.eq.1) then
                 Area = Ro*Ri
             endif
             Volu_B = Volu_B + L*Area
         elseif (etyp.eq.W_eLink) then
             Area = rvr(1)
             Volu_L = Volu_L + L*Area
         endif                 
      enddo      
      call MemFree (PtrNdXYZ)

c -- save to output arguments
      call parNameGetSubs (2,cName1,nLeng,nSubs,Sub(1,1)) 
      call pardef (cName1,0,nSubs,Sub(1,1),Volu_B,i,string)

      call parNameGetSubs (3,cName2,nLeng,nSubs,Sub(1,1)) 
      call pardef (cName2,0,nSubs,Sub(1,1),Volu_L,j,string)

      if (i.gt.0 .and. j.gt.0) then
          call erhandler ('user04',5000,3,'Error with output arguments',
     &                     ZERO,' ') 
          goto 999
      endif 
      
      user04 = 0

 999  continue            
       
      return
      end
