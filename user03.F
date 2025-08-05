      function user03()
c --- primary function: "distribute" fine material
c
#include "impcom.inc"
#include "constants.inc"
c
      PTRFTN  MemAlloc
      
      integer user03, ndgxyz, intinfun, ndinqr, 
     &        elmiqr, elmget, etyiqr
      
      double precision dpinfun, CalDist2P, seed, rnd, dist,
     &                 dens, area
            
      integer i, j, l, node, elem, ndmax, elmax, eType, eReal, eMat,
     &        elData(EL_DIM), nodes(W_NdMax), nnod, etyp

      double precision sRad, mFines, iFines
      
      pointer (PtrNdXYZ,NdXYZ)
      pointer (PtrNdStat,NdStat)
      integer  NdStat(2,*)
      double precision NdXYZ(3,*)
      
      user03 = -999

      sRad = dpinfun(2)
      eType = intinfun(3)
      eReal = intinfun(4)
      eMat = intinfun(5)
      mFines = dpinfun(6)
      
      if (sRad.lt.0.0d0.or.eType.le.0.or.eReal.le.0.or.eMat.lt.0) then
          call erhandler ('user03',5000,3,'Wrong paramteres passed',
     &                     0.0d0,' ') 
          goto 999
      endif
      
      call prope1 (1,eMat,13,0.0d0,dens)
      if (dens.le.0.0d0) then
          call erhandler ('user03',5000,3,
     &                   'Specified material has zero density',
     &                   0.0d0,' ') 
          goto 999
      endif

      call rlgetn(eReal,1,1,area)
      if (area.le.0.0d0) then
          call erhandler ('user03',5000,3,
     &                   'Specified real constant has zero area',
     &                   0.0d0,' ') 
          goto 999
      endif
      
      ndmax = ndinqr(0,DB_MAXDEFINED)
      elmax = elmiqr(0,DB_MAXDEFINED)
      
      if (ndmax.le.1 .or. elmax.le.1) then
          call erhandler ('user03',5000,3,'No elements or nodes',
     &                     0.0d0,' ') 
          goto 999
      endif
      
      PtrNdStat = MemAlloc(ndmax*2,W_LenInt,'PtrNdStat       ')
      PtrNdXYZ = MemAlloc(ndmax*3,W_LenDouble,'PtrNdXYZ        ')
      call v_zero (NdXYZ(1,1),ndmax*3)

c --- get all the coordinates for distance calculations
c     selection status is recorded as well to NdStat(*)
      do node = 1,ndmax
         NdStat(1,node) = ndgxyz (node,NdXYZ(1,node))
         NdStat(2,node) = 0
      enddo
c --- fill in the fiber number for each node
      do elem = 1,elmax
         nnod = elmget(elem,elData(1),nodes(1))
         etyp = etyiqr(elData(EL_TYPE),-JBEAM)
         if (etyp.gt.0) then
             do i = 1,nnod
                node = abs(nodes(i))
                if (node.eq.0) cycle
                if (NdStat(2,node).eq.0) then
                    NdStat(2,node) = elData(EL_REAL)
                elseif (NdStat(2,node).ne.elData(EL_REAL)) then
                    call erhandler ('user03',5000,3,
     &                              'Cannot be used with this model',
     &                               0.0d0,' ') 
                    goto 999
                endif
             enddo
         endif                 
      enddo
      call vi_zero (elData(1),EL_DIM)
      elData(EL_SHAPE) = 2
      elData(EL_DIM) = 3
      elData(EL_TYPE) = eType
      elData(EL_REAL) = eReal
      elData(EL_MAT) = eMat
      elem = elmax
      iFines = 0.0d0
      seed = 0.0d0
      if (mFines.gt.0.0d0) then
100      continue   
             call v_rand(rnd,1,seed)
             i = int(rnd*ndmax)+1
             if (NdStat(1,i).ne.1) goto 100
             nodes(1) = i
             call v_rand(rnd,1,seed)
             j = int(rnd*ndmax) + 1
             do l = 1,ndmax
                j = j + 1
                if (j.gt.ndmax) j = 1
                if (NdStat(1,j).ne.1) cycle
                if (i.eq.j) cycle
c --- are the nodes close enough?
                dist = CalDist2P(NdXYZ(1,i),NdXYZ(1,j))
                if (dist.le.sRad) then
c --- insure that nodes are not on the same fiber
                    if (NdStat(2,i).ne.NdStat(2,j)) then
c --- store the element
                        nodes(2) = j
                        elem = elem + 1
                        call elmput (elem,elData(1),2,nodes(1))
                        iFines = iFines + sRad*dens*area
c --- traget number of fines reached
                        if (iFines.ge.mFines) goto 200
                        goto 100
                    endif
                endif    
             enddo
          goto 100
      else
          do i = 1,ndmax-1
             if (NdStat(1,i).ne.1) cycle
             nodes(1) = i
             do j = i+1,ndmax
                if (NdStat(1,j).ne.1) cycle
c --- are the nodes close enough?
                if (CalDist2P(NdXYZ(1,i),NdXYZ(1,j)).le.sRad) then
c --- insure that nodes are not on the same fiber
                    if (NdStat(2,i).ne.NdStat(2,j)) then
c --- store the element
                        nodes(2) = j
                        elem = elem + 1
                        call elmput (elem,elData(1),2,nodes(1))
                    endif
                endif    
             enddo
          enddo   
      endif

200   continue
      
      call MemFree (PtrNdStat)
      call MemFree (PtrNdXYZ)

      user03 = 0
 999  continue            
       
      return
      end