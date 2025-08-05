      function user07()
c --- primary function: calculate the volume of all fibers and NFCs
c
#include "impcom.inc"
#include "constants.inc"
c
      integer MemIncPercent
      parameter ( MemIncPercent   = 10)

      
      integer ndiv, npoly
      double precision MinLen
      
      PTRFTN  MemAlloc, MemRealloc, MemAllocL, MemReallocL
      
      integer user07, intinfun, lenstr,
     &        ndinqr, elmiqr, rlinqr,
     &        elmget, ndgxyz, rlget, 
     &        GetWriteUnit
      
      integer   SplineCoeffs
      double precision  dpinfun, SegLength, PolyIntersect
      logical   IfContain
      
      pointer (PtrxyPoly,xyPoly)
      pointer (PtrNdXYZ,NdXYZ)
      pointer (PtrElData,ElData)
      pointer (PtrFibEndNd,FibEndNd)
      pointer (PtrFibLen,FibLen)
      pointer (PtrCoeffs,Coeffs)
      pointer (PtrNd2XYZ,Nd2XYZ)
      pointer (PtrSegLen,SegLen)
      
      character*(128)  Command     
      double precision  xyPoly(2,*)
      double precision  NdXYZ(3,*)
      integer ElData(EL_DIM,*)
      integer  FibEndNd(2,*)
      double precision FibLen(*)
      double precision Coeffs(12,*)
      double precision Nd2XYZ(3,*)
      double precision SegLen(*)
      
      integer  ElNd(4)
      integer  ElData0(EL_DIM)
      
      character*(32) string, cName1, cName2
      data string /'    '/
      double precision Sub(5,3), subc(1), dpValue
      integer nLeng, nSubs, kerr, lvl, iott
      
      integer ecount,  fcount, maxElPerFib
      integer count0,  fcount0, nrl, maxf, ncount, ncount0
      integer*8 iL
      integer elem  , node, fib, count, fnode, realfcount
      double precision rvr(6)
      double precision  t, t0, t1, xy(2), xyz(3), ElAvgLen, ElLen,len
      integer div

      integer i,j
c*************************************************************
c*************************************************************
      user07 = -1
c --- read input arguments (2-number of average length divisions,
c                           3-number of polygon sides
c                           4-polygon x's
c                           5-polygon y's
      iott = GetWriteUnit()
      
      write (iott,*) 'Enter the cutting and meshing module'
      
      ndiv   = intinfun(2)
      if (ndiv.eq.0) ndiv = 10
      nPoly  = abs(intinfun(3))
      MinLen = abs(dpinfun(6))
      
      if (nPoly.le.2) goto 900

      call parNameGetSubs (4,cName1,nLeng,nSubs,Sub(1,1)) 
      if (nLeng.le.0) goto 901

      call parNameGetSubs (5,cName2,nLeng,nSubs,Sub(1,1)) 
      if (nLeng.le.0) goto 901

      PtrxyPoly = MemAlloc(nPoly*2,W_LenDouble,'PtrxyPoly       ')

      do i = 1,nPoly
         subc(1) = dble(i)
         call parevl (cName1,1,subc(1),lvl,dpValue,string,kerr) 
         xyPoly(1,i) = dpValue
         call parevl (cName2,1,subc(1),lvl,dpValue,string,kerr) 
         xyPoly(2,i) = dpValue
      enddo   
      
c --- innitialize pointers
      ncount = ndinqr(0,DB_MAXDEFINED)
      ecount = elmiqr(0,DB_MAXDEFINED)
      fcount = rlinqr(0,DB_MAXDEFINED)
      fcount0 = nint(fcount*(1.d2+dble(MemIncPercent))/1.d2)
      ncount0 = nint(ncount*(1.d2+dble(MemIncPercent))/1.d2)
      iL = ncount0*3

      PtrNdXYZ = MemAllocL(iL,W_LenDouble,'PtrNdXYZ        ')
      PtrElData = MemAlloc(fcount0*EL_DIM,W_LenInt,'PtrElData       ')
      PtrFibEndNd = MemAlloc(fcount0*2 ,W_LenInt   ,'PtrFibEndNd     ')
      PtrFibLen = MemAlloc(fcount0 ,W_LenDouble,'PtrFibLen      ')
      
      write (iott,*) 'The memory has been allocated successfully'
      
      do iL = 1,ncount0
         do j = 1,3
            NdXYZ(j,iL) = 0.0d0
         enddo
      enddo
      
      call vi_zero(ElData(1,1)  ,fcount0*EL_DIM)
      call vi_zero(FibEndNd(1,1),fcount0*2)
      call v_zero (FibLen(1)    ,fcount0)
c --- read data from Ansys
      fib   = 0
      count = 0
      node  = 0
      do elem = 1,ecount
         i = elmget (elem,ElData0(1),ElNd(1))
         if (i .eq. 0) cycle
         if (ElData0(EL_REAL) .ne. fib) then
           fib = ElData0(EL_REAL)
           ElData(:,fib) = ElData0(:)
           node = node + 1
           i = ndgxyz (ElNd(1),NdXYZ(1,node))
           FibEndNd(1,fib) = node
         else
           do i=1,EL_DIM
              if (ElData0(i) .ne. ElData(i,fib)) goto 910
           enddo
         endif
         i = ndgxyz (ElNd(1),xyz)
         if ( NdXYZ(1,node) .ne. xyz(1) .or.
     &        NdXYZ(2,node) .ne. xyz(2) .or.
     &        NdXYZ(3,node) .ne. xyz(3) ) goto 920
         node = node+1
         i = ndgxyz (ElNd(3),NdXYZ(1,node))
         node = node+1
         i = ndgxyz (ElNd(2),NdXYZ(1,node))
         FibEndNd(2,fib) = node
      enddo
      
      write (iott,*) 'The nodes and elements are read into memory'
            
      maxElPerFib = maxval(FibEndNd(2,1:fcount)-FibEndNd(1,1:fcount))
      
c --- clear all nodes and elements
      call elsel(0,0)
      call ndsel(0,0)
      write (iott,*) 'All previous nodes and elements are cleared'
      
c      Command = 'edel,all'
c      call mainan (lenstr(Command),Command)
c      Command = 'ndel,all'
c      call mainan (lenstr(Command),Command)
c      Command = 'numcmp,elem'
c      call mainan (lenstr(Command),Command)
c      Command = 'numcmp,node'
c      call mainan (lenstr(Command),Command)
      
c*************************************************************
c   *  CUTTING AND CALCULATING AVERAGE LENGTH                * 
c*************************************************************
c --- innitialize coefficients pointer
      PtrCoeffs=MemAlloc(maxElPerFib*12,W_LenDouble,'PtrCoeffs       ')
c --- for each fiber ...
      maxf = fcount
      realfcount = 0
      fib = 0
      do while (fib .lt. fcount)
        fib = fib + 1
        fnode = FibEndNd(2,fib)-FibEndNd(1,fib)+1
        if (fnode .le. 1) cycle
        if (fib .le. maxf) realfcount = realfcount +1
c --- calculate spline coefficients
        i = SplineCoeffs(NdXYZ(:,FibEndNd(1,fib):FibEndNd(2,fib)),
     &                   fnode,Coeffs)
        if (i .ne. 0) goto 930
c --- calculate length of fiber
        if (ndiv .ge. 0) then
          t0 = 0.0d0
          t1 = 1.0d0
          do i = 1, fnode-1
            FibLen(fib) = FibLen(fib) + SegLength(Coeffs(:,i),t0,t1)
          enddo
        endif
c --- cut the fiber with polygon
c     cut the starting part of fiber which is out of polygon
        elem = 0
        do i = FibEndNd(1,fib), FibEndNd(2,fib)
          xy(1) = NdXYZ(1,i)
          xy(2) = NdXYZ(2,i)
          if (IfContain(nPoly,xyPoly,xy,0)) goto 120
          elem = elem + 1
        enddo
120     node = i
        if ( node .gt. FibEndNd(2,fib) ) then
          FibEndNd(1,fib) = 0
          FibEndNd(2,fib) = 0
        elseif ( node .gt. FibEndNd(1,fib) ) then
          t = PolyIntersect(Coeffs(:,elem),xyPoly,nPoly,xyz)
          FibEndNd(1,fib) = node - 1
          NdXYZ(:,node-1) = xyz
        else
          elem = 1
        endif
c     delete elements which are going out of polygon (until one goes in again)
        do i = FibEndNd(1,fib)+1,FibEndNd(2,fib)
          xy(1) = NdXYZ(1,i)
          xy(2) = NdXYZ(2,i)
          if (.not. (IfContain(nPoly,xyPoly,xy,0)) ) goto 130
          elem = elem + 1
        enddo
130     node = i
c     if fiber goes out of polygon ...
c     replace outgoing element with a new one (ending at the intersection point)
        if (node .le. FibEndNd(2,fib)) then
          fcount = fcount + 1
c         reallocate memory for fibers if needed
          if (fcount .ge. fcount0) then
            count0 = fcount0
            fcount0 = nint(fcount*(1.d2+dble(MemIncPercent))/1.d2)
            
            PtrElData   = MemRealloc(PtrElData  ,fcount0*EL_DIM,
     &                                  W_LenInt   ,'PtrElData       ')
            PtrFibEndNd = MemRealloc(PtrFibEndNd,fcount0*2,
     &                                  W_LenInt   ,'PtrFibEndNd     ')
            PtrFibLen   = MemRealloc(PtrFibLen  ,fcount0,
     &                                  W_LenDouble,'PtrFibLen       ')
      
            call vi_zero(ElData(1,count0+1)  ,(fcount0-count0)*EL_DIM)
            call vi_zero(FibEndNd(1,count0+1),(fcount0-count0)*2)
            call v_zero (FibLen(count0+1)    ,(fcount0-count0))
            write (iott,*) 
     &            'FIBER MEMORY INCREASED FROM',count0,' TO',fcount0
          endif
c         end of reallocate memory for fibers if needed
          FibEndNd(1,fcount) = ncount + 1
c         reallocate memory for nodes if needed
          ncount = ncount + 1 + FibEndNd(2,fib) - node
          if (ncount .ge. ncount0) then
            count0 = ncount0
            ncount0 = nint(ncount*(1.d2+dble(MemIncPercent))/1.d2)
            iL = ncount0*3
            PtrNdXYZ = MemReallocL(PtrNdXYZ,iL,W_LenDouble,
     &                                         'PtrNdXYZ        ')
            do iL = count0,ncount0
               do j = 1,3
                  NdXYZ(j,iL) = 0.0d0
               enddo
            enddo            
            
            write(iott,*) 
     &           'NODE MEMORY INCREASED FROM',count0,' TO',ncount0
          endif
c         end of reallocate memory for nodes if needed
          NdXYZ(:,(ncount- (FibEndNd(2,fib) - node)):ncount) =
     &                          NDXYZ(:,node:FibEndNd(2,fib))
          FibEndNd(2,fcount) = ncount
          t = PolyIntersect(Coeffs(:,elem),xyPoly,nPoly,xyz)
          FibEndNd(2,fib) = node
          NdXYZ(:,node) = xyz
          ElData(:,fcount) = ElData(:,fib)
          ElData(EL_REAL,fcount) = fcount
          nrl = rlget (fib,rvr(1))
          call rlput (fcount,nrl,rvr(1))
        endif
      enddo
c --- ... for each fiber
      if (ndiv .gt. 0) then
          ElAvgLen =(sum(fibLen(1:fcount))/realfcount)/ndiv
          write (iott,*) 'Calculated element size', ElAvgLen
      endif
      
      write (iott,*) 'Proceed to meshing...'
c*************************************************************
c   *  MESHING                                               *
c*************************************************************
c --- innitialize pointer for segment lengths
      ncount0   = fcount * (abs(ndiv)+1) * 2 + 1
      iL = ncount0*3
      PtrNd2XYZ = MemAllocL(iL   ,W_LenDouble,'PtrNd2XYZ       ')
      PtrSegLen = MemAlloc(maxElPerFib ,W_LenDouble,'PtrSegLen       ')
      
      do iL = 1,ncount0
         do j = 1,3
            Nd2XYZ(j,iL) = 0.0d0
         enddo
      enddo      

      ncount = 0
c --- for each fiber ...
      do fib = 1,fcount
        fnode = FibEndNd(2,fib)-FibEndNd(1,fib)+1
        if (fnode .lt. 2) cycle
c --- calculate spline coefficients
        i = SplineCoeffs(NdXYZ(:,FibEndNd(1,fib):FibEndNd(2,fib)),
     &                   fnode,Coeffs)
        if (i .ne. 0) goto 930
c --- calculate length of fiber and half_element_size(ElLen) for the fiber
        call v_zero(SegLen,maxElPerFib)
        t0 = 0.0d0
        t1 = 1.0d0
        do i = 1, fnode-1
           SegLen(i) = SegLength(Coeffs(:,i),t0,t1)
        enddo
        FibLen(fib) = sum(SegLen(1:fnode-1))
        if (FibLen(fib) .lt. MinLen) then
            FibEndNd(:,fib) = 0
            cycle
        endif
        if (ndiv .gt. 0) then
            div = max( nint(FibLen(fib) / ElAvgLen) , 1 )
        else
            div = abs(ndiv)
        endif
        ElLen = 0.5d0*FibLen(fib)/div
c --- start meshing
        ncount = ncount + 1
        Nd2XYZ(:,ncount) = NdXYZ(:,FibEndNd(1,fib))
        FibEndNd(1,fib) = ncount
        len = 0.0d0
        j = 0
        do i = 2, div*2
          do while (len .lt. ElLen)
            j = j+1
            len = len + SegLen(j)
          enddo
          len = len - ElLEn
          call NodeAtDist(SegLen(j)-len, Coeffs(:,j),xyz)
          ncount = ncount + 1
          Nd2XYZ(:,ncount) = xyz
        enddo
        ncount = ncount + 1
        Nd2XYZ(:,ncount) = NdXYZ(:,FibEndNd(2,fib))
        FibEndNd(2,fib) = ncount
        
      enddo
c --- ... for each fiber

      do fib =1,fcount
         if (FibEndNd(1,fib) .lt. FibEndNd(2,fib)) then
             do j = FibEndNd(1,fib),FibEndNd(2,fib)
                call ndpxyz(j,Nd2XYZ(:,j))
             enddo
         endif
      enddo
      
      ElNd(4) = 0
      count = 0
      do fib = 1,fcount
         if (FibEndNd(1,fib) .lt. FibEndNd(2,fib)) then
             do j = FibEndNd(1,fib),FibEndNd(2,fib)-1,2
                count = count + 1
                ElNd(1) = j
                ElNd(2) = j + 2
                ElNd(3) = j + 1
                call elmput(count,ElData(1,fib),4,ElNd)
             enddo
         else
             if (ElData(EL_REAL,fib).gt.0) then
                 call rlsel(ElData(EL_REAL,fib),0)
             endif
         endif
      enddo
      
      call MemFree(PtrxyPoly)
      call MemFree(PtrNdXYZ)
      call MemFree(PtrElData)
      call MemFree(PtrFibEndNd)
      call MemFree(PtrFibLen)
      call MemFree(PtrCoeffs)
      call MemFree(PtrNd2XYZ)
      call MemFree(PtrSegLen)
      
      write (iott,*) 'Meshing and cutting are finished, memory freed'
      
001   user07 = 0
      goto 999
900   call erhandler ('User07', 5000, 3, 'Number of polygon sides is %I
     &                     which is smaller than minimum value of 3.',
     &                     npoly,' ') 
      goto 999
901   call erhandler ('User07', 5000, 3, 'Error in polygon coordinates',
     &                     npoly,' ') 
      goto 999
910   call erhandler ('User07', 5000, 3, 'There is inconsistency in
     &                      element properties of fiber %I. Properties 
     &                      of different elements are not the same.',
     &                      fib,' ') 
      goto 999
920   call erhandler ('User07', 5000, 3, 'There is inconsistency in
     &                      nodes of fiber %I. Elements are not 
     &                      properly connected.',
     &                      fib,' ') 
      goto 999
930   call erhandler ('SplineCoeffs', 5000, 3, 'Less than 2 data points 
     &                      are sent for calculationg spline 
     &                      coefficients of fiber %I.',
     &                      fib,' ') 
      goto 999


999   return
      end