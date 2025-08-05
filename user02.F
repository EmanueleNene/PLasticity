      function user02()
c --- primary function: read in Stefan's files
c
#include "impcom.inc"
#include "constants.inc"
#include "file_sys.inc"
c
      integer user02, Sys_Open, GetWriteUnit, ndgxyz, Sys_Close,
     &        FindIndxMax, FindIndxMin, rlgetn, elmget, rlinqr, etyget
      integer nxyz, kerr, iott, i, j, k, nunit, node, ncount,
     &        elem, ecount, ndpel, nodes(4), elData(EL_DIM), 
     &        elChar(IELCSZ), ireal, rcount, imat, mcount, ET1, ET2,
     &        jac_err, dub_err, intinfun, to_check, jreal, 
     &        cRl, nRl, ishift, elDatab(EL_DIM)

      double precision xyz(3), rvr(6), Ex, Gxy, nuxy, SI_Coef, tol,
     &                 xyz3(3,3), CheckJacBeam, Ex1, Ex2,
     &                 FibMinDist,length1, length2, CalDist2P,
     &                 dpinfun

      character*(B_NameLen)  fname
      character*(B_DirLen)  fdir  
      character*(B_PathLen)  pname
      character*32  ch32infun
      character*28  elcdn
      PTRFTN  MemAlloc
      pointer (PtrElNd,ElNd)
      pointer (PtrNdXYZ,NdXYZ)
      pointer (PtrElBox,ElBox)
      pointer (PtrRl,Rl)
      integer  ElNd(4,*), ElBox(6,*), Rl(*)
      double precision  NdXYZ(3,*)
      logical ElBoxOverlap,IfEqvDP
      parameter (nxyz = 3)

      call Buf_Ini
      
      fname = ch32infun(2)
      call getjobinfo ('CDIR',fdir,i)
      nunit = 9
      ndpel = 3
      iott = GetWriteUnit()
      kerr = 0
     
c --- IMPORT NODES

      call MakeName (fname,'xyz     ',fdir,pname,kerr)
      
      kerr = Sys_Open (nunit,pname,'UN','SERIAL','FORM','READ','BEGIN')
      
      if (kerr .ne. 0) then
          write (iott,*) 'Problem reading ',trim(fname),'.xyz'
          goto 998
      else
          write (iott,*) 'Reading ',trim(fname),'.xyz'
      endif
 
      to_check = intinfun(3)
      
      SI_Coef = dpinfun(4)     
      if (SI_Coef.le.TINY) SI_Coef = 1.0d6
      
      ncount = 0
      
  10  read (nunit,*,err=990,end=990)  node,(xyz(i),i=1,nxyz)
      
      if (node .le. 0) goto 990

      ncount = ncount + 1
c --- convert to microns
      call cv1_v1 (xyz(1),nxyz,SI_Coef)
c --- store the node      
      call ndpxyz (node,xyz(1))

      goto 10

 990  write (iott,2001) ncount
 2001 format (/'Nodes imported:',i8/)
 
      kerr = Sys_Close(nunit,'KEEP')

      if (ncount.eq.0) then
          write (iott,*) 'No nodes found in the file'
          goto 998
      endif

c --- IMPORT FIBERS
      call MakeName (fname,'nod     ',fdir,pname,kerr)
     
      kerr = Sys_Open (nunit,pname,'UN','SERIAL','FORM','READ','BEGIN')

      if (kerr .ne. 0) then
          write (iott,*) 'Problem reading ',trim(fname),'.nod'
          goto 998
      else
          write (iott,*) 'Reading ',trim(fname),'.nod'
      endif

      ecount = 0

      i = etyget(1,elChar(1))
      if (i.eq.0) call vi_zero (elChar(1),IELCSZ)
      
      call uec104 (elcdn,elChar(1),kerr)
      call etyput (1,IELCSZ,elChar(1))
      call vi_zero (elData(1),EL_DIM)
      elData(EL_TYPE) = 1
      elData(EL_SHAPE) = 2
      elData(EL_DIM) = 3
      nodes(4) = 0
      dub_err = 0
      jac_err = 0
      
  20  read (nunit,*,err=991,end=991)  elem,(nodes(i),i=1,ndpel), 
     &                                elData(EL_REAL), elData(EL_MAT)

      if (elem .le. 0) goto 991
      
      if (nodes(3).gt.0) then
         
         do i = 1,3
            j = ndgxyz (nodes(i),xyz3(1,i))
         enddo

         if (to_check.le.0 .or. to_check.gt.2) goto 21
c --- check dublicated fibers
         do i = 1,ncount
            j = ndgxyz (i,xyz(1))
            do k = 1,3
               if (nodes(k).eq.i) cycle
               if (IfEqvDP(xyz3(1,k),xyz(1)) .and. 
     &             IfEqvDP(xyz3(2,k),xyz(2)) .and. 
     &             IfEqvDP(xyz3(3,k),xyz(3))) then
                   if (nodes(k).gt.i) then
                       dub_err = dub_err + 1
                       goto 20
                   endif
               endif
            enddo
         enddo
         if (to_check.eq.1) goto 21
c --- check jacobian
         if (CheckJacBeam(xyz3(1,1)).lt.0.0d0) then
             jac_err = jac_err + 1
             goto 20
         endif

      endif

  21  ecount = ecount + 1
c --- store the element + 1 zero node
      call elmput (ecount,elData(1),ndpel+1,nodes(1))

      goto 20

 991  write (iott,2002) ecount
 2002 format (/'Element imported:',i8/)

      
      if (dub_err.gt.0) 
     &    call erhandler ('user02',5000,2,'Number of dublicated 
     &                     fibers %i',dble(dub_err),' ') 

      if (jac_err.gt.0) 
     &    call erhandler ('user02',5000,2,'Number of distorted 
     &                     fibers %i',dble(jac_err),' ') 

      kerr = Sys_Close(nunit,'KEEP')

      if (ecount.eq.0) then
          write (iott,*) 'No element found'
          goto 998
      endif      

c --- IMPORT DIAMETERS
      call MakeName (fname,'typ     ',fdir,pname,kerr)
     
      kerr = Sys_Open (nunit,pname,'UN','SERIAL','FORM','READ','BEGIN')
      
      if (kerr .ne. 0) then
          write (iott,*) 'Problem reading ',trim(fname),'.typ'
          goto 998
      else
          write (iott,*) 'Reading ',trim(fname),'.typ'
      endif
      
      rcount = 0
  30  read (nunit,*,err=992,end=992)  ireal, rvr(1:4)

      if (ireal .le. 0) goto 992
c --- real number already defined
      if (rlinqr(ireal,DB_SELECTED).ne.0) goto 992
      
      rcount = rcount + 1
c --- scale with a SI coefficient
      call cv1_v1 (rvr(2),3,SI_Coef)
      call rlput (ireal,4,rvr(1))

c --- store the real constants

      goto 30

 992  write (iott,2003) rcount
 2003 format (/'Cross section imported:',i8/)

      kerr = Sys_Close(nunit,'KEEP')

      if (rcount.eq.0) then
          write (iott,*) 'No real constants found'
          goto 998
      endif      
      
c --- IMPORT ELASTIC PROP
      call MakeName (fname,'mat     ',fdir,pname,kerr)
     
      kerr = Sys_Open (nunit,pname,'UN','SERIAL','FORM','READ','BEGIN')

      if (kerr .ne. 0) then
          write (iott,*) 'Problem reading ',trim(fname),'.mat'
          goto 998
      else
          write (iott,*) 'Reading ',trim(fname),'.mat'
      endif
      
      mcount = 0
  40  read (nunit,*,err=993,end=993)  imat, Ex, Gxy

      if (imat .le. 0) goto 993
      mcount = mcount + 1
c --- store the material prop
      Ex = Ex/SI_Coef
      Gxy = Gxy/SI_Coef
      call mpput (imat,1,1,0.0d0,Ex)
c      call mpput (imat,7,1,0.0d0,Gxy)
c      nuxy = 0.5d0*Ex/Gxy - 1 
      nuxy = 0.3d0
      call mpput (imat,4,1,0.0d0,nuxy)
            
      goto 40

 993  write (iott,2004) mcount
 2004 format (/'Mat prop imported:',i8/)

      kerr = Sys_Close(nunit,'KEEP')
      
      if (mcount.eq.0) then
          write (iott,*) 'No materila properties found'
          goto 998
      endif      
      
      if (to_check.eq.3) goto 997
      
      PtrElNd  = MemAlloc(ecount*4,W_LenInt,'PtrElNd         ')
      PtrNdXYZ = MemAlloc(ncount*3,W_LenDouble,'PtrNdXYZ        ')
      PtrElBox = MemAlloc(ecount*6,W_LenInt,'PtrElBox        ')
      PtrRl = MemAlloc(rcount*rcount,W_LenInt,'PtrRl           ')
      call v_zero (NdXYZ(1,1),ncount*3)
      call vi_zero (ElNd(1,1),ecount*4)
      
c --- get elements
      do elem = 1,ecount
         i = elmget (elem,elData(1),ElNd(1,elem))
      enddo

c --- get nodes
      do node = 1,ncount
         i = ndgxyz (node,NdXYZ(1,node))
      enddo

c --- create element boxes
      do i = 1,ecount
          do j = 1,nxyz
             ElBox(j,i)      = FindIndxMax(ElNd(1,i),NdXYZ(1,1),ndpel,j)
             ElBox(j+nxyz,i) = FindIndxMin(ElNd(1,i),NdXYZ(1,1),ndpel,j)
          enddo
      enddo

      i = etyget(3,elChar(1))
      if (i.eq.0) call vi_zero (elChar(1),IELCSZ)

      call uec101 (elcdn,elChar(1),kerr)
      call etyput (3,IELCSZ,elChar(1))

      i = etyget(2,elChar(1))
      if (i.eq.0) call vi_zero (elChar(1),IELCSZ)
     
      call uec102 (elcdn,elChar(1),kerr)
      call etyput (2,IELCSZ,elChar(1))
      
      call vi_zero (elData(1),EL_DIM)
      elData(EL_SHAPE) = 2
      elData(EL_DIM) = 3

      elem = ecount

      call v_zero(rvr(1),6)
      nRl = rcount
      call vi_zero(Rl(1),nRl*nRl)
      do i = 1,ecount
c --- get Ex for the first fiber
         node = elmget (i,elDatab(1),nodes(1))
         ireal = elDatab(EL_REAL)
c --- get radius for the first fiber
         k = rlgetn(ireal,2,1,rvr(1))
         ishift = (ireal-1)*nRl
         call prope1 (i,elDatab(EL_MAT),1,0.0d0,Ex1)
         length1 = CalDist2P(NdXYZ(1,ElNd(1,i)),NdXYZ(1,ElNd(2,i)))
         do j = i+1, ecount
            node = elmget(j,elDatab(1),nodes(1))
            jreal = elDatab(EL_REAL)
c --- do not consider self contact
            if (ireal.eq.jreal) cycle
c --- get radius for the second fiber
            k = rlgetn(jreal,2,1,rvr(2))
            tol = 4.0d0*max(rvr(1),rvr(2)) ! 4 mostly for ZD
c --- find out if elements overlap (super fast)
            if (ElBoxOverlap(ElBox(1,1),NdXYZ(1,1),i,j,tol)) then
              length2 = CalDist2P(NdXYZ(1,ElNd(1,j)),NdXYZ(1,ElNd(2,j)))
              tol = sqrt(tol**2 + max(length1/100.d0,length2/100.d0)**2)
c --- find out if distace exceeds the threshold (fast, but still safe)
              if (FibMinDist(ElNd(1,1),NdXYZ(1,1),i,j).le.tol) then
c --- do we need a new real number or the fibers are already in contact?
                  cRl = Rl(ishift+jreal)
                  if (cRl.eq.0) then
c --- if this is the first contact for this pair of beams
                      rcount = rcount + 1
                      cRl = rcount
c --- get Ex the for second fiber
                      call prope1 (j,elDatab(EL_MAT),1,0.0d0,Ex2)
c --- store real data for this foundation
                      call rlput (rcount, 5,rvr(1))
c --- decide which one should be slave (stiffer, thicker)
                      if (Ex1*rvr(1)**3 .lt. Ex2*rvr(2)**3) then
                          cRl = -cRl
                          ET1 = 2
                          ET2 = 3
                      else
                          ET1 = 3
                          ET2 = 2
                      endif
c --- the desision about the contact pair is defined in the sign of cRl
                      Rl(ishift+jreal) = cRl
                  else
                      if (cRl.lt.0) then
                          ET1 = 2
                          ET2 = 3
                      else
                          ET1 = 3
                          ET2 = 2
                      endif
                  endif

                  elData(EL_MAT) = 0
                  elData(EL_REAL) = abs(cRl)
                  elem = elem + 1
c --- store first segment
                  elData(EL_TYPE) = ET1
                  call elmput (elem,elData(1),ndpel,ElNd(1,i))

c --- store second segment
                  elem = elem + 1
                  elData(EL_TYPE) = ET2
                  call elmput (elem,elData(1),ndpel,ElNd(1,j))
              endif
            endif
         enddo
      enddo

      call MemFree (PtrElNd)
      call MemFree (PtrNdXYZ)
      call MemFree (PtrElBox)
      call MemFree (PtrRl)
      
997   user02 = 0
      goto 999

998   continue
      
      user02 = -999
      call erhandler ('user02', 5000, 3, 'Wrong input data.
     &                 Check the input filename and data format.',
     &                 0.0d0,' ')
      call StopAll (-5)          

999   continue

      return
      end
      
      function FindIndxMax (Vi,Vdp,n,ind1)
c --- primary function: find index of the maximum element
c
#include "impcom.inc"
      
      integer FindIndxMax, n, Vi(n), i, ind1
      double precision CurMax, Vdp(3,*)
      
      FindIndxMax = Vi(1)
      CurMax = Vdp(ind1,Vi(1))
      do i=2,n
         if (Vdp(ind1,Vi(i)).gt.CurMax) then
           FindIndxMax = Vi(i)
           CurMax = Vdp(ind1,Vi(i))
         endif
      enddo

      return
      end

      function FindIndxMin (Vi,Vdp,n,ind1)
c --- primary function: find index of the maximum element
c
#include "impcom.inc"
c      
      integer FindIndxMin, n, Vi(n), i, ind1
      double precision CurMin, Vdp(3,*)
      
      FindIndxMin = Vi(1)
      CurMin = Vdp(ind1,Vi(1))
      do i=2,n
         if (Vdp(ind1,Vi(i)).lt.CurMin) then
           FindIndxMin = Vi(i)
           CurMin = Vdp(ind1,Vi(i))
         endif
      enddo

      return
      end

      function ElBoxOverlap (ElBox,NdXYZ,i1,i2,tol)
c --- primary function: find out if element boxes overlap
c
#include "impcom.inc"
c
      logical ElBoxOverlap
      integer ElBox(6,*), i1, i2
      double precision NdXYZ(3,*), tol, aL, aU

      ElBoxOverlap = .false.
c --- try to build overlap volume
      aU=min(NdXYZ(1,ElBox(1,i1)),NdXYZ(1,ElBox(1,i2)))      
      aL=max(NdXYZ(1,ElBox(4,i1)),NdXYZ(1,ElBox(4,i2)))
 
      if (aL.gt.(aU+tol)) goto 999
 
      aU=min(NdXYZ(2,ElBox(2,i1)),NdXYZ(2,ElBox(2,i2)))      
      aL=max(NdXYZ(2,ElBox(5,i1)),NdXYZ(2,ElBox(5,i2)))

      if (aL.gt.(aU+tol)) goto 999
      
      aU=min(NdXYZ(3,ElBox(3,i1)),NdXYZ(3,ElBox(3,i2)))      
      aL=max(NdXYZ(3,ElBox(6,i1)),NdXYZ(3,ElBox(6,i2)))

      if (aL.gt.(aU+tol)) goto 999
      
      ElBoxOverlap = .true.
 
 999  continue
 
      return
      end
      
      function FibMinDist (ElNd,NdXYZ,i1,i2)
c --- primary function: find minimum distanse between two fibers
c
#include "impcom.inc"
c
      double precision FibMinDist, CalDist2P
      integer ElNd(4,*), i1, i2, k, i, j, nDiv
      parameter (nDiv = 50)  ! 50 should cover everything (was 30)
      double precision NdXYZ(3,*), s, step, mind, c1, c2, c3,
     &                 xi(3,3), xj(3,3), x1(3,nDiv), x2(3,nDiv)

      do i = 1,3    ! nodes
         do j = 1,3     ! coords
            xi(j,i) = NdXYZ(j,ElNd(i,i1))
            xj(j,i) = NdXYZ(j,ElNd(i,i2))
         enddo
      enddo

c --- form coordinate array
      step = 2.0d0/dble(nDiv)
      s = -1.0d0
      do i = 1,nDiv
         c1 = -s+s*s
         c2 =  s+s*s
         c3 =  1-s*s
         do k = 1,3
            x1(k,i) = 0.5d0*(c1*xi(k,1)+c2*xi(k,2)) + c3*xi(k,3)
         enddo
         do k = 1,3
            x2(k,i) = 0.5d0*(c1*xj(k,1)+c2*xj(k,2)) + c3*xj(k,3)
         enddo
         s = s + step
      enddo

      mind = HUGE
c --- go along two fibers
      do i = 1,nDiv
         do j = 1,nDiv
            mind = min(mind,CalDist2P(x1(1,i),x2(1,j)))
         enddo
      enddo  

      FibMinDist = mind
      
      return
      end
