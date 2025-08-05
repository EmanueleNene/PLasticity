      subroutine CalBeamSect(ireal, keyopt, secDataI, secDataR,
     &                       Summary, pResults, nSecRstNodes, SaveVars)
c --- primary function: compute the section integrated quantities
c
c  Output variables:
c
c   secDataI     (int,ar(*),out) - computed integer section data
c   secDataR     (dp,ar(*),out)  - computed real section data
c   NumBNodes    (int,sc,out)    - number of nodes on the boundary of the cross section
c   numNodes     (int,sc,in/out) - number of nodes in the cross section
c   numElems     (int,sc,in/out) - number of elements in the cross section
c   numSectPts   (int,sc,out)    - number of intg. points/section element
c   Summary      (dp,ar(*),out)  - computed section data summary
c   pResults     (dp,ar(*),out)  - solution from CalSect_1 for input section
c
c----------------------------------------------------------------------
c      Beam Section Record
c----------------------------------------------------------------------
c
c      1        section ID
c      2        sectionType
c      3        sectionSubType
c               1 = circular beam
c               2 = circular tube
c      4        variable NumBNodes
c               = start location of flags (0/1) that identify exterior nodes of cross section
c      5        number of Nodes in Section (NumNodes)
c      6        number of Elements in Section (NumElems)
c      7        number of integraion points (NumSectPts)
c                                            NumSectPts= 4*NumElems
c      8 ..27   Cross section Summary
c               8      area
c               9      centroid - y
c               10     centroid - z
c               11     first moment of inertia about y  Iy
c               12     first moment of inertia about z  Iz
c               13     second moment of inerta about y  Iyy
c               14     product moment of inertia        Iyz
c               15     second moment of inertia about z Izz
c               16     torsion Constant                 J
c               17     polar moment of inertia          Ipp
c               18     int(sectoral area)
c               19     int(sectoral area*y)
c               20     int(sectoral area*z)
c               21     int(sectoral area**2) warping coefficient
c               22     shear Center y
c               23     shear Center z
c               24     shear Correction Factor yy
c               25     shear Correction Factor yz
c               26     shear Correction Factor zz
c     & .. 7+NumElems*9+NumNodes integer words
c               1      number of nodes/element == 9
c               2      subsection type
c               3      mesh size in y direction
c               4      mesh size in z direction
c               5      number of nodes in mesh
c               6      number of elements in mesh
c               7      offset to connectivity data
c               NumElems*9  connectivity data
c               NumNodes    boundary flags
c     &  .. 2*NumNodes+NumElems d.p words, real data
c	          NumNodes  y coordinates
c               NumNodes  z coordinates
c               NumElems  matId for each element
c     &  .. 11*NumSectPts d.p words
c               at each section point -
c               1      dA
c               2      function phi
c               3      d(phi)/dy  torsion
c               4      d(phi)/dz  torsion
c               5      y coord
c               6      z coord
c               7      Er*dA
c               8      d(phi)/dy  shear force in z
c               9      d(phi)/dz  shear force in z
c               10     d(phi)/dy  shear force in y
c               11     d(phi)/dz  shear force in y
c     &  .. 19 words
c               1      NumSectRstNodes
c               2      offset y
c               3      offset z
c               4      position flag (1=centroid,2=shear,3=origin)
c               ..     15 double precision input data
c     &  .. 2 dp words
c               2      section name
c
c     Total = 55+NumElems+9*NumElems+3*NumNodes+11*NumSectPts
c
c
#include "impcom.inc"
#include "solu_info.inc"
#include "sec_info.inc"
c
      integer SectionType, numSectPts, NNODI, rlgetn, np_x, np_y,
     &        numNodes, numElems, ireal, j, numEleX, numEleY,
     &        numCircElems, numThruElems, ErrorFlag, loc,
     &        nSecRstNodes, numIntPts, keyopt(12), GetSecChars
c
      double precision secDataR(*),pResults(B_NumResPts,*), HALF,
     &                 rInner, rOuter, Summary(*), rvr(4), thick,
     &                 SaveVars(*), breadth, temp(8), height
c
      parameter (numCircElems=8, numThruElems=2, NNODI=9, numIntPts = 4)

      parameter (HALF=0.5d0)
      
      integer secDataI(*)

      integer uElemConn, NumBNodes
      parameter (uElemConn = 5) ! 4 items ahead
c

#if defined(DEBUG)
      call SubBeg('CalBeamSect')
#endif

      ErrorFlag = 0

      j = rlgetn(ireal,1,4,rvr(1))

      SectionType = nint(abs(rvr(1)))
      rOuter = abs(rvr(2))
      rInner = abs(rvr(3))
      if (rInner.eq.0.0d0 .and. SectionType.ne.3) then
          if (SectionType.eq.1.or.SectionType.eq.2) then
              rInner = rOuter
          elseif (SectionType.eq.4) then
              rInner = rOuter*1.d-3
          endif    
      endif

      if (SectionType.eq.1) then
          numEleX = 2
          numEleY = 2
          breadth = rOuter
          height = rInner
          np_x = 2*numEleX+1
          np_y = 2*numEleY+1
          numElems = numEleX * numEleY
          temp(1) = -breadth*HALF
          temp(2) = -height*HALF
          temp(3) =  breadth*HALF
          temp(4) = -height*HALF
          temp(5) =  breadth*HALF
          temp(6) =  height*HALF
          temp(7) = -breadth*HALF
          temp(8) =  height*HALF
          call BeamRect(temp(1), np_x, np_y, secDataR(1),
     &                  secDataI(7+uElemConn), 
     &                  numNodes, numElems)
      elseif (SectionType.eq.2)then
          thick = abs(rvr(4))
          if (thick.eq.0.0d0) then
              thick = (HALF-1.d-3)*min(rOuter,rInner)
          endif
c --- hollow rectangular (does not work in v12 and above)
      
c NumBNodes       {key = 4}
c NumNodes        {key = 5}
c NumElems        {key = 6}
c NumSectPts      {key = 7}
c NumSectRstNodes (key = 8)

c --- v13
c 4         113
c 5          72
c 6          12 (ok)
c 7          48 (ok)
c 8          24 (ok)
c --- v11
c 4         149
c 5          81
c 6          16
c 7          64
c 8          25

c --- we have to mesh the cross section to eliminate pseudo areas 

          temp(1) = rOuter
          temp(2) = rInner
          temp(3) = thick
          temp(4) = thick
          temp(5) = thick
          temp(6) = thick

          call BeamHollRect (temp(1), temp(2), temp(3), 
     &                       temp(4), temp(5), temp(6),
     &                       secDataR(1), secDataI(7+uElemConn),
     &                       numNodes, numElems)

      elseif (SectionType.eq.3) then
          call BeamSolidCirc (numCircElems, numThruElems, rOuter,
     &                        secDataI(7+uElemConn),secDataR(1), 
     &                        numNodes, numElems)
      elseif (SectionType.eq.4) then
          call BeamPipe (numCircElems, rInner, rOuter,
     &                   secDataR(1), secDataI(7+uElemConn),
     &                   numNodes, numElems)
      endif
      
      NumBNodes = uElemConn+NNODI*numElems
      call GetBoundNodes(secDataI(7+uElemConn), numNodes,
     &                   numElems, secDataI(7+NumBNodes))

      if (B_PtrRealN.ne.PTRFTNNULL .and. B_RealN(ireal).gt.0) then
          numSectPts = GetSecChars (ireal,1)
          nSecRstNodes = GetSecChars (ireal,2)
          loc = B_RealN(ireal)
          call v_move (B_SectData(loc),Summary(1),20)
          call v_move (B_SectData(loc+20),pResults(1,1),
     &                 B_NumResPts*NumSectPts)
      elseif (keyopt(7).eq.1 .and. A_KeyFirstPass.eq.0) then
c          numSectPts = numElems * numIntPts
          numSectPts = GetSecChars (ireal,1)
          nSecRstNodes = GetSecChars (ireal,2)
          call v_move (SaveVars(1),Summary(1),20)
          call v_move (SaveVars(21),pResults(1,1),
     &                 B_NumResPts*NumSectPts)
      else
          call CalSect_11 (secDataR(1), NNODI, numElems, numNodes,
     &                     secDataI(7+uElemConn), numSectPts,
     &                     pResults(1,1), Summary(1), 
     &                     ErrorFlag, nSecRstNodes)
          if (keyopt(7).eq.1) then
c --- save pResults(B_NumResPts=11,NumSectPts=numElems*numIntPts), Summary(20)
              call v_move (Summary(1),SaveVars(1),20)
              call v_move (pResults(1,1),SaveVars(21),
     &                     B_NumResPts*NumSectPts)
          endif
      endif

      secDataI(1) = 1
      secDataI(2) = 1
      secDataI(3) = SectionType
      secDataI(4) = NumBNodes
      secDataI(5) = NumNodes
      secDataI(6) = NumElems
      secDataI(7) = NumSectPts
      secDataI(8) = NNODI
      secDataI(9) = SectionType
      secDataI(10) = NumNodes
      secDataI(11) = NumElems

#if defined(DEBUG)
      call SubEnd('CalBeamSect')
#endif

      return
      end
