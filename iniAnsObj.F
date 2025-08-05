      subroutine iniAnsObj (flag)
c --- primary function:  create/delete ANSYS in-core database objects used only during solution
c
c  input arguments:
c     flag     (int,sc,in)  - flag -> 0 create objects
c                                     1 delete objects
c
#include "impcom.inc"
#include "constants.inc"
c
c --- define common block values for ANSYS
      integer  
     & ceflag, elmflag, nxrefflag, cpflag, etypeflag, mpflag, 
     & nodeflag, realflag, csysflag, dispflag, forceflag,
     & sectflag, tbflag, eslflag, eblflag, nblflag

      common /sdbcom/
     & ceflag, elmflag, nxrefflag, cpflag, etypeflag, mpflag, 
     & nodeflag, realflag, csysflag, dispflag, forceflag,
     & sectflag, tbflag, eslflag, eblflag, nblflag 
c
      integer
     &       initAnsElm, deleteAnsElm,
     &       initAnsNxref, deleteAnsNxref,
     &       initAnsCp, deleteAnsCp,
     &       initAnsEType, deleteAnsEType,
     &       initAnsMp, deleteAnsMp,
     &       initAnsNode, deleteAnsNode,
     &       initAnsReal, deleteAnsReal,
     &       initAnsCsys, deleteAnsCsys,
     &       initAnsDisp, deleteAnsDisp,
     &       initAnsForce, deleteAnsForce,
     &       initAnsTb, deleteAnsTb,
     &       initAnsElmSurfLoad, deleteAnsElmSurfLoad,
     &       initAnsElmBodyLoad, deleteAnsElmBodyLoad,
     &       initAnsNodeBodyLoad,deleteAnsNodeBodyLoad,
     &       initPcgMemManager

      integer  elmiqr, elmget
c
      integer  flag, kerr, elem, numNodes, numElems, i, 
     &         elData(EL_DIM), nodes(W_NdMax), nLen
      
#if defined(DEBUG)
      call SubBeg('iniAnsObj')
#endif

      if (flag .eq. 0) goto 100
      if (flag .eq. 1) goto 500
      goto 999

c --- set flag indicating that objects have not been created
 100  nLen = 0
c --- initialize memory manager (always required for object init)
      kerr = initPcgMemManager(nLen)

      call setSoluFlag(1,-1)

c --- allocate temporary workspace to be used by the objects
      call solveDBwork(0)

c --- create Element type object
      kerr = initAnsEType()
      
c --- create Element object
      kerr = initAnsElm()

c --- create Node object
      kerr = initAnsNode()

c --- create Node-to-Element cross-reference object
      kerr = initAnsNxref()

c --- create CP object
      kerr = initAnsCp()

c --- create Linear material property object
      kerr = initAnsMp()

c --- create Nonlinear material property object
      kerr = initAnsTb()

c --- create Real constant object
      kerr = initAnsReal()

c --- create Coordinate system object
      kerr = initAnsCsys()

c --- create Displacement object
      kerr = initAnsDisp()

c --- create Force object
      kerr = initAnsForce()

c --- create Element surface load object
      kerr = initAnsElmSurfLoad()

c --- create Element body load object
      kerr = initAnsElmBodyLoad()

c --- create Nodal body load object
      kerr = initAnsNodeBodyLoad()

c --- deallocate temporary workspace used by the objects
      call solveDBwork(1)

c --- set flag indicating that objects have been created
      call setSoluFlag(1,1)

      goto 999

 500  continue

c --- UPDATE ELEMENT DATABASE

c --- get the number of defined elements
      numElems = elmiqr(0,DB_MAXDEFINED)
c --- loop over the number of defined elements looking for modified elements
      do elem = 1,numElems
c --- get modified element from solution database object
         numNodes = elmget(elem,elData(1),nodes(1))
         if (numNodes.le.0) cycle
c --- toggle solution database off so this element goes into database
         elmflag = 0
c --- put element into database
         call elmput(elem,elData(1),numNodes,nodes(1))
c --- toggle solution database back on for the next element
         elmflag = 1
      enddo

c --- delete Nodal body load object
      kerr = deleteAnsNodeBodyLoad()
      nblflag = 0

c --- delete Element body load object
      kerr = deleteAnsElmBodyLoad()
      eblflag = 0

c --- delete Element surface load object
      kerr = deleteAnsElmSurfLoad()
      eslflag = 0

c --- delete Force object
      kerr = deleteAnsForce()
      forceflag = 0

c --- delete Displacement object
      kerr = deleteAnsDisp()
      dispflag = 0

c --- delete Coordinate system object
      kerr = deleteAnsCsys()
      csysflag = 0
      
c --- delete Real constant object
      kerr= deleteAnsReal()
      realflag = 0

c --- delete Nonlinear material property object
      kerr = deleteAnsTb()
      tbflag = 0

c --- delete Linear material property object
      kerr = deleteAnsMp()
      mpflag = 0

c --- delete CP object
      kerr = deleteAnsCp()
      cpflag = 0

c --- delete Node-to-Element cross-reference object
      kerr = deleteAnsNxref()
      nxrefflag = 0

c --- delete Node object
      kerr = deleteAnsNode()
      nodeflag = 0

c --- delete Element object
      kerr = deleteAnsElm()
      elmflag = 0

c --- delete Element type object
      kerr = deleteAnsEType()
      etypeflag = 0

c --- set flag indicating that objects have been deleted
      call setSoluFlag(1,0)

c --- cleanup solver memory
      call deletePcgMemManager()
      
  999 continue

#if defined(DEBUG)
      call SubEnd('iniAnsObj')
#endif

      return
      end
