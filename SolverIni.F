      subroutine SolverIni (kerr)
c --- primary function:   subroutine to allocate solver memory space
c --- secondary function: open files
c
c  output arguments:
c     kerr     (int,sc,out)      - error flag:  0 = OK, 1= Error
c
#include "impcom.inc"
#include "constants.inc"
#include "file_sys.inc"
#include "solu_info.inc"
#include "SVR_file.inc"
#include "RST_file.inc"
c
      integer   Bin_Set

      integer   i,kerr

      character*(B_PathLen) pname


#if defined(DEBUG)
      call SubBeg('SolverIni')
#endif

      kerr = 0

c --- create solution database objects
      call iniAnsObj(0)

c --- allocate memory needed for solution
      call AllocSolMem(kerr)

c --- open SVR for solution
      call GetPath (S_Name,pname)
      S_BlockNum = 1
      S_Unit = 3
      i = Bin_Set(S_BlockNum,S_Unit,W_Write,B_NumOfPages,pname)

      if (i .ne. 0) then
          call erhandler ('SolverIni',5000,3,
     &                    'Cannot open SVR file %c.',
     &                     0.0d0,pname)
          S_Unit = 0
          kerr = 999
          goto 999
      endif

      if (A_KeyFirstPass.eq.1) then
          S_FilePos = 0
          call Bin_Head (S_BlockNum,S_Unit,S_FilePos)
      endif

c --- open results file (set read-write key)
      R_BlockNum = 3
      R_Unit = 12
      R_CurSet = 0

      call GetPath(R_Name,pname)
c --- only 4-8 pages for RST file should be enough
      i = Bin_Set(R_BlockNum,R_Unit,W_Write,B_NumOfPages,pname)
      
#if !defined (IEEE_SYS)
c --- store RST file always in IEEE format for ANSYS,
c     other files are saved in native formats
      B_PageFormat(R_BlockNum) = 1
#endif

      if (i .ne. 0) then
          call erhandler ('SolverIni',5000,3,
     &                    'Cannot open the results file %c.',
     &                     0.0d0,pname)
          R_Unit = 0
          kerr = 999
          goto 999
      endif

  999 continue

#if defined(DEBUG)
      call SubEnd('SolverIni')
#endif

      return
      end


