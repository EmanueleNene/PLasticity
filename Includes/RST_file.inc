      character*8  R_Name 
      parameter  (R_Name='rst     ')
      integer    R_OutMask,  R_OutCode
      integer*8  R_FilePos,  R_ResIndxLoc, R_TimIndxLoc, R_LdStepLoc
      integer    R_BlockNum, R_Unit, R_CurSet, R_MaxSet, R_MaxNR,
     &           R_NRes, R_Stop 
      parameter (R_MaxSet = 10000) ! safety reasons
      common /rst_file/ 
     &	R_OutMask,R_OutCode,
     &	R_FilePos, R_ResIndxLoc, R_TimIndxLoc, R_LdStepLoc,
     &	R_BlockNum, R_Unit, R_CurSet, R_MaxNR, R_NRes, R_Stop


c  R_OutMask  - global mask: bit set to 1 if any item to be output
c  R_OutCode  - integer indicating the item specified
c            bit+1 - individual item
c            101   - all result set
c            102   - basic result set

c  R_FilePos    - file position on rst file
c  R_ResIndxLoc - location of index start -1
c  R_LdStepLoc  - location of loadstep/iteration -1
c  R_TimIndxLoc - location of time index start -1

c  R_BlockNum   - block number for res file (usually 6)
c  R_Unit       - file unit for file rst (12)
c  R_CurSet     - current result file load set (one just read/write)
c  R_MaxSet     - maximum number of load sets allowed on file


c  see up-to-date ansys doc for discription of the format
c  we use revison 11.0
c  ALL LOCAL



