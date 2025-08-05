      character*8  S_Name
      parameter   (S_Name='svr     ')
      integer      S_nIndx
      parameter   (S_nIndx=3)
      integer*8    S_FilePos, S_FileLengL, S_RecStart
      integer      S_BlockNum, S_Unit
	  logical      S_InCore
      common /SVR_file/ S_FilePos, S_BlockNum, S_Unit, S_InCore, 
     &                  S_FileLengL, S_RecStart

c	S_FilePos  - file position on file SVR
c	S_BlockNum - block number for file SVR
c	S_Unit     - file unit for file SVR
