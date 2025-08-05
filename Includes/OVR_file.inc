c --- description of element save data file
      character*8  O_Name 
      parameter   (O_Name='ovr     ')

      integer*8         O_FilePos
      integer           O_BlockNum, O_Unit, O_SaveKey
      common /OVR_file/ O_FilePos,  O_BlockNum, O_Unit, O_SaveKey

c	O_SaveKey  - key to save backup file for stepback/restart
c	O_FilePos  - file position on file OVR
c	O_BlockNum - block number for file OVR
c	O_Unit     - file unit for file OVR
