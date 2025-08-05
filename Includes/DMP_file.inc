c --- description of element save data file
      character*8  Z_Name 
      parameter   (Z_Name='dmp     ')

      integer*8       Z_FilePos
      integer         Z_BlockNum, Z_Unit, Z_SaveKey
      common /DMP_file/ Z_FilePos, Z_BlockNum, Z_Unit, Z_SaveKey

c	Z_SaveKey  - key to save backup file for stepback/restart
c	Z_FilePos  - file position on file DMP
c	Z_BlockNum - block number for file DMP
c	Z_Unit     - file unit for file DMP
