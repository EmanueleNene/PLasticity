      integer  D_NumDof,D_DofList,D_NumEl,D_NumNd,
     & D_BacLen,D_FwdLen,D_ULen,D_NumCp,D_LenCp

      common /dof_info/ D_NumDof,D_DofList(7),D_NumEl,D_NumNd,
     & D_BacLen,D_FwdLen,D_ULen,D_NumCp,D_LenCp

c	D_NumDof       - the number of currently active dof per node
c	D_DofList(7)   - the list of currently active dof
c	D_BacLen       - length of M_BackList 
c	D_FwdLen       - length of M_FrwdList  (maximum node number)
c	D_ULen         - length of each of the U vectors
c	D_NumEl        - number of elements in M_ElemOrd vector
c	D_NumNd        - number of nodes in the model
c	D_NumCp        - number of coupled set
c	D_LenCp        - maximum lenght of the coupling record