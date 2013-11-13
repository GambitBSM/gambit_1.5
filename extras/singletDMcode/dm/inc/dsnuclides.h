*         -*- mode: fortran -*-
      integer nnucld,inucld
      parameter (nnucld=3179)
      character*6 nucldsym(nnucld)
      integer nuclda(nnucld),nucldz(nnucld),nucldn(nnucld)
      real*8 nucldm(nnucld),nucldnatab(nnucld),nucldj(nnucld),
     &     nucldmu(nnucld),nucldam(nnucld)
      logical nucldstable(nnucld)
      common /ddnucld/ nucldm,nucldnatab,nucldj,nucldmu,nucldam,
     &     nuclda,nucldz,nucldn,inucld,nucldstable,nucldsym
      save /ddnucld/

