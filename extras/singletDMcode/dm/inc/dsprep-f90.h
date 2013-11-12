!         -*- mode: fortran -*-
!######################################################################*
!                       i n c l u d e     f i l e                      *
!######################################################################*

!***********************************************************************
!**                            dsprep.h                              ***
!**         this piece of code is needed as a separate file          ***
!**              the rest of the code 'includes' dsprep.h            ***
!----------------------------------------------------------------------c
!  author: joakim edsjo, 00-08-16
!...new model - common block for new model switches
      logical :: newmodelanwx,newmodelandwdcosnn,newmodelep, &
       newmodelntdk,newmodelsigmav, &
       dsprepcalled
      common /prep/newmodelanwx,newmodelandwdcosnn,newmodelep, &
       newmodelntdk,newmodelsigmav, &
       dsprepcalled 

!....susypar - susy model dependent parameters used by various routines
      real*8 :: mx,sigmav,sigv,wtot
      common/prep2/mx,sigmav,sigv(30),wtot






