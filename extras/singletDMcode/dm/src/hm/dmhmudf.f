****************************************************************
*** Dark matter halo velocity profile.                       ***
*** This routine gives back u*DF(u) in units of (km/s)^(-2)  ***
***                                                          ***
***    u is the modulus of \vec{u} = \vec{v} - \vec{v}_{MY}  ***
***    with \vec{v} the 3-d velocity of a WIMP in the        ***
***    galactic frame, and \vec{v}_{MY} the projection on    ***
***    the frame you are considering                         ***
***                                                          ***
***    DF(u) = int dOmega DF(\vec{v}), where                 ***
***    DF(\vec{v}) is the halo local velocity distribution   ***
***    function in the galactic frame                        ***
***                                                          ***
*** Note: u*DF(u) is the same as f(u)/u, where f(u) is the   ***
*** one-dimensional distribution function as defined in e.g. ***
*** Gould, ApJ 321 (1987) 571.                               ***
***                                                          ***
*** Note: it is also the same as the one-dimensional         ***
*** distribution function g(u) as defined in, e.g.,          ***
*** Ullio & Kamionkowski, JHEP ....                          ***
***                                                          ***
*** dshmuDF is normalized such that                          ***
*** int_0^\infty u*dshmuDF du = int_0^\infty u^2 DF(u) du =  ***
*** int_0^\infty f(u) du = 1                                 ***
***                                                          ***
*** Input:  u = Speed in km/s                                ***
*** Output: u*DF(u) in (km/s)^(-2)                           ***
***                                                          ***
*** Calls other routines depending of choice of velocity     ***
*** distribution function (as set by veldf in the common     ***
*** blocks in dshmcom.h)                                     ***
*** Author: Joakim Edsjo                                     ***
*** Date: 2004-01-29                                         ***
****************************************************************

      real*8 function dshmuDF(u)
      implicit none
      include 'dshmcom.h'
      real*8 u,dshmuDFgauss,dshmuDFiso,dshmudfnum,dshmudfnumc,
     &  dshmudftab

      if (veldf.eq.'gauss'.or.veldf.eq.'iso') then
        dshmuDF=dshmuDFiso(u)       ! Piero's routine
c        dshmuDF=dshmuDFgauss(u)     ! Joakim's routine, identical results

      elseif (veldf.eq.'num') then
        dshmuDF=dshmuDFnum(u) 

      elseif (veldf.eq.'numc') then
        dshmuDF=dshmuDFnumc(u) 

      elseif (veldf.eq.'user') then
        dshmuDF=dshmudftab(u,1)  ! load file specified by udffile is dshmcom.h
                                 ! set udfload to .true. to force reload 

      else
        write(*,*) 'in function dshmuDF unrecognized veldftype flag'
        write(*,*) 'veldf = ',veldf
        write(*,*) 'program stopped'
        stop

      endif

      return
      end





