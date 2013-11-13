***********************************************************************
*** function dshapbyieldf gives the distributions of antiprotons for    ***
*** basic annihilation channels chi=1-8. parameterizations to the   ***
*** distributions are used.                                         ***
*** input: mx - neutralino mass (gev)                               ***
***        tp - antiproton kinetic energy (gev)                     ***
***        chi - annihilation channel, 1-8, (short version)         ***
*** output: differential distribution, p-bar gev^-1 annihilation^-1 ***
*** author: joakim edsjo, edsjo@physto.se                           ***
*** date: 1998-10-27                                                ***
***********************************************************************

      real*8 function dshapbyieldf(mx,tp,chi)
      implicit none

      real*8 mx,tp,x,p(6),yield
      real*8 a(0:3,1:4,1:8)
      integer chi,i
      data a/
     &  1.79,1.40,0.0,0.0,                       ! p1 c c-bar
     &  3.12,0.04,0.0,0.0,                       ! p2 c c-bar
     &  -2.22,0.0,0.0,0.0,                       ! p3 c c-bar
     &  -0.39,-0.076,0.0,0.0,                    ! p4 c c-bar
     &  1.75,1.40,0.0,0.0,                       ! p1 b b-bar
     &  1.54,0.11,0.0,0.0,                       ! p2 b b-bar
     &  -2.22,0.0,0.0,0.0,                       ! p3 b b-bar
     &  -0.31,-0.052,0.0,0.0,                    ! p4 b b-bar
     &  1.35,1.45,0.0,0.0,                       ! p1 t t-bar
     &  1.18,0.15,0.0,0.0,                       ! p2 t t-bar
     &  -2.22,0.0,0.0,0.0,                       ! p3 t t-bar
     &  -0.21,0.0,0.0,0.0,                       ! p4 t t-bar
     &  16*0.0,                                  ! p1-p4, tau+ tau-
     &  306,0.28,7.2e-4,2.25,                    ! p1 w+ w-
     &  2.32,0.05,0.0,0.0,                       ! p2 w+ w-
     &  -8.5,-0.31,0.0,0.0,                      ! p3 w+ w-
     &  -0.39,-0.17,-2.0e-2,0.23,                ! p4 w+ w-
     &  480.,0.26,9.6e-4,2.27,                   ! p1 z0 z0
     &  2.17,0.05,0.0,0.0,                       ! p2 z0 z0
     &  -8.5,-0.31,0.0,0.0,                      ! p3 z0 z0
     &  -0.33,-0.075,-1.5e-4,0.71,               ! p4 z0 z0
     &  16*0.0,                                  ! p1-p4, mu+ mu-
     &  2.33,1.49,0.0,0.0,                       ! p1 glue glue
     &  3.85,0.06,0.0,0.0,                       ! p2 glue glue
     &  -2.17,0.0,0.0,0.0,                       ! p3 glue glue
     &  -0.312,-0.053,0.0,0.0/                   ! p4 glue glue
      save a

c...start-up
c      x=tp
      x=tp/mx  ! if input is tp
      do i=1,4
         if (a(0,i,chi).ne.0.0d0) then
           p(i)=1.0/(a(0,i,chi)*mx**(a(1,i,chi))+
     &       a(2,i,chi)*mx**(a(3,i,chi)))
         else
           p(i)=0.0d0
         endif
      enddo

      if (p(1).ne.0.0d0.or.p(2).ne.0.0d0) then
        yield=1.0d0/(p(1)*x**(p(3))+p(2)*abs(log10(x))**p(4))
      else
        yield=0.0d0
      endif

c...the yield is now dn/dx

c...change from dn/dx to dn/dt

      if (x.lt.1.0d0) then
c        dshapbyieldf=yield       ! if output is dn/dx
        dshapbyieldf=yield/mx    ! if output is dn/dt
      else
        dshapbyieldf=0.0d0
      endif

      end


