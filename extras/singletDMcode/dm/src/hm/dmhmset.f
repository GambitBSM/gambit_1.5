      subroutine dmhmset(alpha, beta, gamma, rho)
****************************************************************
*** subroutine dshmset:                                      ***
*** initialize the density profile and or the small clump    ***
*** probability distribution                                 ***
*** type of halo:                                            ***
***   hclumpy=1 smooth, hclumpy=2 clumpy                     ***
***                                                          ***
*** a few sample cases are given; specified values of the    ***
*** local halo density 'rho0' and of the length scale        ***
*** parameter 'a' should be considered just indicative       ***
***                                                          ***
*** author: piero ullio (piero@tapir.caltech.edu)            ***
*** date: 00-07-13                                           ***
*** small modif: paolo gondolo 00-07-19                      ***
*** mod: 03-11-19 je, 04-01-13 pu                            ***
****************************************************************
      implicit none
      include 'dshmcom.h'
      include 'dssusy.h'
    
      real*8 alpha, beta, gamma, rho

c...Set default flags
      udfload=.true.  ! load udffile (if chosen by veldf='user') on next
                      ! call to dshmudf.f
      udfearthload=.true. ! load udfearthfile (if chosen by veldfearth='user')
                          ! on next call to dshmudfearth.f

      isodfload=.true.   ! load isodf file on next call to dshmisotrnum.f

c...Three-dimensional velocity dispersion of the WIMPs in the halo
c...Note that in a simple isothermal sphere, the circular speed
c...(approx. solar speed v_sun) is sqrt(2/3)*vd_3d
      vd_3d = 270.0d0     ! WIMP 3D velocity dispersion, km/s
      v_sun = 220.0d0     ! circular velocity at the Sun location, km/s
      v_obs = v_sun       ! Velocity of observer
      vgalesc = 600.d0          ! galactic escape speed in km/s
c...Observer speed w.r.t. the halo, i.e. including Sun + Earth speed,
c...yearly average,
      vobs = 264.d0       ! observer speed w.r.t. the halo in km/s
      v_earth = 29.78d0    ! Earth speed in the solar system in km/s

c a few options for the density profile:

c ... standard parametrization 

      hclumpy = 1
      r_0 = 8.d0          ! sun galactocentric distance (kpc)
      rho0 = rho          ! local halo density (gev cm^-3)
      haloshape = 'spherical' !here you choose between spherical and axisymm.
      halotype = 'albega' !this picks the alpha-beta-gamma profile, set by:
      alphah = alpha        
      betah = beta
      gammah = gamma
      ah = 20.d0          ! length scale (kpc)
      Rref = r_0
      rhoref = rho0       ! normalizing the profile to the local halo density
      rhcut = 1.d-5       ! cut radius (kpc) 
      haloid=' '          ! tag used to identify pbar and dbar files
      veldf='gauss'       ! tag used to identify velocity profile
      veldfearth='sdbest' ! tag to identify earth vel. profile


c...unrescaled local density
      rhox = rho0

      return
      end
