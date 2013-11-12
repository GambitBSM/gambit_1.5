subroutine darkmatter(Flag, InP, InN, DMInOut)

  !----------------------------------------------
  ! interfacing routine DM with MCMC  
  ! authors: Roberto Ruiz, Pat Scott, 
  !          Christoph Weniger
  !---------------------------------------------

  use parameters

  implicit none  

  include'dssusy-f90.h'  
  include'dsidtag-f90.h' 
  include 'dshacom-f90.h'
 
  Type(FLAGS), INTENT(IN) :: Flag
  Type(Input_params), INTENT(IN) :: InP
  Type(Nuisance_params), INTENT(IN) :: InN
  Type(DM), INTENT(INOUT) :: DMInOut

  real (8) :: dmhmj, dmhmjave, dmhrgacontdiff, dmhrgacont, dmhrgacsusy, dmhaloyield ! ID functions
  integer  :: istat                         ! err or war flags
  real (8), save :: jpsi_gamma = 0.d0
  logical, save :: first = .true.

  character(LEN=35), save :: phigam_units, dphigam_units

  real(8) :: dsabsq ! useful function: sqrt(abs(x))
  integer :: i, j, k, l
  real (8) :: step, egam, em, delta, err = 0.d0
  ! Fermi declarations
  integer, parameter :: njsamples = 150
  integer, parameter :: jsamples_power = 4
  double precision :: cospsi0max, Z_gamma_energy
  real (8) :: cospsi0min, jsamples(3,njsamples), gagaRate, gazRate, Emin, Emax
  real (8) :: centralPix(2), centralVal, leftVal, cornerVal, topVal
  real (8), allocatable, save :: GC_js(:,:), GC_susyBits(:)

  !gamma fluxes background 
  real (8) :: phig_back


  DMInOut%PrintOutParams=.false. !reseting


  ! The amount of output onto standard output can be controlled by the
  ! variable prtlevel at deep level and debug for a rough testing. 
  ! Setting prtlevel=0 suppresses all messages,
  ! except fatal errors that terminate the program. Setting prtlevel=1
  ! displays informational messages. Setting prtlevel to higher values
  ! displays more and more information. 
  ! Setting debug = .true. make a diagnosis of what is doing 
  ! at the level of the call to the routines which compute 
  ! the different observables.

  prtlevel = 5 !to get avoid tags
  debug = .false. 

  ! This call initializes the DarkSUSY package. This call should
  ! be the first call in any program using DarkSUSY. This call initial
  ! some global variables and calls various other modules to initializ
  ! them with their default values. Check out src/ini/dsinit.f if you
  ! want to see what it does.

  call dminit(InP, InN, DMInOut) 

  ! All the SUSY spectrum passed to DarkSUSY is computed by Softsusy. 
  ! The  supersymmetric parameters are: the gaugino masses m1, m2, m3; the
  ! Higgs pseudoscalar mass ma; the ratio of Higgs vacuum expectation
  ! values tanbe; the square of the squark and slepton mass parameters
  ! mass2q(i), mass2u(i), mass2d(i), mass2l(i), mass2e(i); the
  ! trilinear soft parameters asofte(i), asoftu(i), asoftd(i). Here
  ! i=1,2,3 is a generation index. They are defined at the 
  ! weak scale. Physical parameters are passed too.    

  !For serious debugging
  ! writes the spectrum in dstest1.tmp
  !        open (unit=30,file='dsspect.tmp')
  !        write(30,*) ' '
  !        call dswspectrum(30)
  !        close(30)

  ! writes the vertex in dsvtx.tmp
  !        open (unit=40,file='dsvtx.tmp')
  !        call dswvertx(40)
  !        close(40)

  ! Now we calculate the ID parameters

  if(Flag%ID_predict) then

    if (debug) write ( * , * ) '    >>> Calculating ID parameters...' 

    if (first) then

      ! Now we initialize the spherically symmetric density profile
      ! provided (see dm/src/hm/dmhmset.f for details).
      !if they are fixed them just call it once
      !for the actual profiles, see dm/src/hm/dmhmsphrho.f
      call dmhmset(InP%alpha, InP%beta, InP%gamma, InP%rho)   

      if(DMInOut%ID_in%gammas%delta_gamma == 0.d0) then
        jpsi_gamma = dmhmj(DMInOut%ID_in%gammas%cospsi0)
        phigam_units = ' (cm^-2 s^-1 sr^-1)'
        dphigam_units = ' (cm^-2 s^-1 MeV^-1 sr^-1)' 
      else
        jpsi_gamma = dmhmjave(DMInOut%ID_in%gammas%cospsi0,DMInOut%ID_in%gammas%delta_gamma)*DMInOut%ID_in%gammas%delta_gamma
        phigam_units = ' (cm^-2 s^-1)'
        dphigam_units = ' (cm^-2 s^-1 MeV^-1)'  
      endif

      if (debug) write(* ,'(A,G10.4)')  '        >>> jpsi for gamma rays: ',jpsi_gamma

    endif

    !if halo params fixed, it is set to the value in the previous interation 
    DMInOut%AddOut%J = jpsi_gamma         

    if (debug) write(*,*) 'setting  j to =',  jpsi_gamma ,  DMInOut%AddOut%J     

    if (Flag%ID_Flags_gamma%GC_region) then

      if(Flag%use_data == 1) then   !fermi data

        Z_gamma_energy = InP%mx * (1.d0 - mass(kz)**2.d0 / (4.d0 * InP%mx**2.d0))
        call dmhrgalsusy(gagaRate, gazRate)

        allocate(GC_susyBits(DMInOut%ID_in%gammas%GCDims(3)))

        if(first) then
        !computes the J factor 

          allocate(GC_js(DMInOut%ID_in%gammas%GCDims(1),DMInOut%ID_in%gammas%GCDims(2)))
      
          GC_js = 0d0
 
          ! Make sure the central region's J factor is correctly integrated over the central peak.
          centralPix = (/DMInOut%ID_in%gammas%GCDims(1)/2,DMInOut%ID_in%gammas%GCDims(2)/2/)
          centralVal = dmhmjave(1.d0,DMInOut%ID_in%gammas%GC_pixArea*4.d0)
      
          leftVal = dmhmjave(DMInOut%ID_in%gammas%GC_angSep(centralPix(1)-1,centralPix(2)),DMInOut%ID_in%gammas%GC_pixArea)
          cornerVal = dmhmjave(DMInOut%ID_in%gammas%GC_angSep(centralPix(1)-1,centralPix(2)-1),DMInOut%ID_in%gammas%GC_pixArea)
          topVal = dmhmjave(DMInOut%ID_in%gammas%GC_angSep(centralPix(1),centralPix(2)-1),DMInOut%ID_in%gammas%GC_pixArea)
          GC_js(centralPix(1):centralPix(1)+1,centralPix(2):centralPix(2)+1) = centralVal
          GC_js(centralPix(1)-1:centralPix(1)+2:3,centralPix(2):centralPix(2)+1) = leftVal
          GC_js(centralPix(1)-1:centralPix(1)+2:3,centralPix(2)-1:centralPix(2)+2:3) = cornerVal
          GC_js(centralPix(1):centralPix(1)+1,centralPix(2)-1:centralPix(2)+2:3) = topVal

     
          ! Fill in the surrounding areas
          cospsi0max = minval(DMInOut%ID_in%gammas%GC_angSep(centralPix(1)-1:centralPix(1)+2,centralPix(2)-1:centralPix(2)+2))
          cospsi0min = minval(DMInOut%ID_in%gammas%GC_angSep)
          forall(j=1:njsamples) jsamples(1,j) = cospsi0max - (dble(njsamples-j)/dble(njsamples-1))**jsamples_power * (cospsi0max-cospsi0min)
          do j=1,njsamples
            jsamples(2,j) = dmhmj(jsamples(1,j))
          enddo

          call dmspline(jsamples(1,:),jsamples(2,:),njsamples,0.d0,0.d0,jsamples(3,:)) 

          !fill left and right of the grid ussing interpolation in angles
          do k=1,DMInOut%ID_in%gammas%GCDims(2)
            do j=1,centralPix(1)-2
              call dmsplint(jsamples(1,:),jsamples(2,:),jsamples(3,:),njsamples,DMInOut%ID_in%gammas%GC_angSep(j,k),GC_js(j,k))
            enddo
            do j=centralPix(1)+3,DMInOut%ID_in%gammas%GCDims(1)
              call dmsplint(jsamples(1,:),jsamples(2,:),jsamples(3,:),njsamples,DMInOut%ID_in%gammas%GC_angSep(j,k),GC_js(j,k))
            enddo
          enddo

          !fill lower and upper of the grid
          do j=centralPix(1)-1,centralPix(1)+2
            do k=1,centralPix(2)-2
              call dmsplint(jsamples(1,:),jsamples(2,:),jsamples(3,:),njsamples,DMInOut%ID_in%gammas%GC_angSep(j,k),GC_js(j,k))
            enddo
            do k=centralPix(2)+3,DMInOut%ID_in%gammas%GCDims(2)
              call dmsplint(jsamples(1,:),jsamples(2,:),jsamples(3,:),njsamples,DMInOut%ID_in%gammas%GC_angSep(j,k),GC_js(j,k))
            enddo
          enddo

          GC_js = GC_js

        endif !end first

        Emin = DMInOut%ID_in%gammas%GC_Ebins(1,1) 
        do j=1,DMInOut%ID_in%gammas%GCDims(3)
          Emax = DMInOut%ID_in%gammas%GC_Ebins(j,2)

          !Fast continuum/IB spectrum
          !GC_susyBits(j) = dmhrgacdiffsusy(DMInOut%ID_in%gammas%GC_Ebins(j,3),istat)

          !Accurate continuum/IB spectrum - *can be rewritten to halve calling time if necessary*
          GC_susyBits(j) = (dmhrgacsusy(DMInOut%ID_in%gammas%GC_Ebins(j,1),istat) - dmhrgacsusy(DMInOut%ID_in%gammas%GC_Ebins(j,2),istat))/(Emax-Emin)

          !Addition of line spectra
          if (Emin .lt. InP%mx .and. Emax .gt. InP%mx) GC_susyBits(j) = GC_susyBits(j) + gagaRate/(Emax-Emin)
          if (Emin .lt. Z_gamma_Energy .and. Emax .gt. Z_gamma_Energy) GC_susyBits(j) = GC_susyBits(j) + gazRate/(Emax-Emin)

          Emin = Emax

        enddo

        !GC_susyBits(l) -> Particle physics
        !DMInOut%ID_in%gammas%GC_model(j,k,l) -> On input, it should contain the background
        !                                     -> On output, background + SUSY signal * GC_J factor *(1+BF)
        forall(j=1:DMInOut%ID_in%gammas%GCDims(1),k=1:DMInOut%ID_in%gammas%GCDims(2),l=1:DMInOut%ID_in%gammas%GCDims(3)) &
         DMInOut%ID_in%gammas%GC_model(j,k,l) = DMInOut%ID_in%gammas%GC_model(j,k,l) + GC_susyBits(l) * GC_js(j,k) * (1.d0 + DMInOut%ID_in%gammas%GCBF)
        !WARNING!! The second term in the expression above needs to be rescaled if the neutralino is not the only CDM
        !component; the scaling factor is (rhox/rho0)**2, where rho0 is the local dark matter density and rhox is the
        !local neutralino density.

        deallocate(GC_susyBits)

        if(Flag%ID_Flags_gamma%halo_fix) then
          first = .false.
        else 
          deallocate(GC_js)
        endif    

        ! PS - The call to flatlib is easier done in the likelihood code, as it requires a continuously-defined (ie interpolated) theory prediction
        !      and a number of constants which are probably easier not to pass over to this routine.
        !      The same is true for the background, as it is easier to to read in the background files once and for all in the main DMBayes
        !      routines than over and over here in the DMFit interface.

      else

        !old non-Fermi data mock data handling (RT)
        step = log10(DMInOut%ID_in%gammas%efluxes_f/DMInOut%ID_in%gammas%efluxes_i)/DMInOut%ID_in%gammas%nbins
 
        do i = 1, DMInOut%ID_in%gammas%nbins !nbins is now set from the data file in fermi_ini
          egam = DMInOut%ID_in%gammas%efluxes_i * 10.d0**((i-1) * step)

          !all the energies are set to be common 
          DMInOut%ID%gammas%fluxgadiff(i) = dmhrgacontdiff(egam, jpsi_gamma, istat) * 1.d-3 !MeV^-1
          DMInOut%ID%gammas%Ekin(i) = egam * 1.d3 !MeV

          DMInOut%ID%gammas%fluxgadiff(i) = DMInOut%ID%gammas%fluxgadiff(i) + phig_back(egam, InN)

          if(haerr /=0) then
            DMInOut%Err%haerr = 1
            return
          endif
        enddo
        ! for getting fake spectra 
        !do i = 1, DMInOut%ID_in%gammas%nbins
        !  em = DMInOut%ID%gammas%Ekin(i)
        !  delta=(DMInOut%ID%gammas%Ekin(i+1)-DMInOut%ID%gammas%Ekin(i))/2.d0
        !  write (* ,'(3x,G10.4,3x,G10.4,3x,G10.4,3x,G10.4,3x,G10.4,3x,G10.4)') em* 1.d3,(em-delta)* 1.d3,(em+delta)* 1.d3, em**2*DMInOut%ID%gammas%fluxgadiff(i)*1.d3, err, err 

        !  print*,(DMInOut%ID%gammas%Ekin(i+1)*1.d6-DMInOut%ID%gammas%Ekin(i)*1.d6)/2.d0

        !enddo
        !stop

        if (debug) write (* ,'(A,G10.4,A)') '       >>> Old data: Max Gamma ray fluxes ', maxval(DMInOut%ID%gammas%fluxgadiff(1:DMInOut%ID_in%gammas%nbins+1)) , dphigam_units

      endif !ends using Fermi or fake data

    endif !ends GC_region

    ! CW: Calculate dN/dE for dwarf limits
    ! PS: ...and for CTA limits 
    if (Flag%ID_Flags_gamma%dwarfs .or. Flag%ID_Flags_gamma%cta_gc) then

      step = log10(DMInOut%ID_in%gammas%efluxes_f/DMInOut%ID_in%gammas%efluxes_i)/(DMInOut%ID_in%gammas%nbins)

      do i = 1, DMInOut%ID_in%gammas%nbins !nbins is set in the calling routine
        egam = DMInOut%ID_in%gammas%efluxes_i * 10.d0**((i-1) * step)

        !all the energies are set to be common 
        DMInOut%ID%gammas%fluxdiff(i) = dmhaloyield(egam, 152, istat) * 1.d-3 !MeV^-1
        DMInOut%ID%gammas%Ekin(i) = egam * 1.d3 !MeV

        if(haerr /=0) then
          DMInOut%Err%haerr = 1
          return
        endif
      enddo
    endif
 
    if (Flag%ID_Flags_gamma%gac) then

      ! ii) gamma-ray flux with continuum energy spectrum integrated above
      ! some given threshold egath (GeV)
       
      DMInOut%ID%gammas%fluxgac = dmhrgacont(DMInOut%ID_in%gammas%egath,jpsi_gamma,istat) !ph cm^-2 s^-1
       
      if (debug)  write (* ,'(A,F8.3,A,G10.4,A)') '         >>> Gamma ray fluxes with threshold enery', &
       DMInOut%ID_in%gammas%egath, ' GeV: ', DMInOut%ID%gammas%fluxgac, dphigam_units 

    endif !gamma cont
     
    if(debug)  write ( *, *)

  endif !ID flag


  first_idmod = 0 !undecouples model dependent tabulation


end subroutine darkmatter

