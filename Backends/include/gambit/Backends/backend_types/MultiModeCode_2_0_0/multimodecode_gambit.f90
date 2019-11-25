MODULE multimodecode_gambit
!Module that includes the function that mimicks the driver in multimodecode 
!so that the inflationary solvers can be called from a gambit source file 
!in the cosmobit (gambit) module.
  use modpkparams
  use potential
  use background_evolution
  use modpk_utils
  use camb_interface
  use access_modpk, only : evolve, potinit ! , reheat_evolve
  use internals
  use modpk_icsampling
  use modpk_io, only : out_opt
  use modpk_deltaN
  use modpk_observables, only : observables
  use csv_file, only : csv_write
  use modpk_errorhandling, only : raise, run_outcome, assert
  use modpk_rng, only : init_random_seed

  implicit none

  public :: multimodecode_gambit_driver

! The structure of the output gambit reads.
  type :: gambit_inflation_observables
    logical :: check_ic_ok
    !Spectra amplitudes
    real(dp) :: As
    real(dp) :: A_iso
    real(dp) :: A_pnad
    real(dp) :: A_ent
    real(dp) :: A_cross_ad_iso
    !Spectral indices
    real(dp) :: ns
    real(dp) :: nt
    real(dp) :: n_iso
    real(dp) :: n_pnad
    real(dp) :: n_ent
    !Tensor-to-scalar
    real(dp) :: r
    !Running, etc
    real(dp) :: alpha_s
    real(dp) :: runofrun
    !Non-Gaussianity
    real(dp) :: f_NL
    real(dp) :: tau_NL
    real(dp) , dimension(100) :: k_array  !<- Added for the FULL POW SPEC
    real(dp) , dimension(100) :: pks_array  !<- Added for the FULL POW SPEC
    real(dp) , dimension(100) :: pks_iso_array  !<- Added for the FULL POW SPEC
    real(dp) , dimension(100) :: pkt_array  !<- Added for the FULL POW SPEC
    integer :: k_size
  end type gambit_inflation_observables


!Define some macros for global use
! #include "modpk_macros.f90"

contains

  subroutine multimodecode_gambit_driver(goutput_inflation_observables,&
                                         ginput_num_inflaton,&
                                         ginput_potential_choice,&
                                         ginput_slowroll_infl_end,&
                                         ginput_instreheat,&
                                         ginput_vparam_rows,&
                                         ginput_use_deltaN_SR,&
                                         ginput_evaluate_modes,&
                                         ginput_use_horiz_cross_approx,&
                                         ginput_get_runningofrunning,&
                                         ginput_ic_sampling,&
                                         ginput_energy_scale,&
                                         ginput_numb_samples,&
                                         ginput_save_iso_N,&
                                         ginput_N_iso_ref,&
                                         ginput_param_sampling,&
                                         ginput_vp_prior_min,&
                                         ginput_vp_prior_max,&
                                         ginput_varying_N_pivot,&
                                         ginput_use_first_priorval,&
                                         ginput_phi_init0,&
                                         ginput_dphi_init0,&
                                         ginput_vparams,&
                                         ginput_N_pivot,&
                                         ginput_k_pivot,&
                                         ginput_dlnk, &
                                         ginput_calc_full_pk, &
                                         ginput_steps, &
                                         ginput_kmin, &
                                         ginput_phi0_priors_min, &
                                         ginput_phi0_priors_max, &
                                         ginput_dphi0_priors_min, &
                                         ginput_dphi0_priors_max, &
                                         ginput_N_pivot_prior_min, &
                                         ginput_N_pivot_prior_max &
                                         )
    ! Gambit output for MultiModeCode results
    type(gambit_inflation_observables), intent(inout) :: goutput_inflation_observables

    ! Gambit interface input parameters
    integer, intent(in) :: ginput_num_inflaton
    integer, intent(in) :: ginput_potential_choice
    logical, intent(in) :: ginput_slowroll_infl_end
    logical, intent(in) :: ginput_instreheat
    integer, intent(in) :: ginput_vparam_rows
    logical, intent(in) :: ginput_use_deltaN_SR
    logical, intent(in) :: ginput_evaluate_modes
    logical, intent(in) :: ginput_use_horiz_cross_approx
    logical, intent(in) :: ginput_get_runningofrunning
    integer, intent(in) :: ginput_ic_sampling
    real(dp), intent(in) :: ginput_energy_scale
    integer, intent(in) :: ginput_numb_samples
    logical, intent(in) :: ginput_save_iso_N
    integer, intent(in) :: ginput_N_iso_ref
    integer, intent(in) :: ginput_param_sampling
    real(dp), dimension(ginput_num_inflaton,&
					ginput_num_inflaton),&
                    intent(in) :: ginput_vp_prior_min
    real(dp), dimension(ginput_num_inflaton,&
					ginput_num_inflaton),&
					intent(in) :: ginput_vp_prior_max
    logical, intent(in) :: ginput_varying_N_pivot
    logical, intent(in) :: ginput_use_first_priorval
    real(dp), dimension(ginput_num_inflaton), &
                    intent(in) :: ginput_phi_init0
    real(dp), dimension(ginput_num_inflaton), &
                    intent(in) :: ginput_dphi_init0
    real(dp), dimension(ginput_vparam_rows, &
                    ginput_num_inflaton), &
                    intent(in) :: ginput_vparams
    real(dp), intent(in) :: ginput_N_pivot
    real(dp), intent(in) :: ginput_k_pivot
    real(dp), intent(in) :: ginput_dlnk
    real(dp), intent(in) :: ginput_kmin
    integer, intent(in) :: ginput_steps
    logical, intent(in) :: ginput_calc_full_pk
    real(dp), dimension(ginput_num_inflaton), intent(in) :: ginput_phi0_priors_min
    real(dp), dimension(ginput_num_inflaton), intent(in) :: ginput_phi0_priors_max
    real(dp), dimension(ginput_num_inflaton), intent(in) :: ginput_dphi0_priors_min
    real(dp), dimension(ginput_num_inflaton), intent(in) :: ginput_dphi0_priors_max
    real(dp), intent(in) :: ginput_N_pivot_prior_min
    real(dp), intent(in) :: ginput_N_pivot_prior_max

    !---------------------------------------------------------------
    !At the moment, the options related to the reheating stuff from the code are
    !not available to be controlled from the gambit function. MultiModeCode at
    !this point does not have reheating option (other than instant-).

    type(observables) :: observs, observs_SR

    !Run-specific input params
    integer :: sample_looper, vparam_rows

    !Parallel variables
    integer :: numtasks, rank

    !Cosmology
    real(dp) :: dlnk

    !Sampling parameters for ICs
    integer :: numb_samples,steps
    integer :: out_adiab
    real(dp) :: energy_scale
    real(dp), dimension(:,:), allocatable :: icpriors_min, icpriors_max


    !Other sampling params
    real(dp) :: N_pivot_prior_min, N_pivot_prior_max, k_min
    integer :: num_inflaton_prior_min, num_inflaton_prior_max
    logical :: varying_N_pivot, varying_num_inflaton
    logical :: more_potential_params
    logical :: get_runningofrunning
    logical :: use_horiz_cross_approx
    logical :: calc_full_pk
    integer :: pfile

    ! Gambit interface giving the values
    ! to the MultiModeCode IC parameters
    num_inflaton = ginput_num_inflaton
    potential_choice = ginput_potential_choice
    slowroll_infl_end = ginput_slowroll_infl_end
    instreheat = ginput_instreheat
    vparam_rows = ginput_vparam_rows
    use_deltaN_SR = ginput_use_deltaN_SR
    evaluate_modes = ginput_evaluate_modes
    use_horiz_cross_approx = ginput_use_horiz_cross_approx
    get_runningofrunning = ginput_get_runningofrunning
    k_min = ginput_kmin
    calc_full_pk = ginput_calc_full_pk

    ! {1=reg_samp, 2=eqen_samp, 3=slowroll_samp, 6=isoN}
    ic_sampling = ginput_ic_sampling
	  energy_scale = ginput_energy_scale
    numb_samples = ginput_numb_samples
    save_iso_N = ginput_save_iso_N
    N_iso_ref = ginput_N_iso_ref

    call deallocate_vars()
    !Make some arrays
    call allocate_vars()

    param_sampling = ginput_param_sampling
    vp_prior_min = ginput_vp_prior_min
    vp_prior_max = ginput_vp_prior_max
    varying_N_pivot = ginput_varying_N_pivot
    use_first_priorval = ginput_use_first_priorval
    phi_init0 = ginput_phi_init0
    dphi_init0 = ginput_dphi_init0
    vparams = ginput_vparams
    N_pivot = ginput_N_pivot
    k_pivot = ginput_k_pivot
    dlnk = ginput_dlnk

    !---------------------------------------------------------------
    ! setting up the printer options at the code-level.
    ! TODO: make this from GAMBIT yaml file. 
    ! NB: this may not be interesting or necessary.
    out_opt%modpkoutput = .false.
    out_opt%output_reduced = .true.
    out_opt%output_badic =.false.

    out_opt%save_traj = .false.
    out_opt%fields_horiz = .false.
    out_opt%fields_end_infl = .false.
    out_opt%spectra = .false.
    out_opt%modes = .false.

	  steps = ginput_steps

    !---------------------------------------------------------------
    ! setting up the integrater options at the code-level.
    ! TODO: make this from GAMBIT yaml file.
    assert%use_assertions = .true.

    tech_opt%accuracy_setting = 1
    tech_opt%use_dvode_integrator = .false.

    !---------------------------------------------------------------
    !---------------------Factory settings--------------------------
    !---------------------------------------------------------------
    tech_opt%rk_accuracy_modes = 1.0e-7
    tech_opt%rk_accuracy_back = 1.0e-6

    tech_opt%dvode_rtol_back = 1.0e-6
    tech_opt%dvode_rtol_modes = 1.0e-6

    tech_opt%dvode_atol_back(1:4) = (/1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14/)

    tech_opt%dvode_atol_modes_real(1:2) = (1e-8, 1e-8)
    tech_opt%dvode_atol_modes_real(3:4) = (1e-6, 1e-6)
    tech_opt%dvode_atol_modes_real(5:8) = (/1e-5, 1e-5, 1e-5, 1e-5/)
    tech_opt%dvode_atol_modes_real(9:12) = (/1e-5, 1e-5, 1e-5, 1e-5/)
    tech_opt%dvode_atol_modes_real(13:14) = (1e-8, 1e-8)
    tech_opt%dvode_atol_modes_real(15:16) = (1e-3, 1e-3)

    tech_opt%dvode_atol_modes_imag(1:2) = (1e0, 1e0)
    tech_opt%dvode_atol_modes_imag(3:4) = (1e0, 1e0)
    tech_opt%dvode_atol_modes_imag(5:8) = (/1e-5, 1e-5, 1e-5, 1e-5/)
    tech_opt%dvode_atol_modes_imag(9:12) = (/1e-5, 1e-5, 1e-5, 1e-5/)
    tech_opt%dvode_atol_modes_imag(13:14) = (1e-8, 1e-8)
    tech_opt%dvode_atol_modes_imag(15:16) = (1e-3, 1e-3)


    !---------------------------------------------------------------
    call observs%set_zero()
    call observs_SR%set_zero()


	!Set random seed
	call init_random_seed()

    call output_initial_data()

    if (ic_sampling==ic_flags%reg_samp) then

      call out_opt%open_files(SR=use_deltaN_SR)

	  call calculate_pk_observables(goutput_inflation_observables,observs,observs_SR,k_pivot,dlnk,calc_full_pk,steps,k_min)

      call out_opt%close_files(SR=use_deltaN_SR)

  else if &
    !Eqen sampling
    (ic_sampling == ic_flags%eqen_samp .or. &
    !Set vels in SR and fields on iso-N surface for N-quad
    ic_sampling == ic_flags%iso_N .or.&
    !Set vels in SR
    ic_sampling == ic_flags%slowroll_samp) then

    call out_opt%open_files(ICs=.true., SR=use_deltaN_SR)

      do sample_looper=1,numb_samples

        if (out_opt%modpkoutput) write(*,*) &
          "---------------------------------------------"
        if (out_opt%modpkoutput) write(*,*) &
          "Sample numb", sample_looper, "of", numb_samples
        if (out_opt%modpkoutput) write(*,*) &
          "---------------------------------------------"

        call calculate_pk_observables(goutput_inflation_observables,observs,observs_SR,k_pivot,dlnk,calc_full_pk,steps,k_min)

      end do

      call out_opt%close_files(ICs=.true., SR=use_deltaN_SR)

    else
      print*, "MODECODE: sampling technique = ",ic_sampling
      call raise%fatal_code(&
        "This sampling technique is not implemented.",&
        __FILE__, __LINE__)

    end if
    ! it might make sense instead to put zero's etc. to the values here instead of just returning
	if (potential_choice == 18) then
      if (.not. observs%is_ic_ok) then
        ! print*, "SMASH: exiting pk_observables "
        RETURN
	  else
		! at the moment, smash potential is only calculated by solving the background dynamics
		! and using slow roll approximation. going beyong would require solving first order 
		! mode equations in perturbation theory.
	    goutput_inflation_observables%As = observs_SR%As
	    goutput_inflation_observables%ns = observs_SR%ns
	    goutput_inflation_observables%nt = observs_SR%nt
	    goutput_inflation_observables%r = observs_SR%r
	    goutput_inflation_observables%f_NL = observs_SR%f_NL
	    goutput_inflation_observables%tau_NL = observs_SR%tau_NL
	    goutput_inflation_observables%alpha_s = observs_SR%alpha_s
      end if
    else
      !-------------setting up the observables--------------------
      goutput_inflation_observables%As = observs%As
      goutput_inflation_observables%A_iso = observs%A_iso
      goutput_inflation_observables%A_pnad = observs%A_pnad
      goutput_inflation_observables%A_ent = observs%A_ent
      goutput_inflation_observables%A_cross_ad_iso = observs%A_cross_ad_iso
      goutput_inflation_observables%ns = observs%ns
      goutput_inflation_observables%nt = observs%nt
      goutput_inflation_observables%n_iso = observs%n_iso
      goutput_inflation_observables%n_pnad = observs%n_pnad
      goutput_inflation_observables%n_ent = observs%n_ent
      goutput_inflation_observables%r = observs%r
      goutput_inflation_observables%alpha_s = observs%alpha_s
      goutput_inflation_observables%runofrun = observs%runofrun
      goutput_inflation_observables%f_NL = observs%f_NL
	    goutput_inflation_observables%tau_NL = observs%tau_NL
      !-------------setting up the observables--------------------

    end if


  contains

    subroutine gambit_get_full_pk(pk_arr,calc_full_pk,ginput_steps,ginput_kmin,ginput_kmax)
      !Find P(k) for all the power spectra from kmin to kmax, as given in the
      !parameters file.
      real(dp), intent(in) :: ginput_kmin, ginput_kmax
      integer, intent(in) :: ginput_steps

      real(dp), dimension(:,:), allocatable, intent(out) :: pk_arr

      real(dp) :: kmin, kmax, incr
      logical, intent(in) :: calc_full_pk
      real(dp) :: p_scalar, p_tensor, p_zeta, p_iso
      real(dp), dimension(:), allocatable :: k_input
      integer :: i, steps, u

      type(power_spectra) :: pk

	  kmin = ginput_kmin
      kmax = ginput_kmax
      steps = ginput_steps

      ! print*,"ModeCode DEBUG: we are inside gambit_get_full_pk"

      ! print*,"calc_full_pk = ",calc_full_pk
      !If don't want full spectrum, return
      if (calc_full_pk) then

        ! print*,"ModeCode DEBUG: we are pass the checkpoint"

        !Make the output arrays
        if (allocated(pk_arr)) deallocate(pk_arr)

       !  print*,"ModeCode DEBUG: steps (inside) = ",steps

        allocate(pk_arr(steps, 9))

        pk_arr=0e0_dp

        !Make the arrays for k values to sample
        allocate(k_input(steps))
        incr=(kmax/kmin)**(1/real(steps-1,kind=dp))
        do i=1,steps
          k_input(i) = kmin*incr**(i-1)
        end do

        do i=1,steps
          call evolve(k_input(i), pk)

          pk_arr(i,:)=(/k_input(i),&
            pk%adiab, &
            pk%isocurv, &
            pk%entropy, &
            pk%pnad, &
            pk%tensor, &
            pk%pressure, &
            pk%press_ad, &
            pk%cross_ad_iso /)

          ! print*,"k_input = ",k_input(i)
          ! print*,"pk%adiab = ",pk%adiab
          ! print*,"pk%isocurv = ",pk%isocurv
          ! print*,"pk%entropy = ",pk%entropy
          ! print*,"pk%pnad = ",pk%pnad
          ! print*,"pk%tensor = ",pk%tensor

	    end do
      end if

	end subroutine gambit_get_full_pk

    subroutine allocate_vars()
!      !Allocate all the necessary arrays that we can with the information given
!      !in the parameters file.

      !Model dependent
      if (potential_choice==8) then
        allocate(vparams(1,4))
      else
        allocate(vparams(vparam_rows,num_inflaton))
      end if

      allocate(icpriors_max(2,num_inflaton))
      allocate(icpriors_min(2,num_inflaton))

      allocate(vp_prior_max(vparam_rows,num_inflaton))
      allocate(vp_prior_min(vparam_rows,num_inflaton))

      allocate(phi_init0(num_inflaton))
      allocate(phi_init(num_inflaton))
      allocate(phidot_sign(num_inflaton))
      allocate(phiarr(num_inflaton, 1:nsteps))
      allocate(dphiarr(num_inflaton, 1:nsteps))
      allocate(param_arr(1:nsteps))
      allocate(phi_infl_end(num_inflaton))
      allocate(phi_pivot(num_inflaton))
      allocate(dphi_pivot(num_inflaton))
      allocate(dphi_init0(num_inflaton))
      allocate(dphi_init(num_inflaton))

    end subroutine allocate_vars

!-------------------DEVELOPEMENT-------------------------------
    subroutine deallocate_vars()

      if (allocated(vparams)) deallocate(vparams)
      if (allocated(icpriors_max)) deallocate(icpriors_max)
      if (allocated(icpriors_min)) deallocate(icpriors_min)
      if (allocated(vp_prior_max)) deallocate(vp_prior_max)
      if (allocated(vp_prior_min)) deallocate(vp_prior_min)

      if (allocated(phi_init0)) deallocate(phi_init0)
      if (allocated(phi_init)) deallocate(phi_init)
      if (allocated(phidot_sign)) deallocate(phidot_sign)
      if (allocated(phiarr)) deallocate(phiarr)
      if (allocated(dphiarr)) deallocate(dphiarr)
      if (allocated(param_arr)) deallocate(param_arr)
      if (allocated(phi_infl_end)) deallocate(phi_infl_end)
      if (allocated(phi_pivot)) deallocate(phi_pivot)
      if (allocated(dphi_pivot)) deallocate(dphi_pivot)
      if (allocated(dphi_init0)) deallocate(dphi_init0)
      if (allocated(dphi_init)) deallocate(dphi_init)

    end subroutine deallocate_vars
!-------------------DEVELOPEMENT-------------------------------

    subroutine output_observables(pk_arr,&
        calc_full_pk, &
        observ_modes, observ_SR)
      !Write the observables to screen/file at the end of a successful run.

      type(observables), intent(in), optional :: observ_modes
      type(observables), intent(in), optional :: observ_SR
      real(dp), dimension(:,:), intent(in) :: pk_arr

      type(observables) :: SR_pred

      logical :: calc_full_pk

      integer :: i

      if (present(observ_SR)) then
        SR_pred = observ_SR
      else
        call SR_pred%set_zero()
      end if

      !Background stats

      if (.not. out_opt%output_reduced) then
        do i=1,size(vparams,1)
          print*, "vparams =", vparams(i,:)
        end do
      end if

      write(*, out_opt%i_fmt) &
        "Number of Inflaton =", num_inflaton
      write(*, out_opt%i_fmt) &
        "Potential Choice =", potential_choice
      write(*, out_opt%e_fmt) &
        "N_pivot =", N_pivot
      if (potential_choice==1) then
        write(*, out_opt%e2_fmt)&
          "N_tot =", N_tot,'(', &
          0.25e0_dp*sum(phi_init0**2) , ')'
      else
        write(*, out_opt%e2_fmt)&
          "N_tot =", N_tot
      end if


      !Mode stats

      if (.not. evaluate_modes) then
        write(*, out_opt%e2_fmt)&
          "Ps =", SR_pred%As
        if (potential_choice==1) then
          write(*, out_opt%e2_fmt)&
            "r = Pt/Ps =", SR_pred%r, '(', 8.0/N_pivot, ')'
        else if (potential_choice==16) then
          write(*, out_opt%e2_fmt)&
            "r = Pt/Ps =", SR_pred%r, '(', 4.0*vparams(2,1)/N_pivot, ')'
        else
          write(*, out_opt%e2_fmt)&
            "r =", SR_pred%r
        end if
        write(*, out_opt%e2_fmt)&
          "n_s =", SR_pred%ns
        write(*, out_opt%e2_fmt)&
          "n_t =", SR_pred%nt
        write(*, out_opt%e2_fmt)&
          "r/n_t =", -(8.0-SR_pred%r/2.0)
        write(*, out_opt%e2_fmt)&
          "alpha_s =", SR_pred%alpha_s
        write(*, out_opt%e2_fmt)&
          "f_NL =", SR_pred%f_NL
        write(*, out_opt%e2_fmt)&
          "tau_NL =", SR_pred%tau_NL, &
          '(>', ((6.0/5.0)*SR_pred%f_NL)**2, ')'

        return
      end if

      write(*, out_opt%e2_fmt)&
        "Ps =", observ_modes%As, '(', SR_pred%As , ')'
      write(*, out_opt%e2_fmt),&
        "Isocurvature P =", observ_modes%A_iso
      write(*, out_opt%e2_fmt),&
        "Pnad P =", observ_modes%A_pnad
      write(*, out_opt%e2_fmt),&
        "Entropy P =", observ_modes%A_ent
      write(*, out_opt%e2_fmt),&
        "Cross Ad-Iso P =", observ_modes%A_cross_ad_iso
      write(*, out_opt%e2_fmt),&
        "Bundle Width =", field_bundle%exp_scalar
      write(*, out_opt%e2_fmt)&
        "r = Pt/Ps =", observ_modes%r, '(', SR_pred%r, ')'

      if (potential_choice==1 .or. potential_choice==15) then
        write(*, out_opt%e2_fmt)&
          "r (m^2 phi^2) =", observ_modes%r, '(', 8.0/N_pivot, ')'
      end if

      write(*, out_opt%e2_fmt)&
        "n_s =", observ_modes%ns, '(', SR_pred%ns,')'
      if (num_inflaton>1) then
        write(*, out_opt%e2_fmt)&
          "n_iso =",  observ_modes%n_iso
        write(*, out_opt%e2_fmt)&
          "n_pnad =", observ_modes%n_pnad
        write(*, out_opt%e2_fmt)&
          "n_ent =",  observ_modes%n_ent
      end if
      write(*, out_opt%e2_fmt)&
        "n_t =", observ_modes%nt, '(', SR_pred%nt , ')'
      write(*, out_opt%e2_fmt)&
        "r/n_t =", observ_modes%r/observ_modes%nt, '(', -(8.0-observ_modes%r/2.0) , ')'
      write(*, out_opt%e2_fmt)&
        "alpha_s =", observ_modes%alpha_s, '(', SR_pred%alpha_s , ')'
      if (get_runningofrunning) then
        write(*, out_opt%e2_fmt)&
          "d^2n_s/dlnk^2 =", observ_modes%runofrun
      end if
      write(*, out_opt%e2_fmt)&
        "Slow-roll f_NL =", observ_modes%f_NL
      write(*, out_opt%e2_fmt)&
        "Slow-roll tau_NL =", observ_modes%tau_NL, &
        '(>', ((6.0/5.0)*observ_modes%f_NL)**2, ')'


      if (calc_full_pk .and. evaluate_modes) then

        !Open output files
        open(newunit=out_adiab,file="out_pk.csv")

        !Write the column header
        call csv_write(&
          out_adiab,&
          (/&
          !Strings have to be same length for gfortran
          'k        ', &
          'P_ad     ', &
          'P_iso    ', &
          'P_ent    ', &
          'P_nad    ', &
          'P_tens   ', &
          'P_press  ', &
          'P_pressad', &
          'P_cross  '/), &
          advance=.true.)

		do i=1,size(pk_arr,1)
		  call csv_write(out_adiab,&
					   pk_arr(i,:),&
					   advance=.true.)
		end do

        if (out_opt%modpkoutput)&
          write(*,*) "Full P(k) written to out_pk.csv"

        !Close output files
        close(out_adiab)


      end if

    end subroutine output_observables

    subroutine output_initial_data()
      !Write the initial data to screen.

      integer :: i

      call out_opt%formatting(num_inflaton)

      do i=1, size(vparams,1)
        if (out_opt%modpkoutput .and. .not. out_opt%output_reduced) &
          write(*, '(A8,I1,A5,100E12.3)'),&
          "vparams(",i,",:) =", vparams(i,:)
      end do

    end subroutine output_initial_data


    !Calculate observables, optionally grab a new IC or a new set of parameters
    !each time this routine is called.
    subroutine calculate_pk_observables(observs_gambit,observs,observs_SR,k_pivot,dlnk,calc_full_pk,steps,kmin)

      real(dp), intent(in) :: k_pivot,dlnk
      logical, intent(in) :: calc_full_pk
      real(dp) :: kmax
      type(observables), intent(inout) :: observs, observs_SR
      type(gambit_inflation_observables), intent(inout) :: observs_gambit
	  integer, intent(in) :: steps
      real(dp), intent(in) :: kmin

      real(dp), dimension(:,:), allocatable :: pk_arr
      logical :: leave

      type(power_spectra) :: pk0, pk1, pk2, pk3, pk4

      real(dp), dimension(:), allocatable :: k_a, pks_a, pkt_a

      character(1024) :: cname
      integer :: ii

      call observs%set_zero()
      call observs_SR%set_zero()

      pk_bad = run_outcome%success
      leave = .false.

!-------------------DEVELOPEMENT-------------------------------
	  ! print*,"we are inside calculate_pk_observables."
!      !Get e-folds after pivot scale leaves horizon

!      !num_inflaton sets many array sizes, have to treat it slightly differently
!      if (varying_num_inflaton) call get_new_num_inflaton()

      !Get vparams
      if (param_sampling == 1) then
	    call get_vparams()
      end if

      !Get ICs
      if (ic_sampling/=ic_flags%reg_samp) then
	    ! print*,"we are inside ic_sampling condition"
        call get_ic(phi_init0, dphi_init0, &
          icpriors_min, icpriors_max, &
          numb_samples,energy_scale)
      end if

      !Load ics
      allocate(observs%ic(2*num_inflaton))
      observs%ic(1:num_inflaton)=phi_init0
      observs%ic(num_inflaton+1:2*num_inflaton)=dphi_init0
      ! print*,"dphi_init0 = ", dphi_init0
      if (use_deltaN_SR) then
        allocate(observs_SR%ic(2*num_inflaton))
        observs_SR%ic = observs%ic
      end if

      ! print*,"N_pivot = ", N_pivot

      !Initialize potential and calc background
      call potinit(observs)

      ! print*,"outside of potinit"
      ! if (potential_choice == 18) then
      !   if (.not. observs%is_ic_ok) then
      !     print*, "ICs for SMASH potential resulted in too long(or short) inflation."
      !     RETURN
      !   end if
      ! end if


      !For outputting field values at horiz crossing
      if (out_opt%fields_horiz) then

        if (out_opt%first_fields_h_out) then

          !First column
          call csv_write(out_opt%fields_h_out,&
            'k',&
            advance=.false.)

          !Next num_inflaton columns
          do ii=1,num_inflaton
            write(cname, "(A3,I4.4)") "phi_piv", ii
            if (ii==num_inflaton) then
              call csv_write(&
                out_opt%fields_h_out,&
                trim(cname), &
                advance=.true.)
            else
              call csv_write(&
                out_opt%fields_h_out,&
                trim(cname), &
                advance=.false.)
            end if
          end do

          out_opt%first_fields_h_out = .false.
        end if

        call csv_write(out_opt%fields_h_out,&
          (/k, phi_pivot/),&
          advance=.true.)
      end if

      call test_bad(pk_bad, observs, leave)
      if (leave) return

      !Calculate SR approximation for sum-separable potentials
      if (use_deltaN_SR) then
        call calculate_SR_observables(observs_SR)

        !Only calculating these in SR
        observs%f_NL = observs_SR%f_NL
        observs%tau_NL = observs_SR%tau_NL
      end if

      ! print*,"observs_SR%ns = ",observs_SR%ns

      !Evaluate the mode functions
      if (evaluate_modes) then

        call evolve(k_pivot, pk0)
          call test_bad(pk_bad, observs, leave)
          if (leave) then
			return
		  end if
!DEBUG
!print*, "Not evaluating second and third evolve routines"
!stop
	     ! print*, "second evolve"
		   ! print*, "k_pivot*exp(-dlnk) = " , k_pivot*exp(-dlnk)
		   ! print*, "dlnk = ", dlnk

        call evolve(k_pivot*exp(-dlnk), pk1)
          call test_bad(pk_bad, observs, leave)
          if (leave) return
		  ! print*, "left second evolve"
		  ! print*, "third evolve"
		  ! print*, "k_pivot*exp(dlnk) = " , k_pivot*exp(dlnk)
		  call evolve(k_pivot*exp(dlnk), pk2)
          call test_bad(pk_bad, observs, leave)
          if (leave) return
		  ! print*, "left third evolve"

        if (get_runningofrunning) then
			print*, "get_runningofrunning = True"
            !Alpha_s from 5-pt stencil
            !or running of running
          call evolve(k_pivot*exp(-2.0e0_dp*dlnk), pk3)
          ! print*, "fourth evolve"
            call test_bad(pk_bad, observs, leave)
            if (leave) return
		  call evolve(k_pivot*exp(2.0e0_dp*dlnk), pk4)
		  ! print*, "fifth evolve"
            call test_bad(pk_bad, observs, leave)
            if (leave) return
        end if

        ! print*, "pk0%adiab = ", pk0%adiab

          !Construct the observables
        if (get_runningofrunning) then
          call observs%set_finite_diff(dlnk, &
            pk0,pk1,pk2,pk3,pk4, &
            field_bundle%exp_scalar)
          ! print*,"here we calculate the observables!"
        else
		  ! print*,"here we calculate the observables!"
          call observs%set_finite_diff(dlnk, &
            pk0,pk1,pk2,&
            bundle_width=field_bundle%exp_scalar)
        end if

        !Get full spectrum for adiab and isocurv at equal intvs in lnk
!		print*,"will call gambit_get_full_pk"
		!steps = 100
		kmax = 1e6
        !kmin = 1e-7

        call gambit_get_full_pk(pk_arr,calc_full_pk,steps,kmin,kmax)

        if (calc_full_pk) then

          observs_gambit%k_array   		 = pk_arr(:steps,1)
          observs_gambit%pks_array 		 = pk_arr(:steps,2)
          observs_gambit%pks_iso_array = pk_arr(:steps,3)
          observs_gambit%pkt_array 		 = pk_arr(:steps,6)

        end if


      end if

      !Write output to stdout
      if (out_opt%modpkoutput) then
        print*,"printing modpkpoutput "
        if (evaluate_modes) then
          call output_observables(pk_arr, &
            calc_full_pk, observs, observs_SR)
        else
          call output_observables(pk_arr, &
            calc_full_pk, observ_SR = observs_SR)
        end if
      end if

      ! print*,"save_iso_N=",save_iso_N

      if (save_iso_N) then
        observs%ic(1:num_inflaton) = phi_iso_N
        observs%ic(num_inflaton+1:2*num_inflaton) = dphi_iso_N
        if (use_deltaN_SR) &
          observs_SR%ic(1:num_inflaton) = phi_iso_N
        if (use_deltaN_SR) &
          observs_SR%ic(num_inflaton+1:2*num_inflaton) = dphi_iso_N

        if (out_opt%output_badic .or. pk_bad==run_outcome%success) then

          if (evaluate_modes) then
            if (out_opt%first_outsamp_N_iso) then
              call observs%print_header(out_opt%outsamp_N_iso)
              out_opt%first_outsamp_N_iso=.false.
            end if
            call observs%printout(out_opt%outsamp_N_iso)
          end if

          if (use_deltaN_SR) then
            if (out_opt%first_outsamp_N_iso_SR) then
              call observs%print_header(out_opt%outsamp_N_iso_SR)
              out_opt%first_outsamp_N_iso_SR=.false.
            end if
            call observs_SR%printout(out_opt%outsamp_N_iso_SR)
          end if

        end if
      end if

   !  print*,"no issue before here!"

    end subroutine calculate_pk_observables

    !Calculate observables for the power spectrum, as well as fNL, using the
    !delta-N formalism in slow-roll
    subroutine calculate_SR_observables(observs_SR)
      use modpk_numerics, only : locate, array_polint, polint
!	  use modpk_reheat, only : reheater
      type(observables), intent(inout) :: observs_SR
      integer :: j, i
      real(dp) :: ah, alpha_ik, dalpha, N_end, del_N, Npiv_renorm
      real(dp), dimension(num_inflaton) :: phi_pivot, phi_end, del_phi

      !Find field values at end of inflation
      !Note that eps=1 perhaps twice, so take the last one.
      call array_polint(epsarr(nactual_bg-4:nactual_bg),&
        phiarr(:,nactual_bg-4:nactual_bg),&
        1.0e0_dp,  phi_end, del_phi)
      call polint(epsarr(nactual_bg-4:nactual_bg),&
        lna(nactual_bg-4:nactual_bg),&
        1.0e0_dp,  N_end, del_N)

      !Find field values at horizon crossing
      Npiv_renorm = N_end - N_pivot

      i= locate(lna(1:nactual_bg), Npiv_renorm)
      j=min(max(i-(4-1)/2,1),nactual_bg+1-4)
      call array_polint(lna(j:j+4), phiarr(:,j:j+4),&
        Npiv_renorm, phi_pivot, del_phi)

      if (use_horiz_cross_approx) then
        HC_approx=.true.
      else
        HC_approx=.false.
      end if

      call observs_SR%set_zero()
      observs_SR%As = PR_SR(phi_pivot,phi_end)
      observs_SR%ns = ns_SR(phi_pivot,phi_end)
      observs_SR%nt = nt_SR(phi_pivot)
      observs_SR%r  = r_SR(phi_pivot,phi_end)
      observs_SR%f_NL  = fnl_SR(phi_pivot,phi_end)
      observs_SR%tau_NL  = taunl_SR(phi_pivot,phi_end)
      observs_SR%alpha_s  = alpha_s_SR(phi_pivot,phi_end)

    end subroutine calculate_SR_observables


    !Check if an IC had some properties where we declare it to be unphysical
    subroutine test_bad(pk_bad,observ,leave)

      integer,  intent(in)     :: pk_bad
      logical,  intent(inout)  :: leave
      type(observables) :: observ

      if (pk_bad /= run_outcome%success) then

        call run_outcome%print_outcome(pk_bad)
        call observ%set_zero()

        !Flag for voiding calculation
        leave = .true.
      end if

    end subroutine test_bad

  end subroutine multimodecode_gambit_driver

end module multimodecode_gambit
