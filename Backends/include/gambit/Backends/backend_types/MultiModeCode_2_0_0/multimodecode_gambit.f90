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
    real(dp) , dimension(1000) :: k_array
    real(dp) , dimension(1000) :: pks_array
    real(dp) , dimension(1000) :: pks_iso_array
    real(dp) , dimension(1000) :: pkt_array
    integer :: k_size
  end type gambit_inflation_observables


!Define some macros for global use
! #include "modpk_macros.f90"

contains

  function multimodecode_gambit_driver(ginput_num_inflaton,&
                                         ginput_potential_choice,&
                                         ginput_evaluate_modes,&
                                         ginput_get_runningofrunning,&
                                         ginput_phi_init0,&
                                         ginput_dphi_init0,&
                                         ginput_vparams,&
                                         ginput_N_pivot,&
                                         ginput_k_pivot,&
                                         ginput_dlnk, &
                                         ginput_steps, &
                                         ginput_kmin, &
                                         ginput_kmax, &
                                         ginput_vparam_rows &
																				 ) result(gambit_obs)

    type(gambit_inflation_observables) :: gambit_obs ! output

    integer :: i

    ! Gambit interface input parameters
    integer, intent(in) :: ginput_num_inflaton
    integer, intent(in) :: ginput_potential_choice
    integer, intent(in) :: ginput_vparam_rows
    logical, optional :: ginput_evaluate_modes
    logical, optional :: ginput_get_runningofrunning
    real(dp), dimension(ginput_num_inflaton), intent(in) :: ginput_phi_init0
    real(dp), dimension(ginput_num_inflaton), intent(in) :: ginput_dphi_init0
    real(dp), dimension(ginput_vparam_rows, ginput_num_inflaton), intent(in) :: ginput_vparams
    real(dp), intent(in) :: ginput_N_pivot
    real(dp), intent(in) :: ginput_k_pivot
    real(dp), intent(in) :: ginput_dlnk
    integer, intent(in) :: ginput_steps
    real(dp), intent(in) :: ginput_kmin
    real(dp), intent(in) :: ginput_kmax

    !---------------------------------------------------------------
    !At the moment, the options related to the reheating stuff from the code are
    !not available to be controlled from the gambit function. MultiModeCode at
    !this point does not have reheating option (other than instant-).

    type(observables) :: observs, observs_SR

    !Run-specific input params
    integer :: sample_looper
    integer :: vparam_rows

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
    real(dp) :: N_pivot_prior_min, N_pivot_prior_max, k_min, k_max
    integer :: num_inflaton_prior_min, num_inflaton_prior_max
    logical :: varying_N_pivot, varying_num_inflaton
    logical :: more_potential_params
    logical :: get_runningofrunning
    logical :: use_horiz_cross_approx
    logical :: calc_full_pk
    integer :: pfile

		call deallocate_vars()

    slowroll_infl_end = .true.
    instreheat = .true.
    use_deltaN_SR = .false.
    use_horiz_cross_approx = .false.
    ic_sampling = 1
    energy_scale = 0.1
    numb_samples = 1
    save_iso_N = .false.
    N_iso_ref = 55
    param_sampling = 1
    varying_N_pivot = .false. ! --
    use_first_priorval = .false. ! --

    ! Gambit interface giving the values
    ! to the MultiModeCode IC parameters
    num_inflaton = ginput_num_inflaton
    potential_choice = ginput_potential_choice
    vparam_rows = ginput_vparam_rows

    !Make some arrays
    call allocate_vars()

    evaluate_modes = ginput_evaluate_modes
    get_runningofrunning = ginput_get_runningofrunning
    k_min = ginput_kmin
    k_max = ginput_kmax
    steps = ginput_steps
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

		call calculate_pk_observables(gambit_obs,observs,observs_SR,k_pivot,dlnk,calc_full_pk,steps,k_min,k_max)

		! it might make sense instead to put zero's etc. to the values here instead of just returning
    if (potential_choice == 18) then
      if (.not. observs%is_ic_ok) then
        RETURN
      else
        ! at the moment, smash potential is only calculated by solving the background dynamics
        ! and using slow roll approximation. going beyong would require solving first order 
        ! mode equations in perturbation theory.
        gambit_obs%As = observs_SR%As
        gambit_obs%ns = observs_SR%ns
        gambit_obs%nt = observs_SR%nt
        gambit_obs%r = observs_SR%r
        gambit_obs%f_NL = observs_SR%f_NL
        gambit_obs%tau_NL = observs_SR%tau_NL
        gambit_obs%alpha_s = observs_SR%alpha_s
      end if
    else
      !-------------setting up the observables--------------------
      gambit_obs%As = observs%As
      gambit_obs%A_iso = observs%A_iso
      gambit_obs%A_pnad = observs%A_pnad
      gambit_obs%A_ent = observs%A_ent
      gambit_obs%A_cross_ad_iso = observs%A_cross_ad_iso
      gambit_obs%ns = observs%ns
      gambit_obs%nt = observs%nt
      gambit_obs%n_iso = observs%n_iso
      gambit_obs%n_pnad = observs%n_pnad
      gambit_obs%n_ent = observs%n_ent
      gambit_obs%r = observs%r
      gambit_obs%alpha_s = observs%alpha_s
      gambit_obs%runofrun = observs%runofrun
      gambit_obs%f_NL = observs%f_NL
      gambit_obs%tau_NL = observs%tau_NL
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

      !If don't want full spectrum, return
      if (calc_full_pk) then

        !Make the output arrays
        if (allocated(pk_arr)) deallocate(pk_arr)

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

    !Calculate observables, optionally grab a new IC or a new set of parameters
    !each time this routine is called.
    subroutine calculate_pk_observables(observs_gambit,observs,observs_SR,k_pivot,dlnk,calc_full_pk,steps,kmin,kmax)

      real(dp), intent(in) :: k_pivot,dlnk
      logical, intent(in) :: calc_full_pk
      type(observables), intent(inout) :: observs, observs_SR
      type(gambit_inflation_observables), intent(inout) :: observs_gambit
      integer, intent(in) :: steps
      real(dp), intent(in) :: kmin, kmax

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

      !Get vparams
      if (param_sampling == 1) then
        call get_vparams()
      end if

      !Load ics
      allocate(observs%ic(2*num_inflaton))
      observs%ic(1:num_inflaton)=phi_init0
      observs%ic(num_inflaton+1:2*num_inflaton)=dphi_init0
      if (use_deltaN_SR) then
        allocate(observs_SR%ic(2*num_inflaton))
        observs_SR%ic = observs%ic
      end if

      !Initialize potential and calc background
      call potinit(observs)

      ! if (potential_choice == 18) then
      !   if (.not. observs%is_ic_ok) then
      !     print*, "ICs for SMASH potential resulted in too long(or short) inflation."
      !     RETURN
      !   end if
      ! end if

      call test_bad(pk_bad, observs, leave)
      if (leave) return

      !Calculate SR approximation for sum-separable potentials
      if (use_deltaN_SR) then
        call calculate_SR_observables(observs_SR)

        !Only calculating these in SR
        observs%f_NL = observs_SR%f_NL
        observs%tau_NL = observs_SR%tau_NL
      end if

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

        call evolve(k_pivot*exp(-dlnk), pk1)
          call test_bad(pk_bad, observs, leave)
          if (leave) return

          call evolve(k_pivot*exp(dlnk), pk2)
          call test_bad(pk_bad, observs, leave)
          if (leave) return

        if (get_runningofrunning) then
            !Alpha_s from 5-pt stencil
            !or running of running
          call evolve(k_pivot*exp(-2.0e0_dp*dlnk), pk3)

            call test_bad(pk_bad, observs, leave)
            if (leave) return
          call evolve(k_pivot*exp(2.0e0_dp*dlnk), pk4)

            call test_bad(pk_bad, observs, leave)
            if (leave) return
        end if

          !Construct the observables
        if (get_runningofrunning) then
          call observs%set_finite_diff(dlnk, &
            pk0,pk1,pk2,pk3,pk4, &
            field_bundle%exp_scalar)
          !print*,"here we calculate the observables!"
        else
          !print*,"here we calculate the observables!"
          call observs%set_finite_diff(dlnk, &
            pk0,pk1,pk2,&
            bundle_width=field_bundle%exp_scalar)
        end if

        !Get full spectrum for adiab and isocurv at equal intvs in lnk

        call gambit_get_full_pk(pk_arr,calc_full_pk,steps,kmin,kmax)

        if (calc_full_pk) then

          observs_gambit%k_array         = pk_arr(:steps,1)
          observs_gambit%pks_array       = pk_arr(:steps,2)
          observs_gambit%pks_iso_array   = pk_arr(:steps,3)
          observs_gambit%pkt_array       = pk_arr(:steps,6)

        end if
      end if

    end subroutine calculate_pk_observables

    !Calculate observables for the power spectrum, as well as fNL, using the
    !delta-N formalism in slow-roll
    subroutine calculate_SR_observables(observs_SR)
      use modpk_numerics, only : locate, array_polint, polint
!     use modpk_reheat, only : reheater
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

  end function multimodecode_gambit_driver

  function multimodecode_parametrised_ps(ginput_num_inflaton,&
                                         ginput_potential_choice,&
                                         ginput_evaluate_modes,&
                                         ginput_get_runningofrunning,&
                                         ginput_phi_init0,&
                                         ginput_dphi_init0,&
                                         ginput_vparams,&
                                         ginput_N_pivot,&
                                         ginput_k_pivot,&
                                         ginput_dlnk, &
                                         ginput_vparam_rows &
                                         ) result(gambit_obs)

    type(gambit_inflation_observables) :: gambit_obs ! output

    integer :: i

    ! Gambit interface input parameters
    integer, intent(in) :: ginput_num_inflaton
    integer, intent(in) :: ginput_potential_choice
    integer, intent(in) :: ginput_vparam_rows
    logical, optional :: ginput_evaluate_modes
    logical, optional :: ginput_get_runningofrunning
    real(dp), dimension(ginput_num_inflaton), intent(in) :: ginput_phi_init0
    real(dp), dimension(ginput_num_inflaton), intent(in) :: ginput_dphi_init0
    real(dp), dimension(ginput_vparam_rows, ginput_num_inflaton), intent(in) :: ginput_vparams
    real(dp), intent(in) :: ginput_N_pivot
    real(dp), intent(in) :: ginput_k_pivot
    real(dp), intent(in) :: ginput_dlnk

    !---------------------------------------------------------------
    !At the moment, the options related to the reheating stuff from the code are
    !not available to be controlled from the gambit function. MultiModeCode at
    !this point does not have reheating option (other than instant-).
    type(observables) :: observs, observs_SR

    !Run-specific input params
    integer :: sample_looper
    integer :: vparam_rows

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
    real(dp) :: N_pivot_prior_min, N_pivot_prior_max
    integer :: num_inflaton_prior_min, num_inflaton_prior_max
    logical :: varying_N_pivot, varying_num_inflaton
    logical :: more_potential_params
    logical :: get_runningofrunning
    logical :: use_horiz_cross_approx
    logical :: calc_full_pk
    integer :: pfile

		call deallocate_vars()

    calc_full_pk = .false.
    slowroll_infl_end = .true.
    instreheat = .true.
    use_deltaN_SR = .false.
    use_horiz_cross_approx = .false.
    ic_sampling = 1
    energy_scale = 0.1
    numb_samples = 1
    save_iso_N = .false.
    N_iso_ref = 55
    param_sampling = 1
    varying_N_pivot = .false.
    use_first_priorval = .false.

    ! Gambit interface giving the values
    ! to the MultiModeCode IC parameters
    num_inflaton = ginput_num_inflaton
    potential_choice = ginput_potential_choice
    vparam_rows = ginput_vparam_rows

    !Make some arrays
    call allocate_vars()

    evaluate_modes = ginput_evaluate_modes
    get_runningofrunning = ginput_get_runningofrunning
    phi_init0 = ginput_phi_init0
    dphi_init0 = ginput_dphi_init0
    vparams = ginput_vparams
    N_pivot = ginput_N_pivot
    k_pivot = ginput_k_pivot
    dlnk = ginput_dlnk

    out_opt%modpkoutput = .false.
    out_opt%output_reduced = .true.
    out_opt%output_badic =.false.
    out_opt%save_traj = .false.
    out_opt%fields_horiz = .false.
    out_opt%fields_end_infl = .false.
    out_opt%spectra = .false.
    out_opt%modes = .false.

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

		call calculate_pk_observables(gambit_obs,observs,observs_SR,k_pivot,dlnk)

		! it might make sense instead to put zero's etc. to the values here instead of just returning
    if (potential_choice == 18) then
      if (.not. observs%is_ic_ok) then
        ! print*, "SMASH: exiting pk_observables "
        RETURN
      else
        ! at the moment, smash potential is only calculated by solving the background dynamics
        ! and using slow roll approximation. going beyong would require solving first order
        ! mode equations in perturbation theory.
        gambit_obs%As = observs_SR%As
        gambit_obs%ns = observs_SR%ns
        gambit_obs%nt = observs_SR%nt
        gambit_obs%r = observs_SR%r
        gambit_obs%f_NL = observs_SR%f_NL
        gambit_obs%tau_NL = observs_SR%tau_NL
        gambit_obs%alpha_s = observs_SR%alpha_s
      end if
    else
      !-------------setting up the observables--------------------
      gambit_obs%As = observs%As
      gambit_obs%A_iso = observs%A_iso
      gambit_obs%A_pnad = observs%A_pnad
      gambit_obs%A_ent = observs%A_ent
      gambit_obs%A_cross_ad_iso = observs%A_cross_ad_iso
      gambit_obs%ns = observs%ns
      gambit_obs%nt = observs%nt
      gambit_obs%n_iso = observs%n_iso
      gambit_obs%n_pnad = observs%n_pnad
      gambit_obs%n_ent = observs%n_ent
      gambit_obs%r = observs%r
      gambit_obs%alpha_s = observs%alpha_s
      gambit_obs%runofrun = observs%runofrun
      gambit_obs%f_NL = observs%f_NL
      gambit_obs%tau_NL = observs%tau_NL
      !-------------setting up the observables--------------------

    end if

  contains

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

    !Calculate observables, optionally grab a new IC or a new set of parameters
    !each time this routine is called.
    subroutine calculate_pk_observables(observs_gambit,observs,observs_SR,k_pivot,dlnk)

      real(dp), intent(in) :: k_pivot,dlnk
      type(observables), intent(inout) :: observs, observs_SR
      type(gambit_inflation_observables), intent(inout) :: observs_gambit

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

      !Get vparams
      if (param_sampling == 1) then
        call get_vparams()
      end if

      !Load ics
      allocate(observs%ic(2*num_inflaton))
      observs%ic(1:num_inflaton)=phi_init0
      observs%ic(num_inflaton+1:2*num_inflaton)=dphi_init0
      if (use_deltaN_SR) then
        allocate(observs_SR%ic(2*num_inflaton))
        observs_SR%ic = observs%ic
      end if

      !Initialize potential and calc background
      call potinit(observs)

      call test_bad(pk_bad, observs, leave)
      if (leave) return

      !Calculate SR approximation for sum-separable potentials
      if (use_deltaN_SR) then
        call calculate_SR_observables(observs_SR)

        !Only calculating these in SR
        observs%f_NL = observs_SR%f_NL
        observs%tau_NL = observs_SR%tau_NL
      end if

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

        call evolve(k_pivot*exp(-dlnk), pk1)
          call test_bad(pk_bad, observs, leave)
          if (leave) return

          call evolve(k_pivot*exp(dlnk), pk2)
          call test_bad(pk_bad, observs, leave)
          if (leave) return

        if (get_runningofrunning) then
            !Alpha_s from 5-pt stencil
            !or running of running
          call evolve(k_pivot*exp(-2.0e0_dp*dlnk), pk3)

            call test_bad(pk_bad, observs, leave)
            if (leave) return
          call evolve(k_pivot*exp(2.0e0_dp*dlnk), pk4)

            call test_bad(pk_bad, observs, leave)
            if (leave) return
        end if

          !Construct the observables
        if (get_runningofrunning) then
          call observs%set_finite_diff(dlnk, &
            pk0,pk1,pk2,pk3,pk4, &
            field_bundle%exp_scalar)
          !print*,"here we calculate the observables!"
        else
          !print*,"here we calculate the observables!"
          call observs%set_finite_diff(dlnk, &
            pk0,pk1,pk2,&
            bundle_width=field_bundle%exp_scalar)
        end if

        !Get full spectrum for adiab and isocurv at equal intvs in lnk
      end if

    end subroutine calculate_pk_observables

    !Calculate observables for the power spectrum, as well as fNL, using the
    !delta-N formalism in slow-roll
    subroutine calculate_SR_observables(observs_SR)
      use modpk_numerics, only : locate, array_polint, polint
!     use modpk_reheat, only : reheater
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

  end function multimodecode_parametrised_ps


end module multimodecode_gambit
