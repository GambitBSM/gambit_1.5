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
  ! use modpk_sampling
  use modpk_io, only : out_opt
  use modpk_deltaN
  use modpk_observables, only : observables
  use csv_file, only : csv_write
  use modpk_errorhandling, only : raise, run_outcome, assert
!  use modpk_reheat, only : reheat_opts, reheater
  use modpk_rng, only : init_random_seed ! , geometric

  implicit none

  public :: multimodecode_gambit_driver

! The structure of the output gambit reads.
  type :: gambit_inflation_observables
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
!                                         ginput_more_potential_params,&
!                                         ginput_varying_num_inflaton,&
!                                         ginput_num_inflaton_prior_min,&
!                                         ginput_num_inflaton_prior_max,&
!                                         ginput_inflaton_sampling,&
                                         ginput_use_deltaN_SR,&
                                         ginput_evaluate_modes,&
                                         ginput_use_horiz_cross_approx,&
                                         ginput_get_runningofrunning,&
                                         ginput_ic_sampling,&
                                         ginput_energy_scale,&
                                         ginput_numb_samples,&
                                         ginput_save_iso_N,&
                                         ginput_N_iso_ref,&
!                                         ginput_set_Nisoref_by_Npivot,&
                                         ginput_param_sampling,&
                                         ginput_vp_prior_min,&
                                         ginput_vp_prior_max,&
                                         ginput_varying_N_pivot,&
                                         ginput_use_first_priorval,&
!                                         ginput_numb_auxparams,&
!                                         ginput_prior_auxparams_min,&
!                                         ginput_prior_auxparams_max,&
                                         ginput_phi_init0,&
                                         ginput_dphi_init0,&
                                         ginput_vparams,&
                                         ginput_N_pivot,&
                                         ginput_k_pivot,&
                                         ginput_dlnk, &
!                                         ginput_effective_V_choice,&
!                                         ginput_turning_choice &
!                                         ginput_number_knots_qsfrandom,&
!                                         ginput_stand_dev_qsfrandom,&
!                                         ginput_knot_range_min,&
!                                         ginput_knot_range_max,&
!                                         ginput_custom_knot_range,&
!                                         ginput_opt_out,&
!                                         ginput_tech_opt,&
!                                         ginput_assert,&
!                                         ginput_reheat_opts &
                                         ginput_calc_full_pk, &
                                         ginput_steps, &
                                         ginput_kmin, &
!                                         ginput_kmax, &
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
!    logical, intent(in) :: ginput_more_potential_params
!    logical, intent(in) :: ginput_varying_num_inflaton
!    integer, intent(in) :: ginput_num_inflaton_prior_min
!    integer, intent(in) :: ginput_num_inflaton_prior_max
!    integer, intent(in) :: ginput_inflaton_sampling
    logical, intent(in) :: ginput_use_deltaN_SR
    logical, intent(in) :: ginput_evaluate_modes
    logical, intent(in) :: ginput_use_horiz_cross_approx
    logical, intent(in) :: ginput_get_runningofrunning
    integer, intent(in) :: ginput_ic_sampling
    real(dp), intent(in) :: ginput_energy_scale
    integer, intent(in) :: ginput_numb_samples
    logical, intent(in) :: ginput_save_iso_N
    integer, intent(in) :: ginput_N_iso_ref
!    logical, intent(in) :: ginput_set_Nisoref_by_Npivot

    integer, intent(in) :: ginput_param_sampling
    real(dp), dimension(ginput_num_inflaton,&
					ginput_num_inflaton),&
                    intent(in) :: ginput_vp_prior_min
    real(dp), dimension(ginput_num_inflaton,&
					ginput_num_inflaton),&
					intent(in) :: ginput_vp_prior_max
    logical, intent(in) :: ginput_varying_N_pivot
    logical, intent(in) :: ginput_use_first_priorval
!    integer, intent(in) :: ginput_numb_auxparams
!    real(dp), dimension(min(100,2*num_inflaton)), &
!                    intent(in) :: ginput_prior_auxparams_min
!    real(dp), dimension(max(100,2*num_inflaton)), &
!                    intent(in) :: ginput_prior_auxparams_max

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

!    integer, intent(in) :: ginput_effective_V_choice
!    integer, intent(in) :: ginput_turning_choice
!    integer, dimension(ginput_num_inflaton-1),&
!             intent(in) :: ginput_number_knots_qsfrandom
!    real(dp), dimension(ginput_num_inflaton-1),&
!             intent(in) :: ginput_stand_dev_qsfrandom
!    real(dp), dimension(ginput_num_inflaton-1),&
!             intent(in) :: ginput_knot_range_min
!    real(dp), dimension(ginput_num_inflaton-1),&
!             intent(in) :: ginput_knot_range_max
!    logical, intent(in) ::ginput_custom_knot_range

    real(dp), intent(in) :: ginput_kmin !, ginput_kmax
    integer, intent(in) :: ginput_steps
    logical, intent(in) :: ginput_calc_full_pk


    real(dp), dimension(ginput_num_inflaton), intent(in) :: ginput_phi0_priors_min
    real(dp), dimension(ginput_num_inflaton), intent(in) :: ginput_phi0_priors_max
    real(dp), dimension(ginput_num_inflaton), intent(in) :: ginput_dphi0_priors_min
    real(dp), dimension(ginput_num_inflaton), intent(in) :: ginput_dphi0_priors_max
    real(dp), intent(in) :: ginput_N_pivot_prior_min
    real(dp), intent(in) :: ginput_N_pivot_prior_max

    !---------------------------------------------------------------
    !At the moment, the options related to the printing out from the code are 
    !not available to be controlled from the gambit function. MultiModeCode at
    !this point does not print out anything.

    ! opt_out options (for GAMBIT) are:

    ! out_opt%modpkoutput = .false.
    ! out_opt%output_reduced = .true.
    ! out_opt%output_badic =.false.

    ! out_opt%save_traj = .false.
    ! out_opt%phi0 = .false.
    ! out_opt%fields_horiz = .false.
    ! out_opt%fields_end_infl = .false.
    ! out_opt%spectra = .false.
    ! out_opt%modes = .false.

    ! out_opt%save_reheat = .false.
    ! out_opt%save_cij = .false.

    ! out_opt%detailed%needmoreopts = .false.
    ! out_opt%detailed%write_Npiv = .false.
    ! out_opt%detailed%write_num_inflaton = .false.
    ! out_opt%detailed%write_vparams = .false.
    ! out_opt%detailed%write_auxparams = .false.

    !---------------------------------------------------------------
    !At the moment, the options related to the technichal stuff (i.e.
    !things related to integration and stepsize etc. in the code are (largely)
    !not available to be controlled from the gambit function.

    !full tech_opts are

    ! assert%use_assertions = .true.

    ! tech_opt%use_integ_with_t= .false.
    ! tech_opt%automate_singlefield_ic = .false.

    ! tech_opt%accuracy_setting = 2
    ! tech_opt%use_dvode_integrator = .true.
    ! tech_opt%use_ode_constraints = .false.

    ! tech_opt%use_analytical_jacobian = .false.

    ! tech_opt%rk_accuracy_modes = 1.0e-7
    ! tech_opt%rk_accuracy_back = 1.0e-9

    ! tech_opt%dvode_rtol_back = 1.0e-9
    ! tech_opt%dvode_rtol_modes = 1.0e-6

    ! tech_opt%dvode_atol_back = 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14

    ! tech_opt%dvode_atol_modes_real(1:2) = 1e-8 1e-8
    ! tech_opt%dvode_atol_modes_real(3:4) = 1e-6 1e-6
    ! tech_opt%dvode_atol_modes_real(5:8) = 1e-5 1e-5 1e-5 1e-5
    ! tech_opt%dvode_atol_modes_real(9:12) = 1e-5 1e-5 1e-5 1e-5
    ! tech_opt%dvode_atol_modes_real(13:14) = 1e-8 1e-8
    ! tech_opt%dvode_atol_modes_real(15:16) = 1e-3 1e-3

    ! tech_opt%dvode_atol_modes_imag(1:2) = 1e0 1e0
    ! tech_opt%dvode_atol_modes_imag(3:4) = 1e0 1e0
    ! tech_opt%dvode_atol_modes_imag(5:8) = 1e-5 1e-5 1e-5 1e-5
    ! tech_opt%dvode_atol_modes_imag(9:12) = 1e-5 1e-5 1e-5 1e-5
    ! tech_opt%dvode_atol_modes_imag(13:14) = 1e-8 1e-8
    ! tech_opt%dvode_atol_modes_imag(15:16) = 1e-3 1e-3

    ! tech_opt%dvode_dN_r = 0.001
    ! tech_opt%dvode_dN_c = 1e-4
    ! tech_opt%dvode_dt = 1e16

    !---------------------------------------------------------------
    !At the moment, the options related to the reheating stuff from the code are
    !not available to be controlled from the gambit function. MultiModeCode at
    !this point does not have reheating option (other than instant-).

    ! reheat_opts options are:

    ! reheat_opts%use_reheat = .true.
    ! reheat_opts%reheat_model = 1

    ! reheat_opts%gamma_sampler = 2

    !For run-time alloc w/out re-compile
    ! namelist /init/ num_inflaton, potential_choice, &
    !   slowroll_infl_end, instreheat, hijack_sampling, vparam_rows, &
    !   more_potential_params

    ! namelist /inflaton_sampling_nml/ varying_num_inflaton, &
    !   num_inflaton_prior_min, num_inflaton_prior_max,&
    !   inflaton_sampling

    ! namelist /analytical/ use_deltaN_SR, evaluate_modes, &
    !   use_horiz_cross_approx, get_runningofrunning

    ! namelist /ic_sampling_nml/ ic_sampling, energy_scale, numb_samples, &
    !   save_iso_N, N_iso_ref, set_Nisoref_by_Npivot

    ! namelist /param_sampling_nml/ param_sampling, vp_prior_min, vp_prior_max, &
    !   varying_N_pivot, use_first_priorval, &
    !   numb_auxparams, prior_auxparams_min, prior_auxparams_max

    ! namelist /params/ phi_init0, dphi_init0, vparams, &
    !   N_pivot, k_pivot, dlnk

    ! namelist /more_params/ effective_V_choice, turning_choice, &
    !   number_knots_qsfrandom, stand_dev_qsfrandom, &
    !   knot_range_min, knot_range_max, custom_knot_range

    ! namelist /print_out/ out_opt

    ! namelist /technical/ tech_opt, assert

    ! namelist /reheat/ reheat_opts
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    !Read initializing params from file
   	  ! open(newunit=pfile, file="parameters_multimodecode.txt", &
      ! status="old", delim = "apostrophe")
    ! read(unit=pfile, nml=init)
    ! read(unit=pfile, nml=analytical)
    ! read(unit=pfile, nml=inflaton_sampling_nml)


    !Instead the idea is to set these inital parameters to be an
    !argument for the function.
    !---------------------------------------------------------------

    type(observables) :: observs, observs_SR

    !Run-specific input params
    integer :: sample_looper, vparam_rows

    !Parallel variables
    integer :: numtasks, rank

    !Cosmology
    real(dp) :: dlnk

    !Sampling parameters for ICs
    integer :: numb_samples
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
!    logical :: hijack_sampling
    integer :: pfile

    ! Gambit interface giving the values
    ! to the MultiModeCode IC parameters
    num_inflaton = ginput_num_inflaton
    potential_choice = ginput_potential_choice
    slowroll_infl_end = ginput_slowroll_infl_end
    instreheat = ginput_instreheat
    vparam_rows = ginput_vparam_rows
!    more_potential_params = ginput_more_potential_params
!    varying_num_inflaton = ginput_varying_num_inflaton
!    num_inflaton_prior_min = ginput_num_inflaton_prior_min
!	 num_inflaton_prior_max = ginput_num_inflaton_prior_max
!    inflaton_sampling = ginput_inflaton_sampling
    use_deltaN_SR = ginput_use_deltaN_SR
    evaluate_modes = ginput_evaluate_modes
    use_horiz_cross_approx = ginput_use_horiz_cross_approx
    get_runningofrunning = ginput_get_runningofrunning

    ! {1=reg_samp, 2=eqen_samp, 3=slowroll_samp, 6=isoN}
    ic_sampling = ginput_ic_sampling
	energy_scale = ginput_energy_scale
    numb_samples = ginput_numb_samples
    save_iso_N = ginput_save_iso_N
    N_iso_ref = ginput_N_iso_ref
!    set_Nisoref_by_Npivot = ginput_set_Nisoref_by_Npivot

    print*,"num_inflaton = ",num_inflaton
    print*,"instreheat = ",instreheat
    print*,"ic_sampling = ",ic_sampling
    print*,"potential_choice = ",potential_choice
    print*,"vparam_rows = ",vparam_rows
    print*,"slowroll_infl_end = ",slowroll_infl_end


    call deallocate_vars()
    !Make some arrays
    call allocate_vars()

    !---------------------------------------------------------------
    !Read other params from file
	 ! read(unit=pfile, nml=ic_sampling_nml)
	 ! read(unit=pfile, nml=param_sampling_nml)
	 ! read(unit=pfile, nml=params)
	 ! read(unit=pfile, nml=print_out)
	 ! read(unit=pfile, nml=reheat)
	 ! read(unit=pfile, nml=technical)
	 ! close(unit=pfile)

    !Instead the idea is to set these inital parameters to be an
    !argument for the function.
    !---------------------------------------------------------------

    param_sampling = ginput_param_sampling
    vp_prior_min = ginput_vp_prior_min
    vp_prior_max = ginput_vp_prior_max
    varying_N_pivot = ginput_varying_N_pivot
    use_first_priorval = ginput_use_first_priorval
!    numb_auxparams = ginput_numb_auxparams
!    prior_auxparams_min = ginput_prior_auxparams_min
!    prior_auxparams_max = ginput_prior_auxparams_max

    phi_init0 = ginput_phi_init0
    dphi_init0 = ginput_dphi_init0
    vparams = ginput_vparams
    N_pivot = ginput_N_pivot
    k_pivot = ginput_k_pivot
    dlnk = ginput_dlnk

    print*, "k_pivot=",k_pivot
	print*, "dlnk=",dlnk
    print*, "ginput_N_pivot=",ginput_N_pivot
    print*, "N_pivot = ",N_pivot
    print*, "ginput_vparams = ",ginput_vparams
    print*, "vparams = ",vparams
    print*, "vparams(1,:)", vparams(1,:)
    print*, "vparams(:,1)", vparams(:,1)

!    effective_V_choice = ginput_effective_V_choice
!    turning_choice = ginput_turning_choice
!    number_knots_qsfrandom = ginput_number_knots_qsfrandom
!    stand_dev_qsfrandom = ginput_stand_dev_qsfrandom
!    knot_range_min = ginput_knot_range_min
!    knot_range_max = ginput_knot_range_max
!    custom_knot_range = ginput_custom_knot_range

    !---------------------------------------------------------------
    ! setting up the printer options at the code-level.
    ! TODO: make this from GAMBIT yaml file. 
    ! NB: this may not be interesting or necessary.
    out_opt%modpkoutput = .false.
    out_opt%output_reduced = .true.
    out_opt%output_badic =.false.

    out_opt%save_traj = .true.
!    out_opt%phi0 = .false.
    out_opt%fields_horiz = .false.
    out_opt%fields_end_infl = .false.
    out_opt%spectra = .false.
    out_opt%modes = .false.

    ! out_opt%save_reheat = .false.
    ! out_opt%save_cij = .false.

!    out_opt%detailed%needmoreopts = .false.
!    out_opt%detailed%write_Npiv = .false.
!    out_opt%detailed%write_num_inflaton = .false.
!    out_opt%detailed%write_vparams = .false.
!    out_opt%detailed%write_auxparams = .false.
    !---------------------------------------------------------------


    !---------------------------------------------------------------
    ! setting up the integrater options at the code-level.
    ! TODO: make this from GAMBIT yaml file.
    assert%use_assertions = .true.

!    tech_opt%use_integ_with_t= .false.
!    tech_opt%automate_singlefield_ic = .false.

    tech_opt%accuracy_setting = 1
    tech_opt%use_dvode_integrator = .false.
!    tech_opt%use_ode_constraints = .false.

!    tech_opt%use_analytical_jacobian = .false.

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

!    tech_opt%dvode_dN_r = 0.001
!    tech_opt%dvode_dN_c = 1e-4
!    tech_opt%dvode_dt = 1e16
!---------------------------------------------------------------
!------------------My authantic settings------------------------
!    tech_opt%rk_accuracy_modes = 1.0e-7
!    tech_opt%rk_accuracy_back = 1.0e-9

!    tech_opt%dvode_rtol_back = 1.0e-9
!    tech_opt%dvode_rtol_modes = 1.0e-6

!    tech_opt%dvode_atol_back = 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14

!    tech_opt%dvode_atol_modes_real(1:2) = 1e-8 1e-8
!    tech_opt%dvode_atol_modes_real(3:4) = 1e-6 1e-6
!    tech_opt%dvode_atol_modes_real(5:8) = 1e-5 1e-5 1e-5 1e-5
!    tech_opt%dvode_atol_modes_real(9:12) = 1e-5 1e-5 1e-5 1e-5
!    tech_opt%dvode_atol_modes_real(13:14) = 1e-8 1e-8
!    tech_opt%dvode_atol_modes_real(15:16) = 1e-3 1e-3

!    tech_opt%dvode_atol_modes_imag(1:2) = 1e0 1e0
!    tech_opt%dvode_atol_modes_imag(3:4) = 1e0 1e0
!    tech_opt%dvode_atol_modes_imag(5:8) = 1e-5 1e-5 1e-5 1e-5
!    tech_opt%dvode_atol_modes_imag(9:12) = 1e-5 1e-5 1e-5 1e-5
!    tech_opt%dvode_atol_modes_imag(13:14) = 1e-8 1e-8
!    tech_opt%dvode_atol_modes_imag(15:16) = 1e-3 1e-3

!    tech_opt%dvode_dN_r = 0.001
!    tech_opt%dvode_dN_c = 1e-4
!    tech_opt%dvode_dt = 1e16
!---------------------------------------------------------------
!---------------------------------------------------------------

    !---------------------------------------------------------------
    call observs%set_zero()
    call observs_SR%set_zero()

	!Set random seed
	call init_random_seed()

    call output_initial_data()

    if (ic_sampling==ic_flags%reg_samp) then

      call out_opt%open_files(SR=use_deltaN_SR)

	  call calculate_pk_observables(observs,observs_SR,k_pivot,dlnk)

      call out_opt%close_files(SR=use_deltaN_SR)

!-------------------DEVELOPEMENT-------------------------------
!    else if (ic_sampling == &
!					 ic_flags%scanning_N_pivots) then
!      call out_opt%open_files(SR=use_deltaN_SR)
!      do sample_looper=1,numb_samples
!        call calculate_pk_observables(k_pivot,dlnk)
!	  end do
!
!      call out_opt%close_files(SR=use_deltaN_SR)
!	end if
!
!    else if &
!      !Eqen sampling
!      (ic_sampling == ic_flags%eqen_samp .or. &
!      !Set vels in SR and fields on iso-N surface for N-quad
!      ic_sampling == ic_flags%iso_N .or.&
!      !Set vels in SR
!      ic_sampling == ic_flags%slowroll_samp .or.&
!      !Grab IC from file
!      ic_sampling == ic_flags%fromfile_samp .or. &
!      !Loop over different vparams for given num_inflaton
!      ic_sampling == ic_flags%parameter_loop_samp .or. &
!      ic_sampling == ic_flags%param_unif_prior .or. &
!      !QSF trajectories
!      ic_sampling == ic_flags%qsf_random .or.  &
!      ic_sampling == ic_flags%qsf_parametric .or. &
!      !Single-field axion
!      ic_sampling == ic_flags%single_axion &
!      ) then

!      call out_opt%open_files(ICs=.true., SR=use_deltaN_SR)

      !Initialize the IC sampler
 !     call init_sampler(icpriors_min, icpriors_max)
!-------------------DEVELOPEMENT-------------------------------

  else if &
    !Eqen sampling
    (ic_sampling == ic_flags%eqen_samp .or. &
    !Set vels in SR and fields on iso-N surface for N-quad
    ic_sampling == ic_flags%iso_N .or.&
    !Set vels in SR
    ic_sampling == ic_flags%slowroll_samp) then

    call out_opt%open_files(ICs=.true., SR=use_deltaN_SR)

!Initialize the IC sampler
!    call init_sampler(icpriors_min, icpriors_max)


      do sample_looper=1,numb_samples

        if (out_opt%modpkoutput) write(*,*) &
          "---------------------------------------------"
        if (out_opt%modpkoutput) write(*,*) &
          "Sample numb", sample_looper, "of", numb_samples
        if (out_opt%modpkoutput) write(*,*) &
          "---------------------------------------------"

        call calculate_pk_observables(observs,observs_SR,k_pivot,dlnk)

      end do

      call out_opt%close_files(ICs=.true., SR=use_deltaN_SR)

!-------------------DEVELOPEMENT-------------------------------
!    else if (ic_sampling == ic_flags%iterate_reheat_only) then
!      !Iterating only over reheating params
!      !requires us to bypass some of the more typical functionality

!      call out_opt%open_files(ICs=.true.)

!      do sample_looper=1,numb_samples

!        if (out_opt%modpkoutput) write(*,*) &
!          "---------------------------------------------"
!        if (out_opt%modpkoutput) write(*,*) &
!          "Sample numb", sample_looper, "of", numb_samples
!        if (out_opt%modpkoutput) write(*,*) &
!          "---------------------------------------------"

!        if (sample_looper==1) then
!          call calculate_pk_reheating(k_pivot,firsttime=.true.)
!        else
!          call calculate_pk_reheating(k_pivot,firsttime=.false.)
!        end if

!      end do

!      call out_opt%close_files(ICs=.true.)
!-------------------DEVELOPEMENT-------------------------------

    else
      print*, "MODECODE: sampling technique=",ic_sampling
      call raise%fatal_code(&
        "This sampling technique is not implemented.",&
        __FILE__, __LINE__)

    end if

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

  contains

    subroutine gambit_get_full_pk(pk_arr,calc_full_pk,ginput_steps,ginput_kmin,ginput_kmax)
      !Find P(k) for all the power spectra from kmin to kmax, as given in the
      !parameters file.
      real(dp), intent(in) :: ginput_kmin, ginput_kmax
      integer, intent(in) :: ginput_steps

      real(dp), dimension(:,:), allocatable, intent(out) :: pk_arr

      real(dp) :: kmin, kmax, incr
      logical, intent(inout) :: calc_full_pk
      real(dp) :: p_scalar, p_tensor, p_zeta, p_iso
      real(dp), dimension(:), allocatable :: k_input
      integer :: i, steps, u

      type(power_spectra) :: pk

      ! namelist /full_pk/ kmin, kmax, steps, calc_full_pk

      ! open(newunit=u, file="parameters_multimodecode.txt", &
      ! status="old", delim = "apostrophe")
      ! read(unit=u, nml=full_pk)
      ! close(u)
	  kmin = ginput_kmin
      kmax = ginput_kmax
      steps = ginput_steps

!      !DEBUG
!      print*, "testing calc_full_pk", calc_full_pk

      !If don't want full spectrum, return
      if (.not. calc_full_pk) return

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

	end subroutine gambit_get_full_pk

    subroutine allocate_vars()
!      !Allocate all the necessary arrays that we can with the information given
!      !in the parameters file.

!      !Prepare extra params if necessary
!      if (more_potential_params) then
!        allocate(turning_choice(num_inflaton-1))
!        allocate(number_knots_qsfrandom(num_inflaton-1))
!        allocate(stand_dev_qsfrandom(num_inflaton-1))
!        allocate(knot_range_min(num_inflaton-1))
!        allocate(knot_range_max(num_inflaton-1))

!        read(unit=pfile, nml=more_params)

!      end if

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

!      allocate(auxparams(max(100,2*num_inflaton)))
!      allocate(prior_auxparams_max(max(100,2*num_inflaton)))
!      allocate(prior_auxparams_min(max(100,2*num_inflaton)))

      allocate(phi_init0(num_inflaton))
      allocate(phi_init(num_inflaton))
      allocate(phidot_sign(num_inflaton))
      allocate(phiarr(num_inflaton, 1:nsteps))
      allocate(dphiarr(num_inflaton, 1:nsteps))
      allocate(param_arr(1:nsteps))
      allocate(phi_infl_end(num_inflaton))
!      allocate(dphi_infl_end(num_inflaton))
      allocate(phi_pivot(num_inflaton))
      allocate(dphi_pivot(num_inflaton))
      allocate(dphi_init0(num_inflaton))
      allocate(dphi_init(num_inflaton))
!      allocate(dphidt_init0(num_inflaton))

    end subroutine allocate_vars

!-------------------DEVELOPEMENT-------------------------------
    subroutine deallocate_vars()
!      !Deallocate all the arrays that we allocate in allocate_vars()

!      if (allocated(turning_choice)) deallocate(turning_choice)
!      if (allocated(number_knots_qsfrandom)) deallocate(number_knots_qsfrandom)
!      if (allocated(stand_dev_qsfrandom)) deallocate(stand_dev_qsfrandom)
!      if (allocated(knot_range_min)) deallocate(knot_range_min)
!      if (allocated(knot_range_max)) deallocate(knot_range_max)
      if (allocated(vparams)) deallocate(vparams)
!      if (allocated(auxparams)) deallocate(auxparams)
      if (allocated(icpriors_max)) deallocate(icpriors_max)
      if (allocated(icpriors_min)) deallocate(icpriors_min)
      if (allocated(vp_prior_max)) deallocate(vp_prior_max)
      if (allocated(vp_prior_min)) deallocate(vp_prior_min)
!      if (allocated(prior_auxparams_max)) deallocate(prior_auxparams_max)
!      if (allocated(prior_auxparams_min)) deallocate(prior_auxparams_min)

      if (allocated(phi_init0)) deallocate(phi_init0)
      if (allocated(phi_init)) deallocate(phi_init)
      if (allocated(phidot_sign)) deallocate(phidot_sign)
      if (allocated(phiarr)) deallocate(phiarr)
      if (allocated(dphiarr)) deallocate(dphiarr)
      if (allocated(param_arr)) deallocate(param_arr)
      if (allocated(phi_infl_end)) deallocate(phi_infl_end)
!      if (allocated(dphi_infl_end)) deallocate(dphi_infl_end)
      if (allocated(phi_pivot)) deallocate(phi_pivot)
      if (allocated(dphi_pivot)) deallocate(dphi_pivot)
      if (allocated(dphi_init0)) deallocate(dphi_init0)
      if (allocated(dphi_init)) deallocate(dphi_init)
!      if (allocated(dphidt_init0)) deallocate(dphidt_init0)

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
!      if (potential_choice == 14 .or. potential_choice==15) then
!        write(*, out_opt%i_fmt) &
!          "Turning Choice =", turning_choice
!      end if
      write(*, out_opt%e_fmt) &
        "N_pivot =", N_pivot
!      if (out_opt%output_reduced) then
!        write(*, out_opt%e_fmt) &
!          "phi_pivot =", phi_pivot(1:min(10,num_inflaton))
!      else
!        write(*, out_opt%e_fmt) &
!          "phi_pivot =", phi_pivot(:)

!      end if
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
!        write(*, out_opt%e2_fmt),&
!          "Bundle Width =", field_bundle%exp_scalar

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
    subroutine calculate_pk_observables(observs,observs_SR,k_pivot,dlnk)

      real(dp), intent(in) :: k_pivot,dlnk

      type(observables), intent(inout) :: observs, observs_SR

      real(dp), dimension(:,:), allocatable :: pk_arr
      logical :: calc_full_pk, leave

      type(power_spectra) :: pk0, pk1, pk2, pk3, pk4

      character(1024) :: cname
      integer :: ii

      call observs%set_zero()
      call observs_SR%set_zero()

      pk_bad = run_outcome%success
      leave = .false.

!-------------------DEVELOPEMENT-------------------------------
	  print*,"we are inside calculate_pk_observables."
!      !Get e-folds after pivot scale leaves horizon
!      if (varying_N_pivot) then
!        save_iso_N = .false.
!        call get_new_N_pivot(N_pivot,&
!          N_pivot_prior_min, N_pivot_prior_max)
!      end if

!      !num_inflaton sets many array sizes, have to treat it slightly differently
!      if (varying_num_inflaton) call get_new_num_inflaton()

      !Get vparams
      if (param_sampling == 1) then
	    call get_vparams()
      end if

      !Get ICs
      if (ic_sampling/=ic_flags%reg_samp) then
	    print*,"we are inside ic_sampling condition"
        call get_ic(phi_init0, dphi_init0, &
          icpriors_min, icpriors_max, &
          numb_samples,energy_scale)
      end if

      !Load ics
      allocate(observs%ic(2*num_inflaton))
      observs%ic(1:num_inflaton)=phi_init0
      observs%ic(num_inflaton+1:2*num_inflaton)=dphi_init0
      print*,"dphi_init0 = ", dphi_init0
      if (use_deltaN_SR) then
        allocate(observs_SR%ic(2*num_inflaton))
        observs_SR%ic = observs%ic
      end if

!      !DEBUG
!      !Estimate # of e-folds vs N_pivot
!      if (potential_choice==1) then
!        if (0.25e0_dp*sum(phi_init0**2)<0.75e0_dp*N_pivot) then
!          pk_bad = 2 !Not enough inflation for pivot scale to leave horizon
!        end if
!      end if

!      call test_bad(pk_bad, observs, leave)
!      if (leave) return

      print*,"N_pivot = ", N_pivot

      !Initialize potential and calc background
      call potinit

      print*,"outside of potinit"


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

      print*,"observs_SR%ns = ",observs_SR%ns

      !Evaluate the mode functions
      if (evaluate_modes) then

        call evolve(k_pivot, pk0)
          call test_bad(pk_bad, observs, leave)
          if (leave) then
			return
		  end if

!        if (reheat_opts%use_reheat) then
!          !Set observables from the sudden decay calculation
!          observs%As = reheater%observs%As
!          observs%ns = reheater%observs%ns
!          observs%r  = reheater%observs%r
!        else
!DEBUG
!print*, "Not evaluating second and third evolve routines"
!stop
	      print*, "second evolve"
		  print*, "k_pivot*exp(-dlnk) = " , k_pivot*exp(-dlnk)
		  print*, "dlnk = ", dlnk

        call evolve(k_pivot*exp(-dlnk), pk1)
          call test_bad(pk_bad, observs, leave)
          if (leave) return
		print*, "left second evolve"
		print*, "third evolve"
		print*, "k_pivot*exp(dlnk) = " , k_pivot*exp(dlnk)
        call evolve(k_pivot*exp(dlnk), pk2)
          call test_bad(pk_bad, observs, leave)
          if (leave) return
		print*, "left third evolve"

        if (get_runningofrunning) then
			print*, "get_runningofrunning = True"
            !Alpha_s from 5-pt stencil
            !or running of running
          call evolve(k_pivot*exp(-2.0e0_dp*dlnk), pk3)
          print*, "fourth evolve"
            call test_bad(pk_bad, observs, leave)
            if (leave) return
          call evolve(k_pivot*exp(2.0e0_dp*dlnk), pk4)
		  print*, "fifth evolve"
            call test_bad(pk_bad, observs, leave)
            if (leave) return
        end if

          !Construct the observables
          if (get_runningofrunning) then
            call observs%set_finite_diff(dlnk, &
              pk0,pk1,pk2,pk3,pk4, &
              field_bundle%exp_scalar)
            print*,"here we calculate the observables!"
          else
			print*,"here we calculate the observables!"
            call observs%set_finite_diff(dlnk, &
              pk0,pk1,pk2,&
              bundle_width=field_bundle%exp_scalar)
          end if


          !Get full spectrum for adiab and isocurv at equal intvs in lnk
!		  print*,"will call gambit_get_full_pk"
          ! call gambit_get_full_pk(pk_arr,calc_full_pk,steps,kmin,kmax)
!		  print*,"endof get_full_pk"

!        end if

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

      !Print output array
      !Only get the SR arrays if use_deltaN_SR
      !Only print observs if evaluated modes
!      if (out_opt%output_badic .or. pk_bad==run_outcome%success) then
!        print*,"printing output array "
!
!		if (evaluate_modes) then
!          if (out_opt%first_outsamp) then
!            call observs%print_header(out_opt%outsamp)
!            out_opt%first_outsamp = .false.
!          end if
!          call observs%printout(out_opt%outsamp)
!        end if
!
!        if (use_deltaN_SR) then
!          if (out_opt%first_outsamp_SR) then
!            call observs%print_header(out_opt%outsamp_SR)
!            out_opt%first_outsamp_SR=.false.
!          end if
!          call observs_SR%printout(out_opt%outsamp_SR)
!        end if
!
!      end if

      print*,"save_iso_N=",save_iso_N

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

     print*,"no issue before here!"

    end subroutine calculate_pk_observables

!-------------------DEVELOPEMENT-------------------------------
!    !Calculate observables when only reheating parameters are being altered
!    !each time this routine is called.  Loads background and mode solutions
!    !from file
!    subroutine calculate_pk_reheating(k_pivot,firsttime)

!      real(dp), intent(in) :: k_pivot
!      real(dp), dimension(:,:), allocatable :: pk_arr
!      logical :: calc_full_pk, leave
!      logical, intent(in) :: firsttime

!      type(power_spectra) :: pk0, pk1, pk2, pk3, pk4
!      type(observables) :: observs, observs_SR

!      character(1024) :: cname
!      integer :: ii

!      !Get previously tabulated data from file
!      call reheater%load_from_file('out_reheaterfile.csv',firsttime)

!      call observs%set_zero()
!      call observs_SR%set_zero()

!      pk_bad = run_outcome%success
!      leave = .false.

!      calc_full_pk = .false.

!      !Don't vary anything else!

!      !Evaluate the post reheating epoch
!      call reheat_evolve()
!        call test_bad(pk_bad, observs, leave)
!        if (leave) return

!      observs%As = reheater%observs%As
!      observs%ns = reheater%observs%ns
!      observs%r  = reheater%observs%r

!      !Write output to stdout
!      if (out_opt%modpkoutput) then
!        if (evaluate_modes) then
!          call output_observables(pk_arr, &
!            calc_full_pk, observs, observs_SR)
!        else
!          call output_observables(pk_arr, &
!            calc_full_pk, observ_SR = observs_SR)
!        end if
!      end if


!      !Print output array
!      !Only get the SR arrays if use_deltaN_SR
!      !Only print observs if evaluated modes
!      if (out_opt%output_badic .or. pk_bad==run_outcome%success) then


!        if (evaluate_modes) then
!          if (out_opt%first_outsamp) then
!            call observs%print_header(out_opt%outsamp)
!            out_opt%first_outsamp = .false.
!          end if
!          call observs%printout(out_opt%outsamp)
!        end if

!        if (use_deltaN_SR) then
!          if (out_opt%first_outsamp_SR) then
!            call observs%print_header(out_opt%outsamp_SR)
!            out_opt%first_outsamp_SR=.false.
!          end if
!          call observs_SR%printout(out_opt%outsamp_SR)
!        end if

!      end if

!    end subroutine calculate_pk_reheating
!-------------------DEVELOPEMENT-------------------------------

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
!      observs_SR%A_bundle = field_bundle%exp_scalar

!	  reheater%phi_pivot = phi_pivot
!	  print*, "phi_pivot given as = " , phi_pivot

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

!    subroutine init_sampler(icpriors_min, icpriors_max,read_nml2)
!    subroutine gambit_init_sampler(icpriors_min, icpriors_max, &
!                 ginput_phi0_priors_min,ginput_dphi0_priors_min, &
!				 ginput_dphi0_priors_max,ginput_N_pivot_prior_min, &
!                ginput_N_pivot_prior_max)
!      !Initialize the parameter or IC sampler if we're using it.
!      !Need to allocate some arrays and set the random seed.
!
!      use modpk_rng, only : init_random_seed
!
!      real(dp), dimension(:), intent(in) :: ginput_phi0_priors_min
!      real(dp), dimension(:), intent(in) :: ginput_phi0_priors_max
!	  real(dp), dimension(:), intent(in) :: ginput_dphi0_priors_min
!      real(dp), dimension(:), intent(in) :: ginput_dphi0_priors_max
!      real(dp), intent(in) :: ginput_N_pivot_prior_min
!	  real(dp), intent(in) :: ginput_N_pivot_prior_max
!
!      real(dp), dimension(:,:), intent(out) :: icpriors_min, &
!        icpriors_max

!      real(dp), dimension(:), allocatable :: phi0_priors_min, &
!        dphi0_priors_min, phi0_priors_max, dphi0_priors_max
!      real(dp) ::  N_pivot_prior_max
!      real(dp) ::  N_pivot_prior_min
!
!      integer :: u, i

!      logical, intent(in), optional :: read_nml2

!      namelist /priors/ phi0_priors_min, phi0_priors_max, &
!        dphi0_priors_min, dphi0_priors_max, &
!        N_pivot_prior_min, N_pivot_prior_max, allow_superplanckian

!      if (allocated(phi0_priors_max)) deallocate(phi0_priors_max)
!      if (allocated(dphi0_priors_max)) deallocate(dphi0_priors_max)
!      if (allocated(phi0_priors_min)) deallocate(phi0_priors_min)
!      if (allocated(dphi0_priors_min)) deallocate(dphi0_priors_min)

!      allocate(phi0_priors_max(num_inflaton))
!      allocate(dphi0_priors_max(num_inflaton))
!      allocate(phi0_priors_min(num_inflaton))
!      allocate(dphi0_priors_min(num_inflaton))
!      phi0_priors_min=0e0_dp
!      phi0_priors_max=0e0_dp
!      dphi0_priors_min=0e0_dp
!      dphi0_priors_max=0e0_dp

!      if (save_iso_N) then
!        if (allocated(phi_iso_N)) deallocate(phi_iso_N)
!        if (allocated(dphi_iso_N)) deallocate(dphi_iso_N)
!        allocate(phi_iso_N(num_inflaton))
!        allocate(dphi_iso_N(num_inflaton))
!      end if

      !When varying num_inflaton this will be useful to stop
      !segfaults when reading from the parameter file
 !     if (present(read_nml2) .and. .not. read_nml2) return

!      !Read phi0 priors from file
!	    open(newunit=u, file="parameters_multimodecode.txt", &
!        status="old", delim = "apostrophe")
!      read(unit=u, nml=priors)
!      close(u)

!      phi0_priors_max = ginput_phi0_priors_max
!      phi0_priors_min = ginput_phi0_priors_min
!      dphi0_priors_max = ginput_dphi0_priors_max
!      dphi0_priors_min = ginput_dphi0_priors_min
!      N_pivot_priors_min = ginput_N_pivot_prior_min
!      N_pivot_priors_max = ginput_N_pivot_prior_max

!      icpriors_max(1,:) = phi0_priors_max
!      icpriors_max(2,:) = dphi0_priors_max
!      icpriors_min(1,:) = phi0_priors_min
!      icpriors_min(2,:) = dphi0_priors_min
!
!    end subroutine gambit_init_sampler

    !*********************
    !LIMITED FUNCTIONALITY
    !*********************
    !Gets a new number of inflatons and resets arrays that have already been allocated.
    !Checks for conflicts with other parameter settings.
!    subroutine get_new_num_inflaton()
!       real(dp) :: rand, p_val

       !Checks
!       if (ic_sampling==ic_flags%reg_samp .or. &
!         !ic_sampling==ic_flags%eqen_samp .or. &
!         ic_sampling==ic_flags%parameter_loop_samp .or. &
!         ic_sampling==ic_flags%param_unif_prior .or. &
!         ic_sampling==ic_flags%qsf_random .or. &
!         ic_sampling==ic_flags%qsf_parametric .or. &
!         ic_sampling==ic_flags%slowroll_samp) then


!         print*, "MODECODE: sampling technique=",ic_sampling
!         call raise%fatal_code(&
!           "This sampling technique is not safe when varying num_inflaton.",&
!           __FILE__, __LINE__)
!       end if

!       !Get rid of arrays you've just allocated
!       call deallocate_vars()

!       !Uniform sampling for num_inflaton
!       if (inflaton_sampling==inflaton_flags%unif) then

!         call random_number(rand)
!         call random_number(rand)
!         rand = rand*(num_inflaton_prior_max-num_inflaton_prior_min) + num_inflaton_prior_min

!         num_inflaton = floor(rand)

       !Uniform sampling for log(num_inflaton)
!       else if (inflaton_sampling==inflaton_flags%logunif) then
!         call random_number(rand)
!         call random_number(rand)
!         rand = rand*&
!           (log10(real(num_inflaton_prior_max,dp))-&
!           log10(real(num_inflaton_prior_min,dp))) &
!           + log10(real(num_inflaton_prior_min,dp))

!         num_inflaton = floor(10.0e0_dp**rand)

!       else if (inflaton_sampling==inflaton_flags%geometric) then
!       !Geometric sampling for num_inflaton, which is max ent if mean is known

!         p_val = 1.0e0_dp/num_inflaton_prior_max !Using n.._prior_min=mean(num_inflaton)
!         num_inflaton = geometric(p_val)

!       else

!         print*, "MODECODE: inflaton_sampling technique=",inflaton_sampling
!         call raise%fatal_code(&
!           "This sampling technique is not implemented.",&
!           __FILE__, __LINE__)

!       end if

!       !Remake the necessary arrays
!       call allocate_vars()

       !You have to reinitialize the sampling priors intelligently.
       !Can't just call init_sampler because param file will likely depend on a
       !set number of fields
       !This can go wrong very easily...
!       if (ic_sampling == ic_flags%iso_N) then
!	       open(newunit=pfile, file="parameters_multimodecode.txt", &
!           status="old", delim = "apostrophe")
!         read(unit=pfile, nml=param_sampling_nml)
!	       close(unit=pfile)

!       else if (ic_sampling == ic_flags%eqen_samp) then

!	       open(newunit=pfile, file="parameters_multimodecode.txt", &
!           status="old", delim = "apostrophe")
!         read(unit=pfile, nml=param_sampling_nml)
!	       close(unit=pfile)

!         call init_sampler(icpriors_min, icpriors_max, read_nml2=.false.)

         !Fixing field ranges by Planck mass
!         icpriors_max(1,:) = 1.0e0_dp
!         icpriors_min(1,:) = -1.0e0_dp

         !Ignoring KE
!         icpriors_min(2,:) = 0e0_dp
!         icpriors_max(2,:) = 0e0_dp


!       else
!         print*, "MODECODE: sampling technique=",ic_sampling
!         call raise%fatal_code(&
!           "This sampling technique is not implemented with variable num_inflaton.",&
!           __FILE__, __LINE__)
!       end if


!       !Get new formatting, etc
!       call output_initial_data()

!    end subroutine get_new_num_inflaton

  end subroutine multimodecode_gambit_driver

end module multimodecode_gambit
