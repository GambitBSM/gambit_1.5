MODULE camb_interface
  !Simple module to tell an external program, such as CAMB or the
  !multimodecode_driver, that a given set of parameters does not give an
  !appropriate inflationary realization.
  IMPLICIT NONE
  INTEGER :: pk_bad
  LOGICAL :: pk_initialized
END MODULE camb_interface

MODULE modpkparams
  !Module defining many "global" variables that various cosmology portions of
  !the code will need access to.  A variable added here can be seen in most
  !places.

  IMPLICIT NONE

  !Double precision.
  integer, parameter :: dp = selected_real_kind(15, 307)

  !Quad precision; remember to compile with -r16
  !INTEGER, parameter :: DP = selected_real_kind(33, 4931)

  LOGICAL :: use_modpk, vnderivs, instreheat

! increase max_vparams to use more potential parameters
  INTEGER*4, parameter :: max_vparams = 9

  INTEGER :: potential_choice

  INTEGER*4 :: nactual_bg, nactual_mode
  INTEGER, PARAMETER :: nsteps=1e5
  real(dp), PARAMETER :: M_Pl=1.0e0_dp
  real(dp), PARAMETER :: Mpc2Mpl=2.6245e-57_dp
  real(dp) :: k_pivot, N_pivot, N_tot, H_pivot
  real(dp) :: a_end, a_pivot
  real(dp) :: a_init
  real(dp) :: h_init, rescale_factor

  real(dp), ALLOCATABLE :: phidot_sign(:)
  real(dp) :: Nefold_max=100000.e0_dp
  real(dp) :: lna(nsteps)
  real(dp) :: hubarr(nsteps), log_aharr(nsteps), epsarr(nsteps), dtheta_dN(nsteps)
  LOGICAL :: slowroll_infl_end
  LOGICAL :: slowroll_start=.false.

  !MULTIFIELD
  integer :: num_inflaton
  real(dp), dimension(:,:), allocatable :: vparams
  real(dp), allocatable :: phi_init0(:), phi_init(:)
  real(dp), allocatable :: dphi_init0(:), dphi_init(:)
  real(dp), allocatable:: phi_pivot(:), dphi_pivot(:), phi_infl_end(:)

  real(dp), allocatable :: phiarr(:,:), dphiarr(:,:) !The first index is the multifield index
  real(dp), allocatable :: param_arr(:)
  real(dp) :: sigma_arr(nsteps)
  real(dp) :: delsigma = M_Pl      !specify total field distance travelled before inflation ends
  !END MULTIFIELD

  real(dp) :: findiffdphi

  real(dp) :: modpk_ns, modpk_nt, modpk_nrun, modpk_As, modpk_r

  !Flags for analytical calculations
  logical :: use_deltaN_SR
  logical :: evaluate_modes

  !Technical options
  type :: tech_options
    integer :: accuracy_setting
    logical :: use_dvode_integrator
    real(dp) :: rk_accuracy_modes, rk_accuracy_back
    real(dp) :: dvode_rtol_modes, dvode_rtol_back
    real(dp), dimension(10000) :: dvode_atol_modes_real, dvode_atol_modes_imag, dvode_atol_back
  end type tech_options

  type(tech_options) :: tech_opt

END MODULE modpkparams


MODULE ode_path
  !Module for control parameters for the numerical integration of either the
  !background equations or the mode equations.
  use modpkparams, only : dp
  implicit none

  INTEGER*4 :: nok,nbad,kount
  LOGICAL, SAVE :: save_steps=.false.
  LOGICAL :: ode_underflow
  LOGICAL :: ode_ps_output
  LOGICAL :: ode_infl_end
  LOGICAL :: infl_ended
  real(dp) :: dxsav
  real(dp), DIMENSION(:), POINTER :: xp
  real(dp), DIMENSION(:), POINTER :: param_p
  real(dp), DIMENSION(:,:), POINTER :: yp
END MODULE ode_path


MODULE internals
  !Module that defines some variables for internal use in the code.
  use modpkparams, only : num_inflaton, dp
  IMPLICIT NONE
  real(dp), PARAMETER :: PI=3.141592653589793238462643383279502884197e0_dp
  real(dp) :: h_ik
  !MULTIFIELD
  integer :: index_ptb_y, index_ptb_vel_y, index_tensor_y, index_uzeta_y
  !END MULTIFIELD
  real(dp) :: k, a_ik


END MODULE internals


module modpk_observables
  !Module that defines the various observables one could calculate around the
  !pivot scale.  Defines objects for observables, power spectra, etc.
  use modpkparams, only : dp, num_inflaton
  use modpk_io, only : out_opt
  use csv_file, only : csv_write
  implicit none

  integer*4 :: ik
  real(dp) :: eval_ps,k_start, useq_ps

  !Power spectrum type, used for one k
  !Simplifies all the various defns for isocurv ptbs
  type :: power_spectra
    !Mode
    real(dp) :: k
    !Ptb spectra
    complex(dp), dimension(:,:), allocatable :: phi_ij
    !Proj ptb spectra onto adiab direction
    real(dp) :: adiab
    !Tensor mode spectrum
    real(dp) :: tensor
    !Approximated zeta spectrum
    real(dp) :: powz
    !Proj ptb spectra onto directions perpend to adiab
    real(dp) :: isocurv
    !Cross-spectrum for isocurv + adiab
    real(dp) :: cross_ad_iso
    !Total non-adiab pressure ptb spectrum
    real(dp) :: pnad
    !Total pressure ptb spectrum
    real(dp) :: pressure
    !Adiab pressure ptb spectrum
    real(dp) :: press_ad
    !Total entropic ptb spectrum, proportional to non-adiab pressure ptb
    real(dp) :: entropy
    !Expansion scalar for field space bundle width
    real(dp) :: bundle_exp_scalar

  end type

  !For temporary calc of spectra in odeint
  type(power_spectra) :: power_internal

  type, public :: KahanSum
    real(dp) :: summand = 0e0_dp
    real(dp) :: remainder = 0e0_dp
    contains
      procedure, public :: add => kahan_summation
      procedure, public :: clear => kahan_clear_mem
  end type

  !Type to save the ICs and observs. Add new observables here
  type :: observables
    real(dp), dimension(:), allocatable :: ic
    !Spectra amplitudes
    real(dp) :: As
    real(dp) :: A_iso
    real(dp) :: A_pnad
    real(dp) :: A_ent
    real(dp) :: A_cross_ad_iso
    !Bundle width from arXiv:1203.2635
    real(dp) :: A_bundle
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
    contains
      procedure, public :: printout => ic_print_observables
      procedure, public :: print_header => ic_print_headers
      procedure, public :: set_zero => set_observs_to_zero
      procedure, public :: set_finite_diff => calculate_observs_finitediff
  end type observables

  private :: ic_print_observables, set_observs_to_zero, &
    calculate_observs_finitediff,kahan_summation,kahan_clear_mem

  contains

    !Procedures for Kahan summation objects

    !Algorithm for Kahan summation
    subroutine kahan_summation(self,val)

      class(KahanSum) :: self
      real(dp), intent(in) :: val
      real(dp) :: ytemp, ttemp

      ytemp = val - self%remainder
      ttemp = self%summand + ytemp
      self%remainder = (ttemp - self%summand) - ytemp
      self%summand = ttemp

    end subroutine kahan_summation

    subroutine kahan_clear_mem(self)

      class(KahanSum) :: self

      self%remainder=0e0_dp
      self%summand=0e0_dp

    end subroutine kahan_clear_mem

    !Write the cosmo observables headers to file
    subroutine ic_print_headers(self, outunit)
      class(observables) :: self
      integer, intent(in) :: outunit

      character(1024) :: cname
      integer :: ii

      !First columns for IC
      do ii=1,size(self%ic)
        write(cname, "(A3,I4.4)") "phi_piv", ii
        call csv_write(&
          outunit,&
          trim(cname), &
          advance=.false.)
      end do

      !Remaining columns
      call csv_write(outunit,&
        (/character(len=10) ::&
          'As', 'ns', 'r', 'nt', 'alpha_s', &
          'A_iso', 'A_pnad', 'A_ent', 'A_bundle', &
          'n_iso', 'n_pnad', 'n_ent', &
          'A_cross', 'f_NL', 'tau_NL'/),&
        advance=.true.)

    end subroutine ic_print_headers


    !Print the cosmo observables to file
    subroutine ic_print_observables(self, outunit)

      class(observables) :: self
      integer, intent(in) :: outunit


      call csv_write(outunit,&
        (/ self%ic(:), &
        self%As, &
        self%ns,&
        self%r, &
        self%nt, &
        self%alpha_s, &
        self%A_iso, &
        self%A_pnad,&
        self%A_ent, &
        self%A_bundle, &
        self%n_iso, &
        self%n_pnad, &
        self%n_ent, &
        self%A_cross_ad_iso, &
        self%f_NL, &
        self%tau_NL /), &
        advance=.true.)

    end subroutine ic_print_observables

    subroutine set_observs_to_zero(self)
      class(observables) :: self

        self%As = 0e0_dp
        self%A_iso = 0e0_dp
        self%A_pnad = 0e0_dp
        self%A_ent = 0e0_dp
        self%A_cross_ad_iso = 0e0_dp
        self%A_bundle = 0e0_dp
        self%ns = 0e0_dp
        self%nt = 0e0_dp
        self%n_iso = 0e0_dp
        self%n_pnad = 0e0_dp
        self%n_ent = 0e0_dp
        self%r = 0e0_dp
        self%alpha_s = 0e0_dp
        self%runofrun = 0e0_dp
        self%f_NL = 0e0_dp
        self%tau_NL = 0e0_dp

    end subroutine set_observs_to_zero

    subroutine calculate_observs_finitediff(self, dlnk, &
        pk0, pklow1, pkhigh1, &
        pklow2, pkhigh2, &
        bundle_width)
      class(observables) :: self
      type(power_spectra), intent(in) :: pk0, pklow1, pkhigh1
      type(power_spectra), intent(in), optional :: pklow2, pkhigh2
      real(dp), intent(in) :: dlnk
      real(dp), intent(in), optional :: bundle_width
      logical :: runofrun

      if (present(pklow2)) then
        runofrun = .true.
      else
        runofrun = .false.
      endif

      !Amplitudes
      self%As = pk0%adiab
      self%A_iso=pk0%isocurv
      self%A_pnad=pk0%pnad
      self%A_ent=pk0%entropy
      self%A_cross_ad_iso = pk0%cross_ad_iso

      !Bundle width
      self%A_bundle=bundle_width

      !Finite difference evaluation of spectral indices
      self%ns = 1.e0_dp+log(pkhigh1%adiab/pklow1%adiab)/dlnk/2.e0_dp
      self%nt = log(pkhigh1%tensor/pklow1%tensor)/dlnk/2.e0_dp
      self%n_iso=log(pkhigh1%isocurv/pklow1%isocurv)/dlnk/2.e0_dp
      self%n_pnad=log(pkhigh1%pnad/pklow1%pnad)/dlnk/2.e0_dp
      self%n_ent=log(pkhigh1%entropy/pklow1%entropy)/dlnk/2.e0_dp

      !Tensor-to-scalar
      self%r = pk0%tensor/pk0%adiab

      if (runofrun) then

        !alpha_s from 5-pt stencil
        self%alpha_s = (1.0e0_dp/12.0e0_dp/dlnk**2)*&
          (-log(pkhigh2%adiab) + 16.0e0_dp*log(pkhigh1%adiab) - &
          30.0e0_dp*log(pk0%adiab) + 16.0e0_dp*log(pklow1%adiab) - &
          log(pklow2%adiab))

        self%runofrun = (1.0e0_dp/2.0e0_dp/dlnk**3)*&
          (log(pkhigh2%adiab) -2.0e0_dp* log(pkhigh1%adiab) &
          + 2.0e0_dp*log(pklow1%adiab) -log(pklow2%adiab))

      else

        self%alpha_s = log(pkhigh1%adiab*pklow1%adiab/pk0%adiab**2)/dlnk**2

        !Default
        self%runofrun = 0e0_dp

      end if

    end subroutine calculate_observs_finitediff

end module modpk_observables
