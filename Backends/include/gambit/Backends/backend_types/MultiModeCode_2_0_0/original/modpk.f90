MODULE access_modpk
  use modpkparams, only : dp
  USE camb_interface
  use modpk_io, only : out_opt
  use modpk_errorhandling, only : raise, run_outcome
  implicit none
  private
  public :: potinit, evolve, total_efold

! Number of k values for computing spline of P(k).
! Set to lower numbers for smoother potentials.
  integer*4, parameter, public :: pkspline_n = 500

  real(dp), parameter, public :: pkspline_kmin = log(1.e-5_dp), pkspline_kmax = log(5.e0_dp)
  real(dp), public :: pkspline_k(pkspline_n), pkspline_p(pkspline_n), &
	pkspline_p2der(pkspline_n), pkspline_pt(pkspline_n), &
	pkspline_pt2der(pkspline_n)

CONTAINS

  SUBROUTINE potinit
    USE modpkparams
    USE modpk_observables
    USE background_evolution, ONLY : backgrnd
    USE potential, ONLY : initialphi
    use modpk_icsampling, only : ic_flags, ic_sampling
    IMPLICIT NONE


    !
    !     Solve the background equations
    !
    pk_bad = run_outcome%success
    phi_init = initialphi(phi_init0)

    !NB: For eqen sampling, dphi_init set in trial_background
    CALL backgrnd


  END SUBROUTINE potinit

  SUBROUTINE total_efold
    USE modpkparams
    USE background_evolution, ONLY : backgrnd_efold
    USE potential, ONLY : initialphi
    IMPLICIT NONE

    pk_bad = run_outcome%success
    phi_init = initialphi(phi_init0)
    CALL backgrnd_efold
    RETURN
  END SUBROUTINE total_efold

  SUBROUTINE evolve(kin, powerspectrum_out)
    USE modpk_odeint
    USE ode_path
    USE modpkparams
    USE internals
    USE modpk_observables
    USE potential, ONLY: pot,powerspectrum, dVdphi, getH, getdHdalpha, field_bundle, getEps, &
      pot, d2Vdphi2
    USE modpk_utils, ONLY : derivs, qderivs, rkqs_c
    use modpk_numerics, only : locate, polint, array_polint
    use modpk_icsampling, only : ic_sampling, ic_flags
    IMPLICIT NONE

    type(power_spectra), intent(out) :: powerspectrum_out

    integer*4 :: i,j
    real(dp) :: accuracy,h1,hmin,x1,x2
    complex(kind=dp), dimension(2*num_inflaton + 2*(num_inflaton**2)+4) :: y
    real(dp), INTENT(IN) :: kin
    real(dp) :: pow_isocurvature
    real(dp) :: pow_adiabatic,  powt, powz
    real(dp), dimension(:,:), allocatable :: power_matrix
    real(dp) :: dum, ah, alpha_ik, dalpha, dh
    real(dp), DIMENSION(num_inflaton) :: p_ik,delphi
    real(dp), DIMENSION(num_inflaton) :: dp_ik

    character(36) :: e2_fmt = '(a25, es12.4, a3, es11.4, a1)'

    ! the psi portion of the y-vector is 1:n=psi_1(1:n) and
    ! 2n:3n=psi_2(1:n), etc.

    !     x = alpha     e-folds
    ! Background
    !     y(1:n) = phi                 dydx(1:n)=dphi/dalpha
    !     y(n+1:2n) = dphi/dalpha      dydx(n+1:2n)=d^2phi/dalpha^2

    ! Mode matrix ptb, psi_IJ
    !     y(2n+1:2n+n**2) = psi               dydx(2n+1:3n)=dpsi/dalpha
    !     y(2n+n**2+1:2n+2n**2) = dpsi/dalpha       dydx(3n+1:4n)=d^2psi/dalpha^2

    ! Tensors
    !     y(2n+2n**2+1) = v                  dydx(4n+1)=dv/dalpha
    !     y(2n+2n**2+2) = dv/dalpha          dydx(4n+2)=d^2v/dalpha^2

    ! --- u_zeta is the adiabatic mode ignoring coupings to other modes, used to compare with the full zeta perturbation
    !! --- full_zeta - u_zeta gives the super-horizon evolution
    !     y(2n+2n**2+3) = u_zeta     dydx(4n+4)=d^2u_zeta/dalpha^2
    !     y(2n+2n**2+4) = du_zeta/dalpha     dydx(4n+4)=d^2u_zeta/dalpha^2

    ! Set aliases for indices for above
    index_ptb_y = 2*num_inflaton+1
    index_ptb_vel_y = 2*num_inflaton+1+num_inflaton**2
    index_tensor_y = 2*num_inflaton+2*num_inflaton**2+1
    index_uzeta_y = index_tensor_y + 2

    ! Make the powerspectrum array.
    if (allocated(powerspectrum_out%phi_ij)) then
      deallocate(powerspectrum_out%phi_ij)
    end if
    if (allocated(power_internal%phi_ij)) then
      deallocate(power_internal%phi_ij)
    end if
    allocate(power_internal%phi_ij(num_inflaton,num_inflaton))
    allocate(powerspectrum_out%phi_ij(num_inflaton,num_inflaton))

    !Evaluation scale
    k=kin*Mpc2Mpl
    powerspectrum_out%k=k

    !When to start evaluating P(k), k<aH/eval_ps
    if (num_inflaton==1) then
      eval_ps = 5.0e2_dp
    else
      eval_ps = 1.0e0_dp
    end if
    useq_ps = 1.0e2_dp !When switch variables to q=\delta \phi (k<aH/useq_ps)

    !How far inside the horizon to set the modes' (Bunch-Davies) IC; k = k_start*aH
    call set_consistent_BD_scale(k_start)

    !! start where k = k_start* aH, deep in the horizon, ah = log(aH)
    ah=LOG(k/k_start)
    i= locate(log_aharr(1:nactual_bg), ah)

    IF(i.eq.0.) THEN
       PRINT*,'MODECODE: The background solution worked,'
       PRINT*,'MODECODE: but the k you requested', k,' is outside'
       PRINT*,'MODECODE: the bounds of the background you solved for.'

       !Override the stop.
       if (ic_sampling/=ic_flags%reg_samp) then
         pk_bad = run_outcome%cant_set_modeIC
         return
       end if

       write(*,e2_fmt) "log(k/k_start):", ah
       write(*,e2_fmt) "log_aharr(1):", log_aharr(1)

       call raise%fatal_cosmo(&
         'Reconsider your phi_init and N_pivot combo.', &
         __FILE__, __LINE__)

    END IF

    ! nactual_bg here is set by the background evolution
    j=MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)

    !MULTIFIELD
    CALL array_polint(log_aharr(j:j+4), phiarr(:,j:j+4), ah, p_ik, delphi)
    CALL array_polint(log_aharr(j:j+4), dphiarr(:,j:j+4), ah, dp_ik, delphi)
    !END MULTIFIELD

    CALL polint(log_aharr(j:j+4), lna(j:j+4), ah,  alpha_ik, dalpha)
    CALL polint(log_aharr(j:j+4), hubarr(j:j+4), ah,  h_ik, dh)

    a_ik=exp(alpha_ik)*a_init

    x1=alpha_ik

    IF(x1.le.0.) THEN

      print*, "MODECODE: alpha=",x1
      print*, "MODECODE: dalpha", dalpha

      call raise%fatal_cosmo(&
         'The phi_init you specified is too small to give &
         sufficient efolds of inflation. We cannot self-consistently &
         solve this for you. Please adjust phi_init and try again.', &
         __FILE__, __LINE__)

    END IF

    IF((sqrt(dot_product(delphi,delphi))/num_inflaton) .GT. 0.1 .OR. dalpha .GT. 0.1 .OR. dh .GT. 0.1) THEN
       if ((sqrt(dot_product(delphi,delphi))/num_inflaton) .GT. 0.1)&
         print*, "MODECODE: Error in dphi interpolation.",&
         (sqrt(dot_product(delphi,delphi))/num_inflaton)
       if (dalpha .GT. 0.1) print*, "MODECODE: Error in alpha interpolation.", dalpha
       if (dh > 0.1) print*, "MODECODE: Error in Hubble interpolation", dh

       call raise%fatal_code(&
          'The interpolation in SUBROUTINE evolve &
          has suspiciously large errors.  &
          Your model smells fishy.',&
          __FILE__, __LINE__)

    ENDIF

    call set_background_and_mode_ic(y)

    power_internal = powerspectrum_out

    !     Call the integrator
    !
    ode_underflow = .false.
    ode_ps_output = .true.
    ode_infl_end = .true.

    save_steps = .false.

    pk_bad = run_outcome%success

    !MULTIFIELD, need to evolve towards the end of inflation
    x2 = lna(nactual_bg) + 5.e0_dp
    !END MULTIFIELD


    !h1=0.1 !guessed start stepsize
    h1=1e-5 !guessed start stepsize

    !Some fast-roll cases need high accuracy; activate conditionally in odeint_c
    if (tech_opt%accuracy_setting==2) then
      !has a big impact on the speed of the code
      accuracy=1.0e-7_dp
    else if (tech_opt%accuracy_setting==1) then
      accuracy=1.0e-6_dp
    else if (tech_opt%accuracy_setting==0) then
      accuracy=1.0e-4_dp
    else if (tech_opt%accuracy_setting==-1) then
      accuracy=tech_opt%rk_accuracy_modes
    else
      print*, "MODECODE: accuracy_setting ==", tech_opt%accuracy_setting

      call raise%fatal_code(&
        "This accuracy_setting is not supported.",&
        __FILE__, __LINE__)

    end if

    hmin=1e-30_dp !minimum stepsize

    CALL odeint(y, x1, x2, accuracy, h1, hmin, derivs, qderivs, rkqs_c)
    nactual_mode = kount  ! update nactual after evolving the modes

    if(.not. ode_underflow) then
      powerspectrum_out = power_internal
      powerspectrum_out%bundle_exp_scalar=field_bundle%exp_scalar
    else
      powerspectrum_out%adiab= 0e0_dp
      powerspectrum_out%isocurv=0e0_dp
      powerspectrum_out%tensor=0e0_dp
      pk_bad = run_outcome%underflow
    endif


    contains

      ! A "packed" identity vector analog of identity matrix.
      subroutine make_identity(identityvector)

        real(dp), dimension(:), intent(out) :: identityvector
        integer :: i, j

        identityvector=0e0_dp

        forall (i=1:num_inflaton)&
          identityvector((i-1)*num_inflaton+i)=1e0_dp

      end subroutine make_identity

      subroutine set_background_and_mode_ic(y)

        ! Note that ah=log(aH) and overall scaled by sqrt(2k)

        complex(dp), dimension(:), intent(out) :: y
        real(dp), dimension(num_inflaton**2) :: identity

        ! Identity vector analog of identity matrix
        call make_identity(identity)

        ! Background - from previous evolution
        y(1:num_inflaton) = cmplx(p_ik,kind=dp)             !phi(x1)
        y(num_inflaton+1:2*num_inflaton) = cmplx(dp_ik,kind=dp) !Not in exact SR

        ! mode matrix - diagonalize, Bunch-Davies
        y(index_ptb_y:index_ptb_vel_y-1) = (1.e0_dp, 0)*identity  !cmplx(1/sqrt(2*k))
        y(index_ptb_vel_y:index_tensor_y-1) = cmplx(0., -k/exp(ah),kind=dp)*identity

        ! tensors
        y(index_tensor_y) = (1.e0_dp, 0) !cmplx(1/sqrt(2*k))
        y(index_tensor_y+1) = cmplx(0., -k/exp(ah),kind=dp)

        ! u_zeta
        y(index_uzeta_y) = (1.e0_dp, 0) !cmplx(1/sqrt(2*k))
        y(index_uzeta_y+1) = cmplx(0., -k/exp(ah),kind=dp)

      end subroutine set_background_and_mode_ic

      !Find the scale at which we can set the Bunch-Davies initial state
      !self-consistently
      !Note that there might be a correction due to massive modes m_heavy>H
      subroutine set_consistent_BD_scale(k_start)
        real(dp), intent(out) :: k_start

        real(dp) :: ah
        real(dp) :: horiz_fract
        integer :: ah_index

        real(dp) :: tol

        real(dp) :: alpha_ik, dalpha
        real(dp) :: h_ik, dh, a_ik
        real(dp) :: eps, V, dV(num_inflaton), d2V(num_inflaton,num_inflaton)
        real(dp), dimension(num_inflaton) :: p_ik, dp_ik, delphi
        real(dp) :: check1
        real(dp), dimension(num_inflaton,num_inflaton) :: check2, check3, check4

        logical :: bd_consistent

        bd_consistent = .false.

        horiz_fract=1e0_dp

        do while (.not. bd_consistent)

          horiz_fract = horiz_fract*2.0e0_dp

          ah=LOG(k/horiz_fract)
          ah_index= locate(log_aharr(1:nactual_bg), ah)

          if (ah_index==0) then
            !The background isn't able to set this IC
            !Set the start scale so far inside, that it fails
            !when it returns from this function.
            k_start = 1e20_dp
            return
          end if

          j=min(max(ah_index-(4-1)/2,1),nactual_bg+1-4)
          call array_polint(log_aharr(j:j+4), phiarr(:,j:j+4), ah, p_ik, delphi)
          call array_polint(log_aharr(j:j+4), dphiarr(:,j:j+4), ah, dp_ik, delphi)

          call polint(log_aharr(j:j+4), lna(j:j+4), ah,  alpha_ik, dalpha)
          call polint(log_aharr(j:j+4), hubarr(j:j+4), ah,  h_ik, dh)

          a_ik = exp(alpha_ik)*a_init
          eps = getEps(p_ik, dp_ik)

          V = pot(p_ik)
          dV = dVdphi(p_ik)
          d2V = d2Vdphi2(p_ik)

          !Check the corrections to the conformal-time mode equations.
          !When k is large, the IC should be a plane wave.
          !If the mass-matrix is diagonal and relevant, then could use Hankel
          !function solution
          check1 = abs((eps-2.0e0_dp)/(horiz_fract**2))
          check2 = abs( d2V/ (horiz_fract**2 * h_ik**2))

          forall (i=1:num_inflaton, j=1:num_inflaton)
            check3(i,j) = abs( (dp_ik(i)*dV(j) + &
              dp_ik(j)*dV(j))/(h_ik**2*horiz_fract**2))
            check4(i,j) = abs( (3.0e0_dp - eps)*dp_ik(i)*dp_ik(j)/horiz_fract**2)
          end forall

          tol = 1e-5_dp
          if   (check1 < tol  .and. &
            all(check2 < tol) .and. &
            all(check3 < tol) .and. &
            all(check4 < tol)) then

            bd_consistent = .true.

          end if

          k_start = horiz_fract

        end do

      end subroutine set_consistent_BD_scale


  end subroutine evolve

end module access_modpk
