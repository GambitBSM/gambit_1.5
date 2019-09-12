MODULE background_evolution
  !Module that evolves the background equation of motion after the initial
  !conditions and parameters are chosen.  Performs many checks to make sure that
  !the results make sense.
  USE camb_interface
  USE modpkparams
  USE modpk_odeint
  USE ode_path
  USE modpk_observables
  USE potential, ONLY : pot,getH, getdHdalpha, dVdphi, getEps, d2Vdphi2, &
    getH_with_t, stability_check_on_H, getEps_with_t
  USE modpk_utils, ONLY : bderivs, rkqs_r, use_t
  use modpk_numerics, only : locate, polint, array_polint
  use modpk_io, only : out_opt
  use modpk_errorhandling, only : raise, run_outcome

  IMPLICIT NONE
  PUBLIC :: backgrnd

CONTAINS

  SUBROUTINE backgrnd
    use modpk_icsampling, only : save_iso_N, N_iso_ref, phi_iso_N, &
      dphi_iso_N, ic_sampling, ic_flags

    INTEGER*4 :: i,j, rescl_count

    real(dp) :: phi_init_trial(size(phi_init))
    real(dp) :: alpha_e,dalpha,V_end,dv,dh
    real(dp) :: a_end_inst
    real(dp) :: V_i, ph, alpha_pivot, bb(size(phi_init))
    real(dp) :: Np_last

    CHARACTER(16) :: fmt = '(a25,es15.7)'

    !MULTIFIELD
    CHARACTER(19) :: array_fmt
    CHARACTER(len=5) :: ci

    real(dp) :: alpha_iso_N

    write(ci, '(I5)'), size(phi_init)
    ci = adjustl(ci)
    array_fmt = '(a25,'//trim(ci)//'es13.5)'
    !END MULTIFIELD

    !     Set initial conditions
    !
    !     y(1)=phi              dydx(1)=dphi/dx
    !     y(2)=dphi/dx          dydx(1)=d^2phi/dx^2

    if (.not. allocated(phi_init)) then
      call raise%fatal_code(&
        "Please initialize phi_init as an array!",&
        __FILE__, __LINE__)
    end if


    if (out_opt%modpkoutput .and. &
      .not. out_opt%output_reduced) then
      write(*, array_fmt) 'phi_init0 =', phi_init
      write(*, array_fmt) 'dphi_init0 =', dphi_init0
    end if

    !Initialize e-folds and background traj
    alpha_e=0e0_dp
    lna = 0e0_dp
    epsarr = 0e0_dp
    phiarr = 0e0_dp
    dphiarr = 0e0_dp

    if (size(phi_init) .eq. 1) then !! for single field

       phidot_sign = -(dVdphi(phi_init))/ABS(dVdphi(phi_init))
       !if < 0 field rolls from large to small
       !if > 0, the field rolls from small to large

       !Check whether phi_infl_end is in the correct direction from phi_init
       IF(.NOT.(slowroll_infl_end)) THEN
          IF (phidot_sign(1) .GT.0 .AND. phi_init(1) .gt. phi_infl_end(1)) THEN

            call raise%fatal_cosmo(&
             "Initial phi is smaller than final phi. &
             Please check your initial conditions.", &
             __FILE__, __LINE__)

          ENDIF
          IF (phidot_sign(1) .LT.0 .AND. phi_init(1) .lt. phi_infl_end(1)) THEN

            call raise%fatal_cosmo(&
             "Initial phi is larger than final phi. &
             Please check your initial conditions.", &
             __FILE__, __LINE__)

          ENDIF

       ENDIF


       !Try the background evolution for a user-defined phi_init and rescale if necessary
       rescl_count=0
       rescale_factor = 0.1e0_dp


       DO
          phi_init_trial=phi_init
          CALL trial_background(phi_init_trial, alpha_e, V_end)
          !! no rescaling implemented for multifield, so the code will exit here for multifield
          IF ((pk_bad/= run_outcome%success) &
            .OR. phi_init_trial(1).EQ.phi_init(1)) EXIT

          rescl_count=rescl_count+1
          IF (rescl_count .EQ. 50) THEN
             pk_bad= run_outcome%infl_didnt_start
             PRINT*,'MODECODE: phi_init rescaling did not work after 50 tries.'
             EXIT
          END IF
          if (out_opt%modpkoutput) then
             PRINT*,'MODECODE: phi_init was inconsistent. Rescaling:'
             WRITE(*,fmt) ' phi_init =', phi_init_trial
          end if
          phi_init = phi_init_trial
       END DO

    !MULTIFIELD
    else

      CALL trial_background(phi_init, alpha_e, V_end)

      if (pk_bad .ne. run_outcome%success) then
        !Override specific errors when doing IC sampling
        if (ic_sampling/=ic_flags%reg_samp .and. &
          pk_bad/=run_outcome%success) then
          return
        else
          call raise%fatal_cosmo(&
           "Bad set of parameters." , &
           __FILE__, __LINE__)
        end if
      end if
    end if
    !END MULTIFIELD

    IF(pk_bad==run_outcome%success) THEN

       V_i=pot(phi_init)

       !MULTIFIELD
       if ((alpha_e - lna(1)) .lt. N_pivot) then
          if (out_opt%modpkoutput) then
            write(*, *), "MODECODE: alpha_e - lna(1) =", alpha_e - lna(1),"< ",N_pivot
            call raise%warning(&
              'Not enough efolds obtained.')
          end if
          pk_bad = run_outcome%pivot_didnt_leaveH
          return

       else if ((alpha_e - lna(1)) .gt. Nefold_max) then

          if (out_opt%modpkoutput) then
            call raise%warning(&
              'Too many efolds obtained.',&
              __FILE__, __LINE__)
          end if
          pk_bad = run_outcome%cant_init_scalefact

       end if
       !MULTIFIELD

       IF (instreheat) THEN

          a_init = EXP(-71.1e0_dp - alpha_e + LOG(V_i/V_end)/4.e0_dp + LOG((M_Pl**4)/V_i)/4.e0_dp)
          Np_last = 0.5e0_dp*N_pivot
          do while (abs(N_pivot-Np_last)>0.01e0_dp)  ! determine N_pivot iteratively
             Np_last = N_pivot
             alpha_pivot = alpha_e-N_pivot

             i=locate(lna(1:nactual_bg),alpha_pivot)
             j=MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)
             CALL polint(lna(j:j+4), hubarr(j:j+4), alpha_pivot, H_pivot, dh)
             CALL array_polint(lna(j:j+4), phiarr(:, j:j+4), alpha_pivot, phi_pivot, bb)
             CALL array_polint(lna(j:j+4), dphiarr(:, j:j+4), alpha_pivot, dphi_pivot, bb)
             N_pivot = -71.1e0_dp-log(V_end/(M_Pl**4))/4.e0_dp-log(k_pivot*Mpc2Mpl/H_pivot)
          end do
          a_end = EXP(alpha_e)*a_init
       ELSE
          alpha_pivot = alpha_e-N_pivot

          if (save_iso_N) then
            alpha_iso_N = alpha_e - N_iso_ref
            if (alpha_iso_N<0e0_dp) then
              if (out_opt%modpkoutput) then
                print*, "MODECODE: alpha_iso_N=",alpha_iso_N,"<0"
              end if
              pk_bad = run_outcome%ref_efold_didnt_leaveH
              return
            end if
          end if


          i=locate(lna(1:nactual_bg),alpha_pivot)
          j=MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)
          CALL polint(lna(j:j+4), hubarr(j:j+4), alpha_pivot, H_pivot, dh)
          CALL array_polint(lna(j:j+4), phiarr(:,j:j+4), alpha_pivot, phi_pivot, bb)
          CALL array_polint(lna(j:j+4), dphiarr(:, j:j+4), alpha_pivot, dphi_pivot, bb)


          if (save_iso_N) then
            if (.not. allocated(phi_iso_N)) then
              call raise%fatal_code(&
                "The array phi_iso_N is not allocated.",&
                __FILE__, __LINE__)
            endif

            i=locate(lna(1:nactual_bg),alpha_iso_N)
            j=MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)
            CALL array_polint(lna(j:j+4), phiarr(:,j:j+4), &
              alpha_iso_N, phi_iso_N, bb)
            CALL array_polint(lna(j:j+4), dphiarr(:, j:j+4), &
              alpha_iso_N, dphi_iso_N, bb)
          end if

          a_init=k_pivot*Mpc2Mpl/H_pivot/EXP(alpha_pivot)


          a_end=EXP(alpha_e)*a_init

          a_end_inst=EXP(-71.1e0_dp+LOG(V_i/V_end)/4.e0_dp+LOG((M_Pl**4)/V_i)/4.e0_dp)


          if (isnan(a_end) .or. a_end<1.0e-100_dp) then
            print*, "MODECODE: a_end=", a_end
            print*, "MODECODE: N_efold=", alpha_pivot
            print*, "MODECODE: likely too many efolds before"
            print*, "MODECODE: pivot scale leaves horizon"
            if (ic_sampling/=ic_flags%reg_samp) then
              pk_bad = run_outcome%cant_init_scalefact
            else
              call raise%fatal_cosmo(&
                "See above.  &
                Try setting your IC closer to horizon crossing.",&
                __FILE__,__LINE__)
            end if

          end if

          IF (a_end .GT. a_end_inst) THEN
             PRINT*,'MODECODE: inflation ends too late with this N_pivot', N_pivot
             pk_bad = run_outcome%bad_reheat
             RETURN
          ENDIF
       END IF

       a_pivot = EXP(alpha_pivot)*a_init

       N_tot = alpha_e - lna(1)

       if (out_opt%modpkoutput) then
          WRITE(*, fmt) ' H_pivot =', H_pivot
          WRITE(*, fmt) ' a_end = ', a_end
          if (.not. out_opt%output_reduced) WRITE(*, array_fmt) ' phi_pivot =', phi_pivot
          WRITE(*, fmt) ' a_pivot =', a_pivot
          WRITE(*, fmt) ' N_pivot =', N_pivot
          WRITE(*, fmt) ' N_tot =', N_tot
       end if

       DO i=1,nactual_bg
          log_aharr(i)=LOG(a_init*EXP(lna(i))*hubarr(i))
       END DO
    ENDIF

    RETURN
  END SUBROUTINE backgrnd

  SUBROUTINE trial_background(phi_init_trial, alpha_e, V_end)
    use modpk_icsampling, only : ic_sampling, ic_flags

    INTEGER*4 :: i,j
    INTEGER*4, PARAMETER :: BNVAR=2

    !! MULTIFIELD
    real(dp), DIMENSION(:), INTENT(INOUT) :: phi_init_trial
    real(dp), DIMENSION(:) :: y(BNVAR*size(phi_init_trial))

    real(dp), DIMENSION(:) :: z_int_with_t(BNVAR*size(phi_init_trial)+1)

    logical :: H_stable, slowroll_init
    !! END MULTIFIELD

    real(dp) :: accuracy, h1, hmin, x1, x2
    real(dp) :: alpha_e, dalpha, V_end, dv, ep
    real(dp) :: ph, alpha_pivot, aa(size(phi_init_trial)), bb(size(phi_init_trial))
    real(dp) :: vv(nsteps) !!epsarr(nsteps),

    real(dp), DIMENSION(:) :: Vp(num_inflaton)
    real(dp) :: Vz, dotphi, thetaN, grad_V

    CHARACTER(16) :: fmt = '(a25,es15.7)'

    !MULTIFIELD
    CHARACTER(19) :: array_fmt
    CHARACTER(len=5) :: ci

    write(ci, '(I5)'), size(phi_init_trial)
    ci = adjustl(ci)
    array_fmt = '(a25,'//trim(ci)//'es13.5)'
    !END MULTIFIELD

    x1=0.0 !starting value
    x2=Nefold_max !ending value

    !MULTIFIELD
    !Set the ICS
    y(1 : size(y)/2) = phi_init_trial  !phi(x1)

    vv = 0e0_dp

    if (ic_sampling==ic_flags%slowroll_samp .or. ic_sampling==ic_flags%iso_N .or.&
      ic_sampling==ic_flags%reg_samp) then

      if( ic_sampling == ic_flags%reg_samp .and. &
        out_opt%modpkoutput) print*, "Setting velocity in slow-roll"

      !dphi/dalpha(x1) slowroll approx
      !MULTIFIELD

      h_init=pot(phi_init_trial)/(6.e0_dp*M_Pl**2) * &
           (1.e0_dp+SQRT(1.e0_dp+2.e0_dp/3.e0_dp* M_Pl**2 *&
           dot_product(dVdphi(phi_init_trial), dVdphi(phi_init_trial)) &
           / pot(phi_init_trial)**2.))

      if (h_init<0.0e0_dp) then
        print*, "MODECODE: h_init=", h_init
        call raise%fatal_code(&
          "The initial value of H is imaginary.  This is &
          likely related to a problem in the definition of V &
          or its derivatives.",&
          __FILE__, __LINE__)
      else
        h_init = sqrt(h_init)
      end if

      !END MULTIFIELD

      y(size(y)/2+1 : (size(y))) = &
        -dVdphi(phi_init_trial)/3.e0_dp/h_init/h_init

      dphi_init0 = y(size(y)/2+1 : (size(y)))

    else

      h_init = getH(phi_init_trial,dphi_init0)

      !Set not necess close to SR
      y(size(y)/2+1 : (size(y))) = dphi_init0

    end if
    !END MULTIFIELD


    !Call the integrator
    ode_underflow = .FALSE.
    ode_ps_output = .FALSE.
    ode_infl_end = .TRUE.
    save_steps = .true.

    pk_bad = run_outcome%success

    !MULTIFIELD
    IF(getEps(y(1:size(y)/2),y(size(y)/2+1:size(y))) .GT. 0.2) THEN
       slowroll_start=.FALSE.
    ELSE
       slowroll_start=.TRUE.
    ENDIF

    !Check to see if starts outside slowroll, then slow-rolls later
    slowroll_init = slowroll_start

    if (out_opt%modpkoutput) write(*,'(a25, L2)') 'slowroll start =', slowroll_start
    !END MULTIFIELD

    !guessed start stepsize
    if (potential_choice.eq.6) then
       h1 = 0.001e0_dp
    end if

    dxsav=1.e-7_dp

    !Check if ICs give instability in H^2=V/(3-eps)
    !If unstable, then integrate in cosmic time t until reach stable region
    !for e-fold integrator
    H_stable = .false.
    call stability_check_on_H(H_stable,y(1:num_inflaton),&
      y(num_inflaton+1:2*num_inflaton),&
      using_t=.false.)

    if (.not. H_stable) then
      !Decrease initial stepsize guess in case near a point where H is approx unstable.
      !print*, "H is UNSTABLE"
      h1 = 1.0e-12_dp
      accuracy = 1.0e-15_dp
      hmin = 0.0e0_dp
    else
      !print*, "H is STABLE"
      h1 = 1.0e-7_dp
      if (tech_opt%accuracy_setting>0) then
        accuracy=1.0e-8
      else
        accuracy=tech_opt%rk_accuracy_back
      end if

      hmin=1.0e-20_dp
    end if

    CALL odeint(y,x1,x2,accuracy,h1,hmin,bderivs,rkqs_r)

    if (ic_sampling/=ic_flags%reg_samp .and. &
      pk_bad/=run_outcome%success) return

    IF(.NOT. ode_underflow) THEN
      if (size(lna) < kount .or. size(xp) < kount) then
        call raise%fatal_code(&
         "kount is too big.  &
         Giving this error instead of &
         a segmentation fault.", &
        __FILE__, __LINE__)
      end if

       lna(1:kount)=xp(1:kount)
       phiarr(:,1:kount)=yp(1:size(y)/2, 1:kount)
       dphiarr(:,1:kount)=yp(size(y)/2+1:size(y),1:kount)

       !MULTIFIELD

       DO i=1,kount

          vv(i) = pot(phiarr(:,i))
          hubarr(i) = getH(phiarr(:, i),dphiarr(:,i))
          epsarr(i) = getEps(phiarr(:,i),dphiarr(:,i))
          sigma_arr(i) = sqrt(dot_product((phiarr(:,i)-phi_init),(phiarr(:,i)-phi_init)))

       END DO
       !END MULTIFIELD

       !
       !     Determine the parameters needed for converting k(Mpc^-1) to K
       !
       nactual_bg=kount
       IF(slowroll_infl_end) THEN
          ep = 1.0e0_dp

          i=locate(epsarr(1:kount),ep)

          !If didn't start in SR, but SR commenced and ended, then there are two
          !points where epsilon=1.0e0_dp.  Need to find the last one
          if (.not. slowroll_init) i=i+1 +locate(epsarr(i+2:kount),ep)

          j=MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)
          CALL polint(epsarr(j:j+4), lna(j:j+4), ep, alpha_e, dalpha)
          CALL polint(epsarr(j:j+4), vv(j:j+4), ep, V_end, dv)
          !MULTIFIELD
          CALL array_polint(epsarr(j:j+4), phiarr(:, j:j+4), ep, phi_infl_end, bb)

          !Check for interpolation errors
          if(dalpha > 0.1 .or. dv > 0.1 .or. bb(1) > 0.1) THEN

             !Check if didn't get enough e-folds
             if (lna(kount) < N_pivot) then
               if (out_opt%modpkoutput) then
                 write(*, *), "MODECODE: lna(kount) - lna(1) =", lna(kount) - lna(1),"< ",N_pivot
                 call raise%warning(&
                    'Not enough efolds obtained.',&
                    __FILE__, __LINE__)
               end if
               pk_bad = run_outcome%pivot_didnt_leaveH
               return
             else

               print*,"MODECODE: dalpha", dalpha
               print*,"MODECODE: dv", dv
               print*,"MODECODE: bb", bb
               print*,"MODECODE: lna", lna(kount-5:kount), "alpha_e", alpha_e
               print*,"MODECODE: epsarr", epsarr(kount-5:kount), "ep", ep

               call raise%fatal_code(&
               'The interpolation in SUBROUTINE trial_background &
               has suspiciously large errors.  &
               Your model smells fishy.' , &
               __FILE__, __LINE__)

             end if
           endif

       ELSE
          if (size(phi_init) .eq. 1) then
             ep = abs(phi_init(1) - phi_infl_end(1))
             i = locate(sigma_arr(1:kount), ep)
             j = MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)
             CALL polint(sigma_arr(j:j+4), lna(j:j+4), ep, alpha_e, dalpha)
             V_end = pot(phi_infl_end)
          else

            !If we have stopped inflation "by hand", then we can auto-set these
            select case(potential_choice)
            case(13)
              dalpha = 0e0_dp
              dv = 0e0_dp
              bb = 0e0_dp

              alpha_e = lna(nactual_bg)
              V_end = vv(nactual_bg)
              phi_infl_end = phiarr(:,nactual_bg)

            case default
              ep = delsigma
              i = locate(sigma_arr(1:kount), ep)
              j = MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)
              CALL polint(sigma_arr(j:j+4), lna(j:j+4), ep, alpha_e, dalpha)
              CALL polint(sigma_arr(j:j+4), vv(j:j+4), ep, V_end, dv)
              CALL array_polint(sigma_arr(j:j+4), phiarr(:, j:j+4), ep, phi_infl_end, bb)
            end select
          end if
          !END MULTIFIELD
       ENDIF

       if (out_opt%modpkoutput .and. .not. out_opt%output_reduced) &
         WRITE(*,array_fmt) ' phi_end =', phi_infl_end

       IF (instreheat) THEN
          !Set a plausible pivot (ignoring logarithmic V_k and V_end terms) for
          !the smallest k likely to be requested by CAMB.
          N_pivot=-LOG10(1.d-6)+60.
       ENDIF

       !MULTIFIELD
       !
       !Rescaling the initial value only implemented for single field cases.
       !For multi-field, one may need to specify whether to rescale a sigle compoent field,
       !or to homogeneously rescale all the fields. It is better to work on a case by case basis.
       !

       IF (size(phi_init) .eq. 1) THEN
          IF(alpha_e.LT.(N_pivot+20.e0_dp)) THEN
             IF ((potential_choice.eq.6).and.(vparams(1,1)<-2.e0_dp)) THEN
                phi_init_trial = phi_init*0.9e0_dp
             ELSE
                phi_init_trial = phi_init+(phi_init-phi_infl_end)*rescale_factor
             ENDIF
             RETURN
          END IF

          IF (alpha_e.GT.N_pivot+55e0_dp) THEN
             ep = alpha_e-(N_pivot+50e0_dp)
             i = locate(lna(1:kount),ep)
             j = MIN(MAX(i-(4-1)/2,1),nactual_bg+1-4)
             CALL array_polint(lna(j:j+4), phiarr(:,j:j+4), ep, aa, bb)
             phi_init_trial = aa
             RETURN
          END IF
       ENDIF
    ELSE
       pk_bad = run_outcome%underflow
    ENDIF

    !END MULTIFIELD

  END SUBROUTINE trial_background

  SUBROUTINE backgrnd_efold

    INTEGER*4, PARAMETER :: BNVAR=2
    real(dp), DIMENSION(:) :: y(BNVAR*num_inflaton)
    real(dp) :: accuracy, h1, hmin, x1, x2


    h_init=SQRT(pot(phi_init)/(6.0e0_dp*M_Pl**2) * &
         (1.0e0_dp+SQRT(1.0e0_dp+2.0e0_dp/3.0e0_dp* M_Pl**2 &
         * dot_product(dVdphi(phi_init), dVdphi(phi_init)) / pot(phi_init)**2.)))

    x1=0.0e0_dp !starting value
    x2=Nefold_max !ending value

    y(1 : size(y)/2) = phi_init  !phi(x1)
    y(size(y)/2+1 : (size(y))) = -dVdphi(phi_init)/3.e0_dp/h_init/h_init !dphi/dalpha(x1) slowroll approx

    !Call the integrator
    ode_underflow = .FALSE.
    ode_ps_output = .FALSE.
    ode_infl_end = .TRUE.
    save_steps = .TRUE.
    pk_bad = run_outcome%success

    IF(getEps(y(1:size(y)/2),y(size(y)/2+1:size(y))) .GT. 1.) THEN
       slowroll_start=.FALSE.
    ELSE
       slowroll_start=.TRUE.
    ENDIF
    if (out_opt%modpkoutput) write(*,'(a25, L2)') 'slowroll start =', slowroll_start

    !guessed start stepsize
    if (potential_choice.eq.6) then
       h1 = 0.001e0_dp
    else
       h1 = 0.1e0_dp
    end if

    dxsav=1.e-7_dp
    accuracy=1.0e-6_dp
    hmin=0.0e0_dp !minimum stepsize
    CALL odeint(y,x1,x2,accuracy,h1,hmin,bderivs,rkqs_r)

    IF(.NOT. ode_underflow) THEN
       lna(1:kount)=xp(1:kount)
       nactual_bg=kount
    ELSE
       pk_bad = run_outcome%underflow
    END IF

    N_tot = lna(nactual_bg) - lna(1)

  END SUBROUTINE backgrnd_efold


end module background_evolution
