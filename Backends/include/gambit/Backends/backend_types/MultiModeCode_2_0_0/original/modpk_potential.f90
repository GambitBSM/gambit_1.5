module potential
  !Module that defines the inflationary potential and its derivatives.
  !Implement your potential by adding a new case here.  Also contains routines
  !for calculating cosmological parameters and power spectra.
  use modpkparams
  use internals, only : pi
  use modpk_numerics
  use modpk_errorhandling, only : raise, assert
  implicit none
  private
  public :: pot, getH, getdHdalpha, getEps, dVdphi, d2Vdphi2, getdepsdalpha, powerspectrum, &
       tensorpower, initialphi, geteta, zpower, getH_with_t, stability_check_on_H, &
       getEps_with_t, d3Vdphi3

  public :: norm
  public :: bundle, field_bundle

  !Bundle width isocurvature measure
  type :: bundle
    real(dp) :: N=0e0_dp
    real(dp) :: exp_scalar=0e0_dp
    real(dp) :: dlogThetadN=0e0_dp
    contains
      procedure, public :: calc_exp_scalar =>bundle_exp_scalar
  end type bundle

  type(bundle) :: field_bundle

contains


  recursive function pot(phi) result(V_potential)
    !
    !     Returns V(phi) given phi, phi can be an array for multifield,
    !     The code implement multifield potential in the form of V = \sum V(phi_i),
    !     More complicated form of potentials can be customized
    !
    real(dp) :: V_potential
    real(dp), intent(in) :: phi(:)

    real(dp) :: m2_V(size(phi)) !! m2_V is the diagonal mass square matrix
    real(dp) :: c1_V(size(phi)) !! some constants for a more flexible potential construction in case 9
    real(dp) :: lambda(size(phi)), finv(size(phi)), mu(size(phi))
    real(dp) :: M2, theta2, c2, s2, mphi1, potsmall, potlarge, phi1shift ! messy parameters for case 10
    real(dp) :: lambda_hybrid, mu_hybrid, nu_hybrid, &
      mass_hybrid

    integer :: phi_light_index
    real(dp) :: lambda4(size(phi)), alpha2(num_inflaton)

    real(dp) :: lambda2
    real(dp), dimension(size(phi),size(phi)) :: m2_matrix
    integer :: i,j, temp_choice

    real(dp), dimension(size(phi)) :: location_phi, step_size, step_slope

    real(dp) :: m_light2, M_heavy2, param0, param_closest, dist,&
      phi_light, phi_light0

    real(dp) :: p_exp

    real(dp) :: V0
    real(dp), dimension(size(phi)) :: A_i
    real(dp), dimension(size(phi),size(phi)) :: B_ij
    integer :: ii, jj

    !Make sure that vparams is set up properly
    if (.not. allocated(vparams)) then
      call raise%fatal_code(&
        "The array vparams that contains the potential parameters is &
        not allocated.",&
        __FILE__,__LINE__)
    end if

    select case(potential_choice)

    case(1)
      ! m_i^2 phi_i^2 --- N-quadratic

      call assert%check(size(vparams,1)>=1,__FILE__,__LINE__)

      m2_V = 10.e0_dp**(vparams(1,:))
      V_potential = 0.5e0_dp*sum(m2_V*phi*phi)

    case(2)
      ! N-flation (axions)

      call assert%check(size(vparams,1)>=2,__FILE__,__LINE__)

      lambda = 10.e0_dp**vparams(1,:)
      finv = 1.e0_dp/(10.e0_dp**vparams(2,:))
      V_potential = sum(lambda**4*(1.e0_dp+cos(finv*phi)))

    case(3)
      ! lambda_i phi_i^4 --- N-quartic

      call assert%check(size(vparams,1)>=1,__FILE__,__LINE__)

      lambda = 10.e0_dp**vparams(1,:)
      V_potential = sum(0.25e0_dp*lambda*phi**4)

    case(4)
      ! lambda_i phi_i --- N-linear

      call assert%check(size(vparams,1)>=1,__FILE__,__LINE__)

      lambda = 10.e0_dp**vparams(1,:)
      V_potential = sum(lambda*phi)

    case(5)
      ! lambda_i phi_i^(2/3)

      call assert%check(size(vparams,1)>=1,__FILE__,__LINE__)

      lambda = 10.e0_dp**vparams(1,:)
      V_potential = sum(lambda*1.5e0_dp*phi**(2.e0_dp/3.e0_dp))

    case(6)
      ! Lambda^4 - mu*phi^4

      call assert%check(size(vparams,1)>=2,__FILE__,__LINE__)

      lambda = 10.e0_dp**vparams(1,:)
      mu = 10.e0_dp**vparams(2,:)
      V_potential = sum(lambda**4 - mu*phi**4/4.e0_dp)

    case(7)
      ! Product of exponentials

      call assert%check(size(vparams,1)>=2,__FILE__,__LINE__)

      lambda = vparams(2, :)
      V_potential = vparams(1,1)*M_Pl**4 * exp(dot_product(lambda, phi/M_Pl))

    case(8)
      !Canonical two-field hybrid
      if (size(phi) /= 2) then
        print*, "MODECODE: Potential_choice =", Potential_choice
        print*, "MODECODE: Number of fields =", size(phi)
        call raise%fatal_cosmo(&
            "This potential requires two fields.  &
            Set num_inflaton=2.", &
            __FILE__, __LINE__)
      end if


      call assert%check(size(vparams,1)>=1,__FILE__,__LINE__)
      call assert%check(size(vparams,2)>=4,__FILE__,__LINE__)

      lambda_hybrid = vparams(1,1)
      mass_hybrid = vparams(1,2)
      mu_hybrid =vparams(1,3)
      nu_hybrid =vparams(1,4)
      V_potential = (lambda_hybrid**4)*((1.0_dp - phi(1)**2/mass_hybrid**2)**2 +&
        phi(2)**2/mu_hybrid**2 +&
        phi(1)**2*phi(2)**2/nu_hybrid**4)

    case(9)
      ! V0 + m_i^2 phi_i^2

      call assert%check(size(vparams,1)>=2,__FILE__,__LINE__)

      m2_V = vparams(1,:)
      c1_V = vparams(2,:)
      V_potential = 0.5e0_dp*sum(m2_V*phi*phi) + sum(c1_V)

    case(10)
      print*, "MODECODE: Potential_choice=", Potential_choice
      call raise%fatal_code(&
        "This potential choice is broken for some reason.",&
        __FILE__, __LINE__)

    case(11)
      !N-quadratic w/one quartic intxn term
      !phi_i^2 + phi_{lightest}^2*phi_i^2

      call assert%check(size(vparams,1)>=2,__FILE__,__LINE__)

      m2_V = 10.0e0_dp**(vparams(1,:))
      lambda4 = 10.0e0_dp**vparams(2,:)
      phi_light_index = minloc(m2_V,1)

      V_potential = 0.5e0_dp*sum(m2_V*phi*phi)+ &
        (1.0e0_dp/24.0e0_dp)*(phi(phi_light_index)**2)*sum(lambda4*phi*phi)

    case(12)
      !Mass matrix with diagonal terms = m_i^2
      !Off-diagonal terms = \eps

      call assert%check(size(vparams,1)>=2,__FILE__,__LINE__)

      m2_V = 10.e0_dp**(vparams(1,:))
      lambda2 = 10.e0_dp**(vparams(2,1))
      V_potential = 0e0_dp

      do i=1, num_inflaton
        do j=1, num_inflaton

          if (i==j) then
            m2_matrix(i,j) = m2_V(i)
            V_potential = V_potential +0.5e0_dp* m2_matrix(i,j)*phi(i)*phi(j)
          else
            m2_matrix(i,j) = lambda2
            V_potential = V_potential + m2_matrix(i,j)*phi(i)*phi(j)
          end if

        end do
      end do

    case(13)

      print*, "MODECODE: Potential_choice=", Potential_choice
      call raise%fatal_code(&
        "This potential choice is broken for some reason.",&
        __FILE__, __LINE__)

    case(14)

      call assert%check(size(vparams,1)>=4,__FILE__,__LINE__)


      !Multifield step potential
      m2_V = 10.e0_dp**(vparams(1,:))
      location_phi = vparams(2,:)
      step_size = vparams(3,:)
      step_slope = vparams(4,:)


      !Check for /0 error
      if (any(abs(step_slope) < 1e-15)) then
        print*, "MODECODE: step_slope=", step_slope
        call raise%fatal_code(&
          "This is a division by zero error.  &
          Set the slope in tanh(phi/step_slope) greater than about 1e-15.",&
          __FILE__, __LINE__)
      end if

      V_potential=0e0_dp

#define ARG ((phi(i)-location_phi(i))/step_slope(i))
      do i=1,size(phi)
        V_potential = V_potential + &
            0.5*m2_V(i)*(phi(i)**2) * &
            (1.0e0_dp + step_size(i)* &
            tanh(ARG ))
      end do

#undef ARG


    case(15)

      print*, "MODECODE: Potential_choice=", Potential_choice
      call raise%fatal_code(&
        "This potential choice is broken for some reason.",&
        __FILE__, __LINE__)

    case(16)
      ! (1/p) lambda_i |phi_i|^p --- N-monomial

      call assert%check(size(vparams,1)>=2,__FILE__,__LINE__)

      p_exp = vparams(2,1)
      m2_V = vparams(1,:)
      V_potential = (1.0e0_dp/p_exp)*sum(m2_V*abs(phi)**p_exp)

    case(17)
      ! Generalized axions

      call assert%check(size(vparams,1)>=1+num_inflaton,__FILE__,__LINE__)

      A_i = vparams(1,:)
      B_ij = vparams(2:num_inflaton+1,:)


      V_potential = 0e0_dp
      do ii=1,size(A_i)
        V_potential = V_potential + A_i(ii)*(1.0e0_dp-cos(2.0e0_dp*pi*phi(ii)))
      end do


      do ii=1,size(A_i); do jj=1, size(A_i)
        V_potential = V_potential + B_ij(ii,jj)*&
          cos(2.0e0_dp*pi*phi(ii) - 2.0e0_dp*pi*phi(jj))
      end do; end do



    case default
      print*, "MODECODE: potential_choice =", potential_choice
      call raise%fatal_code(&
          "Need to set the potential V(phi) for this potential choice.",&
          __FILE__, __LINE__)

    end select

  END FUNCTION pot


  !Needs "result" because array-valued and recursive.
  recursive function dVdphi(phi) result(first_deriv)
    !
    !     Returns dV/dPhi given phi
    !

    real(dp), INTENT(IN) :: phi(:)
    real(dp) :: first_deriv(size(phi))

    real(dp) :: dphi(size(phi)), phiplus(size(phi))
    real(dp) :: m2_V(size(phi)),  c1_V(size(phi)), lambda(size(phi)), finv(size(phi)), mu(size(phi))
    real(dp) :: M2, theta2, c2, s2, mphi1, potsmall, potlarge, phi1shift ! messy parameters for case 10
    integer :: i,j, temp_choice

    real(dp) :: lambda_hybrid, mu_hybrid, nu_hybrid, &
      mass_hybrid

    integer :: phi_light_index
    real(dp) :: lambda4(size(phi)), alpha2(size(phi))

    real(dp) :: lambda2
    real(dp), dimension(size(phi),size(phi)) :: m2_matrix

    real(dp), dimension(size(phi)) :: location_phi, step_size, step_slope

    real(dp), dimension(size(phi)) :: stepsize
    real(dp), dimension(:), allocatable :: numderiv

    real(dp) :: param_closest
    real(dp) :: m_light2, M_heavy2
    real(dp) :: phi_light, dphi_light
    real(dp), dimension(num_inflaton) :: turnfunct, dturnfunct, &
      dparam_dphi

    real(dp) :: p_exp
    integer :: ii

    real(dp) :: V0, twopi
    real(dp), dimension(size(phi)) :: A_i
    real(dp), dimension(size(phi),size(phi)) :: B_ij
    integer :: alpha, beta

    if (vnderivs) then
       ! MULTIFIELD
       do i=1, num_inflaton
          phiplus = phi
          phiplus(i) = phi(i) + 0.5e0_dp*phi(i)*findiffdphi**(1.0e0_dp/3.0e0_dp)
          dphi = phiplus - phi
          first_deriv(i) = (pot(phi+dphi)-pot(phi-dphi))/(2.0e0_dp*dphi(i))
          if (first_deriv(i).eq.0.e0_dp) then

            call raise%fatal_code(&
             'first_deriv(i)=0, possibly signaling a problem with &
             accuracy of numerical derivatives. &
             Try using vnderivs=F if possible.',&
             __FILE__, __LINE__)
          end if
       end do
       ! END MULTIFIELD
    else
       select case(potential_choice)
       !MULTIFIELD
       case(1)
          m2_V = 10.e0_dp**(vparams(1,:))
          first_deriv = m2_V*phi
       case(2)
          lambda = 10.e0_dp**vparams(1,:)
          finv = 1.e0_dp/(10.e0_dp**vparams(2,:))
          first_deriv = -lambda**4*finv*sin(finv*phi)
       case(3)
          lambda = 10.e0_dp**vparams(1,:)
          first_deriv = lambda*phi**3
       case(4)
          lambda = 10.e0_dp**vparams(1,:)
          first_deriv = lambda
       case(5)
          lambda = 10.e0_dp**vparams(1,:)
          first_deriv = lambda*phi**(-1./3.)
       case(6)
          mu = 10.e0_dp**vparams(2,:)
          first_deriv = -mu*phi**3
       case(7)
          lambda = vparams(2,:)
          first_deriv = lambda*vparams(1,1)*M_Pl**3 * exp(dot_product(lambda, phi/M_Pl))
       case(8)
          !Canonical two-field hybrid
          lambda_hybrid = vparams(1,1)
          mass_hybrid = vparams(1,2)
          mu_hybrid =vparams(1,3)
          nu_hybrid =vparams(1,4)

          first_deriv(1) = (lambda_hybrid**4)*(4.0_dp*(phi(1)**2/mass_hybrid**2 - 1.0_dp)&
            *phi(1)/mass_hybrid**2 +&
            2.0_dp*phi(1)*phi(2)**2/nu_hybrid**4)
          first_deriv(2) = (lambda_hybrid**4)*(2.0_dp*phi(2)/mu_hybrid**2 +&
            2.0_dp*phi(1)**2*phi(2)/nu_hybrid**4)
       case(9)!Saddle type things
          m2_V = (vparams(1,:))
          first_deriv = m2_V*phi
       case(10)
         print*, "MODECODE: Potential_choice=", Potential_choice
         call raise%fatal_code(&
           "This potential choice is broken for some reason.",&
           __FILE__, __LINE__)

       case(11)
         m2_V = 10.e0_dp**(vparams(1,:))
         lambda4 = 10.e0_dp**vparams(2,:)
         phi_light_index = minloc(m2_V,1)

         !for i /= lightest
         first_deriv = m2_V*phi + (1.0e0_dp/12.0e0_dp)*(phi(phi_light_index)**2)*lambda4*phi

         !i=lightest
         first_deriv(phi_light_index) = m2_V(phi_light_index)*phi(phi_light_index) + &
           (1.0e0_dp/6.0e0_dp)*(phi(phi_light_index)**3)*lambda4(phi_light_index)

       case(12)

         !Mass matrix with diagonal terms = m_i^2
         !Off-diagonal terms = \eps
         m2_V = 10.e0_dp**(vparams(1,:))
         lambda2 = 10.e0_dp**(vparams(2,1))
         first_deriv = 0e0_dp

         do i=1, num_inflaton
           do j=1, num_inflaton
             if (i==j) then
               m2_matrix(i,j) = m2_V(i)
             else
               m2_matrix(i,j) = lambda2
             end if

             first_deriv(i) = first_deriv(i) + m2_matrix(i,j)*phi(j)
           end do
         end do

       case(13)

         print*, "MODECODE: Potential_choice=", Potential_choice
         call raise%fatal_code(&
           "This potential choice is broken for some reason.",&
           __FILE__, __LINE__)

       case(14)

         !Multifield step potential
         m2_V = 10.e0_dp**(vparams(1,:))
         location_phi = vparams(2,:)
         step_size = vparams(3,:)
         step_slope = vparams(4,:)


#define ARG ((phi(i)-location_phi(i))/step_slope(i))

         do i=1,size(phi)
           first_deriv(i) = m2_V(i)*phi(i)* &
             (1.0e0_dp + step_size(i)*tanh(ARG )) + &
             0.5e0_dp*m2_V(i)*phi(i)**2*step_size(i)* &
             (MySech(ARG))**(2)/&
             step_slope(i)
         end do

#undef ARG

       case(15)

         print*, "MODECODE: Potential_choice=", Potential_choice
         call raise%fatal_code(&
           "This potential choice is broken for some reason.",&
           __FILE__, __LINE__)

       case(16)
         ! (1/p) lambda_i |phi_i|^p --- N-monomial

         p_exp = vparams(2,1)
         m2_V = vparams(1,:)
         first_deriv = m2_V*abs(phi)**(p_exp-1.0e0_dp)*sign(1.0e0_dp,phi)

         !Regularize around phi=0
         do ii=1,size(phi)
           if (abs(phi(ii))<1e-5_dp) &
             first_deriv(ii)=0e0_dp
         end do

       case(17)
         ! Generalized axions

         A_i = vparams(1,:)
         B_ij = vparams(2:num_inflaton+1,:)

         twopi = 2.0e0_dp*pi

         first_deriv = -twopi*A_i*sin(twopi*phi)

         do alpha=1, size(first_deriv)
           do ii=1, size(first_deriv)
             if (alpha==ii) cycle
             first_deriv(alpha) = first_deriv(alpha) - &
               twopi*(B_ij(alpha,ii)+B_ij(ii,alpha))* &
               sin(twopi*phi(alpha) - twopi*phi(ii))
           end do
         end do


       !END MULTIFIELD
       case default

         print*, "MODECODE: potential_choice =", potential_choice
         call raise%fatal_code(&
           "Need to set first derivative for this potential choice &
           or use numerical derivatives (vnderivs=.true.)",&
           __FILE__, __LINE__)

       end select

    end if

  END FUNCTION dVdphi


  !Needs "result" because array-valued and recursive.
  recursive function d2Vdphi2(phi) result(second_deriv)
    !
    !     Returns d^2V/dPhi^2 given phi
    !

    real(dp), INTENT(IN) :: phi(:)
    real(dp) :: second_deriv(size(phi),size(phi))

    real(dp) :: m2_V(size(phi)), c1_V(size(phi)), lambda(size(phi)), finv(size(phi)), mu(size(phi))
    real(dp) :: M2, theta2, c2, s2, mphi1, potsmall, potlarge, phi1shift ! messy parameters for case 10
    integer :: i,j, temp_choice

    real(dp) :: dphi,phiplus
    real(dp) :: lambda_hybrid, mu_hybrid, nu_hybrid, &
      mass_hybrid

    integer :: phi_light_index
    real(dp) :: lambda4(size(phi)), alpha2(size(phi))

    real(dp) :: lambda2
    real(dp), dimension(size(phi),size(phi)) :: m2_matrix

    real(dp), dimension(size(phi)) :: location_phi, step_size, step_slope
    real(dp), dimension(size(phi)) :: stepsize
    real(dp), dimension(:,:), allocatable :: numderiv
    real(dp), dimension(:), allocatable :: numderiv1

    real(dp) :: param_closest
    real(dp) :: m_light2, M_heavy2
    real(dp) :: phi_light, dphi_light, d2phi_light
    real(dp), dimension(num_inflaton) :: turnfunct, dturnfunct, d2turnfunct, &
      dparam_dphi
    real(dp), dimension(num_inflaton, num_inflaton) :: d2param_dphi2, delta
    integer :: ii, jj, ll

    real(dp) :: p_exp

    real(dp) :: V0, twopi, fourpi_sq
    real(dp), dimension(size(phi)) :: A_i
    real(dp), dimension(size(phi),size(phi)) :: B_ij
    integer :: alpha, beta

    if (vnderivs) then
       !MULTIFIELD
       if (size(phi) .ne. 1) then
         call raise%fatal_code(&
          'The 2nd order numerical derivative has not &
          been implemented for more than one field.',&
          __FILE__, __LINE__)
       end if
       phiplus = phi(1) + 0.2e0_dp*phi(1)*findiffdphi**(1.e0_dp/4.e0_dp)
       dphi = phiplus - phi(1)
       second_deriv = (pot(phi+2.e0_dp*dphi)+pot(phi-2.e0_dp*dphi)- &
         2.e0_dp*pot(phi))/(4.e0_dp*dphi*dphi)
    else
       second_deriv(:,:) = 0e0_dp

       select case(potential_choice)
       case(1)
          m2_V = 10.e0_dp**(vparams(1,:))
          forall (i = 1:size(phi)) second_deriv(i,i) = m2_V(i)
       case(2)
          lambda = 10.e0_dp**vparams(1,:)
          finv = 1.e0_dp/(10.e0_dp**vparams(2,:))
          forall (i = 1:size(phi)) second_deriv(i,i) = -lambda(i)**4*finv(i)*finv(i)*cos(finv(i)*phi(i))
       case(3)
          lambda = 10.e0_dp**vparams(1,:)
          forall (i = 1:size(phi)) second_deriv(i,i) = 3.e0_dp*lambda(i)*phi(i)**2
       case(4)
          forall (i = 1:size(phi)) second_deriv(i,i) = 0.e0_dp
       case(5)
          lambda = 10.e0_dp**vparams(1,:)
          forall (i = 1:size(phi)) second_deriv(i,i) = -lambda(i)/3.e0_dp*phi(i)**(-4./3.)
       case(6)
          mu = 10.e0_dp**vparams(2,:)
          forall (i = 1:size(phi)) second_deriv(i,i) = -3.e0_dp*mu(i)*phi(i)**2
       case(7)
          lambda = vparams(2, :)
          forall (i=1:size(phi), j=1:size(phi)) second_deriv(i,j) = &
               lambda(i)*lambda(j)*vparams(1,1)*M_Pl**2 *exp(dot_product(lambda, phi/M_Pl))
       case(8)
          !Canonical two-field hybrid
          lambda_hybrid = vparams(1,1)
          mass_hybrid = vparams(1,2)
          mu_hybrid =vparams(1,3)
          nu_hybrid =vparams(1,4)

          second_deriv(1,1) = (lambda_hybrid**4)*(4.0_dp*(2.0_dp*phi(1)/mass_hybrid**2)*phi(1)/mass_hybrid**2 +&
          4.0_dp*(phi(1)**2/mass_hybrid**2 - 1.0_dp)/mass_hybrid**2 +&
            2.0_dp*phi(1)*phi(2)**2/nu_hybrid**4)

          second_deriv(2,1) = (lambda_hybrid**4)*(4.0_dp*phi(1)*phi(2)/nu_hybrid**4)

          second_deriv(1,2) = second_deriv(2,1)

          second_deriv(2,2) = (lambda_hybrid**4)*(2.0_dp/mu_hybrid**2 +&
            2.0_dp*phi(1)**2/nu_hybrid**4)
       case(9)
          m2_V = vparams(1,:)
          forall (i = 1:size(phi)) second_deriv(i,i) = m2_V(i)
       case(10)
         print*, "MODECODE: Potential_choice=", Potential_choice
         call raise%fatal_code(&
           "This potential choice is broken for some reason.",&
           __FILE__, __LINE__)

       case(11)
         M2_V = 10.e0_dp**(vparams(1,:))
         Lambda4 = 10.e0_dp**vparams(2,:)
         Phi_light_index = minloc(m2_V,1)

         do i=1,num_inflaton
           do j=1, num_inflaton

             if (i/=Phi_light_index .and. j/= Phi_light_index) then
               if (i/=j) then
                 second_deriv(i,j) = 0e0_dp
               else
                 second_deriv(i,j) = m2_V(i) + &
                   (1.0e0_dp/12.0e0_dp)*Lambda4(i)*phi(Phi_light_index)**2
               end if

             else if (i==Phi_light_index .and. j/=Phi_light_index) then

               second_deriv(i,j) = 0e0_dp

             else if (i/=Phi_light_index .and. j==Phi_light_index) then

               second_deriv(i,j) = (1.0e0_dp/6.0e0_dp)*Lambda4(i)* &
                 phi(Phi_light_index)*phi(i)

             else if (i==Phi_light_index .and. j==Phi_light_index) then

               second_deriv(i,j) = m2_V(Phi_light_index) + &
                 0.5e0_dp*Lambda4(Phi_light_index)*phi(Phi_light_index)**2

             end if

           end do

         end do

       case(12)

         !Mass matrix with diagonal terms = m_i^2
         !Off-diagonal terms = \eps
         m2_V = 10.e0_dp**(vparams(1,:))
         lambda2 = 10.e0_dp**(vparams(2,1))
         second_deriv = 0e0_dp

         do i=1, num_inflaton
           do j=1, num_inflaton
             if (i==j) then
               m2_matrix(i,j) = m2_V(i)
             else
               m2_matrix(i,j) = lambda2
             end if

           end do
         end do

         second_deriv = m2_matrix

       case(13)

         print*, "MODECODE: Potential_choice=", Potential_choice
         call raise%fatal_code(&
           "This potential choice is broken for some reason.",&
           __FILE__, __LINE__)

       case(14)

         !Multifield step potential
         m2_V = 10.e0_dp**(vparams(1,:))
         location_phi = vparams(2,:)
         step_size = vparams(3,:)
         step_slope = vparams(4,:)

         second_deriv=0e0_dp

#define ARG ((phi(i)-location_phi(i))/step_slope(i))
         do i=1,size(phi)

           second_deriv(i,i) = m2_V(i)*&
              (1.0e0_dp + step_size(i)* tanh( ARG )) +&
              2.0e0_dp*m2_V(i)*phi(i)*step_size(i)* &
              (MySech( ARG ))**(2)/step_slope(i) - &
              m2_V(i)*phi(i)**2 * (step_size(i)/step_slope(i)**2) *tanh(ARG) *&
              (MySech(ARG))**(2)

         end do
#undef ARG

       case(15)

         print*, "MODECODE: Potential_choice=", Potential_choice
         call raise%fatal_code(&
           "This potential choice is broken for some reason.",&
           __FILE__, __LINE__)

       case(16)
         ! (1/p) lambda_i |phi_i|^p --- N-monomial

         p_exp = vparams(2,1)
         m2_V = vparams(1,:)
         second_deriv = 0e0_dp
         do ii=1,size(phi)
           second_deriv(ii,ii) =(p_exp-1.0e0_dp)*m2_V(ii)*abs(phi(ii))**(p_exp-2.0e0_dp)
         end do

         !Regularize around phi=0
         do ii=1,size(phi)
           if (abs(phi(ii))<1e-5) &
             second_deriv(ii,ii)=0e0_dp
         end do

       case(17)
         ! Generalized axions

         A_i = vparams(1,:)
         B_ij = vparams(2:num_inflaton+1,:)

         twopi = 2.0e0_dp*pi
         fourpi_sq = 4.0e0_dp*pi**2

         !Identity matrix
         delta=0e0_dp
         do ii=1, num_inflaton
           delta(ii,ii) = 1.0e0_dp
         end do

         second_deriv=0e0_dp
         do alpha=1, size(A_i)
           do beta=1, size(A_i)
             second_deriv(alpha, beta) = &
               -fourpi_sq*A_i(alpha)*delta(alpha,beta)*cos(twopi*phi(alpha)) &
               +fourpi_sq*(B_ij(alpha,beta) + B_ij(beta,alpha))*&
                  cos(twopi*phi(alpha)-twopi*phi(beta))*&
                  (1.0e0_dp-delta(alpha,beta))

              if (alpha/=beta) cycle

              do ii=1,size(A_i)
                if (ii==alpha) cycle
                second_deriv(alpha, beta) = second_deriv(alpha, beta)&
                  -fourpi_sq*(B_ij(alpha,ii) + B_ij(ii,alpha))*&
                  cos(twopi*phi(alpha)-twopi*phi(ii))*delta(alpha,beta)
              end do

           end do
         end do


       case default

         print*, "MODECODE: potential_choice =", potential_choice
         call raise%fatal_code(&
           "Need to set second_deriv for this potential choice.",&
            __FILE__, __LINE__)

       end select
       !END MULTIFIELD
    end if

  END FUNCTION d2Vdphi2

#undef HEAVY
#undef PHI_I
#undef DELTAPHI
#undef EXPTERM

  !Only used to get the 3rd order SR parameters for the SR approximation
  !of the scalar running, alpha_s
  function d3Vdphi3(phi) result(third_deriv)
    real(dp), intent(in) :: phi(:)
    real(dp) :: third_deriv(size(phi),size(phi),size(phi))
    real(dp) :: m2_V(size(phi))
    integer :: ii

    real(dp) :: p_exp

    real(dp), dimension(size(phi)) :: location_phi, step_size, step_slope

    third_deriv = 0e0_dp

    select case(potential_choice)
    case(1)
      m2_V = 10.e0_dp**(vparams(1,:))
      third_deriv=0e0_dp
    case(14)

      !Multifield step potential
      m2_V = 10.e0_dp**(vparams(1,:))
      location_phi = vparams(2,:)
      step_size = vparams(3,:)
      step_slope = vparams(4,:)

#define ARG ((phi(ii)-location_phi(ii))/step_slope(ii))
      third_deriv=0e0_dp
      do ii=1,size(phi)

        third_deriv(ii,ii,ii)=&
          (3.0e0_dp*step_size(ii)*m2_V(ii)/step_slope(ii))*&
            (MySech(ARG)**2) &
          - (step_size(ii)*m2_V(ii)*phi(ii)**2/step_size(ii)**3)*&
            (MySech(ARG)**4) &
          -(6.0e0_dp*step_size(ii)*m2_V(ii)*phi(ii)/step_slope(ii)**2)*&
            (MySech(ARG)**2)*tanh(ARG) &
          +(2.0e0_dp*step_size(ii)*m2_V(ii)*phi(ii)**2/step_slope(ii)**3)*&
            (MySech(ARG)**2)*(tanh(ARG)**2)

      end do
#undef ARG

    case(15)

      print*, "MODECODE: Potential_choice=", Potential_choice
      call raise%fatal_code(&
        "This potential choice is broken for some reason.",&
        __FILE__, __LINE__)

    case (16)
      ! (1/p) lambda_i |phi_i|^p --- N-monomial

      p_exp = vparams(2,1)
      m2_V = 10.e0_dp**(vparams(1,:))

      third_deriv = 0e0_dp
      do ii=1,size(phi)
        third_deriv(ii,ii,ii) = (p_exp-2.0e0_dp)*(p_exp-1.0e0_dp)*&
          m2_V(ii)*abs(phi(ii))**(p_exp-3.0e0_dp)
      end do

    case default

      print*, "MODECODE: potential_choice =", potential_choice
      call raise%fatal_code(&
        "Need to set third derivative for this potential choice.",&
        __FILE__, __LINE__)

    end select

  end function d3Vdphi3

  FUNCTION initialphi(phi0)
    !
    !     Sets initial value of phi (depending on potential, may use
    !     either the user-specified value phi0 or something derived
    !     from the potential parameters)
    !
    real(dp), INTENT(IN) :: phi0(:)
    real(dp) :: initialphi(size(phi0))

    real(dp) :: phii(size(phi0))
    real(dp) :: Ninit, finv, lambda, mu, phesq
    real(dp) :: x1, x2


    Ninit = 70.e0_dp

    if (size(phi0) .gt. 1) then
       phii = phi0 ! MULTIFIELD
    else  !SINGLE FIELD
       select case(potential_choice)
       case(1)
          phii = 2.e0_dp*sqrt(Ninit+0.5e0_dp)
       case(2)
          finv = 1.e0_dp/(10.e0_dp**vparams(2,1))
          phii = 2.e0_dp/finv*asin(exp(-0.5e0_dp*Ninit*finv*finv)/ &
               sqrt(1.e0_dp+0.5e0_dp*finv*finv))
       case(3)
          phii = sqrt(8.e0_dp*(Ninit+1.e0_dp))
       case(4)
          phii = sqrt(2.e0_dp*Ninit+0.5e0_dp)
       case(5)
          phii = sqrt(4.e0_dp/3.e0_dp*Ninit+2.e0_dp/9.e0_dp)
       case(6)
          lambda = 10.e0_dp**vparams(1,1)
          mu = 10.e0_dp**vparams(2,1)
          x1 = lambda**4/mu
          phesq = ((sqrt(2.e0_dp)*x1)**(-4./3.)+1.e0_dp/(4.e0_dp*x1))**(-0.5)
          if (vparams(1,1)<-3.e0_dp) then
             phii = sqrt(phesq/(2.e0_dp**1.5*Ninit/sqrt(phesq)+1.e0_dp))
          else
             x2 = 4.e0_dp*Ninit + 2.e0_dp*x1/phesq + 0.5e0_dp*phesq
             phii = sqrt(x2)*sqrt(1.e0_dp-sqrt(1.e0_dp-4.e0_dp*x1/x2/x2))
          end if
       case default
          phii = phi0
       end select
    end if

    initialphi = phii

  END FUNCTION initialphi

  FUNCTION getEps(phi,dphi)
    !
    !     Returns epsilon given phi and dphi/dalpha, alpha=ln(a)
    !     for single field, slowroll parameter epsilon_H = 2 M_pl^2 [(dH/dphi)/H]^2
    !     for canonical multi-field,
    !
    real(dp) :: getEps
    real(dp), INTENT(IN) :: phi(:), dphi(:)

    !MULTIFIELD
    getEps = 0.5e0_dp*(M_Pl)**2 * dot_product(dphi,dphi)

    if (getEps >=3.0e0_dp) then
      print*, "MODECODE: epsilon =", getEps
      print*, "MODECODE: phi =", phi
      print*, "MODECODE: dphi =", dphi

      call raise%fatal_cosmo(&
        "Epsilon is >3.0 in subroutine getEps. &
         This means H is complex.  &
         (This is not a universe I'd like to live in.) &
         This error might arise if there is a large separation &
         in scales (stiff problem) and the integrator walks &
         to a bad position in parameter space. &
         Try reducing the integration stepsize or use the DVODE integrator.", &
        __FILE__, __LINE__)

    end if

    !END MULTIFIELD

  END FUNCTION getEps


  FUNCTION getH(phi,dphi)
    !
    !     Returns H given phi and dphi/dalpha
    !
    real(dp) :: getH
    real(dp), INTENT(IN) :: phi(:), dphi(:)

    getH=pot(phi)/3.0e0_dp/M_Pl**2 / &
      (1.0e0_dp - dot_product(dphi, dphi)/6.0e0_dp/M_Pl**2)

    if (getH < 0.0e0_dp) then
      call raise%fatal_cosmo(&
        "H is complex. &
        Try smaller stepsize in integrator.",&
        __FILE__,__LINE__)
    else
      getH = sqrt(getH)
    end if

  END FUNCTION getH

  !For when using t-integrator
  FUNCTION getH_with_t(phi,dphidt)
    !
    !     Returns H given phi and dphi/dt
    !
    real(dp) :: getH_with_t
    real(dp), INTENT(IN) :: phi(:), dphidt(:)

    ! MULTIFIELD
    getH_with_t= (pot(phi) + 0.5e0_dp*dot_product(dphidt, dphidt))/ &
      3.0e0_dp/M_pl**2

    if (getH_with_t < 0.0e0_dp) then
      call raise%fatal_cosmo(&
        "H is complex.  &
        Try smaller stepsize in integrator.",&
        __FILE__,__LINE__)
    else
      getH_with_t = sqrt(getH_with_t)
    end if

    ! MULTIFIELD

  END FUNCTION getH_with_t

  FUNCTION getEps_with_t(phi,dphi)
    !
    !     Returns epsilon given phi and dphi/dt
    !
    real(dp) :: getEps_with_t, hubble
    real(dp), INTENT(IN) :: phi(:), dphi(:)

    !MULTIFIELD
    hubble =getH_with_t(phi,dphi)
    getEps_with_t = 0.5e0_dp*(M_Pl)**2 * &
      dot_product(dphi,dphi)/hubble**2

    if (getEps_with_t >3.0e0_dp) then
      print*, "MODECODE: epsilon =", getEps_with_t
      print*, "MODECODE: phi =", phi
      print*, "MODECODE: dphi =", dphi

      call raise%fatal_cosmo(&
        "Epsilon is >3.0 in subroutine getEps_with_t. &
        This means H is complex.  &
        (This is not a universe I'd like to live in.) &
        This error might arise if there is a large separation &
        in scales (stiff problem) and the integrator walks &
        to a bad position in parameter space. &
        Try reducing the integration stepsize or use the DVODE integrator.", &
        __FILE__, __LINE__)
    end if
    !END MULTIFIELD


  END FUNCTION getEps_with_t


  FUNCTION getdHdalpha(phi,dphi)
    !
    !     Returns dH/dalpha given phi and dphi/dalpha
    !
    real(dp) :: getdHdalpha
    real(dp), INTENT(IN) :: phi(:), dphi(:)
    ! MULTIFIELD
    getdHdalpha = -dot_product(dphi, dphi) * getH(phi,dphi)/2.0e0_dp/M_Pl**2
    ! END MULTIFIELD
  END FUNCTION getdHdalpha


  FUNCTION getdepsdalpha(phi,dphi)
    !
    !    Returns depsilon/dalpha given phi and dphi/dalpha
    !    Gets this by differentiating Peiris et al Eq A3 (2003)
    !
    real(dp) :: getdepsdalpha, H, dHdalpha, eps
    real(dp), INTENT(IN) :: phi(:), dphi(:)

    H=getH(phi,dphi)
    dHdalpha=getdHdalpha(phi,dphi)
    eps=getEps(phi,dphi)

    ! MULTIFIELD
    getdepsdalpha=6.0e0_dp*dHdalpha/H*(1.0e0_dp-eps/3.0e0_dp) &
      -dot_product(dVdphi(phi), dphi)/(H**2*M_Pl**2)
    ! END MULTIFIELD

  END FUNCTION getdepsdalpha

  FUNCTION geteta(phi, dphi)
    !
    !    Return the eta parameter eta = deps/dalpha / eps
    !
    real(dp), INTENT(IN) :: phi(:), dphi(:)
    real(dp) :: geteta, eps

    eps = getEps(phi, dphi)
    geteta = getdepsdalpha(phi, dphi) / eps

  END FUNCTION geteta

  ! Powerspectrum for psi ptbs mode matrix, outputting full IJ matrix,
  !adiabatic P(k), and isocurv P(k); includes cross correlations
  subroutine powerspectrum(psi, dpsi, phi, dphi, scalefactor, power_spectrum, using_q)
    use internals
    use modpk_observables

    type(power_spectra), intent(inout) :: power_spectrum
    real(dp), dimension(:), intent(in) :: dphi, phi

    logical, optional, intent(in) :: using_q
    logical :: use_q

    ! Hacked matrix to vect
    complex(dp), dimension(:), intent(in) :: psi, dpsi
    real(dp), intent(in) :: scalefactor

    complex(dp), dimension(size(dphi),size(dphi)) :: power_matrix
    complex(dp), dimension(size(dphi),size(dphi)) :: d_power_matrix
    complex(dp), dimension(size(dphi),size(dphi)) :: cross_matrix
    real(dp) :: power_adiab
    real(dp) :: power_isocurv
    real(dp) :: power_pnad
    real(dp) :: power_cross
    real(dp) :: power_entropy

    real(dp) :: AAprod, ABprod, BAprod, BBprod
    type(KahanSum) :: pnad_sumAA, pnad_sumAB, pnad_sumBA, pnad_sumBB

    real(dp) :: H, Pdot

    !Pre-factors for Pnad P(k)
    real(dp), dimension(num_inflaton) :: A_vect, B_vect
    real(dp), dimension(num_inflaton) :: A_vect_adiab, B_vect_adiab

    ! size(dphi)=num_inflaton
    complex(kind=dp), dimension(size(dphi),size(dphi)) :: ptb_matrix,&
      dptb_matrix
    real(dp), dimension(size(dphi)-1) :: pk_iso_vect

    !proj along adiab dir
    real(dp), dimension(size(dphi)) :: omega_z
    !proj along isocurv dirs
    real(dp), dimension(size(dphi)-1,size(dphi)) :: s_iso, d_s_iso
    real(dp), parameter :: div_zero_tol=1e-20_dp

    real(dp) :: zeta2, phi_dot_0_scaled, prefactor
    integer :: numb_infl
    integer :: i, j, ll, kk

    real(dp) :: prod_exponent, phi_adiab, d_phi_adiab, hubble
    real(dp), dimension(size(s_iso)) :: phi_iso, d_phi_iso
    real(dp), dimension(size(phi)) :: d_omega_z, Vprime

    real(dp) :: power_total
    real(dp) :: power_pressure, power_press_cross, power_press_cross2, &
      power_press_adiab

    !Variable passed to powerspectrum
    !Psi_ij =a q_ij is default
    !q_ij optional
    if (present(using_q)) then
      use_q = using_q
    else
      use_q =.false.
    end if

    numb_infl=size(dphi)

    ! Convert hacked vector to matrix
    forall (i=1:numb_infl, j=1:numb_infl)
      ptb_matrix(i,j) = psi((i-1)*numb_infl+j)
      dptb_matrix(i,j) = dpsi((i-1)*numb_infl+j)
    end forall

    ! Make projection vector along adiabatic and isocurv directions
    ! NB: phi_dot_0_scaled = sqrt(2*epsilon) = phi_dot_0/H
    phi_dot_0_scaled = sqrt(dot_product(dphi,dphi))
    omega_z = dphi/phi_dot_0_scaled

    !If there's a major hierarchy in scales for the vector components, then
    !there can be a spurious component of the isocurvature direction
    !oriented along the adiabatic one.  When you're approaching the
    !adiabatic limit this can be the dominant contribution, so we force the
    !smallest term to be zero whenever adding it doesn't affect the value of
    !the norm.

    !call renormalize_remove_smallest(omega_z)


    !Cosmology params to make phi_adiab and d_phi_adiab
    hubble = getH(phi,dphi)
    Vprime = dVdphi(phi)
    !domega/dalpha
    d_omega_z = (1e0_dp/(hubble**2)*phi_dot_0_scaled)* &
      (omega_z*sum(omega_z*Vprime) - Vprime )

    phi_adiab = sum(omega_z*phi)
    !d phi_adiab/dalpha
    d_phi_adiab = sum(omega_z*dphi + d_omega_z*phi)


    !Build power matrices for fields, vels, and cross correls
    !All other quantities are background, so only power matrices dependent on
    !choice of \psi vs q ptb variable

    !NB: These are different for \psi_ij and q_ij variables
    !NB: Default is for \psi_ij, making the q_ij eqns a little more complicated
    !NB: \psi_ij = a*q_ij
    !NB: d\psi_ij = a*(q_ij + dq_ij)

    power_matrix=0e0_dp !in delta_phi
    d_power_matrix=0e0_dp !in psi
    cross_matrix=0e0_dp !in psi

    if (use_q) then
      !Don't divide out by scalefact
      prefactor= (k**3/2.0e0_dp/(pi**2))/(2e0_dp*k)
    else
      prefactor= (k**3/2.0e0_dp/(pi**2)/scalefactor**2)/(2e0_dp*k)
    end if

    do i=1,numb_infl; do j=1, numb_infl; do ll=1,numb_infl

      power_matrix(i,j) =power_matrix(i,j)+ &
        ptb_matrix(i,ll)*conjg(ptb_matrix(j,ll))*prefactor

      if (.not. use_q) then
        cross_matrix(i,j) =cross_matrix(i,j)+ &
          ptb_matrix(i,ll)*conjg(dptb_matrix(j,ll))*prefactor
      else
        cross_matrix(i,j) =cross_matrix(i,j)+ &
          ptb_matrix(i,ll)*conjg(ptb_matrix(j,ll) + dptb_matrix(j,ll))*prefactor
      end if

      if (.not. use_q) then
        d_power_matrix(i,j) =d_power_matrix(i,j)+&
          dptb_matrix(i,ll)*conjg(dptb_matrix(j,ll))*prefactor
      else
        d_power_matrix(i,j) =d_power_matrix(i,j)+&
          (ptb_matrix(i,ll)+dptb_matrix(i,ll))*&
          conjg(ptb_matrix(j,ll)+dptb_matrix(j,ll))*prefactor
      end if

    end do; end do; end do

    !------------------------------------------------------------
    !Adiabatic power spectrum
    power_adiab = dot_product(matmul(power_matrix,omega_z),omega_z)
    power_adiab = (1e0_dp/phi_dot_0_scaled**2)*power_adiab
    !------------------------------------------------------------

    !Isocurvature power spectra
    if (numb_infl>1) then

      !------------------------------------------------------------
      !Build the isocurv unit vects
      call build_isocurv_basis()

      !Vector of iso-spectra
      !Project power_matrix onto directs perpend to adiab direction
      pk_iso_vect = 0e0_dp
      do i=1,size(s_iso,1); do j=1, numb_infl; do ll=1,numb_infl
        pk_iso_vect(i) = pk_iso_vect(i) + &
          s_iso(i,j)*s_iso(i,ll)*power_matrix(j,ll)
      end do; end do; end do

      pk_iso_vect = pk_iso_vect*(1e0_dp/phi_dot_0_scaled**2)

      power_isocurv = 0e0_dp
      power_isocurv = sum(pk_iso_vect)

      !------------------------------------------------------------
      !Cross spectra between adiab and "projected" isocurvature
      power_cross = 0e0_dp
      do i =1, numb_infl; do j=1, numb_infl; do ll=1, size(s_iso,1)
        power_cross = power_cross + &
          omega_z(i)*s_iso(ll,j)*&
          (power_matrix(i,j) + power_matrix(j,i))
      end do; end do; end do
      power_cross = power_cross/(phi_dot_0_scaled**2)

      !------------------------------------------------------------

      !P(k) of total non-adiab pressure ptbs
      !dP_nad_i(k) = (1/a) Sum_j (A_i*Psi_ij + B_i*dPsi_ij)*\hat{a}_j
      !<P_nad P_nad*>
      power_pnad =0e0_dp

      A_vect = get_A_vect(phi,dphi)
      B_vect = get_B_vect(phi,dphi)

      AAprod = 0e0_dp
      ABprod = 0e0_dp
      BAprod = 0e0_dp
      BBprod = 0e0_dp

      do i=1,numb_infl; do j=1,numb_infl

        AAprod = AAprod +A_vect(i)*A_vect(j)*power_matrix(i,j)
        ABprod = ABprod +A_vect(i)*B_vect(j)*cross_matrix(i,j)
        BAprod = BAprod +B_vect(i)*A_vect(j)*conjg(cross_matrix(j,i))
        BBprod = BBprod +B_vect(i)*B_vect(j)*d_power_matrix(i,j)

      end do; end do
      power_pnad = (AAprod + BBprod) + (ABprod + BAprod)

      !Total pressure spectrum
      !<P P*>
      A_vect = get_A_vect_Ptotal(phi,dphi)
      B_vect = get_B_vect_Ptotal(phi,dphi)
      AAprod = 0e0_dp
      ABprod = 0e0_dp
      BAprod = 0e0_dp
      BBprod = 0e0_dp
      do i=1,numb_infl; do j=1,numb_infl
        AAprod = AAprod +A_vect(i)*A_vect(j)*power_matrix(i,j)
        ABprod = ABprod +A_vect(i)*B_vect(j)*cross_matrix(i,j)
        BAprod = BAprod +B_vect(i)*A_vect(j)*conjg(cross_matrix(j,i))
        BBprod = BBprod +B_vect(i)*B_vect(j)*d_power_matrix(i,j)
      end do; end do
      power_pressure = (AAprod + BBprod) + (ABprod + BAprod)

      !Total-adiab pressure cross-spectrum
      !<P_nad P_ad*>
      A_vect = get_A_vect(phi,dphi)
      B_vect = get_B_vect(phi,dphi)
      A_vect_adiab = get_A_vect_Padiab(phi,dphi)
      B_vect_adiab = get_B_vect_Padiab(phi,dphi)
      AAprod = 0e0_dp
      ABprod = 0e0_dp
      BAprod = 0e0_dp
      BBprod = 0e0_dp
      do i=1,numb_infl; do j=1,numb_infl
        AAprod = AAprod +A_vect(i)*A_vect_adiab(j)*power_matrix(i,j)
        ABprod = ABprod +A_vect(i)*B_vect_adiab(j)*cross_matrix(i,j)
        BAprod = BAprod +B_vect(i)*A_vect_adiab(j)*conjg(cross_matrix(j,i))
        BBprod = BBprod +B_vect(i)*B_vect_adiab(j)*d_power_matrix(i,j)
      end do; end do
      power_press_cross = (AAprod + BBprod) + (ABprod + BAprod)

      !<P_ad P_nad*>
      A_vect = get_A_vect(phi,dphi)
      B_vect = get_B_vect(phi,dphi)
      A_vect_adiab = get_A_vect_Padiab(phi,dphi)
      B_vect_adiab = get_B_vect_Padiab(phi,dphi)
      AAprod = 0e0_dp
      ABprod = 0e0_dp
      BAprod = 0e0_dp
      BBprod = 0e0_dp
      do i=1,numb_infl; do j=1,numb_infl
        AAprod = AAprod +A_vect_adiab(i)*A_vect(j)*power_matrix(i,j)
        ABprod = ABprod +A_vect_adiab(i)*B_vect(j)*cross_matrix(i,j)
        BAprod = BAprod +B_vect_adiab(i)*A_vect(j)*conjg(cross_matrix(j,i))
        BBprod = BBprod +B_vect_adiab(i)*B_vect(j)*d_power_matrix(i,j)
      end do; end do
      power_press_cross2 = (AAprod + BBprod) + (ABprod + BAprod)

      !<P_ad P_ad*>
      A_vect_adiab = get_A_vect_Padiab(phi,dphi)
      B_vect_adiab = get_B_vect_Padiab(phi,dphi)
      AAprod = 0e0_dp
      ABprod = 0e0_dp
      BAprod = 0e0_dp
      BBprod = 0e0_dp
      do i=1,numb_infl; do j=1,numb_infl
        AAprod = AAprod +A_vect_adiab(i)*A_vect_adiab(j)*power_matrix(i,j)
        ABprod = ABprod +A_vect_adiab(i)*B_vect_adiab(j)*cross_matrix(i,j)
        BAprod = BAprod +B_vect_adiab(i)*A_vect_adiab(j)*conjg(cross_matrix(j,i))
        BBprod = BBprod +B_vect_adiab(i)*B_vect_adiab(j)*d_power_matrix(i,j)
      end do; end do
      power_press_adiab = (AAprod + BBprod) + (ABprod + BAprod)

      !The values (AA + BB) --> -(AB+BA) as approaches adiab limit.
      !Taking diff of "large" numbs means large error in the difference
      !Check if power_pnad is smaller than DP accuracy and set to zero

      prod_exponent = abs(real(AAprod)) + abs(real(BBprod)) +&
        abs(real(ABprod)) + abs(real(BAprod))

      if (power_pnad>1e-60_dp) then

        prod_exponent = log10(prod_exponent/power_pnad)

        if( prod_exponent >14e0_dp) then
          !Reached DP accuracy, all subseq results are numer error
          power_pnad = 0e0_dp
        end if

      else
        power_pnad=0e0_dp
      endif


      !Rescale so spectrum for S=(H/Pdot)dP_nad - total entropy ptb
      H=getH(phi,dphi)
      Pdot=getPdot(phi,dphi)
      power_entropy = ((H/Pdot)**2)*power_pnad

    else
      power_spectrum%phi_ij  = 0e0_dp
      power_spectrum%adiab   = 0e0_dp
      power_spectrum%isocurv = 0e0_dp
      power_spectrum%pnad    = 0e0_dp
      power_spectrum%entropy = 0e0_dp
      power_spectrum%pressure = 0e0_dp
      power_spectrum%press_ad = 0e0_dp
      power_spectrum%cross_ad_iso = 0e0_dp
    end if

    power_spectrum%phi_ij  =  power_matrix
    power_spectrum%adiab   =  power_adiab
    power_spectrum%isocurv =  power_isocurv
    power_spectrum%pnad    =  power_pnad
    power_spectrum%entropy =  power_entropy
    power_spectrum%pressure =  power_pressure
    power_spectrum%press_ad =  power_press_adiab
    power_spectrum%cross_ad_iso =  power_cross


    contains



    !Basis perpend to adiab direction.  Gram-Schmidt
    subroutine build_isocurv_basis()

      real(dp), dimension(size(dphi),size(dphi)) :: spanning
      !field dirs
      real(dp), dimension(size(dphi),size(dphi)) :: phi_hat
      real(dp), dimension(size(dphi)-1,size(dphi)) :: field_basis
      real(dp), dimension(size(dphi)) :: phi_hat_temp
      real(dp) :: check_dot
      real(dp), dimension(size(dphi)-1) :: normalization, dnormalization
      integer :: i, j, adiab_index

      !Build field space vects
      do i=1,size(phi_hat,1); do j=1,size(phi_hat,2)
        if (i==j) then
          phi_hat(i,j) = 1e0_dp
        else
          phi_hat(i,j) = 0e0_dp
        end if
      end do; end do

      !Find field direction closest aligned to omega_z
      !This check makes sure get linearly indep field space vects
      check_dot=1e-30_dp
      adiab_index=0
      do i=1,size(phi_hat,1)
        if (abs(dot_product(phi_hat(i,:),omega_z(:)))>check_dot) then
          check_dot = abs(dot_product(phi_hat(i,:),omega_z(:)))
          adiab_index=i
        end if
      end do

      if (adiab_index ==0) then
        call raise%fatal_cosmo(&
        "It appears that no field space directions have projection &
        along the adiab direction.",&
        __FILE__, __LINE__)
      end if

      !Set the "most adiab" dir to 1
      if (adiab_index /= 1) then
        phi_hat_temp = phi_hat(1,:)
        phi_hat(1,:) = phi_hat(adiab_index,:)
        phi_hat(adiab_index,:) = phi_hat_temp
      end if

      spanning = 0e0_dp
      spanning(1,:) = omega_z/norm(omega_z)

      s_iso = 0e0_dp

      do i=2,size(spanning,1)

        spanning(i,:) = phi_hat(i,:)
        do j=1, i-1
          spanning(i,:) = spanning(i,:) - &
            projection(spanning(i,:),spanning(j,:))
        end do

        if (norm(spanning(i,:)) > div_zero_tol) then
          s_iso(i-1,:) = spanning(i,:)/norm(spanning(i,:))
        else
          print*, "MODECODE: i=",i
          call raise%fatal_cosmo(&
            "The component spanning(i,:) has zero norm.",&
            __FILE__, __LINE__)
        end if

        !If there's a major hierarchy in scales for the vector components, then
        !there can be a spurious component of the isocurvature direction
        !oriented along the adiabatic one.  When you're approaching the
        !adiabatic limit this can be the dominant contribution, so we force the
        !smallest term to be zero whenever adding it doesn't affect the value of
        !the norm.

        !call renormalize_remove_smallest(s_iso(i-1,:))


        if (abs(dot_product(omega_z,s_iso(i-1,:)))>1e-12 .or.&
          isnan(abs(dot_product(omega_z,s_iso(i-1,:))))) then

          write(*,*), "omega_z.s_iso =",dot_product(omega_z,s_iso(i-1,:))," for i=",i-1
          call raise%fatal_cosmo(&
            "The isocurvature projection has a large adiabatic component.",&
            __FILE__, __LINE__)

        end if

      end do

    end subroutine build_isocurv_basis

    !Takes a vector and sets to zero those elements that are so small they're
    !not affecting the normalization

    !If there's a major hierarchy in scales for the vector components, then
    !there can be a spurious component of the isocurvature direction
    !oriented along the adiabatic one.  When you're approaching the
    !adiabatic limit this can be the dominant contribution, so we force the
    !smallest term to be zero whenever adding it doesn't affect the value of
    !the norm.
    subroutine renormalize_remove_smallest(vect)

      real(dp), dimension(:), intent(inout) :: vect
      real(dp), dimension(size(vect)) :: vect_temp
      real(dp) :: check
      integer :: smallest_index
      integer :: j

      check= 1e10_dp
      vect_temp = vect
      smallest_index= 1
      do j=1,size(vect)
        if (abs(vect(j)) < check) then
          check= abs(vect(j))
          smallest_index= j
        end if
      end do
      vect_temp(smallest_index) = 0e0_dp

      !1e-14 for double-precision round-off error
      if (abs(norm(vect)-norm(vect_temp))/abs(norm(vect)) < 1e-14) then
        vect = vect_temp/norm(vect_temp)
      end if

    end subroutine renormalize_remove_smallest

    !Functions for calculating the total non-adiabatic pressure perturbation
    !matrices

    !dP_nad(k) = (1/a) Sum_j (A_i*Psi_ij + B_i*dPsi_ij)*\hat{a}_j
    function get_A_vect(phi,dphi) result(A)

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp), dimension(size(phi)) :: A
      real(dp) :: H, H2, firstterm, Vdot, c2, eps, gamm
      real(dp), dimension(size(phi)) :: secondterm, thirdterm, Vprime

      H=getH(phi,dphi)
      H2=H**2
      Vprime=dVdphi(phi)
      !c2 = getcs2(phi,dphi)
      Vdot = H*sum(Vprime*dphi)
      eps = geteps(phi,dphi)
      !gamm=(1.0e0_dp+c2)/(1.0e0_dp-c2)

      firstterm = (2.0e0_dp/3.0e0_dp/H2/sum(dphi*dphi))
      secondterm=(-3.0e0_dp*H2*sum(dphi*dphi) - (Vdot/H))*Vprime
      thirdterm = Vdot*H*dphi*(1.0e0_dp + 0.5e0_dp*sum(dphi*dphi))
      A=firstterm*(secondterm+thirdterm)

    end function get_A_vect


    function get_B_vect(phi,dphi) result(B)

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp), dimension(size(phi)) :: B
      real(dp) :: firstterm, H, H2, Vdot, eps, c2
      real(dp), dimension(size(phi)) :: Vprime

      H=getH(phi,dphi)
      H2=H**2
      c2 = getcs2(phi,dphi)

      B = (1.0e0_dp-c2)*H2*dphi

    end function get_B_vect

    function get_A_vect_Ptotal(phi,dphi) result(A)

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp), dimension(size(phi)) :: A
      real(dp) :: H, H2, eps
      real(dp), dimension(size(phi)) :: Vprime

      H=getH(phi,dphi)
      H2=H**2
      eps = geteps(phi,dphi)
      Vprime =dVdphi(phi)

      A = -H2*dphi - eps*H2*dphi - Vprime

    end function get_A_vect_Ptotal

    function get_B_vect_Ptotal(phi,dphi) result(B)

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp), dimension(size(phi)) :: B
      real(dp) :: H

      H=getH(phi,dphi)

      B = (H**2)*dphi

    end function get_B_vect_Ptotal

    function get_A_vect_Padiab(phi,dphi) result(A)

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp), dimension(size(phi)) :: A
      real(dp) :: H, H2, eps, cs2
      real(dp), dimension(size(phi)) :: Vprime

      H=getH(phi,dphi)
      H2=H**2
      eps = geteps(phi,dphi)
      Vprime =dVdphi(phi)
      cs2 = getcs2(phi,dphi)

      A = cs2*(-H2*dphi - eps*H2*dphi + Vprime)

    end function get_A_vect_Padiab

    function get_B_vect_Padiab(phi,dphi) result(B)

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp), dimension(size(phi)) :: B
      real(dp) :: H, cs2

      H=getH(phi,dphi)
      cs2 = getcs2(phi,dphi)

      B = (H**2)*dphi*cs2

    end function get_B_vect_Padiab


    function getPdot(phi,dphi) result(Pdot)

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp) :: Pdot, H, H2
      real(dp), dimension(size(phi)) :: Vprime

      Vprime=dVdphi(phi)

      H=getH(phi,dphi)
      H2=H**2

      Pdot = -H*sum(3.0e0_dp*H2*dphi*dphi + 2e0_dp*Vprime*dphi)

    end function getPdot


    function getrhodot(phi,dphi) result(rhodot)

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp) :: rhodot, H, H2

      H=getH(phi,dphi)

      rhodot = -3.0e0_dp*H**3*sum(dphi*dphi)

    end function getrhodot


    function getcs2(phi,dphi) result(cs2)

      real(dp), dimension(:), intent(in) :: phi, dphi
      real(dp) :: cs2

      cs2 = getpdot(phi,dphi)/getrhodot(phi,dphi)

    end function getcs2


  end subroutine powerspectrum


  !Projection of v orthogonally onto line spanned by u
  pure function projection(v,u) result(proj)

    real(dp), dimension(:), intent(in) :: v, u
    real(dp), dimension(size(v)) :: proj

    proj = u*(dot_product(v,u)/dot_product(u,u))

  end function projection


  !Norm of v
  pure function norm(v)

    real(dp), dimension(:), intent(in) :: v
    real(dp) :: norm

    norm = sqrt(dot_product(v,v))

  end function norm


    !END MULTIFIELD


  pure FUNCTION tensorpower(v, a)
    USE internals
    real(dp) :: tensorpower
    real(dp), INTENT(IN) :: a
    COMPLEX(KIND=DP), intent(in) :: v

    !MULTIFIELD: Calculates P_h(k) given v, a
    tensorpower = abs(v)**2/(2.0e0_dp*k) / a**2 * (k**3)*4.0e0_dp/(PI**2)/(M_Pl**2)
    !END MULTIFIELD

  END FUNCTION tensorpower

  pure FUNCTION zpower(u_zeta, dsigma, a)
    USE internals
    real(dp) :: zpower
    real(dp), INTENT(IN) :: dsigma
    real(dp), INTENT(IN) :: a
    COMPLEX(KIND=DP), INTENT(IN) :: u_zeta

    zpower = abs(u_zeta**2)/dsigma**2/a**2 /(2.0e0_dp*k) * (k**3)/(2.0e0_dp*PI**2)

  END FUNCTION zpower


  function trace_d2logVdphi2(phi) result(trace)
    !
    !     Returns trace of d^2V/dPhi^2 given phi
    !     Used in calculating bundle exp_scalar
    !

    real(dp), intent(in) :: phi(:)
    real(dp) :: trace, V_ab(size(phi),size(phi)), V_a(size(phi)), V
    integer :: i

    V_ab = d2Vdphi2(phi)
    V_a = dVdphi(phi)
    V = pot(phi)

    trace = 0e0_dp
    do i=1,size(phi)
      trace = trace + V_ab(i,i)/V - (V_a(i)/V)**2
    end do

  end function trace_d2logVdphi2

  !Calc bundle exp_scalar by integrating tr(d2Vdphi2) to efold by Riemann sum
  subroutine bundle_exp_scalar(self,phi, efold)

    class(bundle) :: self
    real(dp), intent(in) :: phi(:), efold
    real(dp) :: dN

    dN = efold - self%N

    self%dlogThetadN= self%dlogThetadN - &
      dN*trace_d2logVdphi2(phi)

    self%exp_scalar=exp(self%dlogThetadN)

    self%N=efold

  end subroutine bundle_exp_scalar


  pure FUNCTION MySech(x)
    real(dp), intent(in)  :: x
    real(dp) :: MySech

    IF(ABS(x).GT.40.0e0_dp) THEN
       MySech=0.0e0_dp
    ELSE
       MySech=1.0e0_dp/COSH(x)
    END IF
  END FUNCTION MySech

  !Checks to see if H^2=V/(3-eps) is stable, ie, if H>0 and if V~0, then it will
  !return FALSE for unstable.
  subroutine stability_check_on_H(stable,phi,dphi,using_t)
    logical, intent(inout) :: stable
    real(dp), dimension(:), intent(in) :: phi, dphi
    logical, intent(in) :: using_t
    real(dp) :: eps, V

    V=pot(phi)

    if (using_t) then
      eps = getEps_with_t(phi, dphi)
    else
      eps = getEps(phi, dphi)
    end if

    stable = (3.0e0_dp -eps >1.0e-6)


  end subroutine stability_check_on_H


END MODULE potential
