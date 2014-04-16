!module containing all global types definitons for variables
!This version August 2007 
!SuperBayes Package - by Roberto Ruiz de Austri (rruiz@delta.ft.uam.es) and Roberto Trotta (rxt@astro.ox.ac.uk)

module parameters


 integer, parameter:: lmax = 200

!**********************************************
! NEW TYPES DEFINITIONS
!**********************************************

!--------- Constant SM particle masses (GeV) ----------
 Type particleMasses
   real*8 :: e       = 0.d0
   real*8 :: mu      = 0.d0
   real*8 :: tau     = 0.d0
   real*8 :: nu1     = 0.d0
   real*8 :: nu2     = 0.d0
   real*8 :: nu3     = 0.d0
   real*8 :: u       = 3.d-3
   real*8 :: d       = 6.75d-3
   real*8 :: s       = 0.1175d0
   real*8 :: c       = 1.224d0
   real*8 :: b       = 4.9d0
   real*8 :: W       = 80.392d0
   real*8 :: Z       = 91.1876d0
   real*8 :: h       = 125.d0
   real*8 :: gam     = 0.d0
   real*8 :: g       = 0.d0
   real*8 :: hidden1 = 0.d0
   real*8 :: hidden2 = 0.d0
 end Type particleMasses
 Type (particleMasses) :: MassData

!--------- Input parameters ---------
 Type Input_params
   real*8 :: mx, sigmav, br_uubar, br_ddbar, br_ssbar, br_bbar, br_ttbar, br_ww, br_zz, br_hh, br_ee, br_mumu, br_tautau, br_ccbar, br_gg, br_zgam, br_gamgam, br_invis1, br_invis2
   real*8 :: alpha, beta, gamma, rho, delta, fpn, apn
 end Type Input_params

 !Shifted initialisation of these into paramdef; the same should be done eventually for all names
 character(LEN=lmax), dimension(21):: DMNames
 logical, dimension(21) :: DMUsed
  
 Type Nuisance_params
   real (8) :: a0, e0, delta1, delta2, mtop 
 end Type Nuisance_params

 character(LEN=lmax), dimension(5):: NuisanceNames = (/ 'A_0', 'E_0', '\delta_1', 'delta_2', 'm_t (GeV)'/)
 logical, dimension(5), parameter :: NuisanceUsed = (/ .true., .true., .true., .true., .true./)

 character(LEN=lmax),dimension(26) :: ParamsNames

!Background parameters are split into two groups: 
!- Nonlinear grid parameters
!- Template parameters: linear/power-law parameters that are not included in the grid, but simply rescale the flux maps and thus are implemented analytically.

Type Grid_params
     real*8 :: bg_xs
end Type Grid_params

character(LEN=lmax),dimension(2) :: GridNames = (/ 'bg_xs', 'BG_param2' /)
logical, dimension(2), parameter :: GridUsed = (/ .true., .false. /)

Type Template_params
     real*8 :: X_CO_1,X_CO_2,X_CO_3,i_e,n_e,i_p,n_p,alpha_CR,beta_CR
end Type Template_params

character(LEN=lmax),dimension(9) :: TemplateNames = (/ 'X_CO_1', 'X_CO_2', 'X_CO_3', 'index_e', 'norm_e', 'index_p', 'norm_p', 'CR_alpha', 'CR_beta' /)
logical, dimension(9), parameter :: TemplateUsed = (/ .true., .true., .true., .true., .true., .true., .true., .true., .true. /)

!Even when only scanning over a smaller number of params these arrays can be as big as you like
! character(LEN=lmax), dimension(3) :: BgNames = (/ 'bg_Z', 'bg_S', 'bg_T' /)
! logical, dimension(3), parameter :: BgUsed = (/ .true., .true., .true. /)

Type PS_params
     real*8 :: PS_l,PS_b,N_0,E_0,PS_alpha,PS_beta,Inv_E_c
end Type PS_params

character(LEN=lmax),dimension(7) :: PSNames = (/ 'l(PS)', 'b(PS)', 'N_0', 'E_0', 'alpha_PS', 'beta_PS', 'E_c^{-1}' /)
logical, dimension(7), parameter :: PSUsed = (/ .true., .true., .true., .true., .true., .true., .true. /)

!--------- DM parameters ---------------

 Type AdditionalOut
  real(8) :: J, br_tautau
 End Type AdditionalOut

 character(LEN=lmax), dimension(2) :: AdditionalOutNames= (/ 'BR_{\tau^+\tau^-}' , 'J(\psi = 0)'/)

!Indirect Detection rates

 Type Gammas_In
   logical :: GCPoissonian, replaceGCObsWithBG
   real (8) :: cospsi0, delta_gamma,   &
            egam, egath, thmax, efluxes_i, efluxes_f
   real (8) :: GCBF
   character(len=lmax) :: IRFs
   integer :: nbins, GC_outerPix, GC_corePix
   integer, pointer :: GCDims(:)
   real (8), pointer :: GC_model(:,:,:), GC_angSep(:,:), GC_Ebins(:,:), GC_pixArea  ! PS - added for passing into dminterface.  A similar thing will need to be done
                                       !      with other things read from FITS files, like pixel size / anglular coordinates, etc.
 End Type Gammas_In 


 Type ID_In
    Type(Gammas_In) :: gammas 
 End Type ID_In
 

 Type Gammas_Out
    real (8), dimension(1000) :: fluxgadiff, fluxdiff, Ekin
    real (8) :: fluxgac=0d0
 End Type Gammas_Out 


 Type ID_Out
    Type(Gammas_Out) :: gammas 
 End Type ID_Out


  character(LEN=lmax), dimension(3) :: IDNames= (/ 'd \Phi_\gamma/dE', '\Phi_\gamma', 'E_{fluxes} (GeV)'/)



!---------- Error messages from SoftSusy and Dark Susy about physical 
!and unphysical points ------------
 Type Error_Out
  integer :: haerr = 0
 end Type Error_Out


!All output from DS is collected in one single type
!for practical purposes when manipulating things
 Type DM
    !IF YOU CHANGE THIS REMEMBER TO ADJUST WriteParams FORMAT BELOW !
    Type(AdditionalOut) :: AddOut
    Type(ID_In)  :: ID_in
    Type(ID_Out) :: ID
    Type(Error_Out) :: Err
    logical :: PrintOutParams
 end Type DM


!-------- Some flags to decide what to compute ------------


 Type ID_Flags_gammas
   logical :: halo_fix
   logical :: gadiff, gac, GC_region, dwarfs, cta_gc
 end Type ID_Flags_gammas


 Type Flags
   logical :: Debug
   integer :: use_data
   logical :: ID_predict, ID_predict_gamma
   Type(ID_Flags_gammas) :: ID_Flags_gamma
  End Type Flags

!---------- Format for likelihood data structure ----------
 integer, parameter :: Gaussian = 1, LowerLimit = 2, UpperLimit = 3

 Type LikeDatum
    real(8) :: mu, sigma, tau
    integer :: datum_type
    logical :: tau_percent
 end Type LikeDatum

 type GCR_data
  integer :: n ! number of entries
  real(8), dimension(:), allocatable :: E_low
  real(8), dimension(:), allocatable :: E_high
  real(8), dimension(:), allocatable :: E_mean
  real(8), dimension(:), allocatable :: value
  real(8), dimension(:), allocatable :: err_minus
  real(8), dimension(:), allocatable :: err_plus
  real(8), dimension(:), allocatable :: theory_prediction  
 end type GCR_data

end module parameters


