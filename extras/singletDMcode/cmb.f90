!CMB likelihood module, as per Cline & Scott arXiv:1301.5908
!Pat Scott patscott@physics.mcgill.ca
!April 28 2013

module cmb_poor

use types

implicit none

integer, parameter, public :: WMAP7 = 1, Planck_predicted = 2
integer, parameter, private ::  num_gamma_masses = 12, num_other_masses = 7
double precision, parameter, private :: c1sq = 4.64*4.64

double precision :: sigmav_ref(2), lambda1(2)
data   sigmav_ref   / 3.2d-27, 2.d-27 /
data   lambda1      / 0.279, 3.16 /

double precision, dimension(num_other_masses) :: e_WMAP7, mu_WMAP7, tau_WMAP7, Vtoe_WMAP7, Vtomu_WMAP7, Vtotau_WMAP7, &
 uds_WMAP7, c_WMAP7, b_WMAP7, t_WMAP7, g_WMAP7, W_WMAP7, Z_WMAP7,h_WMAP7 
double precision, dimension(num_gamma_masses) :: gam_WMAP7
data   e_WMAP7      / 0.74, 0.65, 0.59, 0.56, 0.56, 0.55048696818614440, 0.54619728368586662/ &
       mu_WMAP7     / 0.28, 0.25, 0.23, 0.21, 0.21, 0.20760424713083339, 0.20814896859714740/ &
       tau_WMAP7    / 0.22, 0.21, 0.19, 0.19, 0.18, 0.18235699079366086, 0.18265823098169290/ &
       Vtoe_WMAP7   / 0.77, 0.69, 0.62, 0.56, 0.56, 0.56114667803739304, 0.56061051693353070/ &
       Vtomu_WMAP7  / 0.28, 0.26, 0.23, 0.21, 0.20, 0.19841209351474948, 0.19805737291789355/ &
       Vtotau_WMAP7 / 0.22, 0.22, 0.20, 0.18, 0.18, 0.17696854486925812, 0.17566888895227023/ &
       uds_WMAP7    / 0.33, 0.31, 0.30, 0.28, 0.27, 0.26343805954956578, 0.25869417309463771/ &
       c_WMAP7      / 0.33, 0.32, 0.30, 0.28, 0.27, 0.26639797423228173, 0.26130894816932621/ &
       b_WMAP7      / 0.34, 0.32, 0.31, 0.29, 0.27, 0.26501217722533466, 0.25969412848404894/ &
       t_WMAP7      / 0.00, 0.00, 0.27, 0.27, 0.26, 0.25092998079929663, 0.24734726512346480/ &
       gam_WMAP7    / 0.58, 0.50, 0.48, 0.46, 0.44, 0.49, 0.526, 0.59, 0.56, 0.55, 0.53803319487302259, 0.52319415664461100/ &
       g_WMAP7      / 0.33, 0.32, 0.30, 0.28, 0.27, 0.26631720764984612, 0.26154521053382257/ &
       W_WMAP7      / 0.00, 0.26, 0.26, 0.25, 0.24, 0.23665455855080544, 0.24060609325527879/ &
       Z_WMAP7      / 0.00, 0.25, 0.25, 0.23, 0.22, 0.22023499420108308, 0.22006301367229317/ &
       h_WMAP7      / 0.00, 0.00, 0.28, 0.28, 0.26, 0.24938630562999617, 0.24326537673746756/ 

double precision, dimension(num_other_masses) :: e_Planck_predict, mu_Planck_predict, tau_Planck_predict, Vtoe_Planck_predict, &
 Vtomu_Planck_predict, Vtotau_Planck_predict, uds_Planck_predict, c_Planck_predict, b_Planck_predict, &
 t_Planck_predict, g_Planck_predict, W_Planck_predict, Z_Planck_predict, h_Planck_predict
double precision, dimension(num_gamma_masses) :: gam_Planck_predict
data   e_Planck_predict      / 0.79, 0.70, 0.63, 0.59, 0.59, 0.58666309985572440, 0.58209159334603733/ &
       mu_Planck_predict     / 0.30, 0.27, 0.25, 0.23, 0.22, 0.22125910968580018, 0.22183068209842741/ &
       tau_Planck_predict    / 0.24, 0.23, 0.21, 0.20, 0.20, 0.19438420566781126, 0.19466794852476557/ &
       Vtoe_Planck_predict   / 0.83, 0.74, 0.66, 0.60, 0.60, 0.59802765943447433, 0.59744853721891222/ &
       Vtomu_Planck_predict  / 0.30, 0.28, 0.25, 0.23, 0.21, 0.21147986515708125, 0.21107580534956624/ &
       Vtotau_Planck_predict / 0.24, 0.23, 0.21, 0.20, 0.19, 0.18867131375475305, 0.18722956123979526/ &
       uds_Planck_predict    / 0.35, 0.33, 0.32, 0.30, 0.29, 0.28105916819228838, 0.27586933881657616/ &
       c_Planck_predict      / 0.36, 0.34, 0.32, 0.30, 0.29, 0.28421563247681658, 0.27865553159693757/ &
       b_Planck_predict      / 0.36, 0.34, 0.33, 0.31, 0.29, 0.28274864935558808, 0.27693930032626796/ &
       t_Planck_predict      / 0.00, 0.00, 0.29, 0.29, 0.28, 0.26776566877865765, 0.26378659785334652/ &
       gam_Planck_predict    / 0.62, 0.54, 0.51, 0.49, 0.49, 0.52, 0.56, 0.62, 0.60, 0.58, 0.57338700746082460, 0.55758630890557481/ &
       g_Planck_predict      / 0.35, 0.34, 0.32, 0.30, 0.29, 0.28412394660512508, 0.27890516522069225/ &
       W_Planck_predict      / 0.00, 0.28, 0.28, 0.26, 0.25, 0.25241543559128954, 0.25652326884930715/ &
       Z_Planck_predict      / 0.00, 0.27, 0.27, 0.25, 0.24, 0.23492578550548332, 0.23463205177412441/ &
       h_Planck_predict      / 0.00, 0.00, 0.30, 0.30, 0.28, 0.26607808064240362, 0.25938852381532218/

double precision :: normal_masses(num_other_masses), gamma_masses(num_gamma_masses) 
data   normal_masses / 10.d0, 30.d0, 100.d0, 300.d0, 1000.d0, 3.d3, 1.d4 /
data   gamma_masses  / 10.d0, 30.d0, 40.d0, 50.d0, 60.d0, 70.d0, 80.d0, 100.d0, 300.d0, 1000.d0, 3.d3, 1.d4 /

contains


  double precision function cmblike_ann(BFs, sigmav, mass, experiment)

  Type(BFset), intent(IN) :: BFs
  double precision, intent(IN) :: sigmav, mass
  integer, intent(IN) :: experiment
  double precision :: f_eff

  if (mass .gt. normal_masses(num_other_masses)) then
    cmblike_ann = 0.d0
    return
  endif

  !Build f_eff by adding the portion from each channel
  f_eff = 0.d0
  select case (experiment)
  case (WMAP7)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%e,      e_WMAP7)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%mu,     mu_WMAP7)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%tau,    tau_WMAP7)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%Vtoe,   Vtoe_WMAP7)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%Vtomu,  Vtomu_WMAP7)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%Vtotau, Vtotau_WMAP7)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%u,      uds_WMAP7)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%d,      uds_WMAP7)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%s,      uds_WMAP7)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%c,      c_WMAP7)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%b,      b_WMAP7)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%t,      t_WMAP7)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%gamgam, gam_WMAP7, gamma=.true.)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%g,      g_WMAP7)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%W,      W_WMAP7)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%Z,      Z_WMAP7)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%h,      h_WMAP7)
  case (Planck_predicted)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%e,      e_Planck_predict)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%mu,     mu_Planck_predict)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%tau,    tau_Planck_predict)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%Vtoe,   Vtoe_Planck_predict)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%Vtomu,  Vtomu_Planck_predict)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%Vtotau, Vtotau_Planck_predict)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%u,      uds_Planck_predict)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%d,      uds_Planck_predict)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%s,      uds_Planck_predict)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%c,      c_Planck_predict)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%b,      b_Planck_predict)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%t,      t_Planck_predict)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%gamgam, gam_Planck_predict, gamma=.true.)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%g,      g_Planck_predict)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%W,      W_Planck_predict)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%Z,      Z_Planck_predict)
    f_eff = addToFeff(f_eff, sigmav, mass, BFs%h,      h_Planck_predict)
  end select

  cmblike_ann = -0.5 * f_eff*f_eff * c1sq * lambda1(experiment) * (sigmav/sigmav_ref(experiment))**2 / (mass*mass)

  end function cmblike_ann


  double precision function addToFeff(f_eff, sigmav, mass, ri, fvector, gamma)

  double precision, intent(IN) :: f_eff, sigmav, mass, ri, fvector(:)
  logical, optional, intent(IN) :: gamma
  logical :: isgamma = .false.
  double precision, allocatable :: masses(:)
  integer :: masslen, i

  if (abs(ri) .le. epsilon(ri)) then
    addToFeff = f_eff
    return
  endif

  if (present(gamma)) isgamma = gamma

  masslen = merge(num_gamma_masses, num_other_masses, isgamma)
  allocate(masses(masslen))
  if (isgamma) then
    masses = gamma_masses
  else 
    masses = normal_masses
  endif
  
  if (mass .lt. masses(1) .or. mass .gt. masses(masslen)) then
    write(*,*) 'Error: requested mass in addToFeff is out of range.  Quitting...'
  end if

  do i = 1, masslen
    if (masses(i) .ge. mass) exit
  end do

  addToFeff = f_eff + ri*((fvector(i)-fvector(i-1))*(log10(mass)-log10(masses(i-1)))/(log10(masses(i))-log10(masses(i-1))) + fvector(i-1))

  end function addToFeff


end module cmb_poor
