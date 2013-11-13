!Types and constants module for CMB and combined dwarf likelihoods
!Pat Scott patscott@physics.mcgill.ca
!April 28 2013

module types

implicit none

Type BFset
  double precision :: u, d, s, c, b, t, e, mu, tau 
  double precision :: Vtoe, Vtomu, Vtotau, Wto4, Zto4, g, w, z, h 
  double precision :: gamgam, zgam, hgam, invis_1, invis_2
end Type BFset

double precision, parameter :: pi=3.14159265358979323846264338328d0, twopi=2*pi, fourpi=4*pi
double precision, parameter :: root2 = 1.41421356237309504880168872421d0, sqrt2 = root2
double precision, parameter :: log2 = 0.693147180559945309417232121458d0
double precision, parameter :: GeVm2tocm3sm1 = 1.16619266162d-17

double precision, parameter :: mu=2.3d-3, md=4.8d-3, mst=0.095d0, mc=1.275d0, mb=4.18d0, mt=173.d0
double precision, parameter :: me=0.510998928d-3, mmu=0.1056583715d0, mtau=1.77682d0, mh=125.d0, mW=80.385d0, mZ=91.188d0
double precision, parameter :: Gh0=4.07d-3, lh=0.129d0, v0=246.22d0, as = 0.12d0
double precision :: factorisation_maxmass = mt 

integer :: feedback = 2

contains

  subroutine clear(BFs)
  Type(BFset) :: BFs 

  BFs%u=0.d0
  BFs%d=0.d0
  BFs%s=0.d0
  BFs%c=0.d0
  BFs%b=0.d0
  BFs%t=0.d0
  BFs%e=0.d0
  BFs%mu=0.d0
  BFs%tau=0.d0 
  BFs%Vtoe=0.d0
  BFs%Vtomu=0.d0
  BFs%Vtotau=0.d0
  BFs%g=0.d0
  BFs%w=0.d0
  BFs%z=0.d0
  BFs%h=0.d0 
  BFs%gamgam=0.d0
  BFs%zgam=0.d0
  BFs%hgam=0.d0
  BFs%invis_1=0.d0
  BFs%invis_2=0.d0
  BFs%Wto4=0.d0
  BFs%Zto4=0.d0

  end subroutine clear

end module types
