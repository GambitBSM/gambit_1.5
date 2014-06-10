!Driver program for scalar singlet indirect detection calculations
!Pat Scott patscott@physics.mcgill.ca
!April 28 2013

module singletID_helper

use types
use parameters
use cmb_poor
use dwarfs_combined
use cta_gc
use jimlib

implicit none
double precision :: ms, deltalnlike, relicdens_target, ddsensitivity_target, lglam, p
double precision :: coord_topleft(2), coord_bottomleft(2), coord_topright(2)

contains


  double precision function singletCL_current(loglam)

  Type(BFSet) :: BFs  
  double precision :: dwarflnlike, cmblnlike, sigv, loglam, lambda_hs, f_rel

  if (feedback .gt. 2) then
    write(*,*) 'loglam: ', loglam  
  endif
  lambda_hs = 10.d0**(loglam)
  BFs = fractions(lambda_hs,ms)
  sigv = sigmav(lambda_hs,ms) * GeVm2tocm3sm1
  f_rel = relic_f(lambda_hs,ms)
  dwarflnlike = dwarflike(BFs, sigv*f_rel*f_rel, ms, original)
  cmblnlike = cmblike_ann(BFs, sigv*f_rel*f_rel, ms, WMAP7)
  if (feedback .gt. 1) then
    write(*,*) 'Annihilation xsec: ', sigv  
    write(*,*) 'Relic density fraction: ', f_rel
    write(*,*) 'Dwarf log-likelihood:   ', dwarflnlike
    write(*,*) 'CMB log-likelihood:     ', cmblnlike
    write(*,*) 'Branching fractions:  '
    write(*,*) BFs
    write(*,*) 'Total lnlikelihood:  ', dwarflnlike+cmblnlike
    write(*,*) 'Sought lnlikelihood: ', -deltalnlike
  endif

  singletCL_current = dwarflnlike + deltalnlike + cmblnlike

  end function singletCL_current


  double precision function singletCL_future(loglam)

  Type(BFSet) :: BFs  
  double precision :: dwarflnlike, cmblnlike, ctalnlike, sigv, loglam, lambda_hs, f_rel

  if (feedback .gt. 2) then
    write(*,*) 'loglam: ', loglam  
  endif
  lambda_hs = 10.d0**(loglam)
  BFs = fractions(lambda_hs,ms)
  sigv = sigmav(lambda_hs,ms) * GeVm2tocm3sm1
  f_rel = relic_f(lambda_hs,ms)
  dwarflnlike = dwarflike(BFs, sigv*f_rel*f_rel*sqrt(20./3.), ms, original)
  cmblnlike = cmblike_ann(BFs, sigv*f_rel*f_rel, ms, Planck_predicted)
  ctalnlike = ctalike(BFs, sigv*f_rel*f_rel, ms)
  if (feedback .gt. 1) then
    write(*,*) 'Annihilation xsec: ', sigv  
    write(*,*) 'Relic density fraction: ', f_rel
    write(*,*) 'Dwarf log-likelihood:   ', dwarflnlike
    write(*,*) 'CMB log-likelihood:     ', cmblnlike
    write(*,*) 'CTA log-likelihood:     ', ctalnlike
    write(*,*) 'Branching fractions:  '
    write(*,*) BFs
    write(*,*) 'Total lnlikelihood:  ', dwarflnlike+cmblnlike+ctalnlike
    write(*,*) 'Sought lnlikelihood: ', -deltalnlike
  endif

  singletCL_future = dwarflnlike + ctalnlike + cmblnlike + deltalnlike

  end function singletCL_future


  double precision function singletCL_current_flip(mass)

  Type(BFSet) :: BFs  
  double precision :: dwarflnlike, cmblnlike, sigv, lambda_hs, f_rel, mass

  if (feedback .gt. 2) then
    write(*,*) 'ms: ', mass  
  endif
  lambda_hs = 10.d0**(lglam)
  BFs = fractions(lambda_hs,mass)
  sigv = sigmav(lambda_hs,mass) * GeVm2tocm3sm1
  f_rel = relic_f(lambda_hs,mass)
  dwarflnlike = dwarflike(BFs, sigv*f_rel*f_rel, mass, original)
  cmblnlike = cmblike_ann(BFs, sigv*f_rel*f_rel, mass, WMAP7)
  if (feedback .gt. 1) then
    write(*,*) 'Annihilation xsec: ', sigv  
    write(*,*) 'Relic density fraction: ', f_rel
    write(*,*) 'Dwarf log-likelihood:   ', dwarflnlike
    write(*,*) 'CMB log-likelihood:     ', cmblnlike
    write(*,*) 'Branching fractions:  '
    write(*,*) BFs
    write(*,*) 'Total lnlikelihood:  ', dwarflnlike+cmblnlike
    write(*,*) 'Sought lnlikelihood: ', -deltalnlike
  endif

  singletCL_current_flip = dwarflnlike + cmblnlike + deltalnlike

  end function singletCL_current_flip


  double precision function singletCL_future_flip(mass)

  Type(BFSet) :: BFs  
  double precision :: dwarflnlike, cmblnlike, ctalnlike, sigv, lambda_hs, f_rel, mass

  if (feedback .gt. 2) then
    write(*,*) 'ms: ', mass  
  endif
  lambda_hs = 10.d0**(lglam)
  BFs = fractions(lambda_hs,mass)
  sigv = sigmav(lambda_hs,mass) * GeVm2tocm3sm1
  f_rel = relic_f(lambda_hs,mass)
  dwarflnlike = dwarflike(BFs, sigv*f_rel*f_rel*sqrt(20./3.), mass, original)
  cmblnlike = cmblike_ann(BFs, sigv*f_rel*f_rel, mass, Planck_predicted)
  ctalnlike = ctalike(BFs, sigv*f_rel*f_rel, mass)
  if (feedback .gt. 1) then
    write(*,*) 'Annihilation xsec: ', sigv  
    write(*,*) 'Relic density fraction: ', f_rel
    write(*,*) 'Dwarf log-likelihood:   ', dwarflnlike
    write(*,*) 'CMB log-likelihood:     ', cmblnlike
    write(*,*) 'CTA log-likelihood:     ', ctalnlike
    write(*,*) 'Branching fractions:  '
    write(*,*) BFs
    write(*,*) 'Total lnlikelihood:  ', dwarflnlike+cmblnlike+ctalnlike
    write(*,*) 'Sought lnlikelihood: ', -deltalnlike
  endif

  singletCL_future_flip = dwarflnlike + ctalnlike + cmblnlike + deltalnlike

  end function singletCL_future_flip


  double precision function singletCL_CTA(loglam)

  Type(BFSet) :: BFs  
  double precision :: ctalnlike, sigv, loglam, lambda_hs, f_rel

  if (feedback .gt. 2) then
    write(*,*) 'loglam: ', loglam  
  endif
  lambda_hs = 10.d0**(loglam)

  BFs = fractions(lambda_hs,ms)
  sigv = sigmav(lambda_hs,ms) * GeVm2tocm3sm1
  f_rel = relic_f(lambda_hs,ms)
  if (feedback .gt. 1) then
    write(*,*) 'Annihilation xsec: ', sigv  
    write(*,*) 'Relic density fraction: ', f_rel
    write(*,*) 'Branching fractions:  '
    write(*,*) BFs
  endif
  ctalnlike = ctalike(BFs, sigv*f_rel*f_rel, ms)
  if (feedback .gt. 1) then
    write(*,*) 'CTA log-likelihood:   ', ctalnlike
    write(*,*) 'Sought lnlikelihood: ', -deltalnlike
  endif

  singletCL_CTA = ctalnlike + deltalnlike

  end function singletCL_CTA


  double precision function singletCL_CTA_flip(mass)

  Type(BFSet) :: BFs  
  double precision :: ctalnlike, sigv, lambda_hs, f_rel, mass

  if (feedback .gt. 2) then
    write(*,*) 'ms: ', mass  
  endif
  lambda_hs = 10.d0**(lglam)
  BFs = fractions(lambda_hs,mass)
  sigv = sigmav(lambda_hs,mass) * GeVm2tocm3sm1
  f_rel = relic_f(lambda_hs,mass)
  ctalnlike = ctalike(BFs, sigv*f_rel*f_rel, mass)
  if (feedback .gt. 1) then
    write(*,*) 'Annihilation xsec: ', sigv  
    write(*,*) 'Relic density fraction: ', f_rel
    write(*,*) 'CTA log-likelihood:     ', ctalnlike
    write(*,*) 'Branching fractions:  '
    write(*,*) BFs
    write(*,*) 'Sought lnlikelihood: ', -deltalnlike
  endif

  singletCL_CTA_flip = ctalnlike + deltalnlike

  end function singletCL_CTA_flip


  double precision function singletCL_dwarf(loglam)

  Type(BFSet) :: BFs  
  double precision :: dwarflnlike, sigv, loglam, lambda_hs, f_rel

  if (feedback .gt. 2) then
    write(*,*) 'loglam: ', loglam  
  endif
  lambda_hs = 10.d0**(loglam)

  BFs = fractions(lambda_hs,ms)
  sigv = sigmav(lambda_hs,ms) * GeVm2tocm3sm1
  f_rel = relic_f(lambda_hs,ms)
  if (feedback .gt. 1) then
    write(*,*) 'Annihilation xsec: ', sigv  
    write(*,*) 'Relic density fraction: ', f_rel
    write(*,*) 'Branching fractions:  '
    write(*,*) BFs
  endif
  dwarflnlike = dwarflike(BFs, sigv*f_rel*f_rel*sqrt(20./3.), ms, original)
  if (feedback .gt. 1) then
    write(*,*) 'Dwarf log-likelihood:   ', dwarflnlike
    write(*,*) 'Sought lnlikelihood: ', -deltalnlike
  endif

  singletCL_dwarf = dwarflnlike + deltalnlike

  end function singletCL_dwarf


  double precision function singletCL_dwarf_flip(mass)

  Type(BFSet) :: BFs  
  double precision :: dwarflnlike, sigv, lambda_hs, f_rel, mass

  if (feedback .gt. 2) then
    write(*,*) 'ms: ', mass  
  endif
  lambda_hs = 10.d0**(lglam)
  BFs = fractions(lambda_hs,mass)
  sigv = sigmav(lambda_hs,mass) * GeVm2tocm3sm1
  f_rel = relic_f(lambda_hs,mass)
  dwarflnlike = dwarflike(BFs, sigv*f_rel*f_rel*sqrt(20./3.), mass, original)
  if (feedback .gt. 1) then
    write(*,*) 'Annihilation xsec: ', sigv  
    write(*,*) 'Relic density fraction: ', f_rel
    write(*,*) 'Dwarf log-likelihood:     ', dwarflnlike
    write(*,*) 'Branching fractions:  '
    write(*,*) BFs
    write(*,*) 'Sought lnlikelihood: ', -deltalnlike
  endif

  singletCL_dwarf_flip = dwarflnlike + deltalnlike

  end function singletCL_dwarf_flip


  double precision function singletCL_cmb(loglam)

  Type(BFSet) :: BFs  
  double precision :: cmblnlike, sigv, loglam, lambda_hs, f_rel

  if (feedback .gt. 2) then
    write(*,*) 'loglam: ', loglam  
  endif
  lambda_hs = 10.d0**(loglam)

  BFs = fractions(lambda_hs,ms)
  sigv = sigmav(lambda_hs,ms) * GeVm2tocm3sm1
  f_rel = relic_f(lambda_hs,ms)
  if (feedback .gt. 1) then
    write(*,*) 'Annihilation xsec: ', sigv  
    write(*,*) 'Relic density fraction: ', f_rel
    write(*,*) 'Branching fractions:  '
    write(*,*) BFs
  endif
  cmblnlike = cmblike_ann(BFs, sigv*f_rel*f_rel, ms, Planck_predicted)
  if (feedback .gt. 1) then
    write(*,*) 'CMB log-likelihood:   ', cmblnlike
    write(*,*) 'Sought lnlikelihood: ', -deltalnlike
  endif

  singletCL_cmb = cmblnlike + deltalnlike

  end function singletCL_cmb


  double precision function singletCL_cmb_flip(mass)

  Type(BFSet) :: BFs  
  double precision :: cmblnlike, sigv, lambda_hs, f_rel, mass

  if (feedback .gt. 2) then
    write(*,*) 'ms: ', mass  
  endif
  lambda_hs = 10.d0**(lglam)
  BFs = fractions(lambda_hs,mass)
  sigv = sigmav(lambda_hs,mass) * GeVm2tocm3sm1
  f_rel = relic_f(lambda_hs,mass)
  cmblnlike = cmblike_ann(BFs, sigv*f_rel*f_rel, mass, Planck_predicted)
  if (feedback .gt. 1) then
    write(*,*) 'Annihilation xsec: ', sigv  
    write(*,*) 'Relic density fraction: ', f_rel
    write(*,*) 'CMB log-likelihood:     ', cmblnlike
    write(*,*) 'Branching fractions:  '
    write(*,*) BFs
    write(*,*) 'Sought lnlikelihood: ', -deltalnlike
  endif

  singletCL_cmb_flip = cmblnlike + deltalnlike

  end function singletCL_cmb_flip


  double precision function relicdens_x(loglam)

  double precision :: loglam, lambda_hs, f_rel

  if (feedback .gt. 2) then
    write(*,*) 'loglam: ', loglam  
  endif
  lambda_hs = 10.d0**(loglam)
  f_rel = relic_f(lambda_hs,ms)
  if (feedback .gt. 1) then
    write(*,*) 'Relic density fraction: ', f_rel
  endif

  relicdens_x = f_rel - relicdens_target

  end function relicdens_x


  double precision function ddsensitivity_x(t)

  double precision :: t, lambda_hs, ddsens, f_rel

  if (feedback .gt. 2) then
    write(*,*) 'p: ', p  
  endif

  !Make the transformation
  ms = tp2ms(t,p)
  lglam = tp2lglam(t,p)

  lambda_hs = 10.d0**(lglam)
  f_rel = relic_f(lambda_hs,ms)
  ddsens = xenon_limit(lambda_hs,ms)
  if (feedback .gt. 1) then
    write(*,*) 'Required boost to XENON-100 sensitivity: ', ddsens/f_rel
  endif

  ddsensitivity_x = ddsens/f_rel - ddsensitivity_target

  end function ddsensitivity_x


  double precision function ddsensitivity_x_traditional(ms)

  double precision :: ms, lambda_hs, ddsens, f_rel

  if (feedback .gt. 2) then
    write(*,*) 'ms: ', ms  
  endif

  lambda_hs = 10.d0**(lglam)
  f_rel = relic_f(lambda_hs,ms)
  ddsens = xenon_limit(lambda_hs,ms)
  if (feedback .gt. 1) then
    write(*,*) 'Required boost to XENON-100 sensitivity: ', ddsens/f_rel
  endif

  ddsensitivity_x_traditional = ddsens/f_rel - ddsensitivity_target

  end function ddsensitivity_x_traditional


  double precision function ddsensitivity_x_traditional2(loglam)

  double precision :: loglam, lambda_hs, ddsens, f_rel

  if (feedback .gt. 2) then
    write(*,*) 'loglam: ', loglam  
  endif

  lambda_hs = 10.d0**(loglam)
  f_rel = relic_f(lambda_hs,ms)
  ddsens = xenon_limit(lambda_hs,ms)
  if (feedback .gt. 1) then
    write(*,*) 'Required boost to XENON-100 sensitivity: ', ddsens/f_rel
  endif

  ddsensitivity_x_traditional2 = ddsens/f_rel - ddsensitivity_target

  end function ddsensitivity_x_traditional2


  double precision function tp2ms(t,p)
  double precision, intent(IN) :: t, p
  tp2ms = coord_topleft(1) + &
          p * abs (coord_bottomleft(1) - coord_topleft(1)) + &
          t * abs (coord_topright(1) - coord_topleft(1))
  end function


  double precision function tp2lglam(t,p)
  double precision, intent(IN) :: t, p
  tp2lglam = coord_bottomleft(2) + &
             p * abs (coord_topleft(2) - coord_bottomleft(2)) + &
             t * abs (coord_topright(2) - coord_topleft(2))
  end function


  double precision function singletCL_CTA_limits(logsigv)

  Type(BFSet) :: BFs   
  double precision :: logsigv, sigv

  sigv = 10.**logsigv
  call clear(BFs)
  BFs%b=1.d0
  
  singletCL_CTA_limits = ctalike(BFs, sigv, ms) + deltalnlike

  end function singletCL_CTA_limits


end module singletID_helper


program singletID

use singletID_helper
use parameters
use jimlib

implicit none

  double precision :: zbrent, CL, erfinv, lhs, t, width=0.2d0
  double precision :: hitol = 1.d-6, tol=1.d-3, msmin = 45.d0, msmax = 5.d3, lmmin=-4.0d0, lmmax=0.9d0
  integer :: i, nms = 1000, nlm=150, XeCL, Xe_astro
  logical old
  common/oldb/old
  common/XeCLb/XeCL,Xe_astro

  XeCL = 90
  Xe_astro = 0
  old = .false. ! use new or old Xenon100 limit
  feedback = 0


  call read_data
  call init_ctalike

  open(11,file='output/singletID_invis_wid.dat')
  do i = 1, nms
    ms = msmin * 10.d0**(dble(i-1)/dble(nms-1) * (log10(mh*0.5d0-1.d-10) - log10(msmin)))
    lhs = lm_max(ms,0.19d0)
    write(11,*) ms, log10(lhs) 
  enddo
  close(11)

  open(11,file='output/singletID_invis_wid_future.dat')
  do i = 1, nms
    ms = msmin * 10.d0**(dble(i-1)/dble(nms-1) * (log10(mh*0.5d0-1.d-10) - log10(msmin)))
    lhs = lm_max(ms,0.05d0)
    write(11,*) ms, log10(lhs) 
  enddo
  close(11)

stop

  CL = 68.3d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nlm
    lglam = lmmin + dble(i-1)/dble(nlm-1) * (lmmax - lmmin)
    ms = zbrent(singletCL_current_flip,msmin,msmax,tol)
    if (ms .lt. 45.d0) exit
    write(*,*) ms, lglam
    msmin = ms - 3.d0
    msmax = ms + 3.d0
  enddo
  ms = 45.d0
  lhs = zbrent(singletCL_current,lglam-0.1d0,lglam,tol)
  write(*,*) ms, lhs


stop

  CL = 90.0d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nlm
    lglam = lmmin + dble(i-1)/dble(nlm-1) * (lmmax - lmmin)
    ms = zbrent(singletCL_future_flip,msmin,msmax,tol)
    if (ms .lt. 45.d0) exit
    write(*,*) ms, lglam
    msmin = ms - 3.d0
    msmax = ms + 3.d0
  enddo
  ms = 45.d0
  lhs = zbrent(singletCL_future,lglam-0.1d0,lglam,tol)
  write(*,*) ms, lhs

stop

  CL = 68.3d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nlm
    lglam = lmmin + dble(i-1)/dble(nlm-1) * (lmmax - lmmin)
    ms = zbrent(singletCL_future_flip,msmin,msmax,tol)
    if (ms .lt. 45.d0) exit
    write(*,*) ms, lglam
    msmin = ms - 3.d0
    msmax = ms + 3.d0
  enddo
  ms = 45.d0
  lhs = zbrent(singletCL_future,lglam-0.1d0,lglam,tol)
  write(*,*) ms, lhs


stop

  CL = 90.0d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nlm
    lglam = lmmin + dble(i-1)/dble(nlm-1) * (lmmax - lmmin)
    ms = zbrent(singletCL_CTA_flip,msmin,msmax,tol)
    if (ms .lt. 45.d0) exit
    write(*,*) ms, lglam
    msmin = ms - 3.d0
    msmax = ms + 3.d0
  enddo
  ms = 45.d0
  lhs = zbrent(singletCL_CTA,lglam-0.1d0,lglam,tol)
  write(*,*) ms, lhs

stop

  CL = 90.0d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nlm
    lglam = lmmin + dble(i-1)/dble(nlm-1) * (lmmax - lmmin)
    ms = zbrent(singletCL_cmb_flip,msmin,msmax,tol)
    if (ms .lt. 45.d0) exit
    write(*,*) ms, lglam
    msmin = ms - 3.d0
    msmax = ms + 3.d0
  enddo
  ms = 45.d0
  lhs = zbrent(singletCL_cmb,lglam-0.1d0,lglam,tol)
  write(*,*) ms, lhs

stop

  CL = 90.0d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nlm
    lglam = lmmin + dble(i-1)/dble(nlm-1) * (lmmax - lmmin)
    ms = zbrent(singletCL_dwarf_flip,msmin,msmax,tol)
    if (ms .lt. 45.d0) exit
    write(*,*) ms, lglam
    msmin = ms - 3.d0
    msmax = ms + 3.d0
  enddo
  ms = 45.d0
  lhs = zbrent(singletCL_dwarf,lglam-0.1d0,lglam,tol)
  write(*,*) ms, lhs


stop

  open(11,file='output/singletID_1sigma_current.dat')
  CL = 68.3d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nms
    ms = msmin * 10.d0**(dble(i-1)/dble(nms-1) * (log10(msmax) - log10(msmin)))
    lhs = zbrent(singletCL_current,lmmin,lmmax,tol)
    write(11,*) ms, lhs
    write(*,*) ms, lhs
    lmmin = lhs - width
    lmmax = lhs + width
  enddo
  close(11)

  open(11,file='output/singletID_90_future.dat')
  CL = 90.0d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nms
    ms = msmin * 10.d0**(dble(i-1)/dble(nms-1) * (log10(msmax) - log10(msmin)))
    lhs = zbrent(singletCL_future,lmmin,lmmax,tol)
    write(11,*) ms, lhs 
    write(*,*) ms, lhs
   lmmin = lhs - width
    lmmax = lhs + width
  enddo
  close(11)

  open(11,file='output/singletID_1sigma_future.dat')
  CL = 68.3d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nms
    ms = msmin * 10.d0**(dble(i-1)/dble(nms-1) * (log10(msmax) - log10(msmin)))
    lhs = zbrent(singletCL_future,lmmin,lmmax,tol)
    write(11,*) ms, lhs
    write(*,*) ms, lhs
    lmmin = lhs - width
    lmmax = lhs + width
  enddo
  close(11)

  open(11,file='output/singletID_1sigma_current.dat')
  CL = 68.3d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nlm
    lglam = lmmin + dble(i-1)/dble(nlm-1) * (lmmax - lmmin)
    ms = zbrent(singletCL_current_flip,msmin,msmax,tol)
    write(11,*) ms, lglam 
    write(*,*) ms, lglam
    msmin = ms - 3.d0
    msmax = ms + 3.d0
  enddo
  close(11)

  open(11,file='output/singletID_90_future.dat')
  CL = 90.0d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nlm
    lglam = lmmin + dble(i-1)/dble(nlm-1) * (lmmax - lmmin)
    ms = zbrent(singletCL_future_flip,msmin,msmax,tol)
    write(11,*) ms, lglam 
    write(*,*) ms, lglam
    msmin = ms - 3.d0
    msmax = ms + 3.d0
  enddo
  close(11)

  open(11,file='output/singletID_1sigma_future.dat')
  CL = 68.3d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nlm
    lglam = lmmin + dble(i-1)/dble(nlm-1) * (lmmax - lmmin)
    ms = zbrent(singletCL_future_flip,msmin,msmax,tol)
    write(11,*) ms, lglam 
    write(*,*) ms, lglam
    msmin = ms - 3.d0
    msmax = ms + 3.d0
  enddo
  close(11)

stop

  open(11,file='output/singletID_90_CTA.dat')
  CL = 90.0d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nms
    ms = msmin * 10.d0**(dble(i-1)/dble(nms-1) * (log10(msmax) - log10(msmin)))
      lhs = zbrent(singletCL_CTA,lmmin,lmmax,tol)
    write(11,*) ms, lhs 
    write(*,*) ms, lhs 
    lmmin = lhs - width
    lmmax = lhs + width
  enddo
  close(11)

  open(11,file='output/singletID_90_CTA.dat')
  CL = 90.0d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nlm
    lglam = lmmin + dble(i-1)/dble(nlm-1) * (lmmax - lmmin)
    ms = zbrent(singletCL_CTA_flip,msmin,msmax,tol)
    write(11,*) ms, lglam 
    write(*,*) ms, lglam
    msmin = ms - 3.d0
    msmax = ms + 3.d0
  enddo
  close(11)

  open(11,file='output/singletID_90_cmb.dat')
  CL = 90.0d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nms
    ms = msmin * 10.d0**(dble(i-1)/dble(nms-1) * (log10(msmax) - log10(msmin)))
      lhs = zbrent(singletCL_cmb,lmmin,lmmax,tol)
    write(11,*) ms, lhs 
    write(*,*) ms, lhs 
    lmmin = lhs - width
    lmmax = lhs + width
  enddo
  close(11)

  open(11,file='output/singletID_90_cmb.dat')
  CL = 90.0d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nlm
    lglam = lmmin + dble(i-1)/dble(nlm-1) * (lmmax - lmmin)
    ms = zbrent(singletCL_cmb_flip,msmin,msmax,tol)
    write(11,*) ms, lglam 
    write(*,*) ms, lglam
    msmin = ms - 3.d0
    msmax = ms + 3.d0
  enddo
  close(11)

  open(11,file='output/singletID_90_dwarf.dat')
  CL = 90.0d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nms
    ms = msmin * 10.d0**(dble(i-1)/dble(nms-1) * (log10(msmax) - log10(msmin)))
      lhs = zbrent(singletCL_dwarf,lmmin,lmmax,tol)
    write(11,*) ms, lhs 
    write(*,*) ms, lhs 
    lmmin = lhs - width
    lmmax = lhs + width
  enddo
  close(11)

  open(11,file='output/singletID_90_dwarf.dat')
  CL = 90.0d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nlm
    lglam = lmmin + dble(i-1)/dble(nlm-1) * (lmmax - lmmin)
    ms = zbrent(singletCL_dwarf_flip,msmin,msmax,tol)
    write(11,*) ms, lglam 
    write(*,*) ms, lglam
    msmin = ms - 3.d0
    msmax = ms + 3.d0
  enddo
  close(11)

stop

  open(11,file='output/singletID_dd_100A.dat')
  ddsensitivity_target = 100.d0
  do i = 1, nlm
    lglam = lmmin + dble(i-1)/dble(nlm-1) * (lmmax - lmmin)
    ms = zbrent(ddsensitivity_x_traditional,msmin,msmax,tol)
    write(*,*) ms, lglam
    write(11,*) ms, lglam 
    msmin = ms - 1.d0
    msmax = ms + 1.d0
  enddo
  close(11)

  open(11,file='output/singletID_dd_100B.dat')
  ddsensitivity_target = 100.d0
  do i = 1, nlm
    lglam = lmmin + dble(i-1)/dble(nlm-1) * (lmmax - lmmin)
    ms = zbrent(ddsensitivity_x_traditional,msmin,msmax,tol)
    write(*,*) ms, lglam
    write(11,*) ms, lglam 
  enddo
  close(11)

  !open(11,file='output/singletID_dd_100C.dat')
  ddsensitivity_target = 100.d0
  do i = nms, 1, -1
    ms = msmin * 10.d0**(dble(i-1)/dble(nms-1) * (log10(msmax) - log10(msmin)))
    lhs = zbrent(ddsensitivity_x_traditional2,lmmin,lmmax,tol)
    write(*,*) ms, lhs
    lmmin = lhs - width
    lmmax = lhs + width
  enddo

stop

  open(11,file='output/singletID_dd_1C.dat')
  ddsensitivity_target = 1.d0
  coord_topleft = (/70.d0, lmmax/)
  coord_bottomleft = (/70.0d0, lmmin/)
  coord_topright = (/200.d0, lmmax/)
  do i = 1, nlm
    p = dble(i-1)/dble(nlm-1)
    t = zbrent(ddsensitivity_x,0.d0,1.d0,hitol)
    write(11,*) tp2ms(t,p), tp2lglam(t,p) 
    write(*,*) tp2ms(t,p), tp2lglam(t,p) 
  enddo
  close(11)

  open(11,file='output/singletID_dd_5C.dat')
  ddsensitivity_target = 5.d0
  coord_topleft = (/70.d0, lmmax/)
  coord_bottomleft = (/100.d0, lmmin-1.5d0/)
  coord_topright = (/370.d0, lmmax+1.5d0/)
  do i = 1, nlm
    p = dble(i-1)/dble(nlm-1)
    t = zbrent(ddsensitivity_x,0.d0,1.d0,hitol)
    write(11,*) tp2ms(t,p), tp2lglam(t,p) 
    write(*,*) tp2ms(t,p), tp2lglam(t,p) 
  enddo
  close(11)

  open(11,file='output/singletID_dd_20C.dat')
  ddsensitivity_target = 20.d0
  coord_topleft = (/65.d0, lmmax/)
  coord_bottomleft = (/125.d0, lmmin-2.5d0/)
  coord_topright = (/2500.d0, lmmax+4.0d0/)
  do i = 40, 216
    p = 0.1 + 0.9*(dble(i-1)/dble(nlm-1))
    t = zbrent(ddsensitivity_x,0.d0,1.d0,hitol)
    write(11,*) tp2ms(t,p), tp2lglam(t,p) 
    write(*,*) tp2ms(t,p), tp2lglam(t,p) 
  enddo
  coord_topleft = (/65.d0, lmmax/)
  coord_bottomleft = (/125.d0, lmmin-6.5d0/)
  coord_topright = (/1200.d0, lmmax+6.0d0/)
  do i = 130*2, nlm*2
    p = 0.1 + 0.9*(dble(i-1)/dble(499))
    t = zbrent(ddsensitivity_x,0.d0,1.d0,hitol)
    write(11,*) tp2ms(t,p), tp2lglam(t,p) 
    write(*,*) tp2ms(t,p), tp2lglam(t,p) 
  enddo
  close(11)

stop

  open(11,file='output/singletID_dd_1A.dat')
  ddsensitivity_target = 1.d0
  coord_topleft = (/45.d0, lmmax/)
  coord_bottomleft = (/45.d0, lmmin/)
  coord_topright = (/62.5d0, lmmax/)
  do i = 1, nlm
    p = dble(i-1)/dble(nlm-1)
    t = zbrent(ddsensitivity_x,0.d0,1.d0,tol)
    write(11,*) tp2ms(t,p), tp2lglam(t,p) 
    write(*,*) tp2ms(t,p), tp2lglam(t,p) 
  enddo
  close(11)


  open(11,file='output/singletID_dd_5A.dat')
  ddsensitivity_target = 5.d0
  coord_topleft = (/45.d0, lmmax/)
  coord_bottomleft = (/45.d0, lmmin/)
  coord_topright = (/62.5d0, lmmax/)
  do i = 1, nlm
    p = dble(i-1)/dble(nlm-1)
    t = zbrent(ddsensitivity_x,0.d0,1.d0,tol)
    write(11,*) tp2ms(t,p), tp2lglam(t,p)
    write(*,*) tp2ms(t,p), tp2lglam(t,p)  
  enddo
  close(11)

  open(11,file='output/singletID_dd_20A.dat')
  ddsensitivity_target = 20.d0
  coord_topleft = (/59.0d0, lmmax/)
  coord_bottomleft = (/59.0d0, lmmin/)
  coord_topright = (/62.d0, lmmax/)
  do i = 1, nlm
    p = dble(i-1)/dble(nlm-1)
    t = zbrent(ddsensitivity_x,0.d0,1.d0,tol)
    write(11,*) tp2ms(t,p), tp2lglam(t,p) 
    write(*,*) tp2ms(t,p), tp2lglam(t,p) 
  enddo
  close(11)

stop

  open(11,file='output/singletID_dd_1B.dat')
  ddsensitivity_target = 1.d0
  coord_topleft = (/64.d0, lmmax/)
  coord_bottomleft = (/64.d0, lmmin/)
  coord_topright = (/66.d0, lmmax/)
  do i = 1, nlm
    p = dble(i-1)/dble(nlm-1)
    t = zbrent(ddsensitivity_x,0.d0,1.d0,hitol)
    write(11,*) tp2ms(t,p), tp2lglam(t,p) 
    write(*,*) tp2ms(t,p), tp2lglam(t,p) 
  enddo
  close(11)

  open(11,file='output/singletID_dd_5B.dat')
  ddsensitivity_target = 5.d0
  coord_topleft = (/63.0d0, lmmax/)
  coord_bottomleft = (/63.d0, lmmin/)
  coord_topright = (/64.d0, lmmax/)
  do i = 1, nlm
    p = dble(i-1)/dble(nlm-1)
    t = zbrent(ddsensitivity_x,0.d0,1.d0,hitol)
    write(11,*) tp2ms(t,p), tp2lglam(t,p) 
    write(*,*) tp2ms(t,p), tp2lglam(t,p) 
  enddo
  close(11)

  open(11,file='output/singletID_dd_20B.dat')
  ddsensitivity_target = 20.d0
  coord_topleft = (/62.5d0, lmmax/)
  coord_bottomleft = (/62.5d0, lmmin/)
  coord_topright = (/64.d0, lmmax/)
  do i = 1, nlm
    p = dble(i-1)/dble(nlm-1)
    t = zbrent(ddsensitivity_x,0.d0,1.d0,hitol)
    write(11,*) tp2ms(t,p), tp2lglam(t,p) 
    write(*,*) tp2ms(t,p), tp2lglam(t,p) 
  enddo
  close(11)

stop


  open(11,file='output/singletID_relic_1.dat')
  relicdens_target = 1.d0
  do i = 1, nms
    ms = msmin * 10.d0**(dble(i-1)/dble(nms-1) * (log10(msmax) - log10(msmin)))
    lhs = zbrent(relicdens_x,lmmin,lmmax,tol)
    write(11,*) ms, lhs 
    write(*,*) ms, lhs 
    lmmin = lhs - width
    lmmax = lhs + width
  enddo
  close(11)

  open(11,file='output/singletID_relic_0.1.dat')
  relicdens_target = 1.d-1
  do i = 1, nms
    ms = msmin * 10.d0**(dble(i-1)/dble(nms-1) * (log10(msmax) - log10(msmin)))
    lhs = zbrent(relicdens_x,lmmin,lmmax,tol)
    write(11,*) ms, lhs 
    write(*,*) ms, lhs 
    lmmin = lhs - width
    lmmax = lhs + width
  enddo
  close(11)

  open(11,file='output/singletID_relic_0.01.dat')
  relicdens_target = 1.d-2
  do i = 1, nms
    ms = msmin * 10.d0**(dble(i-1)/dble(nms-1) * (log10(msmax) - log10(msmin)))
    lhs = zbrent(relicdens_x,lmmin,lmmax,tol)
    write(11,*) ms, lhs 
    write(*,*) ms, lhs 
    lmmin = lhs - width
    lmmax = lhs + width
  enddo
  close(11)

  open(11,file='output/singletID_1sigma_current.dat')
  CL = 68.3d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nms
    ms = msmin * 10.d0**(dble(i-1)/dble(nms-1) * (log10(msmax) - log10(msmin)))
    lhs = zbrent(singletCL_current,lmmin,lmmax,tol)
    write(11,*) ms, lhs
    write(*,*) ms, lhs
    lmmin = lhs - width
    lmmax = lhs + width
  enddo
  close(11)

  open(11,file='output/singletID_90_future.dat')
  CL = 90.0d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nms
    ms = msmin * 10.d0**(dble(i-1)/dble(nms-1) * (log10(msmax) - log10(msmin)))
    lhs = zbrent(singletCL_future,lmmin,lmmax,tol)
    write(11,*) ms, lhs 
    write(*,*) ms, lhs
    lmmin = lhs - width
    lmmax = lhs + width
  enddo
  close(11)

  open(11,file='output/singletID_1sigma_future.dat')
  CL = 68.3d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nms
    ms = msmin * 10.d0**(dble(i-1)/dble(nms-1) * (log10(msmax) - log10(msmin)))
    lhs = zbrent(singletCL_future,lmmin,lmmax,tol)
    write(11,*) ms, lhs
    write(*,*) ms, lhs
    lmmin = lhs - width
    lmmax = lhs + width
  enddo
  close(11)

stop

  !For reproducing published CTA annihilation x-sec limits
  open(11,file='90CL_CTAlimits_bb.dat')
  CL = 90.0d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nms
    ms = msmin * 10.d0**(dble(i-1)/dble(nms-1) * (log10(msmax) - log10(msmin)))
    lhs = zbrent(singletCL_CTA_limits,-35.d0,-10.d0,tol)
    write(11,*) ms, lhs 
  enddo
  close(11)

stop

  !For obtaining yield curves (requires hacky print statements in dwarflike)
  !do i = 1, nms
  !  ms = mw * 10.d0**(dble(i-1)/dble(nms-1) * (log10(mw) - log10(msmin)))
  !  BFs%w = 1.d0
  !  sigv = 3.d-26
  !  dwarflnlike = dwarflike(BFs, sigv, ms, original)
  !enddo

end program

