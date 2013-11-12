SUBROUTINE sigip_exp(mchi, sigip_bound, experiment)

!----------------------------------------------
! useful routine in case the spin-dependent neutralino-nucleon 
! wants to be used as a constraint. The "exp" flag regards 
! the experimental data to be applied. 
! exp: 1: CDMS, 2: ZEPLIN-I, 3: EDELWEISS-I, 4: XENON 
! author: Roberto Ruiz
!-----------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER :: DIM_CDMS = 41, DIM_ZEPL = 29, DIM_EDELW = 29, &
                      DIM_XENON = 61 
INTEGER, INTENT(IN) :: experiment
REAL(8), INTENT(IN) :: mchi
REAL(8), INTENT(OUT) :: sigip_bound

!local

REAL(8) :: lmchi, lsigip_bound
REAL(8), DIMENSION(DIM_CDMS) :: mchia_CDMS, sigipa_CDMS, y2a_CDMS, &
                                lmchia_CDMS, lsigipa_CDMS
REAL(8), DIMENSION(DIM_ZEPL) :: mchia_ZEPL, sigipa_ZEPL, y2a_ZEPL,& 
                                lmchia_ZEPL, lsigipa_ZEPL
REAL(8), DIMENSION(DIM_EDELW) :: mchia_EDELW, sigipa_EDELW, y2a_EDELW, &
                                lmchia_EDELW, lsigipa_EDELW
REAL(8), DIMENSION(DIM_XENON) :: mchia_XENON, sigipa_XENON, y2a_XENON, &
                                lmchia_XENON, lsigipa_XENON

lmchi = log10(mchi)

!data from CDMS

if(experiment .eq. 1) then


 mchia_CDMS = (/6.72950096, 7.07106781, 7.42997145, 7.80709182, 8.20335356, &
           8.61972821, 9.05723664, 9.51695153, 10., 11.2201845, 12.5892541, &
           14.1253754, 15.8489319, 19.9526231, 25.1188643, 31.6227766, &
           39.8107171, 50.1187234, 63.0957344, 79.4328235, 100., 125.892541, &
           158.489319, 199.526231, 251.188643, 316.227766, 398.107171, & 
           501.187234, 630.957344, 794.328235, 1000., 1258.92541, 1584.89319, &
           1995.26231, 2511.88643, 3162.27766, 3981.07171, 5011.87234,  &
           6309.57344 , 7943.28235, 10000./)

! MaxGap: 90% 
 sigipa_CDMS = (/17385.3047, 8859.76397, 2079.71203, 846.218975, 457.768092, &
           350.374944, 272.562811, 223.870091, 166.865074, 38.7318626, &
           14.2240198, 6.66039491, 3.70161383, 1.59303347, 0.90501953, &
           0.66081004, 0.56256708, 0.462352483, 0.456081245, &
           0.512564783, 0.562567352, 0.562570538, 0.72647646, &
           0.812385259, 0.98512227, 1.23417611, 1.51151667, 1.83032758, &
           2.30015202, 2.85160797, 3.57134514, 4.45606314, 5.56788263, &
           6.99458566, 8.78428656, 11.0041227, 13.8139187, 17.3752908, &
           21.8460124, 27.4902768, 34.5361013 /)


! MaxGap: 95% 
! sigipa_CDMS = (/22618.8051, 10927.2903, 2648.20503, 1088.7403, 589.659919, &
!           445.14942, 341.943202, 276.27667, 203.403212,  48.8448156, &
!           18.2466503, 8.63031045, 4.80840143, 2.0813092, 1.18757246, &
!           0.859941278, 0.710023916, 0.593828958, 0.593828649, &
!           0.643829319, 0.669681578, 0.807967841, 0.923571711, & 
!           1.07119499, 1.28328038, 1.58494204, 1.96100958, 2.42507033, &
!           2.97587081, 3.69830569, 4.64680087, 5.81475626, 7.25306741, &
!           9.10682081, 11.3901857, 14.3312164, 17.9694274, 22.6091451, &
!           28.4416259, 35.7354209, 44.9467802 /)



 if(mchi .lt. 6.72950096) then
   lsigip_bound = -2.7598d0 !-2.6455d0 95%
 else if(mchi .gt. 1.d4) then !extrapol.
   write(*,*) 'warning: mchi out of tabulation range'
   lsigip_bound = -5.4617d0 !-5.3473d0 95%
 else !interpl.
   lmchia_CDMS = log10(mchia_CDMS)
   lsigipa_CDMS = log10(sigipa_CDMS*1.e-7)
   call spline(lmchia_CDMS, lsigipa_CDMS, DIM_CDMS, 1.0e30, 1.0e30, &
               y2a_CDMS)
   call splint(lmchia_CDMS, lsigipa_CDMS, y2a_CDMS, DIM_CDMS, &
               lmchi, lsigip_bound)
endif 

!data from Fig. 4 in astro-ph/0509259 for ZEPLIN-I

else if(experiment .eq. 2) then

 mchia_ZEPL = (/9.85, 10.3, 11.1, 11.8, 12.8,& 
          14.1, 15.3, 17.2, 19.4, 23.6,& 
          26.6, 31.9, 39.7, 50.3, 61.4,& 
          79.6, 103., 123., 153., 173.,& 
          215., 273., 391., 582., 1.17e+03,& 
          2.21e+03, 2.8e+03, 5.29e+03, 1e+04/)

 sigipa_ZEPL = (/1050., 805., 540., 368., 224.,&
           145., 106., 68.6, 43.7, 28.9,& 
           22.1, 17.2, 13.4, 11.6, 10.8,& 
           10.8, 11.9, 13.4, 15.1, 16.7, &
           19.7, 24.0, 33.5, 50., 97.3,& 
           177., 231., 428., 805./)

 if(mchi .lt. 9.85) then
   lsigip_bound =  -3.98d0
 else if(mchi .gt. 1.d4) then !extrapol.
   write(*,*) 'warning: mchi out of tabulation range'
   lsigip_bound = -4.094d0
 else !interpl.
   lmchia_ZEPL = log10(mchia_ZEPL)
   lsigipa_ZEPL = log10(sigipa_ZEPL*1.e-7)
 call spline(lmchia_ZEPL, lsigipa_ZEPL, DIM_ZEPL, 1.0e30, 1.0e30, &
             y2a_ZEPL)
 call splint(lmchia_ZEPL, lsigipa_ZEPL, y2a_ZEPL, DIM_ZEPL, &
             lmchi, lsigip_bound)
 endif

!data from Fig. 4 in astro-ph/0509259 for EDELWEISS-I

else if(experiment .eq. 3) then


 mchia_EDELW =(/15.71, 18.51, 22.17, 26.13, 31.3, &
          36.88, 44.19, 52.07, 64.46, 75.96, &
          98.78, 120.3, 146.5, 178.4, 217.2,& 
          264.5, 322.1, 392.3, 477.7, 581.7,&
          708.4, 862.6, 1050., 1279., 1558.,& 
          1897., 3e+03, 5e+03, 1e+04/)

 sigipa_EDELW = (/980.2, 449.6, 206.3, 111., 51.95,& 
           32.16, 22., 17.66, 14.46, 13.62, &
           13.62, 14.46, 15.35, 17.31, 19.51,& 
           22.44, 26.34, 31.52, 36.99, 44.28,& 
           53., 63.44, 77.48, 92.74, 113.3,& 
           135.6, 200., 370., 750./)

 if(mchi .lt. 15.71) then
   lsigip_bound = -4.01d0 
 else if(mchi .gt. 1.d4) then !extrapol.
   write(*,*) 'warning: mchi out of tabulation range'
   lsigip_bound = -4.125d0
 else !interpl.
   lmchia_EDELW = log10(mchia_EDELW)
   lsigipa_EDELW = log10(sigipa_EDELW*1.e-7)
   call spline(lmchia_EDELW, lsigipa_EDELW, DIM_EDELW, 1.0e30, 1.0e30, & 
               y2a_EDELW)
   call splint(lmchia_EDELW, lsigipa_EDELW, y2a_EDELW, DIM_EDELW, &
               lmchi, lsigip_bound)
 endif

!data from XENON

else if(experiment .eq. 4) then

 mchia_XENON = (/10., 11.220185, 12.589254, 14.125375, 15.848932, &
           17.782794, 19.952623, 22.387211, 25.118864, 28.183829, &
           31.622777, 35.481339, 39.810717, 44.668359, 50.118723, &
           56.234133, 63.095734, 70.794578, 79.432823, 89.125094, &
           100., 112.20185, 125.89254, 141.25375, 158.48932, 177.82794, &
           199.52623, 223.87211, 251.18864, 281.83829, 316.22777, &
           354.81339, 398.10717, 446.68359, 501.18723, 562.34133, &
           630.95734, 707.94578, 794.32823, 891.25094, 1000., &
           1122.0185, 1258.9254, 1412.5375, 1584.8932, 1778.2794, &
           1995.2623, 2238.7211, 2511.8864, 2818.3829, 3162.2777, &
           3548.1339, 3981.0717, 4466.8359, 5011.8723, 5623.4133, &   
           6309.5734, 7079.4578, 7943.2823, 8912.5094, 10000. /)


! 90% CL 
 sigipa_XENON = (/7.5052544, 3.0852955, 1.5251386, 0.87321348, &
           0.56371985, 0.40269452, 0.31393190, 0.26395642, 0.23666854, &
           0.22373859, 0.22060354, 0.22465645, 0.23438603, 0.24893239, &
           0.26785401, 0.29098760, 0.31837310, 0.35019735, 0.38676584, &
           0.42848607, 0.47586034, 0.52948302, 0.59004236, 0.65832598, &  
           0.73522896, 0.82176408, 0.91907384, 1.0284464, 1.1513299, &
           1.2893530, 1.4443459, 1.6183641, 1.8137160, 2.0329932, &
           2.2791050, 2.5553168, 2.8652936, 3.2131483, 3.6034967, &
           4.0415184, 4.5330255, 5.0845400, 5.7033800, 6.3977573, &
           7.1768858, 8.0511039, 9.0320120, 10.132626, 11.367551, &
           12.753173, 14.307878, 16.052297, 18.009577, 20.205690, &
           22.669776, 25.434534, 28.536649, 32.017285, 35.922627, &
           40.304497, 45.221041 /)    


 if(mchi .lt. 10.) then
   lsigip_bound = -6.1246d0
 else if(mchi .gt. 1.d4) then !extrapol.
   write(*,*) 'warning: mchi out of tabulation range'
   lsigip_bound = -5.3447d0
 else !interpl.
   lmchia_XENON = log10(mchia_XENON)
   lsigipa_XENON = log10(sigipa_XENON*1.e-7)
   call spline(lmchia_XENON, lsigipa_XENON, DIM_XENON, 1.0e30, 1.0e30, &
               y2a_XENON)
   call splint(lmchia_XENON, lsigipa_XENON, y2a_XENON, DIM_XENON, &
               lmchi, lsigip_bound)
endif 


else
  print*,'Experiment option in DD exp. bound not recognized'
  sigip_bound = 0.d0
  return
endif 


sigip_bound = 10.**lsigip_bound

END SUBROUTINE sigip_exp








