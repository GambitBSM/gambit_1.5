function phig_back(e, InN)

 use parameters

 implicit none

 Type(Nuisance_params), INTENT(IN) :: InN

 real*8 :: phig_back
 real*8 :: e, delta


 if(e <= InN%e0) then
  delta = InN%delta1
 else
  delta = InN%delta2
 endif
   

 phig_back = InN%a0 * (e/InN%e0)**delta


end function phig_back
