(* ::Package:: *)

(* ::Section:: *)
(* Dummy Mathematica package for testing *)


(* \author Tomas Gonzalo *)
(*         t.e.gonzalo@fys.uio.no *)
(* \date 2016 Sept *)



CalculateSquare[var_] := var^2;


CalculateSum[var1_,var2_]:=var1+var2


Var=42;


Var2=2.5;


PrintVar[]:=Var;


PrintVarorVar2[check_]:=If[check,Var,Var2];
