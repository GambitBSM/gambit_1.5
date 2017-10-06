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


VarEqualVar2[]:=Var==Var2


StringTest[var_]:="This is not a test, "<>ToString[var];


VoidTest[]:=Var=Var+1;


ExtractElement[list_,i_]:=list[[i]]


SquareList[list_]:=Table[list[[i]]^2,{i,Length[list]}]


\[Delta]var = 75;


\[CapitalGamma]func[var_] := var + 2;
