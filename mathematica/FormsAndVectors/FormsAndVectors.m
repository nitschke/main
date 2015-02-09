(* ::Package:: *)

BeginPackage["FormsAndVectors`"]

OneForm2D::usage = 
 	"OneForm[f1,f2] construct a 1-Form f1dx^1+f2dx^2 as  {f1,f2}"
TwoForm2D::usage =
 	"TwoForm2D[w12] construct a 2-form w12dx^1\[And]dx^2 as {{w12}}"

IsNullForm::usage = "IsNullForm[alpha] is True if alpha is a 0-form"
IsOneForm::usage = "IsOneForm[alpha] is True if alpha is a 1-form"
IsTwoForm::usage = "IsTwoForm[alpha] is True if alpha is a 2-form"

Wedge::usage = "Wedge[alpha,beta] = alpha^beta"

Begin["Private`"]

OneForm2D[f1_, f2_] := {f1, f2}
TwoForm2D[w12_] := {{w12}}

IsNullForm[alpha_] := ArrayDepth[alpha] == 0
IsOneForm[alpha_] := ArrayDepth[alpha] == 1
IsTwoForm[alpha_] := ArrayDepth[alpha] == 2

Wedge[alpha_, beta_] := 
 	Which[
  		IsOneForm[alpha],
  			Which[
    				IsOneForm[beta],
    					{{alpha[[1]] beta[[2]] - alpha[[2]] beta[[1]]}},
    				True,
						Print["Wedge::Cant determine"];
						$Failed
    					(*"Wedge::Cant determine"*)
    			],
   		True,
  			"Wedge::Cant determine"
  	]

End[]

EndPackage[]





"FormsAndVectors`"


"OneForm[f1,f2] construct a 1-Form f1dx^1+f2dx^2 as  {f1,f2}"


"TwoForm2D[w12] construct a 2-form w12dx^1\[And]dx^2 as {{w12}}"


"OneForm[alpha] is True if alpha is a 1-form"


"Wedge[alpha,beta] = alpha^beta"


"Private`"


"Private`"
