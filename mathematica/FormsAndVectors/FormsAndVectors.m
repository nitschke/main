(* ::Package:: *)

BeginPackage["FormsAndVectors`"]

Wedge11::usage = "Wedge product of two 1-forms";

ExD0::usage = "Exterior derivativ of a 0-form";
ExD1::usage = "Exterior derivativ of a 1-form";

Sharp1::usage = "Rise the indices of a 1-form";
Flat1::usage = "Lower the indices of a Vector";

Hodge0::usage = "Hodge-Star of a 0-form";
Hodge1::usage = "Hodge-Star of a 1-form";
Hodge2::usage = "Hodge-Star of a 2-form";

Begin["Private`"]

Wedge11[alpha_,beta_] := {{alpha[[1]]*beta[[2]] - alpha[[2]]*beta[[1]]}}

ExD0[f_,x_,y_] :=  {D[f,x],D[f,y]}
ExD1[alpha_,x_,y_] := {{D[alpha[[2]],x] - D[alpha[[1]],y]}}

Sharp1[alpha_,g_] := Inverse[g].alpha
Flat1[v_,g_] := g.v

muVal[g_]:= Sqrt[Det[g]] (*Sqrt[Abs[Det[g]]]*)
Hodge0[f_,g_] := {{muVal[g]f}}
Hodge1[alpha_,g_] := muVal[g]{{0,-1},{1,0}}.(Transpose[Inverse[g]].alpha)
Hodge2[omega_,g_] := omega[[1,1]]/muVal[g]

End[]

EndPackage[]



