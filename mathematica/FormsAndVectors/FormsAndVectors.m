(* ::Package:: *)

BeginPackage["FormsAndVectors`"]

MetricFromPara::usage = "Metric from parametricmap {x[u,v],y[u,v],z[u,v]}";

GlobalVecFromPara::usage = "Global (x,y,z)-Vector of a given local (u,v)-Tangentialvector from parametrisation";
LocalVecFromPara::usage = "Local (u,v)-Vector (tangential part) of a given (x,y,z)-Vector from parametrisation";

Wedge11::usage = "Wedge product of two 1-forms";

ExD0::usage = "Exterior derivativ of a 0-form";
ExD1::usage = "Exterior derivativ of a 1-form";

Sharp1::usage = "Rise the indices of a 1-form";
Flat1::usage = "Lower the indices of a Vector";

Hodge0::usage = "Hodge-Star of a 0-form";
Hodge1::usage = "Hodge-Star of a 1-form";
Hodge2::usage = "Hodge-Star of a 2-form";

Dot1::usage = "Inner product of two 1-forms";

Rot0::usage = "Rotation (*d) of a 0-form";
Rot1::usage = "Rotation (*d) of a 1-form";

ExCoD1::usage = "Exterior Coderivative (\[Delta]=-*d*=-Div) of a 1-form";
ExCoD2::usage = "Exterior Coderivative (\[Delta]=-*d*) of a 2-form";

Grad0::usage = "Gradient (#\[SmallCircle]d) of a 0-form";

LBeltrami0::usage = "Laplace-Beltrami (-\[Delta]d=*d*d) of a 0-form";
LBeltrami1::usage = "Laplace-Beltrami (-\[Delta]d=*d*d) of a 1-form";
LBeltrami1Vec::usage = "Laplace-Beltrami (-\[Delta]d=*d*d) of a vector (with pre-flat and post-sharp)";

LCoBeltrami1::usage = "Laplace-Co-Beltrami (-d\[Delta]=d*d*) of a 1-form";
LCoBeltrami1Vec::usage = "Laplace-Co-Beltrami (-d\[Delta]=d*d*) of a vector (with pre-flat and post-sharp)";
LCoBeltrami2::usage = "Laplace-Co-Beltrami (-d\[Delta]=d*d*) of a 2-form";

LDeRham0::usage = "Laplace-De-Rham (\[Delta]d) of a 0-form";
LDeRham1::usage = "Laplace-De-Rham (\[Delta]d+d\[Delta]) of a 1-form";
LDeRham1Vec::usage = "Laplace-De-Rham (\[Delta]d+d\[Delta]) of a vector (with pre-flat and post-sharp)";
LDeRham2::usage = "aplace-De-Rham (d\[Delta]) of a 2-form";

LieD0::usage = "Lie-derivative of a 0-form";
LieD1::usage = "Lie-derivative of a 1-form";
LieD2::usage = "Lie-derivative of a 2-form";
LieDT02::usage = "Lie-derivative of a (0,2)-Tensor";

Inner1::usage = "Contraction of a 1-form";
Inner2::usage = "Contraction of a 2-form";

DotForm1::usage = "Dot product of 1-forms";

L2Prod0::usage = "L2 Product of 0-forms"
L2Prod1::usage = "L2 Product of 1-forms"
L2Prod2::usage = "L2 Product of 2-forms"

DoubleDotFormForm11::usage = "Double dot product (:) of 1-forms of 1-forms";

Begin["Private`"]

MetricFromPara[paraMap_,x_,y_] := Module[{Dx=D[paraMap,x], Dy=D[paraMap,y]}, Outer[Dot,{Dx,Dy},{Dx,Dy},1]]

GlobalVecFromPara[locVec_,paraMap_,x_,y_] := Module[{Dx=D[paraMap,x], Dy=D[paraMap,y]}, locVec[[1]]Dx + locVec[[2]]Dy]
LocalVecFromPara[globVec_,paraMap_,x_,y_] := Module[{Dx=D[paraMap,x], Dy=D[paraMap,y]}, 
										          Sharp1[{globVec.Dx, globVec.Dy},Outer[Dot,{Dx,Dy},{Dx,Dy},1]]]

Wedge11[alpha_,beta_] := {{alpha[[1]]*beta[[2]] - alpha[[2]]*beta[[1]]}}

ExD0[f_,x_,y_] :=  {D[f,x],D[f,y]}
ExD1[alpha_,x_,y_] := {{D[alpha[[2]],x] - D[alpha[[1]],y]}}

Sharp1[alpha_,g_] := Inverse[g].alpha
Flat1[v_,g_] := g.v

muVal[g_]:= Sqrt[Det[g]] (*Sqrt[Abs[Det[g]]]*)
Hodge0[f_,g_] := {{muVal[g]f}}
Hodge1[alpha_,g_] := muVal[g]{{0,-1},{1,0}}.Sharp1[alpha,g]
Hodge2[omega_,g_] := omega[[1,1]]/muVal[g]

Dot1[alpha_,beta_,g_] := alpha.Sharp1[beta,g]

Rot0[f_,x_,y_,g_] := Hodge1[ExD0[f,x,y],g]
Rot1[alpha_,x_,y_,g_] := Hodge2[ExD1[alpha,x,y],g]

ExCoD1[alpha_,x_,y_,g_] := -Hodge2[ExD1[Hodge1[alpha,g],x,y],g]
ExCoD2[omega_,x_,y_,g_] := -Hodge1[ExD0[Hodge2[omega,g],x,y],g]

Grad0[f_,x_,y_,g_] := Sharp1[ExD0[f,x,y],g]

LBeltrami0[f_,x_,y_,g_] := -ExCoD1[ExD0[f,x,y],x,y,g]
LBeltrami1[alpha_,x_,y_,g_] := -ExCoD2[ExD1[alpha,x,y],x,y,g]
LBeltrami1Vec[vec_,x_,y_,g_] := Sharp1[LBeltrami1[Flat1[vec,g],x,y,g],g]

LCoBeltrami1[alpha_,x_,y_,g_] := -ExD0[ExCoD1[alpha,x,y,g],x,y]
LCoBeltrami1Vec[vec_,x_,y_,g_] := Sharp1[LCoBeltrami1[Flat1[vec,g],x,y,g],g]
LCoBeltrami2[omega_,x_,y_,g_] := -ExD1[ExCoD2[omega,x,y,g],x,y]

LDeRham0[f_,x_,y_,g_] := -LBeltrami0[f,x,y,g]
LDeRham1[alpha_,x_,y_,g_] := -LBeltrami1[alpha,x,y,g] - LCoBeltrami1[alpha,x,y,g]
LDeRham1Vec[vec_,x_,y_,g_] := Sharp1[LDeRham1[Flat1[vec,g],x,y,g],g]
LDeRham2[omega_,x_,y_,g_] := -LCoBeltrami2[omega,x,y,g]

LieD0[vec_,f_,x_,y_] := vec.ExD0[f,x,y]
(*LieD1[vec_,alpha_,x_,y_] := {vec[[1]]D[alpha[[1]],x] + vec[[2]]D[alpha[[1]],y] + alpha[[1]]D[vec[[1]],x] + alpha[[2]]D[vec[[2]],x],
							  vec[[1]]D[alpha[[2]],x] + vec[[2]]D[alpha[[2]],y] + alpha[[1]]D[vec[[1]],y] + alpha[[2]]D[vec[[2]],y]}*)
LieD1[vec_,alpha_,x_,y_] := Module[{xv={x,y}},Table[Sum[vec[[j]]D[alpha[[i]],xv[[j]]] + alpha[[j]]D[vec[[j]],xv[[i]]],{j,1,2}],{i,1,2}]]
LieD2[vec_,omega_,x_,y_] := {{D[omega[[1,1]]vec[1],x] + D[omega[[1,1]]vec[2],y]}}
LieDT02[vec_,sigma_,x_,y_] := Module[{xv={x,y}},Table[Sum[vec[[k]]D[sigma[[i,j]],xv[[k]]]+sigma[[k,j]]D[vec[[k]],xv[[i]]]+sigma[[i,k]]D[vec[[k]],xv[[j]]],{k,1,2}],{i,1,2},{j,1,2}]]

Inner1[vec_,alpha_] := vec.alpha
Inner2[vec_,omega_] := omega[[1,1]]{-vec[[2]],vec[[1]]}

DotForm1[alpha_,beta_,g_] := alpha.Sharp1[beta,g]

L2Prod0[f1_,f2_,ivalx_,ivaly_,g_] := Integrate[f1*Hodge0[f2,g][[1,1]],ivalx,ivaly]
L2Prod1[alpha_,beta_,ivalx_,ivaly_,g_] := Integrate[Wedge11[alpha,Hodge1[beta,g]][[1,1]],ivalx,ivaly]
L2Prod2[omega1_,omega2_,ivalx_,ivaly_,g_] := Integrate[Hodge2[omega1,g]*omega2-[[1,1]],ivalx,ivaly]

DoubleDotFormForm11[sigma_,tau_,g_] := Module[{gInv=Inverse[g]}, Sum[sigma[[i,k]]gInv[[i,l]]gInv[[k,s]]tau[[l,s]],{i,1,2},{k,1,2},{l,1,2},{s,1,2}]]

End[]

EndPackage[]





















