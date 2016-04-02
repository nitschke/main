(* ::Package:: *)

BeginPackage[ "SurfaceTensorAnalysis`"]
(*required: 
	g :: metric tensor, eg g = {{1, 0}, {0, Sin[u]^2}}
	vars :: varibales, eg vars = {u, v}

  form:
	{array,{typ,typ,...,typ}}, eg t={{{1,2},{3,4}},{flat,sharp}} result in a flat-sharp 2-tensor (first index covariant, second contravariant)
*)

Unprotect[g,gInv,vars,GammaT];

Init::usage "Init[metric, variables]";

gInv::usage "inverse metric tensor (Full Sharp)";

SubTMSwap::usage "SubTMSwap[tensorMatrix, dimension, index] gives tensorMatrix[[All,...,index,...,All]] (->Trans_(1,k)) with index in the dimensionth place.
					index=All is equivalent to a (1,dimension)-Transpose, 
					e.g., if k=dimension, then index list (i_1,...,i_k,...i_n) becomes (i_k,i_2,...,i_(k-1),i_1,i_(k+1),...,i_n).";
SubTMPush::usage "SubTMPush[tensorMatrix, dimension, index] gives tensorMatrix[[All,...,index,...,All]] (->Trans_(k->1)) with index in the dimensionth place.
					index=All is equivalent to the Transpose, which set the dimension on the first place and push the rest back, 
					e.g., if k=dimension, then index list (i_1,...,i_k,...i_n) becomes (i_k,i_1,...,i_(k-1),i_(k+1),...,i_n).";
SubTMPushInv::usage "SubTMPushInv[tensorMatrix, dimension, index] (->Trans_(1->k)) inverse (in tensorMatrix) of SubTMPush.";

TransWithFirstTM::usage = "TransWithFirstTM[tensorMatrix, k] -> Trans_(1,k)";
TransToFirstTM::usage = "TransToFirstTM[tensorMatrix, k] -> Trans_(k->1)";
TransFromFirstTM::usage = "TransFromFirstTM[tensorMatrix, k] -> Trans_(1->k)";
TransFromToTM::usage = "TransFromToTM[tensorMatrix, k, i] -> Trans_(k->i) = Trans_(1->i) after Trans_(k->1)";

TransformT:usage "TransformT[tensor, {typ,typ,...,typ}]";

OuterT::usage "OuterT[tensor_1, tensor_2] gives the outer (free) product.";

ContractT::usage "ContractT[t, dim1, dim2]";

CoDT::usage "CoDT[tensor] gives the covariant derivative (Levi-Cevita-Connection) of the tensor.";
flat=0;
sharp=1;

TensorAnalysis::TensorOrderMismatch = "Tensor order mismatch; `1` tensor has order `3` and `2` tensor has order `4`.";

Begin[ "Private`"]

Init[metric_,variables_]:=Module[{},Unprotect[g,gInv,vars,GammaT];
						g = metric;
						$Assumptions = $Assumptions && Det[g]>0;
						(*Todo: check that vars is not quantifier and protect them*)
						vars = variables;
						gInv = Inverse[metric]//Simplify;
						GammaT = {Table[Sum[gInv[[k,l]]*(D[g[[j,l]],vars[[i]]]+D[g[[i,l]],vars[[j]]]-D[g[[i,j]],vars[[l]]]),{l,2}],{i,2},{j,2},{k,2}]/2//Simplify,{0,0,1}};
						Protect[g,gInv,vars,GammaT];];

SubTMSwap[tM_, k_, i_]:=Transpose[tM,Range[Length[Dimensions[tM]]] /. {1 -> k, k -> 1}][[i]];

SubTMPushInv[tM_, k_, i_]:=Transpose[tM,Prepend[Select[Range[Length[Dimensions[tM]]],#!=k&],k]][[i]];
SubTMPush[tM_, k_, i_]:=Transpose[tM,Insert[Drop[Range[Length[Dimensions[tM]]],1],1,k]][[i]];

TransWithFirstTM[tM_, k_]:= SubTMSwap[tM, k, All];
TransToFirstTM[tM_, k_]:= SubTMPush[tM, k, All];
TransFromFirstTM[tM_, k_]:= SubTMPushInv[tM, k, All];
TransFromToTM[tM_, k_, i_]:= TransFromFirstTM[TransToFirstTM[tM,k],i];

TransformT[t_,types_]:=Module[{len=Length[types],qM=t[[1]]},
								If[len!=Length[t[[2]]],Message[TensorAnalysis::TensorOrderMismatch,"input", "desired", Length[t[[2]]],len]];
								For[k=1,k<len+1,k++,
									If[types[[k]] > t[[2,k]],
										(*qM=SubTM[Table[Sum[gInv[[l,i]]*SubTM[qM,k,i],{i,2}],{l,2}],k,All]//Simplify,*)
										qM=SubTMSwap[gInv.SubTMSwap[qM,k,All],k,All]//Simplify,
										If[types[[k]] < t[[2,k]],
											(*qM=SubTM[Table[Sum[g[[l,i]]*SubTM[qM,k,i],{i,2}],{l,2}],k,All]//Simplify,*)
											qM=SubTMSwap[g.SubTMSwap[qM,k,All],k,All]//Simplify,
											qM=qM
										]
									]
								];
								{qM,types}
							];

ContractT[t_,d1_,d2_]:=Module[{tt, types=t[[2]]},
								If[types[[d1]]==types[[d2]],
									types[[d2]]=Mod[types[[d2]]+1,2];
									tt=TransformT[t,types],
									tt=t
									];
								{Sum[SubTMPush[SubTMPush[tt[[1]],Max[d1,d2],i],Min[d1,d2],i],{i,2}]//Simplify,Delete[types,{{d1},{d2}}]}
							]

OuterT[t1_,t2_]:={Outer[Times,t1[[1]],t2[[1]]]//Simplify,Flatten[Append[t1[[2]],t2[[2]]]]}

CoDT[t_]:=Module[{len=Length[t[[2]]],tM=t[[1]],types=t[[2]],CDtM,tGT=OuterT[t,GammaT]},
				CDtM = TransFromFirstTM[Table[D[tM,x],{x,vars}],len+1] //Simplify;(* t_(i_ 2,...,i_(n+1)|i_ 1) -> t_(i_ 1,...,i_n|i_(n+1))*)
				For[d=1,d<len+1,d++,
					CDtM -= (-1)^types[[d]] * TransFromToTM[ContractT[tGT, d, len+3-types[[d]]][[1]], len+types[[d]], d]//Simplify;
				];
				{CDtM,Append[types,0]}
			];

(* Defaults *)
Init[IdentityMatrix[2],{x,y}]; 

End[]

Protect[g,gInv,vars,GammaT];

EndPackage[]








