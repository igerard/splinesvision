(* ::Package:: *)

BeginPackage["SplineCurve`","BSpline`"]

Sp::usage=
  "Sp[k,n,poles] defini une courbe spline par l'ordre, 
   le vecteur noeud et les coordonnees"
   
SpN::usage=
  "SpN[k,n,poles] normalise le vecteur noeud"
   
Draw::usage=
  "Draw[Sp[k,n,poles], opts] dessine la spline Sp selon la
   dimension des coordonnees, passe les opts a Plot ou ParametricPlot"

openKnot::usage=
  "openKnot[k,n]"

knotsUsed::usage=
  "knotsUsed[k,n]"
  
normalizeKnot::usage=
  "normalizeKnot[k,n]"

Begin["`Private`"] 

Sp::usage=
  "Sp[k,n,coords]";

Sp[k_Integer, n_List, pol_List][u_] := 
  Sum[pol[[i]] Bs[i,k,n][u], {i,1,Length[n]-k}]

openKnot[k_Integer,n_?VectorQ] :=
  Join[{n[[1]]-1},
        Table[n[[1]],{k-2}],
       n,
       Table[n[[-1]],{k-2}],
       {n[[-1]]+1}]

suppressDouble[n_List] :=
  Module[{i},
         For[i=2,i<Length[n],i++,
             If[n[[i]] == n[[i-1]],
                Delete[n,i]]]]

knotsUsed[k_Integer,n_?VectorQ] := 
  n[[Range[k,Length[n]-k+1]]]

normalizeKnot[k_Integer,n_?VectorQ] :=
  Module[{a = n[[k]],b = n[[-k]],f},
          f = (#-a)/(b-a)&;
       Map[f,n]]
       
SpN[Sp[k_,n_,poles_]] :=
  Module[{nn},
         nn = normalizeKnot[k,n];
         Sp[k,nn,poles]]
  
deriveKnot[k_,knots_] :=
  Block[{nn = Drop[Drop[knots,1],-1]},
        If[Equal[nn[[Range[1,k-1] ]] ],
           nn = MapAt[#-1&,nn,{1}]];
        If[Equal[nn[[-Range[1,k-1] ]] ],
           nn = MapAt[#+1&,nn,{-1}]];
        nn]

Sp /:
  Derivative[1][Sp[k_Integer, n_List, pol_List]] :=
    Block[{i,nn={},npol={}},
          nn = deriveKnot[k,n];
          For[i=2,i<=Length[pol],i++,
              AppendTo[npol,
                       (pol[[i]]-pol[[i-1]]) ~dz~
                       (n[[i+k-1]]-n[[i]])]  ];
          Sp[k-1,nn,(k-1)npol]]

Sp /:
  Derivative[l_][Sp[k_Integer, n_List, pol_List]] :=
    Derivative[l-1][Sp[k,n,pol]'] /; l>1
    
Sp[k_Integer,
   h_Integer,
   n_List,
   m_List,
   coords_][u_Symbol,v_?NumberQ] :=
  Module[{coordU={},i,j},
         For[i=1,i<=Length[n]-k,i++,
             AppendTo[coordU, 
                    Sum[coords[[i,j]] Bs[j,h,m][v],
                        {j,1,Length[m]-h}]]];
         Sp[k,n,coordU]]
         
Sp[k_Integer,
   h_Integer,
   n_List,
   m_List,
   coords_][u_?NumberQ,v_Symbol] :=
  Module[{coordV={},i,j},
         For[j=1,j<=Length[m]-h,j++,
             AppendTo[coordV, 
                    Sum[coords[[i,j]] Bs[i,k,n][u],
                        {i,1,Length[n]-k}]]];
         Sp[h,m,coordV]]
         
Sp[k_Integer,
   h_Integer,
   n_List,
   m_List,
   coords_][u_,v_] :=
  Sum[Sum[coords[[i,j]] Bs[i,k,n][u] Bs[j,h,m][v],
          {j,1,Length[m]-h}],
          {i,1,Length[n]-k}]
          
Derivative[1,0][Sp[k_,h_,n_,m_,coords_]] :=
    Sum[Sum[coords[[i,j]] Bs[i,k,n]'[#1] Bs[j,h,m][#2],
            {j,1,Length[m]-h}],
        {i,1,Length[n]-k}]&
          
Derivative[0,1][Sp[k_,h_,n_,m_,coords_]] :=
    Sum[Sum[coords[[i,j]] Bs[i,k,n][#1] Bs[j,h,m]'[#2],
            {j,1,Length[m]-h}],
       {i,1,Length[n]-k}]&

Draw[Sp[k_Integer,n_?VectorQ,pols_?VectorQ],
     opts___] :=
  Plot[Evaluate[Sp[k,n,pols][u]],{u,n[[k]],n[[-k]]},
       opts,
       Compiled -> True]
      
Draw[Sp[k_Integer,n_?VectorQ,pols_?MatrixQ],
     opts___] :=
  Module[{u,fcts},
         fcts = Map[Sp[k,n,#][u]&, Transpose[pols]];
         ParametricPlot[Evaluate[fcts],{u,n[[k]],n[[-k]]},
                        opts,
                        Compiled -> True]] /; 
    Dimensions[pols][[2]] == 2

Draw[Sp[k_Integer,n_?VectorQ,pols_?MatrixQ],opts___] :=
  Module[{u,fcts},
         fcts = Map[Sp[k,n,#][u]&, Transpose[pols]];
         ParametricPlot3D[Evaluate[fcts],{u,n[[k]],n[[-k]]},
                          opts,
                          Compiled -> True]] /; 
    Dimensions[pols][[2]] == 3

End[]
EndPackage[]
