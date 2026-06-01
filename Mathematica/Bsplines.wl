(* ::Package:: *)

(* ::Package:: *)
(* Copyright 1992 Visual Corporation. *)
(*   1992 - 2026   *)

(* Bspline Definitions *)

(* by G\[EAcute]rard Iglesias *)

BeginPackage["BSpline`"]

Bs::usage=
  "Bs[i,k,n][u] defini la i eme base spline d'ordre k..."

dz::usage=
  "dz[a,b] effectue la division en renvoyant 0 si b = 0"

Begin["`Private`"]

ZeroQ[a_?NumberQ] := a == 0
NonZeroQ[a_?NumberQ] := a != 0

a_?ZeroQ ~dz~ b_?ZeroQ := 0
a_ ~dz~ b_?ZeroQ := 0
a_ ~dz~ b_?NonZeroQ := a/b
Derivative[1,0][dz] := (1 ~dz~ #2)&
Derivative[0,1][dz] := (-#1 ~dz~ #2^2)&
SetAttributes[dz, Listable]

Bs[i_Integer, 1, n_List][u_] := Bs[i, n][u]
 
Bs[i_Integer, n_List][u_?NumberQ] := 
  If[u >= n[[i]] && u < n[[i + 1]], 1, 0]
 
Bs[i_Integer, k_Integer, n_List][u_] := 
 ((u-n[[i]]) ~dz~ (n[[i+k-1]]-n[[i]]))Bs[i,k-1,n][u] + 
 ((n[[i+k]]-u) ~dz~ (n[[i+k]]-n[[i+1]]))Bs[i+1,k-1,n][u]
 
Bs[i_Integer, k_Integer, n_List][u_Symbol] := 
Simplify[Expand[
   ((u-n[[i]]) ~dz~ (n[[i+k-1]]-n[[i]])) *
                               Bs[i,k-1,n][u] + 
 ((n[[i+k]]-u) ~dz~ (n[[i+k]]-n[[i+1]])) *
                            Bs[i+1,k-1,n][u]]]
 
Bs[i_Integer, n_List]' := 0&
      
Derivative[1][Bs[i_Integer, k_Integer, n_List]] := 
 ((k - 1) *
  ((Bs[i,k-1,n][#] ~dz~ (n[[i+k-1]]-n[[i]])) - 
   (Bs[i+1,k-1,n][#] ~dz~ (n[[i+k]]-n[[i+1]])))) &

Derivative[l_][Bs[i_Integer, k_Integer, n_List]] := 
  Derivative[l-1][Bs[i, k, n]']

End[]
EndPackage[]
