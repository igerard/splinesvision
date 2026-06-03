(* ::Package:: *)

BeginPackage["SplineInterpolation`","SplineCurve`"]

condsSuppV0::usage=
  "condsSuppV0 function de calcul des equation rendant les vitesses
   nulle au extremites d'un spline"
   
condsSuppA0::usage=
  "condsSuppA0"
  
condsSuppQuad::usage=
  "condsSuppQuad"
  
condsSupThirdEq::usage=
  "condsSupThirdEq"
  
cordLengthKnots::usage=
  "cordLengthKnots"
  
centripeteKnots::usage=
  "centripeteKnots"
  
SplineFit::usage=
  "SplineFit[pts,k,parFunc,condSuppFunc]"
  
Begin["`Private`"]

Dist[a_?NumberQ, b_?NumberQ] := Abs[b - a]

Dist[a_?VectorQ, b_?VectorQ] := 
  Module[{l = (b-a)}, Sqrt[l . l]]
  
unknown[k_,ipts_] := 
  Module[{i,
          l = Length[ipts]+k-2},
    Table[Q[i],
          {i,1,l}]]
          
extendParameters[pars_?VectorQ,
                 k_Integer] :=
   Module[{before,after,i},
          before = Join[{pars[[1]]-1},
                        Table[pars[[1]],{i,k-2}]];
          after = Join[Table[pars[[-1]],{i,k-2}],
                       {pars[[-1]]+1}];
          Join[before,pars,after]]

interpolsPtsEqs[k_, knots_, ipts_, unknowns_] := 
  Module[{i},
         Table[Sp[k,knots,unknowns] [knots[[i+k-1]]] == 
               ipts[[i]],
               {i,1,Length[ipts]}]]
               
condsSuppV0[k_, knots_, unknowns_] := 
  {} /; k <= 2
   
condsSuppA0[k_, knots_, unknowns_] := 
  {} /; k <= 2
   
condsSuppV0[3, knots_, unknowns_] := 
  {Derivative[1][Sp[3,knots,unknowns]][knots[[3]]] == 0}
   
condsSuppA0[3, knots_, unknowns_] := 
  {Derivative[2][Sp[3,knots,unknowns]][knots[[3]]] == 0} 
     
condsSuppV0[k_, knots_, unknowns_] := 
  {Derivative[1][Sp[k,knots,unknowns]][knots[[k]]] == 0,
   Derivative[1][Sp[k,knots,unknowns]][knots[[-k]]] == 0}[[Range[1,k-2]]]
   
condsSuppA0[k_, knots_, unknowns_] := 
  {Derivative[2][Sp[k,knots,unknowns]][knots[[k]]] == 0,
   Derivative[2][Sp[k,knots,unknowns]][knots[[-k]]] == 0}[[Range[1,k-2]]]
   
condsSuppQuad[k_, knots_, unknowns_] := 
  {Derivative[2][Sp[k,knots,unknowns]][knots[[k]]] == 
    Derivative[2][Sp[k,knots,unknowns]][knots[[k+1]]],
   Derivative[2][Sp[k,knots,unknowns]][knots[[-k]]] == 
     Derivative[2][Sp[k,knots,unknowns]][knots[[-k-1]]]}[[Range[1,k-2]]]
   
condsSupThirdEq[4,knots_,uq_] :=
  Module[{d2,poles},
         d2 = Derivative[2][Sp[4,knots,uq]];
           poles = d2[[3]];
         {(poles[[2]] - poles[[1]])/(knots[[4+1]] - knots[[4]])
           == (poles[[3]] - poles[[2]])/(knots[[4+2]] - knots[[4+1]]),
          (poles[[-1]] - poles[[-2]])/(knots[[-4]] - knots[[-4-1]])
           == (poles[[-2]] - poles[[-3]])/(knots[[-5]] - knots[[-6]])}
           ]
           
cordLengthKnots[k_,ipts_] :=
  Module[{knots = {0},i},
           For[i = 2, i <= Length[ipts], i++,
               AppendTo[knots, 
                        knots[[i-1]] +
                        Dist[ipts[[i-1]], ipts[[i]]]]];
           Join[Table[-k+i,{i,1,k-1}],
                         knots,
                Table[knots[[-1]]+i, {i,1,k-1}]]]
                
centripeteKnots[k_, ipts_] :=
  Module[{knots = {0},i},
           For[i = 2, i <= Length[ipts], i++,
               AppendTo[knots, 
                        knots[[i-1]] +
                        Sqrt[Dist[ipts[[i-1]], 
                             ipts[[i]]]]]];
           Join[Table[-k+i,{i,1,k-1}],
                         knots,
                Table[knots[[-1]]+i, {i,1,k-1}]]]
                     
SplineFit[ipts_?VectorQ, 
          k_Integer, 
          parFunc_Symbol, 
          condSuppFunc_Symbol] :=
  Module[{knots, 
          unknowns, 
          ptsEqs, 
          condsSupp, 
          pars, 
          sols,
          fcts,
          u},
         knots = parFunc[k,ipts] // N;
         unknowns = unknown[k,ipts];
         ptsEqs = interpolsPtsEqs[
                     k,knots, ipts, unknowns];
         condsSupp = condSuppFunc[k,knots,unknowns];
         sols = Solve[Join[ptsEqs,condsSupp],unknowns];
         fcts = Sp[k,knots,unknowns][u] /. sols[[1]];
         pars = knots[[Range[k,Length[knots]-k+1]]];
         Sp[k,knots,unknowns /. sols[[1]]]]

SplineFit[pars_?VectorQ, 
          ipts_?VectorQ, 
          k_Integer, 
          condSuppFunc_Symbol] :=
  Module[{knots,
          unknowns, 
          ptsEqs, 
          condsSupp,
          sols,
          fcts,
          u},
         knots = extendParameters[pars,k];
         unknowns = unknown[k,ipts];
         ptsEqs = interpolsPtsEqs[
                     k,knots, ipts, unknowns];
         condsSupp = condSuppFunc[k,knots,unknowns];
         sols = Solve[Join[ptsEqs,condsSupp],unknowns];
         fcts = Sp[k,knots,unknowns][u] /. sols[[1]];
         Sp[k,knots,unknowns /. sols[[1]]]]       

SplineFit[iptsND_?MatrixQ, 
          k_Integer, 
          parFunc_Symbol, 
          condSuppFunc_Symbol] :=
  Module[{ipts = Transpose[iptsND],
          knots, 
          unknowns, 
          ptsEqs, 
          condsSupp, 
          sols,
          fcts,
          u},
         knots = parFunc[k,iptsND] // N;
         unknowns = unknown[k,iptsND];
         ptsEqs = 
           Map[interpolsPtsEqs[
                 k,knots, #, unknowns]&, ipts];
         condsSupp = condSuppFunc[k,knots,unknowns];
         sols = 
           Map[Solve[Join[#,condsSupp],unknowns]&,ptsEqs];
         fcts = Map[Sp[k,knots,unknowns][u] /. #[[1]]&,sols];
         Sp[k,knots,
            Transpose[Map[unknowns /. #[[1]]&,sols]]
           ]
       ]

End[]

EndPackage[]
