(* ::Package:: *)

(* ::Section::Closed:: *)
(*Theoretical summary*)


(* ::Text:: *)
(*Let (\[ScriptCapitalM], Subscript[g, ab]) be a space-time and consider a (normalized) vector field *)
(*\!\(\*OverscriptBox[\(n\), \(\[RightVector]\)]\) defined on \[ScriptCapitalM]. The vector field *)
(*\!\(\*OverscriptBox[\(n\), \(\[RightVector]\)]\) can be taken as the basic element to find the orthogonal splitting of the main geometric quantities of (\[ScriptCapitalM], Subscript[g, ab]):*)


(* ::Item:: *)
(*The orthogonal splitting of the metric tensor:*)


(* ::Text:: *)
(*	Subscript[g, ab]=Subscript[h, ab]+(Subscript[n, a] Subscript[n, b])/g( *)
(*\!\(\*OverscriptBox[\(n\), \(\[RightVector]\)]\),*)
(*\!\(\*OverscriptBox[\(n\), \(\[RightVector]\)]\)),*)


(* ::Text:: *)
(*where Subscript[h, ab] is the first fundamental form.*)


(* ::Item:: *)
(*The orthogonal splitting of the Levi-Civita covariant derivative of  *)
(*\!\(\*OverscriptBox[\(n\), \(\[RightVector]\)]\):*)


(* ::Text:: *)
(*	\!\( *)
(*\*SubscriptBox[\(\[Del]\), \(a\)]*)
(*\*SubscriptBox[\(n\), \(b\)]\)=Subscript[K, ab]+Subscript[A, b] Subscript[n, a],*)


(* ::Text:: *)
(*where Subscript[K, ab] is the second fundamental form. Subscript[K, ab] is symmetric if *)
(*\!\(\*OverscriptBox[\(n\), \(\[RightVector]\)]\) is hypersurface forming.*)


(* ::Item:: *)
(*The induced covariant derivative Subscript[D, a]: *)


(* ::Text:: *)
(*	Subscript[D, a] Subscript[T, Subscript[b, 1]\[Ellipsis] Subscript[b, p]]=Subscript[h, a]^p Subscript[h, Subscript[b, 1]]^Subscript[c, 1] \[CenterEllipsis] Subscript[h, Subscript[b, p]]^Subscript[c, p] \!\( *)
(*\*SubscriptBox[\(\[Del]\), \(p\)]*)
(*\*SubscriptBox[\(T\), \( *)
(*\*SubscriptBox[\(c\), \(1\)] \[Ellipsis]\ *)
(*\*SubscriptBox[\(c\), \(p\)]\)]\) ,*)


(* ::Text:: *)
(*where \[Del] is the Levi-Civita covariant derivative compatible with Subscript[g, ab].*)


(* ::Item:: *)
(*Orthogonal splitting of the Weyl, Riemann and Einstein tensors.*)


(* ::Item:: *)
(*Gau\[SZ] equation: metrich[-a, e] metrich[-b, f] metrich[-c, h] metrich[-d, j] RiemannCD[-e, -f, -h, -j]==Riemanncd[-a, -b, -c, -d]+(ExtrinsicKmetrich[-a, -d] ExtrinsicKmetrich[-b, -c])/xAct`xTensor`Scalar[n[-a] n[a]]-(ExtrinsicKmetrich[-a, -c] ExtrinsicKmetrich[-b, -d])/xAct`xTensor`Scalar[n[-a] n[a]]*)


(* ::Item:: *)
(*Codazzi equation: metrich[-a, e] metrich[-b, f] metrich[-c, h] n[d] RiemannCD[-d, -h, -e, -f]=-(cd[-a][ExtrinsicKmetrich[-b, -c]])+cd[-b][ExtrinsicKmetrich[-a, -c]]*)


(* ::Item:: *)
(*Ricci equation: metrich[-a, e] metrich[-b, f] n[c] n[d] RiemannCD[-c, -e, -d, -f]=ExtrinsicKmetrich[-a, c] ExtrinsicKmetrich[-b, -c]-Accelerationn[-a] Accelerationn[-b] xAct`xTensor`Scalar[n[-a] n[a]]+xAct`xTensor`Scalar[n[-a] n[a]] (cd[-b][Accelerationn[-a]])-xAct`xTensor`LieD[n[a]][ExtrinsicKmetrich[-c, -d]]*)


(* ::Item:: *)
(*ADM formalism (dimension 4). *)


(* ::Text:: *)
(*	Evolution equations: *)
(*	xAct`xTensor`LieD[n[a]][ExtrinsicKmetrich[-b, -d]]==1/6 (-6 Riccicd[-b, -d]+3 metrich[-b, -d] S[a, -a]+3 S[-b, -d]+(6 ExtrinsicKmetrich[a, -a] ExtrinsicKmetrich[-b, -d])/xAct`xTensor`Scalar[n[-a] n[a]]+(6 ExtrinsicKmetrich[-b, a] ExtrinsicKmetrich[-d, -a] (-1+xAct`xTensor`Scalar[n[-a] n[a]]))/xAct`xTensor`Scalar[n[-a] n[a]]-6 Accelerationn[-b] Accelerationn[-d] xAct`xTensor`Scalar[n[-a] n[a]]-3 S[-b, -d] xAct`xTensor`Scalar[n[-a] n[a]]-4 metrich[-b, -d] xAct`xTensor`Scalar[S[a, -a]]+2 metrich[-b, -d] xAct`xTensor`Scalar[n[-a] n[a]] xAct`xTensor`Scalar[S[a, -a]]-metrich[-b, -d] \[CapitalXi][]-(4 metrich[-b, -d] \[CapitalXi][])/xAct`xTensor`Scalar[n[-a] n[a]]+6 xAct`xTensor`Scalar[n[-a] n[a]] (cd[-d][Accelerationn[-b]]))*)
(*	xAct`xTensor`LieD[n[s]][metrich[-a, -b]]=2 ExtrinsicKmetrich[-a, -b]*)
(**)
(*	Constraint equations:*)
(**)
(*	Hamiltonian constraint: \[CapitalXi][]==1/2 (xAct`xTensor`Scalar[ExtrinsicKmetrich[a, -a]]^2-xAct`xTensor`Scalar[ExtrinsicKmetrich[-a, -b] ExtrinsicKmetrich[a, b]]-RicciScalarcd[] xAct`xTensor`Scalar[n[-a] n[a]])*)
(*	Momentum constraint: J[-d]==cd[-c][ExtrinsicKmetrich[-d, c]]-cd[-d][ExtrinsicKmetrich[c, -c]]*)


(* ::Text:: *)
(*A more advanced case study involving the analysis of a general timelike congruence and its kinematical quantities can be found here:*)
(*https://github.com/xAct-contrib/examples/blob/master/TimelikeCongruence.nb*)


(* ::Text:: *)
(**)
(**)


(* ::Section::Closed:: *)
(*Initialization*)


(* ::Text:: *)
(*Package loading*)


(* ::Input:: *)
(*(*Exit*)*)


(* ::Input:: *)
(*<<xAct`xTensor`*)


(* ::Text:: *)
(*Some preliminary arrangements.*)


(* ::Input:: *)
(*<<xAct`ShowTime1`*)


(* ::Input:: *)
(*$ShowTimeThreshold*)


(* ::Input:: *)
(*$PrePrint=ScreenDollarIndices;*)


(* ::Input:: *)
(*SetOptions[ContractMetric,AllowUpperDerivatives->True]*)


(* ::Input:: *)
(*ApplyRule[expr_]:=MakeRule[Evaluate[List@@expr],MetricOn->All,ContractMetrics->True]*)


(* ::Input:: *)
(*ApplyRuleN[expr_]:=MakeRule[Evaluate[List@@expr]]*)


(* ::Input:: *)
(*$Equations={};*)


(* ::Input:: *)
(*SetAttributes[AddEquation,HoldAll]*)
(*AddEquation[expr_Symbol]:=AppendTo[$Equations,Hold[expr]]*)


(* ::Input:: *)
(**)
(**)
(**)


(* ::Section::Closed:: *)
(*Space-time definition*)


(* ::Text:: *)
(*We define a 4-dimensional manifold (spacetime manifold).*)


(* ::Input:: *)
(*DefManifold[M4,4,{a,b,c,d,e,p,q,r,l,f,h,m,x,y,z,j,k,s}]*)


(* ::Text:: *)
(*We define the spacetime metric. *)


(* ::Input:: *)
(*DefMetric[-1,g[-a,-b],CD,{";","\[EmptyDownTriangle]"}]*)


(* ::Text:: *)
(*Output printing options for some quantities*)


(* ::Input:: *)
(*{PrintAs[epsilong]^="\[Eta]",PrintAs[RiemannCD]^="R",PrintAs[WeylCD]^="W",PrintAs[RicciCD]^="R",PrintAs[RicciScalarCD]^="R",PrintAs[EinsteinCD]^="G"};*)


(* ::Text:: *)
(*We define a normal vector:*)


(* ::Input:: *)
(*DefTensor[n[a],M4]*)


(* ::Text:: *)
(*We assume that the normal vector is constant:*)


(* ::Input:: *)
(*DefConstantSymbol[\[CurlyPhi]]*)


(* ::Input:: *)
(*DefInertHead[ABS]*)


(* ::Text:: *)
(*ABS represents absolute value:*)


(* ::Input:: *)
(*\[CurlyPhi]==Sqrt[Scalar@ABS[n[-a]n[a]]]*)


(* ::Input:: *)
(*Define\[CurlyPhi]=%;*)


(* ::Input:: *)
(*AddEquation@Define\[CurlyPhi];*)


(* ::Input:: *)
(*Power[#,2]&/@Define\[CurlyPhi]*)


(* ::Input:: *)
(*Solve[%,Scalar@ABS[n[-a]n[a]]]//Flatten*)


(* ::Input:: *)
(*Equal@@First@%*)


(* ::Input:: *)
(*NoScalar/@%*)


(* ::Input:: *)
(*NormSquareTo\[CurlyPhi]=%;*)


(* ::Input:: *)
(*AddEquation@NormSquareTo\[CurlyPhi];*)


(* ::Text:: *)
(*Introduce a constant symbol representing \[PlusMinus]1:*)


(* ::Input:: *)
(*DefConstantSymbol[\[Epsilon]]*)


(* ::Text:: *)
(*Basic properties:*)


(* ::Input:: *)
(*\[Epsilon]/:Power[\[Epsilon],n_?EvenQ]:=1*)


(* ::Input:: *)
(*\[Epsilon]/:Power[\[Epsilon],n_?OddQ]:=\[Epsilon]*)


(* ::Text:: *)
(*By definition:*)


(* ::Input:: *)
(*n[-a]n[a]==\[Epsilon] \[CurlyPhi]^2*)


(* ::Input:: *)
(*PutScalar/@%*)


(* ::Input:: *)
(*NormSquareTo\[CurlyPhi]NoABS=%;*)


(* ::Input:: *)
(*AddEquation@NormSquareTo\[CurlyPhi]NoABS;*)


(* ::Input:: *)
(*AutomaticRules[n,ApplyRule@NormSquareTo\[CurlyPhi]NoABS]*)


(* ::Input:: *)
(*$Rules*)


(* ::Section::Closed:: *)
(*Definition of the spatial metric*)


(* ::Text:: *)
(*Sign convention to define the extrinsic curvature and the acceleration. These are global variables which can be modified by the user.*)


(* ::Input:: *)
(*$ExtrinsicKSign=1;$AccelerationSign=1;*)


(* ::Text:: *)
(*Next we introduce the induced metric. This is a metric called metrich[-a,-b] and its Levi-Civita covariat derivative is cd (spatial covariant derivative).*)


(* ::Input:: *)
(*DefMetric[1,metrich[-a,-b],cd,{"|","D"},InducedFrom->{g,n},PrintAs->"h"]*)


(* ::Text:: *)
(*Along with the spatial metric and the spatial covariant derivative many new quantities are introduced. Among them we find the acceleration of the unit normal n[a], the extrinsic curvature, the spatial volume element,  the curvature tensors of the spatial covariant derivative and the orthogonal projector operator. *)


(* ::Input:: *)
(*ServantsOf@metrich*)


(* ::Input:: *)
(*ServantsOf@cd*)


(* ::Text:: *)
(* We set the output printing for some of the objects which have been defined together with the induced metric.*)


(* ::Input:: *)
(*{PrintAs[ExtrinsicKmetrich]^="K",PrintAs[Riemanncd]^="\!\(\*OverscriptBox[\(R\), \(_\)]\)",PrintAs[Riccicd]^="\!\(\*OverscriptBox[\(R\), \(_\)]\)",PrintAs[RicciScalarcd]^="\!\(\*OverscriptBox[\(R\), \(_\)]\)",PrintAs[Accelerationn]^="A",PrintAs[epsilonmetrich]^="\[CurlyEpsilon]"};*)


(* ::Text:: *)
(*Output printing.*)


(* ::Input:: *)
(*{Accelerationn[-a],ExtrinsicKmetrich[-a,-b],Riemanncd[-a,-b,-c,d],Riccicd[-a,-b]}*)


(* ::Text:: *)
(*All these quantities are spatial. This means that their contraction with n[-a] vanishes.*)


(* ::Input:: *)
(*Times[#,n[a]]&/@%*)


(* ::Text:: *)
(*Also they are invariant with respect to the action of the spacetime metric*)


(* ::Input:: *)
(*Times[#,metrich[a,f]]&/@%%*)


(* ::Text:: *)
(*Some important predefined relations.*)


(* ::Input:: *)
(*$ExtrinsicKSign=1;$AccelerationSign=1;*)


(* ::Input:: *)
(*CD[-a]@n[-b]==(CD[-a]@n[-b]//GradNormalToExtrinsicK)*)


(* ::Input:: *)
(*CDN=%;*)


(* ::Input:: *)
(*AddEquation@CDN;*)


(* ::Input:: *)
(*$Equations*)


(* ::Input:: *)
(*%//ReleaseHold*)


(* ::Text:: *)
(*From here we get the expression of the acceleration and the second fundamental form.*)


(* ::Input:: *)
(*Times[#,n[a]]&/@CDN*)


(* ::Input:: *)
(*ToCanonical/@(ContractMetric/@%)*)


(* ::Input:: *)
(*Expand/@%*)


(* ::Input:: *)
(*PutScalar/@%*)


(* ::Input:: *)
(*Simplification@%*)


(* ::Text:: *)
(*The projection of any expression in all its free indices is carried out with ProjectWith:*)


(* ::Input:: *)
(*?ProjectWith*)


(* ::Text:: *)
(*Example:*)


(* ::Input:: *)
(*ProjectWith[metrich]/@CDN*)


(* ::Input:: *)
(*Solve[%,ExtrinsicKmetrich[-b,-a]]//Flatten*)


(* ::Input:: *)
(*%[[1]]/.Rule->Equal*)


(* ::Input:: *)
(*ToCanonical/@%*)


(* ::Input:: *)
(*ExtrinsicToGradNormal=%;*)


(* ::Input:: *)
(*AddEquation@ExtrinsicToGradNormal;*)


(* ::Input:: *)
(*$Equations*)


(* ::Input:: *)
(*Projectormetrich[RiemannCD[-a,-b,-c,-d]]==ProjectWith[metrich][RiemannCD[-a,-b,-c,-d]]*)


(* ::Text:: *)
(*The orthogonal projection does not change quantities which are already spatial.  We add these properties here  for the Lie and covariant derivatives of spatial quantities:*)


(* ::Input:: *)
(*HoldPattern@Projectormetrich@LieD[n[s_]]@expr_:=LieD[n[s]]@expr/;OrthogonalToVectorQ[n]@expr;*)


(* ::Input:: *)
(*HoldPattern@Projectormetrich@CD[a_]@expr_:=cd[a]@expr/;OrthogonalToVectorQ[n]@expr;*)


(* ::Input:: *)
(*LieD[n[a]]@Accelerationn[-b]*)


(* ::Input:: *)
(*n[b]%*)


(* ::Text:: *)
(*Predicate function to check whether a tensor or tensorial expressions are orthogonal or not:*)


(* ::Input:: *)
(*?OrthogonalToVectorQ*)


(* ::Input:: *)
(**)
(**)
(**)


(* ::Section::Closed:: *)
(*Orthogonal decomposition of the 4-dimensional volume element*)


(* ::Text:: *)
(*Our convention to define the spatial volume element is given by:*)


(* ::Input:: *)
(*epsilong[-a,-b,-c,-d]n[a]==epsilonmetrich[-b,-c,-d]/\[CurlyPhi]*)


(* ::Input:: *)
(*Epsilon4ToEpsilon3=%*)


(* ::Input:: *)
(*AddEquation@Epsilon4ToEpsilon3;*)


(* ::Text:: *)
(*where epsilong[-a, -b, -c, -d] is the 4-dimensional volume element.  Therefore*)


(* ::Input:: *)
(*Epsilon4ToEpsilon3[[2]]==Epsilon4ToEpsilon3[[1]]*)


(* ::Input:: *)
(*Times[#,\[CurlyPhi]]&/@%*)


(* ::Input:: *)
(*Epsilon3ToEpsilon4=%*)


(* ::Input:: *)
(*AddEquation@Epsilon3ToEpsilon4;*)


(* ::Text:: *)
(*We wish to find the orthogonal splitting of epsilong[-a, -b, -c, -d]. To that end we use the command InducedDecomposition. This command enables us to compute the orthogonal splitting of any tensor or tensorial expression:*)


(* ::Input:: *)
(*?InducedDecomposition*)


(* ::Input:: *)
(*epsilong[-a,-b,-c,-d]==InducedDecomposition[epsilong[-a,-b,-c,-d],{metrich,n}]*)


(* ::Text:: *)
(*Some of these terms are obviously zero. This fact is taken care of during the canonicalisation.*)


(* ::Input:: *)
(*ToCanonical/@%*)


(* ::Text:: *)
(*The spatial volume element appears explicitly in other terms. To see it we apply the rule epsilong[-a, -b, -c, -d] n[a]==epsilonmetrich[-b, -c, -d]*)


(* ::Input:: *)
(*%/.ApplyRule@Epsilon4ToEpsilon3*)


(* ::Text:: *)
(*Finally Projectormetrich[epsilong[-a, -b, -c, -d]] is essentially the four dimensional determinant of the spatial metric and hence it vanishes.*)


(* ::Input:: *)
(*%/._Projectormetrich:>0*)


(* ::Input:: *)
(*%/.ApplyRule@NormSquareTo\[CurlyPhi]NoABS*)


(* ::Text:: *)
(*We keep this important rule for later use:*)


(* ::Input:: *)
(*DecomposeEpsilong=%*)


(* ::Input:: *)
(*AddEquation@DecomposeEpsilong;*)


(* ::Text:: *)
(*Recall that Projectormetrich is an inert head representing the orthogonal projection operation acting on any tensorial expression:*)


(* ::Input:: *)
(*MasterOf@Projectormetrich*)


(* ::Input:: *)
(**)
(**)
(**)


(* ::Section::Closed:: *)
(*Decomposition of a spacetime covariant derivative: orthogonal decomposition of CD[-b][Accelerationn[-a]]*)


(* ::Text:: *)
(*This illustrates how to find the orthogonal splitting of the space-time covariant derivative of a spatial quantity.*)


(* ::Text:: *)
(*We start from the formula:*)


(* ::Input:: *)
(*LieD[n[s]]@Accelerationn[-a]==LieDToCovD[LieD[n[s]]@Accelerationn[-a],CD]*)


(* ::Input:: *)
(*ScreenDollarIndices/@%*)


(* ::Input:: *)
(*MapAt[Hold,%,Position[%,n[b]CD[-b]@Accelerationn[-a]]]*)


(* ::Input:: *)
(*Solve[%,Hold[n[b]CD[-b]@Accelerationn[-a]]]//Flatten//ReleaseHold*)


(* ::Input:: *)
(*%[[1]]/.Rule->Equal*)


(* ::Input:: *)
(*%//GradNormalToExtrinsicK*)


(* ::Input:: *)
(*ToCanonical/@%*)


(* ::Input:: *)
(*NCDAcceleration=%;*)


(* ::Input:: *)
(*AddEquation@NCDAcceleration;*)


(* ::Text:: *)
(*We are ready to obtain the orthogonal splitting of CD[-a][Accelerationn[-b]]*)


(* ::Input:: *)
(*CD[-a]@Accelerationn[-b]==InducedDecomposition[CD[-a]@Accelerationn[-b],{metrich,n}]*)


(* ::Input:: *)
(*%/.ApplyRule@NCDAcceleration*)


(* ::Input:: *)
(*Expand//@%*)


(* ::Text:: *)
(*We want that everything is left in terms of xAct`xTensor`LieD[n[a]][Accelerationn[-b]] which is spatial:*)


(* ::Input:: *)
(*%/.MakeRule[{n[e]CD[-a]@Accelerationn[-e],-CD[-a]@n[e]Accelerationn[-e]},MetricOn->All,ContractMetrics->True]*)


(* ::Input:: *)
(*%//GradNormalToExtrinsicK*)


(* ::Input:: *)
(*Expand//@%*)


(* ::Input:: *)
(*%/.Projectormetrich->ProjectWith[metrich]*)


(* ::Input:: *)
(*ToCanonical/@%*)


(* ::Input:: *)
(*NoScalar/@%*)


(* ::Input:: *)
(*ToCanonical/@%*)


(* ::Input:: *)
(*PutScalar/@%*)


(* ::Input:: *)
(*DecomposeCDAcceleration=%*)


(* ::Input:: *)
(*AddEquation@DecomposeCDAcceleration;*)


(* ::Input:: *)
(**)
(**)


(* ::Section::Closed:: *)
(*Orthogonal decomposition of Subscript[\[EmptyDownTriangle], c] *)
(*\!\(\*SubsuperscriptBox[\(K\), \(ab\), \(\ \ \)]\)*)


(* ::Text:: *)
(*The computation is similar to the previous case with the acceleration:*)


(* ::Text:: *)
(*We start from the formula:*)


(* ::Input:: *)
(*LieD[n[s]]@ExtrinsicKmetrich[-a,-b]==LieDToCovD[LieD[n[s]]@ExtrinsicKmetrich[-a,-b],CD]*)


(* ::Input:: *)
(*ScreenDollarIndices/@%*)


(* ::Input:: *)
(*MapAt[Hold,%,Position[%,n[c]CD[-c]@ExtrinsicKmetrich[-a,-b]]]*)


(* ::Input:: *)
(*Solve[%,Hold[n[c]CD[-c]@ExtrinsicKmetrich[-a,-b]]]//Flatten//ReleaseHold*)


(* ::Input:: *)
(*%[[1]]/.Rule->Equal*)


(* ::Input:: *)
(*%//GradNormalToExtrinsicK*)


(* ::Input:: *)
(*ToCanonical/@%*)


(* ::Input:: *)
(*NCDExtrinsicK=%;*)


(* ::Input:: *)
(*AddEquation@NCDExtrinsicK;*)


(* ::Text:: *)
(*We are ready to obtain the orthogonal splitting of CD[-a][ExtrinsicKmetrich[-b, -c]]*)


(* ::Input:: *)
(*CD[-a]@ExtrinsicKmetrich[-b,-c]==InducedDecomposition[CD[-a]@ExtrinsicKmetrich[-b,-c],{metrich,n}]*)


(* ::Input:: *)
(*%/.ApplyRule@NCDExtrinsicK*)


(* ::Input:: *)
(*Expand//@%*)


(* ::Input:: *)
(*%/.MakeRule[{n[e]CD[-a]@ExtrinsicKmetrich[-b,-e],-CD[-a]@n[e]ExtrinsicKmetrich[-b,-e]},MetricOn->All,ContractMetrics->True]*)


(* ::Input:: *)
(*%//GradNormalToExtrinsicK*)


(* ::Input:: *)
(*Expand//@%*)


(* ::Input:: *)
(*%/.Projectormetrich->ProjectWith[metrich]*)


(* ::Input:: *)
(*ToCanonical/@%*)


(* ::Text:: *)
(*This is the final form of the orthogonal splitting: note how al the terms are products of the normal and spatial quantities:*)


(* ::Input:: *)
(*%/.$Rules*)


(* ::Input:: *)
(*DecomposeCDExtrinsicK=%*)


(* ::Input:: *)
(*AddEquation@DecomposeCDExtrinsicK;*)


(* ::Input:: *)
(**)
(**)


(* ::Section::Closed:: *)
(*Weyl tensor electric and magnetic parts*)


(* ::Text:: *)
(*The Weyl tensor electric part (this illustrates  how to define a tensor which is orthogonal to the normal vector field):*)


(* ::Input:: *)
(*DefTensor[Ee[-a,-b],M4,Symmetric[{-a,-b}],OrthogonalTo->{n[a]},ProjectedWith->{metrich[a,-c]},PrintAs->"E"]*)


(* ::Text:: *)
(*The Weyl tensor magnetic part:*)


(* ::Input:: *)
(*DefTensor[B[-a,-b],M4,Symmetric[{-a,-b}],OrthogonalTo->{n[a]},ProjectedWith->{metrich[a,-c]},PrintAs->"B"]*)


(* ::Text:: *)
(*For later use we define the Weyl tensor right dual:*)


(* ::Input:: *)
(*DefTensor[WCD01[-a,-b,-c,-d],M4,GenSet[-Cycles[{-a,-b}],-Cycles[{-c,-d}]],PrintAs->"\!\(\*SuperscriptBox[\(W\), \(*\)]\)"]*)


(* ::Text:: *)
(*The Weyl tensor electric and magnetic parts are traceless (this illustrates the use of AutomaticRules):*)


(* ::Input:: *)
(*AutomaticRules[Ee,MakeRule[{Ee[a,-a],0},MetricOn->All,ContractMetrics->True]]*)


(* ::Input:: *)
(*AutomaticRules[B,MakeRule[{B[a,-a],0},MetricOn->All,ContractMetrics->True]]*)


(* ::Text:: *)
(*Our conventions for the definition of the electric and magnetic parts of Weyl tensor:*)


(* ::Input:: *)
(*WeylCD[-a,-b,-c,-d]n[b]n[d]==Ee[-a,-c]*)


(* ::Input:: *)
(*WCD01[-a,-b,-c,-d]n[b]n[d]==B[-a,-c]*)


(* ::Input:: *)
(*AutomaticRules[WeylCD,MakeRule[{WeylCD[-a,-b,-c,-d]n[b]n[d],Ee[-a,-c]},MetricOn->All,ContractMetrics->True]]*)


(* ::Input:: *)
(*AutomaticRules[WCD01,MakeRule[{WCD01[-a,-b,-c,-d]n[b]n[d],B[-a,-c]},MetricOn->All,ContractMetrics->True]]*)


(* ::Input:: *)
(**)
(**)
(**)


(* ::Section::Closed:: *)
(*Orthogonal decomposition of the Weyl tensor*)


(* ::Text:: *)
(*We wish to find the orthogonal decomposition of the Weyl tensor. To that end we use again the command  InducedDecomposition. *)


(* ::Input:: *)
(*WeylCD[-a,-b,-c,-d]==InducedDecomposition[WeylCD[-a,-b,-c,-d],{metrich,n}]*)


(* ::Input:: *)
(*ToCanonical/@%*)


(* ::Input:: *)
(*Stage1=%*)


(* ::Text:: *)
(*We need to study separately some terms of this expression. These are Projectormetrich[n[e] WeylCD[-b, -e, -c, -d]] and Projectormetrich[WeylCD[-a, -b, -c, -d]]*)


(* ::Input:: *)
(**)
(**)
(**)
(**)


(* ::Section::Closed:: *)
(*Study of Projectormetrich[n[e] WeylCD[-b, -e, -c, -d]]*)


(* ::Text:: *)
(*We start from the identity*)


(* ::Input:: *)
(*Projectormetrich[n[e]*WeylCD[-b,-e,-c,-d]]==1/2n[e]metrich[-b,a]Hold@epsilonmetrich[-c,-d,-p]epsilonmetrich[f,h,p]WeylCD[-a,-e,-f,-h]*)


(* ::Text:: *)
(*In this last expression \!\( *)
(*TagBox[*)
(*SubsuperscriptBox["\[CurlyEpsilon]", "\<\"cdj\"\>", "\<\"   \"\>"],*)
(*Tensor]\) has been wrapped with Hold to prevent the automatic expansion of the product \!\(\**)
(*TagBox[*)
(*SubsuperscriptBox["\[CurlyEpsilon]", "\<\"   \"\>", "\<\"jfh\"\>"],*)
(*Tensor]\)\!\(\**)
(*TagBox[*)
(*SubsuperscriptBox["\[CurlyEpsilon]", "\<\"cdj\"\>", "\<\"   \"\>"],*)
(*Tensor]\). We replace one of the spatial volume elements by its definition in terms of the spacetime volume element*)


(* ::Input:: *)
(*MapAt[ReplaceAll[#,ApplyRule@Epsilon3ToEpsilon4]&,%,Position[%,epsilonmetrich[f,h,p]]]*)


(* ::Input:: *)
(*%//ReleaseHold*)


(* ::Text:: *)
(*Finally we replace the product {*)
(* {\[Eta], {*)
(*   { , f, h, p},*)
(*   {j,  ,  ,  }*)
(*  }}*)
(*} WeylCD[-a, -e, -f, -h]  by the right dual of the Weyl tensor. This yields an expression in terms of the Weyl tensor magnetic part.*)


(* ::Input:: *)
(*%/.MakeRule[{epsilong[-c,-d,-p,-q]WeylCD[-a,-b,p,q],2WCD01[-a,-b,-c,-d]},MetricOn->All,ContractMetrics->True]*)


(* ::Input:: *)
(*PmetrichNWeylCD=%*)


(* ::Input:: *)
(*AddEquation@PmetrichNWeylCD;*)


(* ::Input:: *)
(**)
(**)
(**)


(* ::Section::Closed:: *)
(*Study of Projectormetrich[WeylCD[-a, -b, -c, -d]]*)


(* ::Text:: *)
(*We start from the identity*)


(* ::Input:: *)
(*Projectormetrich[WeylCD[-a,-b,-c,-d]]==1/4Hold@epsilonmetrich[-a,-b,-p]Hold@epsilonmetrich[e,f,p]Hold@epsilonmetrich[-c,-d,-q]Hold@epsilonmetrich[h,j,q]WeylCD[-e,-f,-h,-j]*)


(* ::Input:: *)
(*MapAt[ReplaceAll[#,ApplyRule@Epsilon3ToEpsilon4]&,%,Position[%,Hold[epsilonmetrich[e,f,p]]]]*)


(* ::Input:: *)
(*MapAt[ReplaceAll[#,ApplyRule@Epsilon3ToEpsilon4]&,%,Position[%,Hold[epsilonmetrich[h,j,q]]]]*)


(* ::Input:: *)
(*%//ReleaseHold*)


(* ::Input:: *)
(*ToCanonical/@(ContractMetric/@%)*)


(* ::Input:: *)
(*ProjectorWeyl=%;*)


(* ::Input:: *)
(*AddEquation@ProjectorWeyl;*)


(* ::Input:: *)
(**)
(**)
(**)


(* ::Section::Closed:: *)
(*Final expression for the orthogonal decomposition of the Weyl tensor*)


(* ::Text:: *)
(*So we insert the results found in previous slides in Stage1*)


(* ::Input:: *)
(*Stage1*)


(* ::Input:: *)
(*%/.ApplyRule@PmetrichNWeylCD*)


(* ::Input:: *)
(*%/.ApplyRule@ProjectorWeyl*)


(* ::Input:: *)
(*%/.ApplyRule@NormSquareTo\[CurlyPhi]NoABS*)


(* ::Input:: *)
(*ToCanonical/@%*)


(* ::Input:: *)
(*%/.\[Epsilon]^(-2)->1*)


(* ::Text:: *)
(*This is the final form for the orthogonal decomposition of Riemann tensor*)


(* ::Input:: *)
(*DecomposeWeyl=%*)


(* ::Input:: *)
(**)
(**)


(* ::Section::Closed:: *)
(*Orthogonal splitting of the Einstein tensor*)


(* ::Text:: *)
(*First of all, we define the different quantities which appear in the orthogonal splitting of the Einstein tensor.*)


(* ::Text:: *)
(*Energy-momentum density:*)


(* ::Input:: *)
(*DefTensor[\[CapitalXi][],M4]*)


(* ::Text:: *)
(*By definition*)


(* ::Input:: *)
(*EinsteinCD[-a,-b]n[a]n[b]==\[CapitalXi][]*)


(* ::Input:: *)
(*EnergyDensity=%;*)


(* ::Input:: *)
(*AddEquation@EnergyDensity;*)


(* ::Text:: *)
(*Momentum flux*)


(* ::Input:: *)
(*DefTensor[J[-a],M4,OrthogonalTo->{n[a]},ProjectedWith->{metrich[a,-d]}]*)


(* ::Text:: *)
(*By definition*)


(* ::Input:: *)
(*Projectormetrich[EinsteinCD[-a,-c]n[c]]==J[-a]*)


(* ::Input:: *)
(*MomentumFlux=%;*)


(* ::Input:: *)
(*AddEquation@MomentumFlux;*)


(* ::Text:: *)
(*Stress tensor:*)


(* ::Input:: *)
(*DefTensor[S[-a,-b],M4,Symmetric[{-a,-b}],OrthogonalTo->{n[b]},ProjectedWith->{metrich[b,-d]}]*)


(* ::Text:: *)
(*By definition*)


(* ::Input:: *)
(*Projectormetrich[EinsteinCD[-a,-b]]==S[-a,-b]*)


(* ::Input:: *)
(*StressTensor=%;*)


(* ::Input:: *)
(*AddEquation@StressTensor;*)


(* ::Text:: *)
(*We perform the 1+3 decomposition of  Einstein tensor in terms of these quantities*)


(* ::Input:: *)
(*EinsteinCD[-a,-b]==InducedDecomposition[EinsteinCD[-a,-b],{metrich,n}]*)


(* ::Input:: *)
(*%/.ApplyRule@EnergyDensity/.ApplyRule@MomentumFlux/.ApplyRule@StressTensor*)


(* ::Input:: *)
(*%/.$Rules*)


(* ::Input:: *)
(*DecomposeEinstein=%;*)


(* ::Text:: *)
(*From here we deduce the decomposition of the Ricci tensor*)


(* ::Input:: *)
(*%//EinsteinToRicci*)


(* ::Input:: *)
(*%//MetricToProjector*)


(* ::Input:: *)
(*Solve[%,RicciCD[-a,-b]]//Flatten*)


(* ::Input:: *)
(*%[[1]]/.Rule->Equal*)


(* ::Input:: *)
(*%/.$Rules*)


(* ::Input:: *)
(*DecomposeRicci=%;*)


(* ::Input:: *)
(*AddEquation@DecomposeRicci;*)


(* ::Input:: *)
(*ReplaceIndex[Evaluate@DecomposeRicci,{-a->b}]*)


(* ::Input:: *)
(*PutScalar/@%*)


(* ::Input:: *)
(*%/.$Rules*)


(* ::Input:: *)
(*IndexSolve[%,RicciScalarCD[]]//Flatten//First*)


(* ::Input:: *)
(*DecomposeRicciScalar=%;*)


(* ::Input:: *)
(*AddEquation@DecomposeRicciScalar;*)


(* ::Input:: *)
(*DecomposeRicci/.DecomposeRicciScalar*)


(* ::Input:: *)
(*PutScalar/@%*)


(* ::Input:: *)
(*ToCanonical@%*)


(* ::Input:: *)
(*Collect[#,{n[-a_]n[-b_],n[-a_]},Simplify]&/@%*)


(* ::Input:: *)
(*DecomposeRicci=%;*)


(* ::Input:: *)
(*AddEquation@DecomposeRicci;*)


(* ::Section::Closed:: *)
(*Orthogonal splitting of the Riemann tensor*)


(* ::Text:: *)
(*We have now all the information to carry out the orthogonal splitting of the Riemann tensor:*)


(* ::Input:: *)
(*RiemannCD[-a,-b,-c,-d]==(RiemannCD[-a,-b,-c,-d]//RiemannToWeyl)*)


(* ::Input:: *)
(*%/.ApplyRule@DecomposeRicci*)


(* ::Input:: *)
(*%/.ApplyRule@DecomposeWeyl*)


(* ::Input:: *)
(*%/.DecomposeRicciScalar*)


(* ::Input:: *)
(*PutScalar/@%;*)


(* ::Input:: *)
(*%//MetricToProjector;*)


(* ::Input:: *)
(*ToCanonical@%;*)


(* ::Input:: *)
(*%/.$Rules;*)


(* ::Input:: *)
(*Collect[#,n[-a_]n[-b_],Simplification]&/@%*)


(* ::Input:: *)
(*%/.ApplyRule@NormSquareTo\[CurlyPhi]NoABS*)


(* ::Input:: *)
(*Collect[#,n[-a_]n[-b_],Simplification]&/@%*)


(* ::Input:: *)
(*DecomposeRiemann=%;*)


(* ::Input:: *)
(*AddEquation@DecomposeRiemann;*)


(* ::Section::Closed:: *)
(*Gau\[SZ], Codazzi and Ricci equations*)


(* ::Text:: *)
(*The following command generalizes the classical Gau\[SZ], Codazzi and Ricci equations:*)


(* ::Input:: *)
(*?GaussCodazzi*)


(* ::Text:: *)
(*Complete decomposition of the Riemann tensor:*)


(* ::Input:: *)
(*RiemannCD[-a,-b,-c,-d]==GaussCodazzi[RiemannCD[-a,-b,-c,-d],metrich]*)


(* ::Text:: *)
(*From the complete decomposition we can obtain the classical Gau\[SZ] relation:*)


(* ::Input:: *)
(*Projectormetrich/@%*)


(* ::Input:: *)
(*%/.Projectormetrich->ProjectWith[metrich]*)


(* ::Input:: *)
(*%/.$Rules*)


(* ::Input:: *)
(*GaussEquation=%;*)


(* ::Input:: *)
(*AddEquation@GaussEquation;*)


(* ::Text:: *)
(*From the complete decomposition we can obtain the classical Codazzi relation :*)


(* ::Input:: *)
(*RiemannCD[-d,-h,-e,-f]==GaussCodazzi[RiemannCD[-d,-h,-e,-f],metrich];*)


(* ::Input:: *)
(*Times[metrich[-a,e]metrich[-b,f]metrich[-c,h]n[d],#]&/@%;*)


(* ::Input:: *)
(*PutScalar/@%*)


(* ::Input:: *)
(*ToCanonical@%*)


(* ::Input:: *)
(*CodazziEquation=%;*)


(* ::Input:: *)
(*AddEquation@CodazziEquation;*)


(* ::Text:: *)
(*From the complete decomposition we can obtain the classical Ricci relation:*)


(* ::Input:: *)
(*RiemannCD[-c,-e,-d,-f]==GaussCodazzi[RiemannCD[-c,-e,-d,-f],metrich];*)


(* ::Input:: *)
(*Times[metrich[-a,e]metrich[-b,f]n[c]n[d],#]&/@%*)


(* ::Input:: *)
(*Expand/@%*)


(* ::Input:: *)
(*%/.ApplyRule@NCDExtrinsicK*)


(* ::Input:: *)
(*Expand/@%*)


(* ::Input:: *)
(*%/.Projectormetrich->ProjectWith[metrich]*)


(* ::Input:: *)
(*PutScalar/@%*)


(* ::Input:: *)
(*ToCanonical@%*)


(* ::Input:: *)
(*%/.$Rules*)


(* ::Input:: *)
(*RicciEquation=%;*)


(* ::Input:: *)
(**)
(**)
(**)


(* ::Section::Closed:: *)
(*ADM formalism*)


(* ::Text:: *)
(*Start with the generalized Gau\[SZ]-Codazzi decomposition of the Riemann tensor:*)


(* ::Input:: *)
(*RiemannCD[-a,-b,-c,-d]==GaussCodazzi[RiemannCD[-a,-b,-c,-d],metrich]*)


(* ::Text:: *)
(*Write this in terms of the Lie derivative of the second fundamental form:*)


(* ::Input:: *)
(*%/.ApplyRule@NCDExtrinsicK*)


(* ::Input:: *)
(*%/.Projectormetrich->ProjectWith@metrich*)


(* ::Input:: *)
(*%/.$Rules*)


(* ::Text:: *)
(*Use the decomposition of the Riemann tensor obtained before:*)


(* ::Input:: *)
(*%/.ApplyRule@DecomposeRiemann*)


(* ::Input:: *)
(*First@%-Last@%==0;*)


(* ::Input:: *)
(*Collect[#,{n[b_]n[a_],n[c_]},ToCanonical]&/@%*)


(* ::Input:: *)
(*MasterEquation=%;*)


(* ::Input:: *)
(*AddEquation@MasterEquation;*)


(* ::Text:: *)
(*The coefficients of the master equation with respect to the unit normal give us the ADM equations:*)


(* ::Input:: *)
(*Coefficient[MasterEquation,n[-a]n[-c]]*)


(* ::Input:: *)
(*Solve[%,LieD[n[a]]@ExtrinsicKmetrich[-b,-d]]//Flatten*)


(* ::Input:: *)
(*Equal@@First@%*)


(* ::Input:: *)
(*LieDExtrinsicKmetrich=%;*)


(* ::Input:: *)
(*AddEquation@LieDExtrinsicKmetrich;*)


(* ::Input:: *)
(*Coefficient[MasterEquation,n[-a]]/.n[-d]->0/.n[-c]->0*)


(* ::Text:: *)
(*We isolate from here the Weyl tensor magnetic part:*)


(* ::Input:: *)
(*Times[#,epsilonmetrich[c,d,-f]]&/@%*)


(* ::Input:: *)
(*ContractMetric@%*)


(* ::Input:: *)
(*Solve[%,B[-f,-b]]//Flatten*)


(* ::Input:: *)
(*Equal@@First@%*)


(* ::Input:: *)
(*MagneticPart=%;*)


(* ::Input:: *)
(*AddEquation@MagneticPart;*)


(* ::Input:: *)
(*ProjectWith[metrich]@MasterEquation*)


(* ::Input:: *)
(*%/.-c->a*)


(* ::Input:: *)
(*ToCanonical@ContractMetric@%*)


(* ::Input:: *)
(*Solve[%,Ee[-b,-d]]//Flatten*)


(* ::Input:: *)
(*Equal@@First@%*)


(* ::Input:: *)
(*ElectricPart=%;*)


(* ::Input:: *)
(*AddEquation@ElectricPart;*)


(* ::Input:: *)
(*LieDExtrinsicKmetrich/.ApplyRule@ElectricPart*)


(* ::Input:: *)
(*Simplification@%*)


(* ::Input:: *)
(*EvolutionExtrinsicADM=%;*)


(* ::Input:: *)
(*AddEquation@EvolutionExtrinsicADM;*)


(* ::Text:: *)
(*These are the constraint eqs in the ADM formalism:*)


(* ::Input:: *)
(*ElectricPart/.-b->d*)


(* ::Input:: *)
(*PutScalar/@%*)


(* ::Input:: *)
(*ToCanonical@%*)


(* ::Input:: *)
(*Solve[%,\[CapitalXi][]]//Flatten*)


(* ::Input:: *)
(*Equal@@First@%*)


(* ::Input:: *)
(*HamiltonianConstraint=%;*)


(* ::Input:: *)
(*AddEquation@HamiltonianConstraint;*)


(* ::Text:: *)
(*We compute next the momentum constraint:*)


(* ::Input:: *)
(*Coefficient[MasterEquation,n[-a]]/.n[-d]->0/.n[-c]->0*)


(* ::Input:: *)
(*%/.-b->c*)


(* ::Input:: *)
(*ToCanonical@%*)


(* ::Input:: *)
(*Solve[%,J[-d]]//Flatten*)


(* ::Input:: *)
(*Equal@@First@%*)


(* ::Input:: *)
(*MomentumConstraint=%;*)


(* ::Input:: *)
(*AddEquation@MomentumConstraint;*)


(* ::Input:: *)
(*{EvolutionExtrinsicADM,HamiltonianConstraint,MomentumConstraint}*)


(* ::Input:: *)
(*LieD[n[b]]@metrich[-a,-b]*)


(* ::Input:: *)
(*%//ToCanonical*)


(* ::Input:: *)
(*ElectricPart*)


(* ::Input:: *)
(*MagneticPart*)


(* ::Input:: *)
(**)
(**)
(**)
(**)


(* ::Title:: *)
(*Evaluating the generalized Gauss-Codazzi and Ricci equations in quadratic gravity*)


(* ::Text:: *)
(*Now we can define the quadratic-gravity equivalent Subscript[G^(quad), \[Mu]\[Nu]] of the Einstein-Tensor.*)


(* ::Code:: *)
(*DefTensor[\[DoubleStruckCapitalA][],M4]*)
(*DefTensor[\[DoubleStruckCapitalB][],M4]*)
(*DefTensor[\[DoubleStruckCapitalC][],M4]*)
(**)
(*DefTensor[Gquad[-a,-b],M4]*)
(*Gquad[a_,b_] := ( *)
(*	+ (2/3)(\[DoubleStruckCapitalA]-3\[DoubleStruckCapitalB])CD[a][CD[b][RicciScalarCD[]]] - 2\[DoubleStruckCapitalA] CD[-c][CD[c][RicciCD[a,b]]] + 1/3 (\[DoubleStruckCapitalA]+6\[DoubleStruckCapitalB])g[a,b]CD[-c][CD[c][RicciScalarCD[]]]*)
(*	- 4\[DoubleStruckCapitalA] RicciCD[c,d]RiemannCD[a,-c,b,-d] + 2(\[DoubleStruckCapitalB]+2\[DoubleStruckCapitalA]/3)RicciScalarCD[]RicciCD[a,b] + 1/2 g[a,b](2\[DoubleStruckCapitalA] RicciCD[c,d]RicciCD[-c,-d] - (\[DoubleStruckCapitalB]+2\[DoubleStruckCapitalA]/3)RicciScalarCD[]RicciScalarCD[])*)
(*);*)


(* ::Input:: *)
(*DefTensor[H[-a,-b],M4]*)
(*H[a_,b_]:=\[DoubleStruckCapitalC] EinsteinCD[a,b] + Gquad[a,b]*)


(* ::Input:: *)
(*H[-a,-b]*)


(* ::Input:: *)
(*GaussCodazzi[H[-a,-b]];*)
(*Projectormetrich/@%;*)
(*%/.Projectormetrich->ProjectWith[metrich];*)
(*%//ScreenDollarIndices*)



