import Examples.MainTheorems.Internal

/-!
# Proof aliases for the public theorem statements

These term macros hide the definitional transport from the internal proof
bundle to the deliberately duplicated reader-facing definitions.
-/

macro "theorem1_2_from_internal" "(" degreeBound:term ")" : term => do
  let publicGroup := Lean.mkIdent
    `SaxlCounterexamples.MainTheorems.FinitePermutationGroup
  let publicBaseSize := Lean.mkIdent
    `SaxlCounterexamples.MainTheorems.FinitePermutationGroup.baseSize
  let publicIsBase := Lean.mkIdent
    `SaxlCounterexamples.MainTheorems.FinitePermutationGroup.IsBase
  let publicSaxlGraph := Lean.mkIdent
    `SaxlCounterexamples.MainTheorems.FinitePermutationGroup.saxlGraph
  `(by
    obtain ⟨Q, hQ⟩ :=
      SaxlCounterexamples.MainTheorems.Internal.theorem1_2 $degreeBound
    let P : $publicGroup := {
      G := Q.G
      Point := Q.Point
      group := Q.group
      action := Q.action
      finiteGroup := Q.finiteGroup
      finitePoint := Q.finitePoint
      faithful := Q.faithful }
    refine ⟨P, ?_⟩
    simpa only [P,
      $(publicBaseSize):term, $(publicIsBase):term,
      $(publicSaxlGraph):term,
      SaxlCounterexamples.MainTheorems.Core.FinitePermutationGroup.baseSize,
      SaxlCounterexamples.MainTheorems.Core.FinitePermutationGroup.IsBase,
      SaxlCounterexamples.MainTheorems.Core.FinitePermutationGroup.saxlGraph]
      using hQ)

macro "theorem1_3_from_internal" "(" B:term "," degreeBound:term ","
    hB:term ")" : term => do
  let publicGroup := Lean.mkIdent
    `SaxlCounterexamples.MainTheorems.FinitePermutationGroup
  let publicBaseSize := Lean.mkIdent
    `SaxlCounterexamples.MainTheorems.FinitePermutationGroup.baseSize
  let publicIsBase := Lean.mkIdent
    `SaxlCounterexamples.MainTheorems.FinitePermutationGroup.IsBase
  let publicSaxlGraph := Lean.mkIdent
    `SaxlCounterexamples.MainTheorems.FinitePermutationGroup.saxlGraph
  `(by
    obtain ⟨Q, hQ⟩ :=
      SaxlCounterexamples.MainTheorems.Internal.theorem1_3
        $B $degreeBound $hB
    let P : $publicGroup := {
      G := Q.G
      Point := Q.Point
      group := Q.group
      action := Q.action
      finiteGroup := Q.finiteGroup
      finitePoint := Q.finitePoint
      faithful := Q.faithful }
    refine ⟨P, ?_⟩
    simpa only [P,
      $(publicBaseSize):term, $(publicIsBase):term,
      $(publicSaxlGraph):term,
      SaxlCounterexamples.MainTheorems.Core.FinitePermutationGroup.baseSize,
      SaxlCounterexamples.MainTheorems.Core.FinitePermutationGroup.IsBase,
      SaxlCounterexamples.MainTheorems.Core.FinitePermutationGroup.saxlGraph]
      using hQ)
