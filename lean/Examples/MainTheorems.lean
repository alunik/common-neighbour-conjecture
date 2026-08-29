import Examples.MainTheorems.ProofAliases

/-!
# Theorems 1.2 and 1.3

The complete reader-facing interface: the four non-Mathlib definitions needed
to read the results, followed by the two results. Proofs are linked from
`Examples.MainTheorems.ProofAliases`.
-/

namespace SaxlCounterexamples.MainTheorems

/-- A finite permutation group: a finite group acting faithfully on a finite
set. -/
structure FinitePermutationGroup where
  G : Type
  Point : Type
  [group : Group G]
  [action : MulAction G Point]
  [finiteGroup : Finite G]
  [finitePoint : Finite Point]
  [faithful : FaithfulSMul G Point]

attribute [instance] FinitePermutationGroup.group
  FinitePermutationGroup.action FinitePermutationGroup.finiteGroup
  FinitePermutationGroup.finitePoint FinitePermutationGroup.faithful

namespace FinitePermutationGroup

/-- A set of points is a base if only the identity fixes every point in it. -/
def IsBase (P : FinitePermutationGroup) (B : Set P.Point) : Prop :=
  ∀ g : P.G, (∀ x ∈ B, g • x = x) → g = 1

/-- The least size of a base. An injective map from `Fin n` represents an
`n`-element base; existential quantification makes its enumeration irrelevant. -/
noncomputable def baseSize (P : FinitePermutationGroup) : Nat :=
  sInf {n : Nat | ∃ b : Fin n → P.Point,
    Function.Injective b ∧ P.IsBase (Set.range b)}

/-- The Saxl graph of a finite permutation group. Two distinct points are
adjacent exactly when they lie together in a base of minimum size. Thus this
is the ordinary Saxl graph at base size two and the generalized Saxl graph at
larger base sizes. -/
noncomputable def saxlGraph (P : FinitePermutationGroup) :
    SimpleGraph P.Point where
  Adj x y :=
    x ≠ y ∧ ∃ b : Fin P.baseSize → P.Point,
      Function.Injective b ∧ P.IsBase (Set.range b) ∧
        x ∈ Set.range b ∧ y ∈ Set.range b
  symm := ⟨by
    rintro x y ⟨hxy, b, hinj, hbase, hx, hy⟩
    exact ⟨hxy.symm, b, hinj, hbase, hy, hx⟩⟩
  loopless := ⟨by
    rintro x ⟨hxx, -⟩
    exact hxx rfl⟩

end FinitePermutationGroup

/-- **Theorem 1.2.** For every prescribed bound there is a primitive
permutation group of degree at least that bound and base size two whose Saxl
graph has two nonadjacent vertices with no common neighbour. -/
theorem theorem1_2 (degreeBound : Nat) :
    ∃ P : FinitePermutationGroup,
      degreeBound ≤ Nat.card P.Point ∧
        MulAction.IsPreprimitive P.G P.Point ∧
        P.baseSize = 2 ∧
        ∃ x y : P.Point, x ≠ y ∧
          ¬ P.saxlGraph.Adj x y ∧
          P.saxlGraph.commonNeighbors x y = ∅ :=
  theorem1_2_from_internal(degreeBound)

/-- **Theorem 1.3.** For every `B ≥ 3` and every prescribed bound there is a
primitive permutation group of degree at least that bound and base size `B`
whose Saxl graph has two distinct vertices with no common neighbour. -/
theorem theorem1_3 (B degreeBound : Nat) (hB : 3 ≤ B) :
    ∃ P : FinitePermutationGroup,
      degreeBound ≤ Nat.card P.Point ∧
        MulAction.IsPreprimitive P.G P.Point ∧
        P.baseSize = B ∧
        ∃ x y : P.Point, x ≠ y ∧
          P.saxlGraph.commonNeighbors x y = ∅ :=
  theorem1_3_from_internal(B, degreeBound, hB)

end SaxlCounterexamples.MainTheorems
