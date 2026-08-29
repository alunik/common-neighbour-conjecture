import Examples.EveryBase.Main
import Examples.MainTheorems.Definitions
import Saxl.Generalized

/-!
# Internal implementation of Theorems 1.2 and 1.3

Construction parameters, bridge lemmas, and proof machinery used by the
minimal public module `Examples.MainTheorems`.
-/

noncomputable section

namespace SaxlCounterexamples.MainTheorems.Internal

open Core Core.FinitePermutationGroup

private theorem isBase_range_iff_isBaseTuple
    (P : FinitePermutationGroup) {n : Nat} (b : Fin n → P.Point) :
    P.IsBase (Set.range b) ↔ Saxl.IsBaseTuple P.G P.Point b := by
  constructor
  · intro h g hg
    apply h g
    intro x hx
    obtain ⟨i, rfl⟩ := hx
    exact hg i
  · intro h g hg
    apply h g
    intro i
    exact hg (b i) ⟨i, rfl⟩

private theorem baseSize_eq_of_exactBaseSize
    (P : FinitePermutationGroup) {n : Nat}
    (h : Saxl.ExactBaseSize P.G P.Point n) : P.baseSize = n := by
  let sizes : Set Nat :=
    {m : Nat | ∃ b : Fin m → P.Point,
      Function.Injective b ∧ P.IsBase (Set.range b)}
  have hn : n ∈ sizes := by
    obtain ⟨b, hinj, hbase⟩ := h.1
    exact ⟨b, hinj, (isBase_range_iff_isBaseTuple P b).2 hbase⟩
  change sInf sizes = n
  apply le_antisymm (Nat.sInf_le hn)
  by_contra hle
  have hlt : sInf sizes < n := Nat.lt_of_not_ge hle
  obtain ⟨b, hinj, hbase⟩ := Nat.sInf_mem (⟨n, hn⟩ : sizes.Nonempty)
  exact h.2 (sInf sizes) hlt
    ⟨b, hinj, (isBase_range_iff_isBaseTuple P b).1 hbase⟩

private theorem faithfulSMul_of_exactBaseSize
    {G Ω : Type} [Group G] [MulAction G Ω] {n : Nat}
    (h : Saxl.ExactBaseSize G Ω n) : FaithfulSMul G Ω := by
  rw [faithfulSMul_iff]
  intro g hfix
  obtain ⟨b, _hinj, hbase⟩ := h.1
  exact hbase g (fun i ↦ hfix (b i))

private theorem saxlGraph_adj_iff_generalizedAdjacent
    (P : FinitePermutationGroup) (tail : Nat)
    (hsize : P.baseSize = tail + 2) (x y : P.Point) :
    P.saxlGraph.Adj x y ↔
      Saxl.GeneralizedAdjacent P.G P.Point tail x y := by
  change
    (x ≠ y ∧ ∃ b : Fin P.baseSize → P.Point,
      Function.Injective b ∧ P.IsBase (Set.range b) ∧
        x ∈ Set.range b ∧ y ∈ Set.range b) ↔ _
  rw [hsize, Saxl.generalizedAdjacent_iff_mem_base]
  constructor
  · rintro ⟨hxy, b, hinj, hbase, hx, hy⟩
    exact ⟨hxy, b, ⟨hinj, (isBase_range_iff_isBaseTuple P b).1 hbase⟩,
      hx, hy⟩
  · rintro ⟨hxy, b, ⟨hinj, hbase⟩, hx, hy⟩
    exact ⟨hxy, b, hinj, (isBase_range_iff_isBaseTuple P b).2 hbase,
      hx, hy⟩

private theorem pointDegrees_unbounded (tail degreeBound : Nat) :
    ∃ (d : Nat) (hd : Odd d) (hd3 : 3 ≤ d),
      degreeBound <
        Nat.card (EveryBase.HqProductModule d hd hd3 tail) := by
  let d := 2 * degreeBound + 3
  have hd : Odd d := by
    refine ⟨degreeBound + 1, ?_⟩
    simp [d]
    omega
  have hd3 : 3 ≤ d := by simp [d]
  have hbound : degreeBound ≤ 3 ^ d - 1 := by
    have hp : 1 + d * (3 - 1) ≤ 3 ^ d := by
      simpa using
        (one_add_le_pow_of_two_add_nonneg
          (R := Nat) (a := 2) (by omega) d)
    omega
  have hpow : 2 ^ degreeBound ≤ 2 ^ (3 ^ d - 1) :=
    Nat.pow_le_pow_right (by omega) hbound
  refine ⟨d, hd, hd3, degreeBound.lt_two_pow_self.trans_le (hpow.trans ?_)⟩
  rw [← EveryBase.vq_card d hd]
  exact EveryBase.vq_card_le_productModule_card d hd hd3 tail

private theorem hqBadVector_ne_zero
    (d : Nat) (hd : Odd d) (hd3 : 3 ≤ d) (tail : Nat) :
    EveryBase.hqBadVector d hd hd3 tail ≠ 0 := by
  classical
  let c0 : EveryBase.HqColourType d hd hd3 tail :=
    EveryBase.hqFirstColour d hd hd3 tail
  have hbad_at :
      EveryBase.hqBadVector d hd hd3 tail c0 =
        EveryBase.hqBadSeedVector d hd := by
    change
      EveryBase.badVector (EveryBase.hqSeed d hd hd3)
          (EveryBase.hqColours d hd hd3 tail) c0 c0 =
        EveryBase.hqBadSeedVector d hd
    rw [EveryBase.badVector_at]
    rfl
  have hseed_ne_zero : EveryBase.hqBadSeedVector d hd ≠ 0 := by
    change (EveryBase.hqSeed d hd hd3).u ≠ 0
    simpa using
      (EveryBase.u_ne_add_smul (EveryBase.hqSeed d hd hd3)
        (1 : (EveryBase.hqSeed d hd hd3).H)
        (0 : (EveryBase.hqSeed d hd hd3).V))
  intro hbad
  apply hseed_ne_zero
  have hat : EveryBase.hqBadVector d hd hd3 tail c0 =
      (0 : EveryBase.Vq d) := by
    simpa using congrFun hbad c0
  rw [hbad_at] at hat
  exact hat

private theorem familyFacts
    (d : Nat) (hd : Odd d) (hd3 : 3 ≤ d) (tail : Nat) :
    ∃ P : FinitePermutationGroup,
      Nat.card P.Point =
          Nat.card (EveryBase.HqProductModule d hd hd3 tail) ∧
        MulAction.IsPreprimitive P.G P.Point ∧
        P.baseSize = tail + 2 ∧
        ∃ x y : P.Point, x ≠ y ∧
          (tail = 0 → ¬ P.saxlGraph.Adj x y) ∧
          P.saxlGraph.commonNeighbors x y = ∅ := by
  have hcore := EveryBase.hq_coreFamily d hd hd3 tail
  let _ : FaithfulSMul
      (EveryBase.GBd d hd hd3 tail)
      (EveryBase.HqProductModule d hd hd3 tail) :=
    faithfulSMul_of_exactBaseSize hcore.1
  let _ : Finite (EveryBase.GBd d hd hd3 tail) :=
    Finite.of_injective
      (MulAction.toPerm : EveryBase.GBd d hd hd3 tail →
        Equiv.Perm (EveryBase.HqProductModule d hd hd3 tail))
      MulAction.toPerm_injective
  let P : FinitePermutationGroup := {
    G := EveryBase.GBd d hd hd3 tail
    Point := EveryBase.HqProductModule d hd hd3 tail }
  have hsize : P.baseSize = tail + 2 :=
    baseSize_eq_of_exactBaseSize P (by simpa [P] using hcore.1)
  have hcommon :
      P.saxlGraph.commonNeighbors
          0 (EveryBase.hqBadVector d hd hd3 tail) = ∅ := by
    apply Set.eq_empty_iff_forall_notMem.2
    intro z hz
    rw [SimpleGraph.mem_commonNeighbors] at hz
    apply hcore.2.2
    refine ⟨z, ?_, ?_⟩
    · exact (saxlGraph_adj_iff_generalizedAdjacent
        P tail hsize 0 z).1 hz.1
    · exact (saxlGraph_adj_iff_generalizedAdjacent
        P tail hsize z (EveryBase.hqBadVector d hd hd3 tail)).1
          hz.2.symm
  have hnonadj : tail = 0 →
      ¬ P.saxlGraph.Adj 0 (EveryBase.hqBadVector d hd hd3 tail) := by
    intro htail hadj
    subst tail
    apply EveryBase.hq_zero_not_generalizedAdjacent_badVector d hd hd3
    exact (saxlGraph_adj_iff_generalizedAdjacent
      P 0 hsize 0 (EveryBase.hqBadVector d hd hd3 0)).1 hadj
  refine ⟨P, rfl, ?_, hsize,
    0, EveryBase.hqBadVector d hd hd3 tail, ?_, hnonadj, hcommon⟩
  · simpa [P] using hcore.2.1
  · simpa [P] using (hqBadVector_ne_zero d hd hd3 tail).symm

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
          P.saxlGraph.commonNeighbors x y = ∅ := by
  obtain ⟨d, hd, hd3, hdegree⟩ := pointDegrees_unbounded 0 degreeBound
  obtain ⟨P, hcard, hprimitive, hsize, x, y, hxy, hnonadj, hcommon⟩ :=
    familyFacts d hd hd3 0
  refine ⟨P, ?_, hprimitive, ?_, x, y, hxy, hnonadj rfl, hcommon⟩
  · simpa [hcard] using hdegree.le
  · simpa using hsize

/-- **Theorem 1.3.** For every `B ≥ 3` and every prescribed bound there is a
primitive permutation group of degree at least that bound and base size `B`
whose Saxl graph has two distinct vertices with no common neighbour. -/
theorem theorem1_3 (B degreeBound : Nat) (hB : 3 ≤ B) :
    ∃ P : FinitePermutationGroup,
      degreeBound ≤ Nat.card P.Point ∧
        MulAction.IsPreprimitive P.G P.Point ∧
        P.baseSize = B ∧
        ∃ x y : P.Point, x ≠ y ∧
          P.saxlGraph.commonNeighbors x y = ∅ := by
  obtain ⟨d, hd, hd3, hdegree⟩ :=
    pointDegrees_unbounded (B - 2) degreeBound
  obtain ⟨P, hcard, hprimitive, hsize, x, y, hxy, _hnonadj, hcommon⟩ :=
    familyFacts d hd hd3 (B - 2)
  have htail : B - 2 + 2 = B := by omega
  refine ⟨P, ?_, hprimitive, hsize.trans htail, x, y, hxy, hcommon⟩
  simpa [hcard] using hdegree.le

end SaxlCounterexamples.MainTheorems.Internal
