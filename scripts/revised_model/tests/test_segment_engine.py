from __future__ import annotations

import sys
import unittest
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from segment_engine import (  # noqa: E402
    ChromosomeTopology,
    GenomeState,
    LoxInstance,
    PairSampler,
    SegmentInstance,
)


def toy_state() -> GenomeState:
    parts = []
    for index in range(5):
        parts.append(
            SegmentInstance(
                template_id=f"S{index}",
                copy_id=f"S{index}@0",
                orientation=1,
                length_bp=100,
                orf_ids=(f"G{index}",),
                essential_orf_ids=(f"G{index}",) if index == 2 else (),
            )
        )
        if index < 4:
            parts.append(LoxInstance(template_id=f"L{index}", copy_id=f"L{index}@0"))
    return GenomeState({"Toy": ChromosomeTopology("Toy", parts)})


class SegmentEngineTests(unittest.TestCase):
    def test_deletion_keeps_one_boundary_and_specific_instances(self) -> None:
        state = toy_state()
        outcome = state.apply_event("deletion", "Toy", 0, 2)
        chromosome = state.chromosomes["Toy"]
        self.assertEqual([lox.template_id for lox in chromosome.lox()], ["L0", "L3"])
        self.assertEqual([seg.template_id for seg in chromosome.segments()], ["S0", "S3", "S4"])
        self.assertEqual(outcome.affected_orf_ids, ("G1", "G2"))
        self.assertEqual(state.gate_missing({"G2"}), ["G2"])

    def test_inversion_changes_order_and_orientation_without_copy_change(self) -> None:
        state = toy_state()
        before = state.orf_copy_numbers()
        state.apply_event("inversion", "Toy", 0, 3)
        chromosome = state.chromosomes["Toy"]
        self.assertEqual([seg.template_id for seg in chromosome.segments()], ["S0", "S3", "S2", "S1", "S4"])
        self.assertEqual([seg.orientation for seg in chromosome.segments()], [1, -1, -1, -1, 1])
        self.assertEqual(state.orf_copy_numbers(), before)

    def test_duplication_creates_segment_and_lox_instances(self) -> None:
        state = toy_state()
        state.apply_event("duplication", "Toy", 0, 2)
        chromosome = state.chromosomes["Toy"]
        self.assertEqual(len(chromosome.segments()), 7)
        self.assertEqual(len(chromosome.lox()), 6)
        counts = state.orf_copy_numbers()
        self.assertEqual(counts["G1"], 2)
        self.assertEqual(counts["G2"], 2)
        self.assertEqual(len({segment.copy_id for segment in chromosome.segments()}), 7)

    def test_deleting_one_duplicated_instance_does_not_zero_all_copies(self) -> None:
        state = toy_state()
        state.apply_event("duplication", "Toy", 0, 2)
        chromosome = state.chromosomes["Toy"]
        duplicated_indices = [
            index for index, lox in enumerate(chromosome.lox()) if "@e1" in lox.copy_id
        ]
        self.assertTrue(duplicated_indices)
        state.apply_event("deletion", "Toy", 2, 4)
        self.assertGreaterEqual(state.orf_copy_numbers()["G1"], 1)
        state.validate()

    def test_signatures_distinguish_inversion_from_copy_number(self) -> None:
        reference = toy_state()
        inverted = toy_state()
        inverted.apply_event("inversion", "Toy", 0, 3)
        self.assertEqual(reference.orf_copy_number_signature()[0], inverted.orf_copy_number_signature()[0])
        self.assertNotEqual(reference.structural_signature()[0], inverted.structural_signature()[0])

    def test_multistep_topology_remains_reconstructable(self) -> None:
        state = toy_state()
        state.apply_event("duplication", "Toy", 0, 2)
        state.apply_event("inversion", "Toy", 1, 4)
        state.apply_event("deletion", "Toy", 0, 2)
        state.validate()
        structural_hash, canonical = state.structural_signature()
        self.assertEqual(len(structural_hash), 32)
        self.assertIn("JUNCTIONS:", canonical)
        self.assertIn("ACTIVE_LOX:", canonical)

    def test_partial_hic_sampler_uses_additive_component_masses(self) -> None:
        direct = pd.DataFrame(
            [{"chromosome": "Toy", "lox_id_1": "L0", "lox_id_2": "L2", "direct_contact_weight": 1.0}]
        )
        sampler = PairSampler("partial_hic_fallback", 1.0, direct)
        rng = np.random.default_rng(1701)
        draws = 20_000
        direct_draws = sum(sampler.sample(toy_state(), rng)[3] == "direct_hic" for _ in range(draws))
        fallback_total = 3 / 2 + 2 / 3 + 1 / 4
        expected = 1.0 / (1.0 + fallback_total)
        self.assertAlmostEqual(direct_draws / draws, expected, delta=0.02)

        deleted = toy_state()
        deleted.apply_event("deletion", "Toy", 0, 2)
        self.assertTrue(all(sampler.sample(deleted, rng)[3] == "distance_fallback" for _ in range(200)))


if __name__ == "__main__":
    unittest.main(verbosity=2)
