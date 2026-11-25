# SPDX-License-Identifier: LGPL-3.0-or-later
"""Synthetic, external-program-free tests for optimizer report analysis."""

from __future__ import annotations

import contextlib
import csv
import io
import json
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

try:
    from . import analyze
except ImportError:
    import analyze


def run_row(**changes):
    row = {
        "case": "water",
        "set": "small",
        "start": "distorted",
        "potential": "gfn2",
        "config": "anc",
        "atoms": 3,
        "status": "converged",
        "steps": 3,
        "evaluations": 4,
        "hessian_evaluations": 1,
        "potential_seconds": 2.0,
        "hessian_seconds": 1.0,
        "optimizer_seconds": 0.04,
        "wall_seconds": 3.04,
        "scf_iterations_total": 12,
        "energy": -5.0,
        "max_gradient": 1e-5,
        "final_xyz": "water.xyz",
        "trace": [[1, 0.1, -4.9, 0.01, "normal"], [4, 3.0, -5.0, 1e-5, "tight"]],
    }
    row.update(changes)
    return row


class AnalysisTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        (self.directory / "water.xyz").write_text("3\nwater\nO 0 0 0\nH 0.7 0.5 0\nH -0.7 0.5 0\n")

    def jsonl(self, name, rows):
        path = self.directory / name
        path.write_text("".join(json.dumps(row) + "\n" for row in rows))
        return path

    def reference(self, **changes):
        row = {name: run_row()[name] for name in analyze.IDENTITY}
        row.update(energy=-5.0, xyz="water.xyz")
        row.update(changes)
        return row

    def load(self, *rows):
        return analyze.read_jsonl([self.jsonl("runs.jsonl", rows)])

    def test_relative_geometry_and_paired_metrics(self):
        rows = self.load(
            run_row(gradient_seconds=0.8, gradient_evaluations=4),
            run_row(config="lbfgs", evaluations=8, wall_seconds=6.08, energy=-4.9),
            run_row(config="xtb-opt-normal", evaluations=2, wall_seconds=1.52),
        )
        references = analyze.read_jsonl([self.jsonl("refs.jsonl", [self.reference()])], reference=True)
        analyze.enrich_runs(rows, references, "anc")
        anc, lbfgs, xtb = rows
        self.assertAlmostEqual(anc["optimizer_ms_per_evaluation"], 10.0)
        self.assertAlmostEqual(anc["hessian_gradient_equivalents"], 5.0)
        self.assertAlmostEqual(lbfgs["hessian_gradient_equivalents"], 4.0)
        self.assertAlmostEqual(lbfgs["delta_energy_hartree"], 0.1)
        self.assertTrue(lbfgs["same_minimum"])  # Geometry arm of the OR rule.
        self.assertEqual(lbfgs["xtb_evaluations_ratio"], 4.0)
        self.assertEqual(lbfgs["baseline_evaluations_ratio"], 2.0)
        self.assertAlmostEqual(lbfgs["xtb_wall_seconds_ratio"], 4.0)
        self.assertEqual(xtb["xtb_evaluations_ratio"], 1.0)
        self.assertEqual(analyze.geometry_path(anc, "final_xyz"), (self.directory / "water.xyz").resolve())

    def test_energy_or_rmsd_and_failure_pairing(self):
        (self.directory / "ref_heavy.xyz").write_text("3\nreference\nC 0 0 0\nN 1 0 0\nO 0 1 0\n")
        (self.directory / "shift.xyz").write_text("3\ndifferent\nC 0 0 0\nN 4 0 0\nO 0 4 0\n")
        rows = self.load(
            run_row(energy=-5 + 0.5e-6, final_xyz="shift.xyz", hessian_gradient_equivalents=7.5),
            run_row(config="different", energy=-4.9, final_xyz="shift.xyz"),
            run_row(config="failed", status="max_steps", evaluations=8),
            run_row(config="xtb-opt-normal", status="scf_failed", evaluations=None),
        )
        refs = analyze.read_jsonl([self.jsonl("refs.jsonl", [self.reference(xyz="ref_heavy.xyz")])], reference=True)
        rows[2]["final_xyz"] = None
        rows[3]["final_xyz"] = None
        analyze.enrich_runs(rows, refs, "anc")
        self.assertGreater(rows[0]["heavy_atom_rmsd_angstrom"], 0.1)
        self.assertTrue(rows[0]["same_minimum"])  # Energy arm of the OR rule.
        self.assertFalse(rows[1]["same_minimum"])
        self.assertEqual(rows[0]["hessian_gradient_equivalents"], 7.5)
        self.assertIsNone(rows[0]["xtb_evaluations_ratio"])
        self.assertIsNone(rows[2]["baseline_evaluations_ratio"])
        stats = analyze.ratio_statistics([{"ratio": 1}, {"ratio": 4}, {"ratio": None}], "ratio")
        self.assertEqual(stats, (2.5, 2.0, 2))

    def test_kabsch_rotation_translation_heavy_atoms_and_reflection(self):
        symbols = ["C", "N", "O", "S", "H"]
        a = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0], [1.0, 1.0, 1.0]])
        rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        b = a @ rotation + [3, -2, 5]
        b[-1] += 50  # Hydrogen must not change the heavy-atom RMSD.
        self.assertLess(analyze.heavy_atom_rmsd((symbols, a), (symbols, b)), 1e-12)
        b = a.copy()
        b[:, 0] *= -1
        self.assertGreater(analyze.heavy_atom_rmsd((symbols, a), (symbols, b)), 0.1)
        with self.assertRaisesRegex(ValueError, "ordered element"):
            analyze.heavy_atom_rmsd((symbols, a), (list(reversed(symbols)), a))

    def test_nulls_preserved_and_zero_overhead(self):
        row = run_row(final_xyz=None, optimizer_seconds=0, hessian_seconds=None, evaluations=4)
        rows = self.load(row)
        analyze.enrich_runs(rows, [], "missing")
        self.assertEqual(rows[0]["optimizer_ms_per_evaluation"], 0)
        self.assertIsNone(rows[0]["hessian_gradient_equivalents"])
        self.assertIsNone(rows[0]["same_minimum"])
        row = run_row()
        for name in analyze.MEASUREMENTS:
            row[name] = None
        rows = self.load(row)
        analyze.enrich_runs(rows, [], "anc")
        self.assertIsNone(rows[0]["optimizer_ms_per_evaluation"])
        self.assertIsNone(rows[0]["baseline_evaluations_ratio"])

    def test_xtb_criterion_evaluations(self):
        trace = [
            [1, 0, 0, 0.0001, "coarse", 0.0005],
            [2, 0, -5e-6, 0.0001, "coarse", 0.0005],  # Energy exactly at threshold.
            [4, 0, -6e-6, 0.0001, "normal", 0.0005],  # Accuracy changed.
            [5, 0, -7e-6, 0.0001, "normal", 0.001],  # Norm exactly at threshold.
            [7, 0, -8e-6, 0.0001, "normal"],  # Legacy sample without a norm.
            [8, 0, -9e-6, 0.0001, "normal", None],
            [9, 0, -10e-6, 0.0001, "normal", 0.0005],
            [11, 0, -11e-6, 0.0001, "normal", 0.0004],
        ]
        rows = self.load(
            run_row(config="hits", steps=7, evaluations=11, trace=trace),
            run_row(config="misses", steps=5, evaluations=8, trace=trace[:6]),
            run_row(config="legacy"),
        )
        analyze.enrich_runs(rows, [], "anc")
        self.assertEqual([row["xtb_criterion_evaluations"] for row in rows], [9, None, None])
        self.assertIsNone(analyze.xtb_criterion_evaluations(trace[:1]))

    def test_invalid_jsonl_and_schema(self):
        for change in (
            {"status": "Converged"},
            {"evaluations": -1},
            {"evaluations": True},
            {"atoms": 0},
            {"evaluations": 1.5},
            {"energy": math.inf},
            {"trace": [[2, 1, -5, 0.1, 1], [1, 2, -5, 0.01, 1]]},
            {"trace": [[1, 0, -5, None, 1]]},
            {"trace": [[1, 0, -5, 0.01]]},
            {"trace": [[1, 0, -5, None, "normal", -0.01]]},
            {"gradient_seconds": -1},
        ):
            with self.subTest(change=change), self.assertRaises(ValueError):
                self.load(run_row(**change))
        missing = run_row()
        del missing["evaluations"]
        with self.assertRaisesRegex(ValueError, "missing measurement"):
            self.load(missing)
        with self.assertRaisesRegex(ValueError, "duplicate run"):
            self.load(run_row(), run_row())
        path = self.directory / "bad.jsonl"
        for text in ('{"foo": }\n', "[]\n", '{"extra": 1e999}\n'):
            path.write_text(text)
            with self.subTest(text=text), self.assertRaises(ValueError):
                analyze.read_jsonl([path])

    def test_invalid_references_and_xyz(self):
        path = self.jsonl("refs.jsonl", [self.reference(), self.reference()])
        with self.assertRaisesRegex(ValueError, "duplicate reference"):
            analyze.read_jsonl([path], reference=True)
        with self.assertRaisesRegex(ValueError, "energy"):
            analyze.read_jsonl([self.jsonl("refs.jsonl", [self.reference(energy=None)])], reference=True)
        (self.directory / "bad.xyz").write_text("2\nshort\nH 0 0 0\n")
        with self.assertRaisesRegex(ValueError, "invalid XYZ"):
            analyze.read_xyz(self.directory / "bad.xyz")
        rows = self.load(run_row(atoms=4))
        with self.assertRaisesRegex(ValueError, "atoms disagrees"):
            analyze.enrich_runs(rows, [], "anc")

    def test_full_cli_outputs(self):
        results = self.jsonl(
            "results.jsonl",
            [
                run_row(
                    preopt={"evaluations": 2},
                    trace=[[2, 0.1, -4.999999, 0.001, "normal", 0.002], [3, 3.0, -5, 0.0001, "normal", 0.0005]],
                ),
                run_row(config="lbfgs", evaluations=8, wall_seconds=6.08),
                run_row(case="failed", config="anc", status="max_steps", energy=None, final_xyz=None, trace=[]),
            ],
        )
        xtb = self.jsonl(
            "xtb.jsonl",
            [
                run_row(
                    config="xtb-opt-normal",
                    evaluations=2,
                    wall_seconds=1.52,
                    trace=[[1, 0, -4.999999, None, "normal", 0.002], [2, 0, -5, None, "normal", 0.0005]],
                )
            ],
        )
        refs = self.jsonl("refs.jsonl", [self.reference()])
        out = self.directory / "report"
        with contextlib.redirect_stdout(io.StringIO()):
            code = analyze.main(
                [
                    str(results),
                    "--refs",
                    str(refs),
                    "--xtb",
                    str(xtb),
                    "--out",
                    str(out),
                    "--configs",
                    "anc,lbfgs",
                    "--molecules",
                    "water",
                ]
            )
        self.assertEqual(code, 0)
        with (out / "per_molecule.csv").open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 4)
        anc = next(row for row in rows if row["case"] == "water" and row["config"] == "anc")
        self.assertEqual(anc["preopt"], '{"evaluations":2}')
        self.assertEqual(float(anc["xtb_evaluations_ratio"]), 2.0)
        self.assertEqual(anc["xtb_criterion_evaluations"], "3")
        self.assertEqual(float(anc["xtb_criterion_evaluations_ratio"]), 1.5)
        self.assertEqual(next(row for row in rows if row["config"] == "lbfgs")["xtb_criterion_evaluations"], "")
        self.assertEqual(next(row for row in rows if row["case"] == "failed")["energy"], "")
        summary = (out / "summary.md").read_text()
        self.assertIn("xTB criterion evals median / mean", summary)
        self.assertIn("Criterion evals ratio vs xTB median / GM (n)", summary)
        self.assertIn("| gfn2 | small | anc | 2 | 1 | 1 | 4 / 4 | 3 / 3 |", summary)
        self.assertIn("relative to `anc`", summary)
        self.assertIn("| 2 / 2 (1) | 1.5 / 1.5 (1) |", summary)
        for name in ("evaluations.png", "convergence.png", "evaluation_cdf.png"):
            data = (out / name).read_bytes()
            self.assertEqual(data[:8], b"\x89PNG\r\n\x1a\n")
            self.assertGreater(len(data), 5000)

    def test_cli_fails_before_writing_for_malformed_data(self):
        results = self.jsonl("results.jsonl", [run_row(status="BAD")])
        out = self.directory / "report"
        with contextlib.redirect_stderr(io.StringIO()) as stderr:
            self.assertEqual(analyze.main([str(results), "--out", str(out)]), 2)
        self.assertIn("lower snake_case", stderr.getvalue())
        self.assertFalse(out.exists())

    def test_energy_gap_trace_uses_reference_and_floors_zero(self):
        row = run_row()
        row["_reference_energy"] = -5.0
        first, final = analyze.energy_gap_trace(row)
        self.assertAlmostEqual(first[2], 0.1)
        self.assertEqual(final, (4, 3.0, 1e-12))
        row["_reference_energy"] = None
        self.assertEqual(analyze.energy_gap_trace(row), [])

    def test_selection_and_empty_traces(self):
        rows = [run_row(case=f"case{index}", atoms=index + 1) for index in range(7)]
        self.assertEqual([key[0] for key in analyze.choose_molecules(rows, None)], ["case0", "case2", "case4", "case6"])
        self.assertEqual(analyze.choose_molecules([run_row(trace=[])], None), [])
        self.assertEqual(
            analyze.choose_configs([run_row(), run_row(config="xtb-opt-normal")], "anc"), ["anc", "xtb-opt-normal"]
        )
        with self.assertRaisesRegex(ValueError, "unknown configurations"):
            analyze.choose_configs(rows, "missing")
        with self.assertRaisesRegex(ValueError, "no trace data"):
            analyze.choose_molecules(rows, "missing")


if __name__ == "__main__":
    unittest.main()
