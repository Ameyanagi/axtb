# SPDX-License-Identifier: LGPL-3.0-or-later
"""Small synthetic-output tests; an installed official xtb is not required."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import xtb_opt


def scf_table(iterations: int) -> str:
    return (
        " iter      E             dE          RMSdq      gap      omega  full diag\n"
        + "\n".join(
            f"{i:4d}     -1.0000000 -0.100000E-01 0.100E-03    7.20       0.0  T" for i in range(1, iterations + 1)
        )
        + "\n\n SCC iter. ... 0 min, 0.001 sec\n"
    )


OUTPUT = (
    scf_table(3)
    + ".............................. CYCLE    1 ..............................\n"
    + scf_table(2)
    + " * total energy  : -1.1000000 Eh     change -1.0000000E-01 Eh\n"
    + "   gradient norm : 0.0001000 Eh/α   predicted -1.0000000E-01\n"
    + ".............................. CYCLE    2 ..............................\n"
    + scf_table(1)
    + " * total energy  : -1.1200000 Eh     change -2.0000000E-02 Eh\n"
    + "   gradient norm : 0.0000800 Eh/α   predicted -2.0000000E-02\n"
    + " *** GEOMETRY OPTIMIZATION CONVERGED AFTER 7 ITERATIONS ***\n"
    + " ANCopt (total)                 0 d,  0 h,  1 min,  2.000 sec\n"
    + "         single point calculation ... 1 min, 0.250 sec (97.177%)\n"
    + "                    model hessian ... 0 min, 0.250 sec (0.403%)\n"
    + " Final Singlepoint\n"
    + scf_table(2)
    + " | TOTAL ENERGY               -1.123456789012D+00 Eh   |\n"
    + " | GRADIENT NORM               0.000123456789 Eh/α |\n"
    + "normal termination of xtb\n"
)


class ParserTests(unittest.TestCase):
    def test_separate_cycle_count_reported_iterations_and_scf_tables(self):
        result = xtb_opt.parse_output(OUTPUT, potential="gfn2")
        self.assertEqual(result["status"], "converged")
        self.assertEqual(result["steps"], 7)
        self.assertEqual(result["reported_iterations"], 7)
        self.assertEqual(result["cycle_blocks"], 2)
        self.assertEqual(result["evaluations"], 2)
        self.assertEqual(result["scf_iterations_total"], 8)
        self.assertEqual(result["potential_seconds"], 60.25)
        self.assertEqual(result["model_hessian_seconds"], 0.25)
        self.assertEqual(result["optimizer_seconds"], 1.75)
        self.assertEqual(result["energy"], -1.123456789012)
        self.assertEqual(result["gradient_norm"], 0.000123456789)
        self.assertIsNone(result["max_gradient"])
        self.assertIsNone(result["hessian_evaluations"])
        self.assertIsNone(result["hessian_seconds"])
        self.assertEqual(
            result["trace"],
            [[1, 0, -1.1, None, "normal", 0.0001], [2, 0, -1.12, None, "normal", 0.00008]],
        )

    def test_failure_and_missing_measurements(self):
        text = ".............................. CYCLE    1 ..............................\n"
        text += "*** FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN 1 ITERATIONS ***\nnormal termination of xtb\n"
        result = xtb_opt.parse_output(text, potential="gfnff")
        self.assertEqual(result["status"], "max_steps")
        self.assertEqual(result["steps"], 1)
        self.assertEqual(result["evaluations"], 1)
        self.assertEqual(result["trace"], [])
        self.assertIsNone(result["energy"])
        self.assertIsNone(result["scf_iterations_total"])
        self.assertEqual(xtb_opt.parse_output(text, potential="gfnff", returncode=1)["status"], "max_steps")
        self.assertEqual(xtb_opt.parse_output(OUTPUT, potential="gfn2", returncode=1)["status"], "error")
        self.assertEqual(
            xtb_opt.parse_output("SCF not converged, aborting...", potential="gfn2", returncode=1)["status"],
            "scf_failed",
        )

    def test_partial_scf_output_does_not_claim_total(self):
        result = xtb_opt.parse_output(OUTPUT.replace(scf_table(1), ""), potential="gfn2")
        self.assertEqual(result["evaluations"], 2)
        self.assertIsNone(result["scf_iterations_total"])

    def test_gfnff_counts_only_cycles(self):
        text = OUTPUT.replace("AFTER 7 ITERATIONS", "AFTER 2 ITERATIONS")
        for iterations in (1, 2, 3):
            text = text.replace(scf_table(iterations), "")
        result = xtb_opt.parse_output(text, potential="gfnff")
        self.assertEqual(result["evaluations"], 2)
        self.assertEqual(result["cycle_blocks"], 2)
        self.assertEqual(result["reported_iterations"], 2)
        self.assertIsNone(result["scf_iterations_total"])

    def test_real_xtb_cycle_summary_format(self):
        text = "".join(
            f".............................. CYCLE    {cycle} ..............................\n"
            " * total energy  :   -29.5846846 Eh     change       -0.1307968E-05 Eh\n"
            "   gradient norm :     0.0000538 Eh/α   predicted    -0.0000000E+00\n"
            for cycle in range(1, 10)
        )
        text += "*** GEOMETRY OPTIMIZATION CONVERGED AFTER 9 ITERATIONS ***\nnormal termination of xtb\n"
        result = xtb_opt.parse_output(text, potential="gfn2")
        self.assertEqual(result["evaluations"], 9)
        self.assertEqual(result["steps"], 9)
        self.assertEqual(result["reported_iterations"], 9)
        self.assertEqual(
            result["trace"],
            [[cycle, 0, -29.5846846, None, "normal", 0.0000538] for cycle in range(1, 10)],
        )

    def test_cycle_trace_excludes_initial_and_final_reports(self):
        initial = " * total energy : -1.0 Eh\n   gradient norm : 1.0 Eh/α\n"
        text = initial
        text += ".............................. CYCLE    3 ..............................\n"
        text += " * total energy : -2.0D+00 Eh\n   gradient norm : 2.0D-03 Eh/α\n"
        text += ".............................. CYCLE    8 ..............................\n"
        text += " * total energy : -3.0 Eh\n"
        final = " * total energy : -4.0 Eh\n   gradient norm : 4.0 Eh/α\n"
        for boundary in (
            "*** GEOMETRY OPTIMIZATION CONVERGED AFTER 8 ITERATIONS ***",
            "*** FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN 8 ITERATIONS ***",
            " ANCopt (total) 0 d, 0 h, 0 min, 0.001 sec",
            " Final Singlepoint",
        ):
            with self.subTest(boundary=boundary):
                result = xtb_opt.parse_output(text + boundary + "\n" + final, potential="gfnff")
                self.assertEqual(result["evaluations"], 2)
                self.assertEqual(result["cycle_blocks"], 2)
                self.assertEqual(result["trace"], [[3, 0, -2.0, None, "normal", 0.002]])

    def test_no_cycles_has_zero_evaluations_and_no_trace(self):
        result = xtb_opt.parse_output(" * total energy : -1.0 Eh\n   gradient norm : 1.0 Eh/α\n", potential="gfn2")
        self.assertEqual(result["evaluations"], 0)
        self.assertEqual(result["cycle_blocks"], 0)
        self.assertEqual(result["trace"], [])
        self.assertIsNone(result["steps"])
        self.assertIsNone(result["reported_iterations"])


class RunnerTests(unittest.TestCase):
    def test_fake_executable_isolated_jobs_and_geometry_order(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            starts = root / "starts" / "set1"
            starts.mkdir(parents=True)
            source = "2\ncharge=-1 uhf=1\nH 0 0 0\nO 0 0 1\n"
            (starts / "valid.given.xyz").write_text(source)
            (starts / "bad.given.xyz").write_text(source)
            binary = root / "fake-xtb"
            binary.write_text(
                f"#!{sys.executable}\n"
                "import os, pathlib, sys\n"
                "assert sys.argv[1:] == ['input.xyz', '--gfn', '2', '--opt', 'normal', '--chrg', '-1', '--uhf', '1']\n"
                "assert os.environ['OMP_NUM_THREADS'] == '1' and os.environ['OMP_STACKSIZE'] == '1G'\n"
                "assert os.environ['XTBHOME'] == str(pathlib.Path.cwd()) and 'XTBPATH' not in os.environ\n"
                "assert sorted(p.name for p in pathlib.Path('.').iterdir()) == ['input.xyz']\n"
                "assert not pathlib.Path('shared-marker').exists()\n"
                "pathlib.Path('shared-marker').touch()\n"
                "pathlib.Path('xtbopt.xyz').write_text(pathlib.Path('input.xyz').read_text())\n"
                f"print({OUTPUT!r})\n"
            )
            binary.chmod(0o755)
            output = root / "results" / "out.jsonl"
            self.assertEqual(
                xtb_opt.main(
                    [
                        "--starts",
                        str(root / "starts"),
                        "--potential",
                        "gfn2",
                        "--xtb",
                        str(binary),
                        "--jobs",
                        "2",
                        "--out",
                        str(output),
                    ]
                ),
                0,
            )
            rows = [json.loads(line) for line in output.read_text().splitlines()]
            self.assertEqual([row["case"] for row in rows], ["bad", "valid"])
            for row in rows:
                self.assertEqual(row["status"], "converged")
                self.assertEqual((output.parent / row["final_xyz"]).read_text(), source)
                self.assertTrue((output.parent / row["log"]).is_file())
                self.assertGreater(row["wall_seconds"], 0)
            binary.write_text(
                binary.read_text().replace(
                    "pathlib.Path('xtbopt.xyz').write_text(pathlib.Path('input.xyz').read_text())",
                    "pathlib.Path('xtbopt.xyz').write_text('2\\n\\nO 0 0 1\\nH 0 0 0\\n')",
                )
            )
            result = xtb_opt.run_start(
                xtb_opt.discover_starts(root / "starts", cases=["valid"])[0],
                xtb=str(binary),
                potential="gfn2",
                level="normal",
                output_parent=output.parent,
                timeout=10,
            )
            self.assertEqual(result["status"], "error")
            self.assertIsNone(result["final_xyz"])
            self.assertIn("atom identities/order", result["error"])


if __name__ == "__main__":
    unittest.main()
