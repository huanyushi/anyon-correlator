import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest


class BudgetTests(unittest.TestCase):
    def test_campaign_timeout_preserves_nested_recovery(self):
        spec = importlib.util.spec_from_file_location("supervise", Path(__file__).with_name("supervise.py"))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "campaign"
            progress = {"attempts": [{"seed": 19}], "best": {"energy": -0.7},
                        "active_directory": str(out / "attempt02-seed23")}
            code = "import pathlib,json,time; p=pathlib.Path(%r); q=p/'attempt02-seed23'; q.mkdir(parents=True); (q/'step0002.json').write_text('{}'); (p/'progress.json').write_text(%r); time.sleep(30)" % (str(out), json.dumps(progress))
            result = module.supervise([sys.executable, "-c", code], out, budget=1)
            self.assertEqual(result["last_checkpoint"], "attempt02-seed23/step0002.json")
            self.assertEqual(result["attempts"], progress["attempts"])
            self.assertEqual(result["best"], progress["best"])
            self.assertFalse(result["accepted"])

    def test_timeout_during_evaluation_and_audit(self):
        spec = importlib.util.spec_from_file_location("supervise", Path(__file__).with_name("supervise.py"))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        for phase in ("evaluation", "audit"):
            with self.subTest(phase=phase), tempfile.TemporaryDirectory() as tmp:
                out = Path(tmp) / "run"
                # A real child process writes a checkpoint, then blocks inside a computation.
                code = "import pathlib,json,time; p=pathlib.Path(%r); p.mkdir(); (p/'step0001.json').write_text(json.dumps({'metadata':{'accepted':False}})); (p/'phase').write_text(%r); time.sleep(30)" % (str(out), phase)
                result = module.supervise([sys.executable, "-c", code], out, budget=1)
                self.assertEqual(result["status"], "budget")
                self.assertFalse(result["accepted"])
                self.assertEqual(result["last_checkpoint"], "step0001.json")
                self.assertEqual(json.loads((out / "summary.json").read_text())["status"], "budget")
                self.assertEqual((out / "phase").read_text(), phase)


if __name__ == "__main__":
    unittest.main()
