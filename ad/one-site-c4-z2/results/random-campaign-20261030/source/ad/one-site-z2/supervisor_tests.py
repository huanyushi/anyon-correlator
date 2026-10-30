import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest


class BudgetTests(unittest.TestCase):
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
