"""Fast, self-contained repository hygiene and core utility unit tests.

Run with either:
    python3 -m unittest discover -s analysis_novel/tests -p "test_unit_*.py" -v
or:
    pytest analysis_novel/tests/test_unit_conifer.py
"""

from __future__ import annotations

import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


class TestRepositoryReproducibility(unittest.TestCase):
    def test_required_metadata_files_present(self) -> None:
        for name in ("README.md", "LICENSE", "CITATION.cff", "SETUP.md", ".gitignore"):
            self.assertTrue((REPO_ROOT / name).is_file(), f"Missing top-level {name}")

    def test_portable_example_config_has_no_personal_paths(self) -> None:
        example = REPO_ROOT / "configs" / "paths.example.yaml"
        self.assertTrue(example.is_file(), "Missing configs/paths.example.yaml")
        text = example.read_text(encoding="utf-8")
        self.assertNotIn("fitzpatria", text)
        self.assertNotIn("/scratch/", text)

    def test_shell_scripts_have_no_hardcoded_personal_email_or_home(self) -> None:
        for sh in REPO_ROOT.rglob("*.sh"):
            text = sh.read_text(encoding="utf-8", errors="ignore")
            self.assertNotIn("amy.fitzpatrick@ucd.ie", text, f"Hardcoded email in {sh}")
            self.assertNotIn("/home/people/fitzpatria", text, f"Hardcoded home path in {sh}")
            self.assertNotIn("/scratch/fitzpatria/", text, f"Hardcoded scratch path in {sh}")

    def test_no_dangling_references_symlink_tracked(self) -> None:
        ref_path = REPO_ROOT / "references"
        self.assertFalse(ref_path.is_symlink(), "Dangling personal symlink 'references' should not exist")

    def test_all_python_scripts_compile_cleanly(self) -> None:
        import py_compile
        py_files = [
            p for p in REPO_ROOT.rglob("*.py")
            if ".git" not in p.parts and "__pycache__" not in p.parts
        ]
        self.assertGreater(len(py_files), 10)
        for p in py_files:
            py_compile.compile(str(p), doraise=True)


if __name__ == "__main__":
    unittest.main()
