---
name: pyinstaller-build-needs-biopython
description: "PyInstaller exe must be built from an interpreter that has requirements.txt installed, or Bio is silently omitted and the app crashes on launch"
metadata: 
  node_type: memory
  type: project
  originSessionId: 8d20ec62-3d92-4b2f-827c-41ca1ed79ffd
---

Before rebuilding the COBLAST+ standalone with `build_standalone_exe.py`, the
**build interpreter must have `requirements.txt` installed** (Flask + biopython,
and biopython pulls numpy). The script runs `sys.executable -m PyInstaller`, so
PyInstaller analyzes imports against whatever interpreter launched it.

**Why:** On 2026-07-01 the global Python312 (where PyInstaller was installed) had
Flask but not biopython, and no `.venv` existed. `--collect-submodules Bio.SeqIO`
found nothing, so the exe shipped **without `Bio`**. It launched, ran the smoke
test, then died with `ModuleNotFoundError: No module named 'Bio'` (raised via
`require_biopython()` in blast_runner.py) — a console double-click just flashed
and closed. Tell: broken build was ~54 MB vs a healthy ~63–66 MB.

**How to apply:** (1) `python -m pip install -r requirements.txt` into the build
interpreter first; verify `python -c "from Bio import SeqIO"`. (2) Build with
`--blast-bin "C:/Program Files/NCBI/blast-2.17.0+/bin"` (the script's default
`../ncbi-blast-2.17.0+/bin` does not exist on this machine). (3) Verify the exe
with `./dist/COBLAST.exe --check-only --no-browser` — this runs the smoke test
(keep smoke ON; `--skip-smoke` hides exactly this class of bug). Relates to
[[coblast-verification-anchor]].
