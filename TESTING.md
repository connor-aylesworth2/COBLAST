# External Testing Guide

These steps are intended for people testing the local BLAST+ Flask prototype on
their own machine.

## 1. Install prerequisites

Install:

- Python 3.11 or newer
- Git
- NCBI BLAST+ command-line tools

On Windows, download BLAST+ from NCBI and extract it somewhere stable, for
example:

```text
C:\Tools\ncbi-blast-2.17.0+\bin
```

The BLAST+ `bin` folder should contain `blastn`, `blastp`, `blastx`, `tblastn`,
`makeblastdb`, and `blastdbcmd`.

## 2. Get the prototype

```powershell
git clone https://github.com/connor-aylesworth2/blast-flask-app.git
cd blast-flask-app
```

If you are testing a feature branch before it is merged:

```powershell
git checkout codex/localhost-remote-safety
```

## 3. Run the one-step launcher

If you were given the standalone Windows executable, download and run:

```powershell
.\COBLAST.exe
```

The standalone executable bundles the prototype interface, Python dependencies,
toy sample data, and required BLAST+ executables. It stores runtime files beside
the executable in `COBLAST_data`.

If you are testing from a source-code checkout instead, use the Python launcher.

Recommended setup and launch:

```powershell
python run_COBLAST.py
```

On Windows, if `python` does not launch Python 3.11 or newer, try:

```powershell
py -3.11 run_COBLAST.py
```

The launcher checks BLAST+, creates `.venv` if needed, installs
`requirements.txt`, runs `smoke_test.py`, starts Flask on `127.0.0.1`, and opens:

```text
http://127.0.0.1:5000
```

If another local process is already using port `5000`, the launcher or
standalone `.exe` will print a different `127.0.0.1` address. Use the address
shown in the newest launcher window.

If BLAST+ is installed somewhere the launcher cannot find, pass the BLAST+ `bin`
directory explicitly:

```powershell
python run_COBLAST.py --blast-bin 'C:\Tools\ncbi-blast-2.17.0+\bin'
```

Useful launcher options:

```powershell
python run_COBLAST.py --check-only
python run_COBLAST.py --skip-smoke
python run_COBLAST.py --no-browser
python run_COBLAST.py --port 5050
```

## 4. Confirm the app is local

Open:

```text
http://127.0.0.1:5000
```

If the launcher printed a different port, open that address instead.

The Flask server should be reachable from the same machine only. Do not launch
the app with a public host binding for external testing.

Remote BLAST is disabled. Query data and BLAST databases should remain local
during prototype testing.

## 5. Manual troubleshooting

Use this section only if the one-step launcher fails.

Check Python:

```powershell
python --version
py -0p
```

If needed, install Python 3.11 or newer and rerun:

```powershell
py -3.11 run_COBLAST.py
```

Check BLAST+ manually:

```powershell
$env:BLAST_BIN = 'C:\Tools\ncbi-blast-2.17.0+\bin'
& "$env:BLAST_BIN\blastn.exe" -version
```

If the browser reports that BLAST+ is missing from `C:\Program Files\NCBI` after
running the standalone `.exe`, close any older COBLAST/Flask windows and reopen
the newest address printed by the `.exe`. That message usually comes from an
older local server still occupying port `5000`.

Create the virtual environment manually:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Run the backend smoke test manually:

```powershell
python smoke_test.py
```

Expected result: the smoke test creates toy nucleotide and protein databases,
runs the supported BLAST programs, parses structured results, and exits without
an error.

Start the app manually:

```powershell
python app.py
```

## 6. Suggested tester workflow

1. Paste or upload a small FASTA query.
2. Choose a compatible BLAST program.
3. Choose a compatible registered database.
4. Try the standard, sensitive, and fast presets.
5. Run BLAST.
6. Confirm the results page includes the BLAST command, return code, runtime,
   stdout, stderr, and a structured results table.
7. Download results as CSV and TSV.
8. Visit `/databases` and confirm registered databases can be checked, added,
   created from FASTA, and removed from the registry without deleting files.

## 7. Reporting feedback

Please include:

- operating system and version
- Python version
- BLAST+ version
- exact command or web action that failed
- traceback, browser error, or terminal output
- whether the problem used toy data or your own local data

Do not include patient-identifiable, confidential, or unpublished biological
sequence data in GitHub issues or screenshots.
