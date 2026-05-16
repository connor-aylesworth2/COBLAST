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

## 2. Get the prototype

```powershell
git clone https://github.com/connor-aylesworth2/blast-flask-app.git
cd blast-flask-app
```

If you are testing a feature branch before it is merged:

```powershell
git checkout codex/localhost-remote-safety
```

## 3. Create a virtual environment

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## 4. Point the app at BLAST+

Set `BLAST_BIN` to the folder containing `blastn`, `blastp`, `blastx`,
`tblastn`, `makeblastdb`, and `blastdbcmd`.

Windows PowerShell example:

```powershell
$env:BLAST_BIN = 'C:\Tools\ncbi-blast-2.17.0+\bin'
```

Confirm BLAST+ is reachable:

```powershell
& "$env:BLAST_BIN\blastn.exe" -version
```

## 5. Run the backend smoke test

```powershell
python smoke_test.py
```

Expected result: the smoke test creates toy nucleotide and protein databases,
runs the supported BLAST programs, parses structured results, and exits without
an error.

## 6. Run the web interface

```powershell
python app.py
```

Open:

```text
http://127.0.0.1:5000
```

The Flask server should be reachable from the same machine only. Do not launch
the app with a public host binding for external testing.

## 7. Suggested tester workflow

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

## 8. Reporting feedback

Please include:

- operating system and version
- Python version
- BLAST+ version
- exact command or web action that failed
- traceback, browser error, or terminal output
- whether the problem used toy data or your own local data

Do not include patient-identifiable, confidential, or unpublished biological
sequence data in GitHub issues or screenshots.
