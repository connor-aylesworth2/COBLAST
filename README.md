# Local Flask + BLAST+ Prototype

This is a minimal local Flask wrapper around the installed NCBI BLAST+ executables.
It validates FASTA input with Biopython, runs allowlisted BLAST+ programs through
`subprocess`, and parses BLAST tabular/XML output with `Bio.SearchIO`.

## 1. Verify BLAST+

```powershell
& '..\ncbi-blast-2.17.0+\bin\blastn.exe' -version
```

Expected output includes `blastn: 2.17.0+`.

## 2. Create and activate a virtual environment

Install Python first if `python --version` does not work in your own PowerShell.

```powershell
cd 'C:\Users\cjohn\OneDrive\Desktop\School Shit\Edinburgh Stuff\eToL-V Dissertation\blast_flask_app'
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

`requirements.txt` installs Flask plus Biopython, which is used through
`Bio.SearchIO` to parse BLAST tabular and XML output.

## 3. Run the backend smoke test

```powershell
python smoke_test.py
```

This creates tiny toy nucleotide and protein databases with `makeblastdb`, then
runs `blastn`, `blastp`, `blastx`, and `tblastn` through the shared backend.
The toy databases are written under your temp directory, usually:

```text
C:\Users\Connor\AppData\Local\Temp\blast_flask_demo\db\toy_nt
C:\Users\Connor\AppData\Local\Temp\blast_flask_demo\db\toy_protein
```

## 4. Run the Flask app

```powershell
python app.py
```

Open:

```text
http://127.0.0.1:5000
```

The default database field points at the toy nucleotide database created by
`smoke_test.py`. For `blastp` or `blastx`, use the toy protein database prefix
instead.

## Current supported BLAST programs

- `blastn`: nucleotide query vs nucleotide database
- `blastp`: protein query vs protein database
- `blastx`: translated nucleotide query vs protein database
- `tblastn`: protein query vs translated nucleotide database

The backend validates query sequence type before running BLAST. Nucleotide
queries accept IUPAC nucleotide ambiguity codes and convert `U` to `T`; protein
queries accept common IUPAC amino acid codes plus `*`.

## Useful note

If BLAST+ is installed somewhere else, set `BLAST_BIN` before running:

```powershell
$env:BLAST_BIN = 'C:\Users\cjohn\OneDrive\Desktop\School Shit\Edinburgh Stuff\eToL-V Dissertation\ncbi-blast-2.17.0+\bin'
```
