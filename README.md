# Local Flask + BLAST+ Prototype

This is a minimal local Flask wrapper around the installed NCBI BLAST+ executables.
It validates FASTA input with Biopython, runs allowlisted BLAST+ programs through
`subprocess`, and parses BLAST tabular/XML output with `Bio.SearchIO`.
The browser interface uses an NCBI-inspired colour palette and keeps routine
clinician-facing controls separate from advanced BLAST parameters.

## 1. Verify BLAST+

```powershell
& '..\ncbi-blast-2.17.0+\bin\blastn.exe' -version
```

Expected output includes `blastn: 2.17.0+`.

## 2. Create and activate a virtual environment

Install Python first if `python --version` does not work in your own PowerShell.

```powershell
cd 'C:\Projects\blast_flask_app'
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
The smoke-test databases are written under your temp directory, usually:

```text
C:\Users\<your-username>\AppData\Local\Temp\blast_flask_demo\db\toy_nt
C:\Users\<your-username>\AppData\Local\Temp\blast_flask_demo\db\toy_protein
```

## 4. Run the Flask app

```powershell
python app.py
```

Open:

```text
http://127.0.0.1:5000
```

The bundled Flask development server is explicitly bound to `127.0.0.1`, so it
listens only on the local machine when launched with `python app.py`. Do not run
the app with `--host 0.0.0.0` unless the project has been deliberately reviewed
for network exposure.

Remote BLAST is disabled. The backend does not expose NCBI BLAST+'s `-remote`
option and rejects any generated BLAST command that contains it, keeping query
sequences and databases local by default.

For tester-facing installation and validation steps, see `TESTING.md`. For
privacy, security, and clinical-use boundaries, see `PRIVACY_SECURITY.md`.

The app also seeds a local SQLite registry and managed toy BLAST databases under:

```text
instance\database_registry.sqlite
instance\databases\
```

The main page lets the user choose a BLAST program, choose a compatible
registered database from a dropdown, review the database description, and run
BLAST. The raw BLAST database prefix remains visible in the advanced database
details panel.

The routine search form includes:

- query paste/upload
- search type
- compatible database
- sensitivity preset
- run button

Advanced settings are collapsed by default and include E-value, maximum target
sequences, word size, BLASTN task, minimum percent identity, output parser
format, timeout, and raw database details.

The results table reports the top BLAST hits with query ID, subject ID, subject
title, percent identity, alignment length, query coverage, E-value, and bit
score.
Each completed run is saved locally under `instance\results\` so the displayed
result table can be downloaded as CSV or TSV without rerunning BLAST.

The database-management page is available at:

```text
http://127.0.0.1:5000/databases
```

From that page, users can view registered databases, check availability, see
database type, add an existing BLAST database, create a new BLAST database from
FASTA with `makeblastdb`, and remove a database from the registry without
deleting BLAST files.

## Current supported BLAST programs

- `blastn`: nucleotide query vs nucleotide database
- `blastp`: protein query vs protein database
- `blastx`: translated nucleotide query vs protein database
- `tblastn`: protein query vs translated nucleotide database

The interface enforces these compatibility rules by filtering registered
databases by `nucl` or `prot` type after the BLAST program is selected.

The backend validates query sequence type before running BLAST. Nucleotide
queries accept IUPAC nucleotide ambiguity codes and convert `U` to `T`; protein
queries accept common IUPAC amino acid codes plus `*`.

On Windows, BLAST-facing database paths may be stored using short path segments
such as `SCHOOL~1`. This avoids BLAST+ parsing issues with spaces in local
folder names while preserving the source FASTA path in the registry.

## Useful note

If BLAST+ is installed somewhere else, set `BLAST_BIN` before running:

```powershell
$env:BLAST_BIN = 'C:\Tools\ncbi-blast-2.17.0+\bin'
```

## License

This repository currently includes a placeholder license notice in `LICENSE`.
Choose a formal license before wider distribution.
