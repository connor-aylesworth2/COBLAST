# COBLAST+ Prototype

## Tester Quick Start: Download the Release Executable First

Most Windows testers should start with the prebuilt executable from the latest
GitHub Release (the **Releases** section of the repository page), not the
source-code setup below. The executable is no longer stored in the repository
tree; it is attached to each Release as:

```text
COBLAST.exe
README.md
```

Download `COBLAST.exe` and its `README.md` into a local folder, read the
`README.md`, then run:

```powershell
.\COBLAST.exe
```

The standalone Windows executable bundles the local Flask interface, Python
dependencies, templates, static assets, toy sample data, and the required NCBI
BLAST+ executables. It starts a local server bound to `127.0.0.1`, opens the
browser, and stores runtime files in a stable per-user data folder at
`%LOCALAPPDATA%\COBLAST_data` (for example
`C:\Users\<you>\AppData\Local\COBLAST_data`).

Because this is an unsigned research prototype executable, Windows SmartScreen
or antivirus software may warn before first launch. The `.exe` is intended for
agreed Windows prototype testing only. It is not the correct entry point for
Linux or macOS testing; those environments should use the source-code workflow
below with a local BLAST+ installation.

If Windows reports that it cannot access the specified device, path, or file,
extract the release folder fully, move it to a simple local folder such as
`C:\COBLAST`, right-click `COBLAST.exe` > Properties > Unblock if available,
and run this diagnostic command from PowerShell:

```powershell
.\COBLAST.exe --check-only --skip-smoke --no-browser
```

If that diagnostic succeeds but the browser launch still fails, run:

```powershell
.\COBLAST.exe --skip-smoke --no-browser
```

Then open the printed `http://127.0.0.1:...` address manually.

## What COBLAST Does

COBLAST+ is a general-purpose, clinician-oriented graphical front end for NCBI
BLAST+. Its core job is ordinary BLAST: paste or upload a query, choose a
compatible local database, and run `blastn`, `blastp`, `blastx`, or `tblastn` on
the local machine. COBLAST+ wraps the stock BLAST+ executables and does not
change how BLAST searches, so a routine run returns the same hits a clinician
would get from the command-line tools. The optional preset analyses described
below (APOE genotyping and the eToL microbiome panel) are conveniences layered on
top of that core, not the main workflow.

Mechanically, COBLAST is a minimal local Flask wrapper for NCBI BLAST+. It
validates FASTA input with Biopython, runs allowlisted BLAST+ programs through
`subprocess`, and parses BLAST tabular/XML output into structured result tables.
The browser interface uses an NCBI-inspired colour palette and keeps routine
clinician-facing controls separate from advanced BLAST parameters.

For SRA-style batch exploration, COBLAST can also run one query across multiple
registered local databases. Two optional presets build on this batch path: an
APOE exact-match preset that searches four stored APOE probes and summarizes
per-sample probe counts in a visual table plus CSV/TSV exports, and eToL presets
that search an electronic Tree of Life probe panel — the full 1,017-probe
microbial set, a one-probe-per-species quick set, or the viral eToL-V panel — and
summarize the species/viruses detected per sample with per-probe and per-species
count exports.

The Windows release executable bundles the required BLAST+ executables. When
running from source, COBLAST expects compatible NCBI BLAST+ command-line tools
to be installed locally or supplied with `BLAST_BIN`.

## Source-code Quick Start

Use this section if you are developing COBLAST, testing from source, or running
on Linux/macOS.

The easiest way to set up and launch the local interface is:

```powershell
cd 'C:\path\to\COBLAST-'
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

If another local process is already using port `5000`, the launcher chooses the
next available `127.0.0.1` port and prints the exact address to open.

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

## Windows .exe Launcher

`run_COBLAST.py` can also be packaged as a standalone Windows executable with
PyInstaller. The standalone executable bundles the Flask interface, Python
dependencies, templates, static assets, toy sample data, and the required BLAST+
executables. When launched, it extracts those bundled files to a temporary
runtime folder, stores persistent app data in a stable per-user data folder at
`%LOCALAPPDATA%\COBLAST_data`, starts Flask on `127.0.0.1`, and opens the
browser.
If port `5000` is already occupied by another local Flask/COBLAST session, the
`.exe` moves to the next available local port and prints the address in its
terminal window.

The standalone `.exe` does not bundle large user-created BLAST databases or
clinical datasets. Those remain local files chosen or created by the user.

For the current prototype, a prebuilt Windows executable is attached to the
latest GitHub Release (the **Releases** section of the repository page) rather
than committed to the repository tree. Download `COBLAST.exe` from there and run
it directly on Windows. Because it is an unsigned research prototype executable,
Windows SmartScreen or antivirus software may warn before first launch.

Build the executable from a clean checkout on Windows:

```powershell
cd 'C:\path\to\COBLAST-'
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip install pyinstaller
python build_standalone_exe.py --blast-bin 'C:\Tools\ncbi-blast-2.17.0+\bin'
```

The build helper writes the executable to:

```text
dist\COBLAST.exe
```

Test the executable before sharing it:

```powershell
.\dist\COBLAST.exe --check-only --skip-smoke --no-browser
```

Then run the app through the executable:

```powershell
.\dist\COBLAST.exe
```

By default the standalone executable stores app data in a stable per-user
folder at `%LOCALAPPDATA%\COBLAST_data`, so the same data is reused no matter
where the `.exe` lives. To override that location, set `COBLAST_DATA_DIR`:

```powershell
$env:COBLAST_DATA_DIR = 'C:\COBLAST_data'
.\dist\COBLAST.exe
```

For GitHub distribution, prefer attaching `COBLAST.exe` to a GitHub Release.
PyInstaller build folders and generated executables are ignored by Git by
default, because the bundled executable may exceed GitHub's normal per-file
repository size limit.

## Updating or Removing the Windows Test Build

The current Windows test build does not use a formal installer. A tester's
local COBLAST installation is simply the folder containing `COBLAST.exe`.
Persistent app data lives separately in a stable per-user folder:

```text
%LOCALAPPDATA%\COBLAST_data
```

(for example `C:\Users\<you>\AppData\Local\COBLAST_data`). Because this folder
is independent of where `COBLAST.exe` sits, registered databases and saved
results are preserved automatically across version updates — there is no longer
any need to copy a data folder next to the new executable.

To install a newer test version while keeping old registered databases and
saved results:

1. Close the COBLAST browser tab and the COBLAST terminal window.
2. Download and fully extract the new `release` folder. Replace the old
   `COBLAST.exe`, or run the new one from any location — its folder no longer
   holds your data.
3. (Optional) Back up the per-user data folder first, in case you want to roll
   back:

```powershell
Copy-Item -Recurse "$env:LOCALAPPDATA\COBLAST_data" "$env:LOCALAPPDATA\COBLAST_data_backup_2026-06-10"
```

4. Run a quick diagnostic from the new folder:

```powershell
.\COBLAST.exe --check-only --skip-smoke --no-browser
```

5. Start the new version:

```powershell
.\COBLAST.exe
```

The new version reads the same `%LOCALAPPDATA%\COBLAST_data` folder, so the
registry, generated BLAST database files, and saved result exports carry over
with no copy step. If databases were registered from external locations, those
external FASTA or BLAST database files must remain at the same paths, or the
database registry may show them as `missing`. Use the database-management page
to verify databases after updating.

If a tester used `COBLAST_DATA_DIR` to store data somewhere else, set
`COBLAST_DATA_DIR` again before launching the new version.

To remove a test version after confirming a newer one works:

1. Delete the folder containing the old `COBLAST.exe`.
2. Delete the per-user data folder `%LOCALAPPDATA%\COBLAST_data` only if its
   databases and results are no longer needed. Back it up first if in doubt.

COBLAST does not install Windows services, browser extensions, or system-wide
BLAST settings. Uninstalling the test build is therefore folder removal unless
the tester manually created shortcuts or environment variables.

## Manual Setup

Use these steps if you need to debug the installation manually.

### 1. Verify BLAST+

```powershell
& '..\ncbi-blast-2.17.0+\bin\blastn.exe' -version
```

Expected output includes `blastn: 2.17.0+`.

### 2. Create and activate a virtual environment

Install Python first if `python --version` does not work in your own PowerShell.

```powershell
cd 'C:\path\to\COBLAST-'
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

`requirements.txt` installs Flask plus Biopython, which is used for FASTA
validation and XML output parsing.

### 3. Run the backend smoke test

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

### 4. Run the Flask app

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

## Troubleshooting Python Setup

If `python run_COBLAST.py` fails with a Python version error, check your current
Python version:

```powershell
python --version
py -0p
```

Install Python 3.11 or newer if needed, then rerun:

```powershell
py -3.11 run_COBLAST.py
```

If virtual environment creation fails, confirm that the Python `venv` module is
available:

```powershell
py -3.11 -m venv .venv
```

If dependency installation fails, check internet access and rerun the launcher.
If BLAST+ is not found, either set `BLAST_BIN` or pass `--blast-bin`:

```powershell
$env:BLAST_BIN = 'C:\Tools\ncbi-blast-2.17.0+\bin'
python run_COBLAST.py
```

If the browser shows an error looking for BLAST+ under `C:\Program Files\NCBI`
after launching the standalone `.exe`, close older COBLAST/Flask terminal
windows and open the address printed by the newest `.exe` window. That usually
means the browser was still pointed at a stale local server on port `5000`.

If Windows says it cannot access the specified device, path, or file when the
`.exe` starts, the most common causes are:

- the `release` folder is still inside a ZIP preview
- the executable is blocked by Windows because it came from the internet
- antivirus or Windows SmartScreen has quarantined the unsigned prototype
- an organization-managed device blocks unsigned apps or bundled executables
- the app is being run from OneDrive, SharePoint, Teams, Outlook, a network
  share, or another protected location

Move the extracted `release` folder to `C:\COBLAST`, unblock the executable
from Properties if that option appears, and run:

```powershell
.\COBLAST.exe --check-only --skip-smoke --no-browser
```

If the diagnostic names a blocked path under a temporary `_MEI...` folder, then
Windows is likely blocking PyInstaller's bundled BLAST+ executables after
COBLAST extracts them. That usually requires whitelisting the prototype,
running on a less restricted test machine, or moving to a signed installer or
one-folder distribution for later testing.

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
- run button

Advanced settings are collapsed by default and include E-value, maximum target
sequences, word size, BLASTN task, minimum percent identity, output parser
format, timeout, and raw database details.

Searches use BLAST+'s own defaults (for example e-value 10 and 500 maximum
target sequences) unless a value is entered in Advanced settings, so the routine
form stays simple. Each search has a default wall-clock timeout of 3,600 seconds
(1 hour); the advanced timeout field can override it, up to a maximum of 3,600
seconds.

When a BLAST search, database verification, or database creation job is running,
the interface shows a waiting screen with elapsed time and status messages. For
large FASTA files, database creation may still take many minutes because the
current prototype waits for `makeblastdb` to finish before returning to the
database page.

The results page shows the BLAST command and parameters used for the run. The
results table reports the top BLAST hits with query ID, subject ID, subject
title, percent identity, alignment length, query coverage, E-value, and bit
score.
Each completed run is saved locally under `instance\results\` so the displayed
result table can be downloaded as CSV or TSV without rerunning BLAST.

The database-management page is available at:

```text
http://127.0.0.1:5000/databases
```

If the launcher printed a different local address, use that same address with
`/databases` added to the end.

From that page, users can view registered databases, check availability, see
database type, add an existing BLAST database, create a new BLAST database from
FASTA with `makeblastdb`, and remove a database from the registry without
deleting BLAST files.

## SRA Pilot and Batch Workflow

The prototype now includes two stepwise tools for exploring patient SRA-scale
questions locally before choosing a larger storage or compute strategy.

Open the SRA workbench at:

```text
http://127.0.0.1:5000/sra
```

The workbench scans local SRA project folders from:

- `COBLAST_SRA_DIR`, if set
- `SRA_DATA_DIR`, if set
- the `sra` subfolder of the per-user data directory in standalone/runtime-data
  mode — for the `.exe` that is `%LOCALAPPDATA%\COBLAST_data\sra` (for example
  `C:\Users\<you>\AppData\Local\COBLAST_data\sra`)
- a sibling `SRA_data` folder beside the source checkout, when present

For a clinician-supplied SRA, the current local-first prototype assumes the SRA
or derived FASTA stays in a stable local folder or controlled local data drive.
COBLAST does not upload patient libraries to NCBI, Galaxy, RNASEQ.COM, or a
commercial server. Those remain separate deployment decisions.

The SRA workbench can:

- fetch full SRA runs by accession: enter one or more run accessions
  (SRR/ERR/DRR) and COBLAST opens a terminal that runs `prefetch` then
  `fastq-dump` to download and convert each run to a complete FASTA in the
  scanned `sra` folder, when the SRA Toolkit is available
- list local `.sra`, FASTA, and BLAST database artifacts
- register an existing SRA-derived BLAST database prefix
- register all discovered SRA-derived BLAST database prefixes, or only selected
  prefixes, in bulk
- create a small pilot BLAST database from the first N records of an existing
  FASTA file
- use SRA Toolkit `fastq-dump` to convert a limited number of spots from a local
  `.sra` file into a pilot FASTA, when `SRA_TOOLKIT_BIN` is set or a sibling
  `sratoolkit` folder is present

### Getting SRA data onto your machine

The workbench can now fetch runs for you, or you can place data on disk
yourself. Every path lands data in a per-project subfolder of the scanned `sra`
folder (one folder per accession/patient, e.g.
`%LOCALAPPDATA%\COBLAST_data\sra\SRRxxxxxxxx\`).

**Path 0 — SRA-mart: fetch SRA runs from the workbench (recommended, needs the toolkit).**
In the SRA workbench's **SRA-mart** box, enter one or more run accessions
(SRR/ERR/DRR, separated by commas, spaces, or newlines) and press **Build &
Send Fetch Query**. COBLAST opens a new terminal that runs `prefetch` then a
full `fastq-dump --fasta` for each accession, downloading the run and converting
it to a complete FASTA under `...\sra\<ACCESSION>\`. This is Path 2 below,
automated. The run appears in the SRA Projects table as `fasta-ready` when the
terminal finishes; build the full database from `/databases` → Create from
FASTA. Notes:

- Only *run* accessions are accepted. Study, experiment, and sample accessions
  (SRP/SRX/SRS/PRJ...) are rejected so a whole multi-terabyte study is never
  queued by accident.
- Requires the SRA Toolkit. Set `SRA_TOOLKIT_BIN` to its `bin` folder or drop
  the unzipped `sratoolkit.<ver>` folder next to the `.exe`; the Fetch box is
  hidden until the toolkit is found.
- Downloads are kept on disk and re-discovered on every visit, so a run only
  downloads once per machine. A full run can be tens of GB — `prefetch` pulls a
  compressed `.sra` (a few GB) and `fastq-dump` expands it to the larger FASTA.

**Path 1 — browser only (small pilot runs, no toolkit).** In the NCBI SRA Run
Browser (`https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRRxxxxxxxx`),
open the **Reads** tab, filter to a manageable read count, and download a
**FASTA** subset. Save it as `...\sra\SRRxxxxxxxx\SRRxxxxxxxx.fasta`. The project
shows as `fasta-ready`; use **Create Pilot DB**. (This export is capped/slow — it
is a sample, not a whole run.) You can also download the raw `.sra` object from
the Run Browser's **Data access** tab, but converting it still needs Path 2.

**Path 2 — SRA Toolkit by hand (what Path 0 automates).** Browser FASTA export
does not scale to a whole run. Path 0 runs these steps for you; do them manually
only if you want to fetch outside the workbench. Install the NCBI SRA Toolkit,
then:

```powershell
prefetch SRRxxxxxxxx
fasterq-dump SRRxxxxxxxx --fasta --outdir "$env:LOCALAPPDATA\COBLAST_data\sra\SRRxxxxxxxx"
```

For a *full* eToL run build the database from the Databases page
(`/databases` → Create from FASTA, type `nucl`) rather than the pilot sampler,
which only takes the first N reads. To let the workbench's **Create Pilot FASTA**
button run `fastq-dump` on a local `.sra` for you, point COBLAST at the toolkit:
set `SRA_TOOLKIT_BIN` to its `bin` folder, or drop the unzipped
`sratoolkit.<ver>` folder next to the `.exe`.

The intended first simulation is deliberately small:

1. Put one SRA project under a scanned folder such as the following (or use
   **SRA-mart** to download one by accession):

   ```text
   %LOCALAPPDATA%\COBLAST_data\sra\patient_001\patient_001.sra
   ```

2. If a full FASTA already exists, use `Create Pilot DB` with a small record
   count such as `1000`.

3. If only `.sra` exists, use `Create Pilot FASTA` first, then create a pilot
   database from the generated FASTA after the page refreshes.

4. Run a normal BLAST query against that pilot database.

For many patients, open the batch BLAST page at:

```text
http://127.0.0.1:5000/batch-blast
```

Batch BLAST runs one query against multiple selected registered databases
concurrently — the databases are searched in parallel, auto-sized to the
machine's cores, with the number searched at once overridable via the Advanced
"Max databases at once" field. This is the prototype path for testing the
"100 patients" problem:
prepare or register each patient as a local nucleotide BLAST database, select
the compatible databases, and run the batch. The batch page includes compatible
database filtering and select-all/deselect-all controls. Completed batch runs
save individual per-database results, expose per-database CSV links, and export
the aggregate raw hit table as CSV or TSV.

The batch page also includes an APOE exact-match probe preset. When selected,
COBLAST uses the four stored APOE probe sequences, runs BLASTN against the
selected nucleotide databases, and saves only hits with 100% identity and 100%
query coverage. APOE batch results include an `APOE Probe Summary` table with
one row per selected sample/database. The summary counts exact matches for
`AE4=C`, `AE4=T`, `AE2=C`, and `AE2=T`, then reports `% C<->T` as:

```text
(AE4=T hits + AE2=T hits) / total exact APOE probe hits * 100
```

When an SRX, ERX, or DRX accession appears in the database display name or
database path, COBLAST uses that accession as the sample label; otherwise it
falls back to the database display name. APOE summary tables can be downloaded
as CSV or TSV, and the underlying raw exact-hit table remains available for
troubleshooting.

### eToL probe presets (optional)

The batch page also includes optional **eToL probe presets** for the electronic
Tree of Life microbiome workflow of Hu, Haas & Lathe (*BMC Microbiology*
2022;22:317). Each preset BLASTNs a stored probe panel — the full 1,017-probe
microbial set (**eToL Full**), a one-probe-per-species quick set (**eToL Quick**),
or the viral structural-protein panel (**eToL-V**) — against the selected
nucleotide databases using BLAST's permissive default-megablast *net* (E-value <
0.01, no identity or coverage filter), then counts matched reads per probe and
species with CSV/TSV exports. Each panel is searched together with its
housekeeping control probes for host-cell normalization; the controls are
appended automatically and are not a separately selectable preset. The presets
offer an optional secondary human filter that removes host-derived reads. Only
one preset can be active at a time.

The eToL presets can also assemble matched reads into contigs, identify
them, and optionally **re-probe** (a second pass that uses contigs as probes to
recover reads the net missed) — leave re-probing off for routine and cohort runs
and characterise it once per data type before relying on it; see
**[docs/eToL.md](docs/eToL.md)** for when to use it.

Contig assembly, identification, and re-probing use the **CAP3** assembler. For
licensing reasons CAP3 is **not** distributed with COBLAST+; instead install
[Unipro UGENE](https://ugene.net/download-all.html) in its default location —
UGENE bundles CAP3, and COBLAST+ auto-detects it under
`…\Unipro UGENE\tools\cap3` with no further setup. So the full workflow is just:
install UGENE, then run COBLAST+. CAP3 is optional — the BLAST net, per-probe and
per-species counts, and all exports work without it, and a run that requests
assembly with no CAP3 present continues and reports that contigs were skipped. If
CAP3 lives somewhere non-standard, point COBLAST+ at its folder with `CAP3_BIN`.

These presets are an optional analysis layered on the general batch workflow
above, not a required part of using COBLAST+. For the full details — the net, the
two-pass megablast/blastn-short split, cross-probe de-duplication, host-cell
normalization, the secondary human filter, and the export formats — see
**[docs/eToL.md](docs/eToL.md)**.

## Adding Clinician Databases

In the intended local-use workflow, the clinician's sequencing data can be made
into a local BLAST database. The clinician can then paste or upload one query
sequence at a time, choose the compatible database, and run BLAST locally.

Before adding data, decide whether the database is nucleotide or protein:

- use `nucl` for DNA/RNA sequencing FASTA files and run `blastn` or `tblastn`
- use `prot` for protein FASTA files and run `blastp` or `blastx`

The input database file must be FASTA. If the sequencing data are still in
FASTQ format, convert or export them to FASTA before using the current
prototype.

Useful field meanings:

- `Display name`: the readable name shown in the database dropdown
- `Database type`: `nucl` or `prot`; this controls which BLAST programs can use it
- `Category`: broad grouping such as `human`, `viral`, `eToL-V`, `toy`, or `custom`
- `Source FASTA path`: the local path to the FASTA file used to create the database
- `BLAST database prefix path`: the path prefix passed to BLAST with `-db`
- `Description`: short text shown to help users choose the right database

The BLAST database prefix is not usually a single visible file. For example, if
`makeblastdb` creates files named:

```text
C:\COBLAST_data\databases\patient_001_reads.nhr
C:\COBLAST_data\databases\patient_001_reads.nin
C:\COBLAST_data\databases\patient_001_reads.nsq
```

then the database prefix path is:

```text
C:\COBLAST_data\databases\patient_001_reads
```

### Create a New Database From a FASTA File

Use this when the clinician has sequencing data as a FASTA file and wants this
prototype to run `makeblastdb` for them.

1. Save the FASTA file somewhere local and stable, for example:

   ```text
   C:\COBLAST_data\input\patient_001_reads.fasta
   ```

2. Start COBLAST and open the local browser address printed by the launcher.

3. Open the database-management page:

   ```text
   http://127.0.0.1:5000/databases
   ```

4. Expand `Create BLAST Database From FASTA`.

5. Enter a clear `Display name`, for example:

   ```text
   Patient 001 sequencing reads
   ```

6. Choose the correct `Database type`:

   ```text
   nucl
   ```

   Use `nucl` for DNA/RNA sequencing reads. Use `prot` only if the FASTA
   contains protein sequences.

7. Choose a `Category`, usually `human` or `custom` for clinician-provided
   sequencing data.

8. Paste the full FASTA path into `Source FASTA path`, for example:

   ```text
   C:\COBLAST_data\input\patient_001_reads.fasta
   ```

9. Leave `Output database prefix path` blank unless you want to choose the exact
   database location. If left blank, COBLAST stores the database under its local
   managed database folder.

10. Add a short `Description`, for example:

    ```text
    Local nucleotide database created from patient 001 sequencing reads.
    ```

11. Select `Create Database`.

12. Confirm that the new row appears under `Registered Databases` with status
    `available`. If needed, select `Verify` for that database or `Verify All`.

13. Return to `Run BLAST`, choose a compatible search type, choose the new
    database from the dropdown, paste or upload the query FASTA, and select
    `Run BLAST`.

### Add an Existing BLAST Database

Use this when a database has already been created outside COBLAST with
`makeblastdb` or another BLAST+ database-preparation workflow.

1. Locate the BLAST database files on the local machine.

   Nucleotide databases usually include files such as:

   ```text
   .nhr
   .nin
   .nsq
   ```

   Protein databases usually include files such as:

   ```text
   .phr
   .pin
   .psq
   ```

2. Identify the database prefix path by removing those BLAST file extensions.
   For example, these files:

   ```text
   D:\BLAST_Databases\viral_panel.nhr
   D:\BLAST_Databases\viral_panel.nin
   D:\BLAST_Databases\viral_panel.nsq
   ```

   should be registered with this prefix:

   ```text
   D:\BLAST_Databases\viral_panel
   ```

3. Start COBLAST and open the database-management page:

   ```text
   http://127.0.0.1:5000/databases
   ```

4. Expand `Add Existing BLAST Database`.

5. Enter a `Display name`, choose the correct `Database type`, and choose a
   `Category`.

6. Paste the prefix into `BLAST database prefix path`. Do not paste the `.nin`,
   `.nsq`, `.pin`, or `.psq` file path itself.

7. If known, paste the original FASTA path into `Source FASTA path`. This field
   is helpful for record keeping but is optional when registering an existing
   database.

8. Add a short `Description` that will help the clinician choose the database
   later.

9. Select `Add Database`.

10. Confirm that the database status is `available`. If it is `missing` or
    `invalid`, check that the prefix path is correct and that the database type
    matches the files being registered.

Removing a database from the database-management page removes it from the
COBLAST registry only. It does not delete the FASTA file or the BLAST database
files from disk. Removed toy databases stay removed instead of being seeded
again on the next page load. To restore any removed database, add its existing
BLAST prefix again through `Add Existing BLAST Database`. The database page's
`Remove All Missing` action removes every entry currently marked `missing`
while leaving `available` and `invalid` entries registered.

## Current supported BLAST programs

- `blastn`: nucleotide query vs nucleotide database
- `blastp`: protein query vs protein database
- `blastx`: translated nucleotide query vs protein database
- `tblastn`: protein query vs translated nucleotide database

The interface enforces these compatibility rules by filtering registered
databases by `nucl` or `prot` type after the BLAST program is selected.

By default, `blastn` runs with the `megablast` task — the same default as
command-line BLAST+ — so a routine search matches what the clinician would get
from the NCBI `blastn` executable. For short queries, the `blastn-short` task can
be selected under Advanced settings (the exact-match probe presets choose it
automatically).

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

No license has been chosen yet, so this repository reserves all rights by
default. Add a `LICENSE` file with a formal license before wider distribution.
