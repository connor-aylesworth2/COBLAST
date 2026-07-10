# COBLAST Prototype Executable

`COBLAST.exe` is a standalone Windows prototype build. It bundles the local
Flask interface, Python dependencies, toy sample data, and the required NCBI
BLAST+ executables.

Run:

```powershell
.\COBLAST.exe
```

The app starts a local Flask server bound to `127.0.0.1` and opens the browser.
Runtime files, generated toy BLAST databases, registered databases, and result
exports are stored in a stable per-user data folder at `%LOCALAPPDATA%\COBLAST_data`
(for example `C:\Users\<you>\AppData\Local\COBLAST_data`).

This executable is unsigned and intended for agreed prototype testing only.

## Choosing Where COBLAST+ Stores Data

COBLAST+ writes everything it builds — BLAST databases, SRA downloads, search
results, the database registry, and temporary scratch files — under a single
data folder. Databases and SRA runs can reach tens or hundreds of GB, so you can
put that folder on any drive with room to spare instead of your system (`C:`)
drive.

- **First launch:** a folder picker asks where to store COBLAST+ data. Choose a
  folder on the drive you want (for example `D:\COBLAST_data`); COBLAST+
  remembers it for every future launch. If you skip the picker, COBLAST+ uses
  the default `%LOCALAPPDATA%\COBLAST_data` on `C:`.
- **Change it later:** open the **Settings** page from the top navigation bar
  (or `http://127.0.0.1:5000/settings`), enter a new folder, and save. The
  change takes effect the next time you start COBLAST+.
- **No spaces in the path.** BLAST+ cannot build databases under a folder whose
  path contains spaces (for example `D:\My Data`). Choose a space-free folder
  such as `D:\COBLAST_data`; COBLAST+ rejects a spaced path and asks again.
- **Existing data is not moved.** Switching folders starts a fresh database
  registry at the new location. Databases you built at the old location stay on
  disk and keep working, but they are not listed until you re-add them on the
  Databases page (Add Existing BLAST Database). Their files are never deleted.
- **Temporary files follow the data folder too**, so large `fastq-dump`/CAP3/
  BLAST scratch during big runs no longer fills your system drive.

Advanced/scripted use: launch `COBLAST.exe --data-dir D:\COBLAST_data` to set
the location without the picker, or `COBLAST.exe --pick-data-dir` to force the
picker to reappear. Setting the `COBLAST_DATA_DIR` environment variable
overrides the saved location for one session.

## Optional: Contig Assembly (eToL re-probing)

The eToL presets can assemble matched reads into contigs and re-probe
with them. Those steps use the **CAP3** assembler, which is not shipped with
COBLAST for licensing reasons — you get it by installing Unipro UGENE, which
bundles CAP3. Everything else (BLAST searches, per-species counts, exports) works
without CAP3; a run that requests assembly without it simply reports that contigs
were skipped, so UGENE is only needed if you want contig assembly or re-probing.

To enable it:

1. Download the Windows **installer** (not the portable ZIP) from
   <https://ugene.net/download-all.html> and run it. Accept the default install
   location. COBLAST only auto-detects UGENE under `C:\Program Files`, so a
   portable ZIP extracted elsewhere (Downloads, Desktop) will not be found.
2. Confirm CAP3 landed at
   `C:\Program Files\Unipro UGENE\tools\cap3\cap3.exe` (both 64-bit and 32-bit
   `Program Files` are checked).
3. Start `COBLAST.exe` — no configuration needed; COBLAST finds CAP3 there
   automatically.
4. To verify detection, run a microbial eToL preset with assembly enabled. The
   results page shows a `Contigs assembled (CAP3)` row when CAP3 was found, or
   `Contig assembly: CAP3 not found - assembly skipped` when it was not.

If you installed UGENE somewhere non-standard, or use the portable ZIP, point
COBLAST at the folder that holds `cap3.exe` before launching:

```powershell
$env:CAP3_BIN = 'C:\path\to\Unipro UGENE\tools\cap3'
.\COBLAST.exe
```

## Updating from an Older Test Version

COBLAST does not use a formal installer. The test installation is the folder
containing `COBLAST.exe`. Persistent app data lives separately in a stable
per-user folder at `%LOCALAPPDATA%\COBLAST_data` (for example
`C:\Users\<you>\AppData\Local\COBLAST_data`). Because that folder is independent
of where `COBLAST.exe` sits, databases and results carry over automatically
across version updates. If you chose a custom data folder (the first-run picker
or the Settings page), the new version reads that same saved location
automatically — there is nothing to re-select.

To install a newer version and keep old databases/results:

1. Close COBLAST completely.
2. Download and fully extract the new `release` folder. Replace the old
   `COBLAST.exe`, or run the new one from any location.
3. (Optional) Back up the per-user data folder first, in case you want to roll
   back:

```powershell
Copy-Item -Recurse "$env:LOCALAPPDATA\COBLAST_data" "$env:LOCALAPPDATA\COBLAST_data_backup_2026-06-10"
```

4. Run this diagnostic from the new folder:

```powershell
.\COBLAST.exe --check-only --skip-smoke --no-browser
```

5. Start the new version:

```powershell
.\COBLAST.exe
```

The new version reads the same `%LOCALAPPDATA%\COBLAST_data` folder, so
COBLAST-managed databases, the database registry, and saved result exports are
preserved with no copy step. If a database was registered from an external
folder, keep that external FASTA or BLAST database path unchanged, then use the
database-management page to verify it after updating.

To remove a version, delete the folder containing `COBLAST.exe`. Delete the
per-user `%LOCALAPPDATA%\COBLAST_data` folder only if its databases and results
are no longer needed.

If Windows says it cannot access the specified device, path, or file:

1. Make sure the `release` folder has been extracted. Do not run `COBLAST.exe`
   from inside the ZIP preview.
2. Move the extracted folder to a simple local path such as `C:\COBLAST`.
3. Right-click `COBLAST.exe`, choose Properties, and select Unblock if Windows
   shows that checkbox.
4. Try the diagnostic command from PowerShell:

```powershell
.\COBLAST.exe --check-only --skip-smoke --no-browser
```

If the diagnostic succeeds but the browser still closes, run:

```powershell
.\COBLAST.exe --skip-smoke --no-browser
```

Then open the printed `http://127.0.0.1:...` address manually in the browser.
