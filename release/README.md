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

## Optional: Contig Assembly (eToL re-probing)

The microbial eToL presets can assemble matched reads into contigs and re-probe
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
across version updates.

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
