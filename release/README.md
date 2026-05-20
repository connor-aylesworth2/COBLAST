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
exports are stored beside the executable in `COBLAST_data`.

This executable is unsigned and intended for agreed prototype testing only.

## Updating from an Older Test Version

COBLAST does not use a formal installer. The test installation is the folder
containing `COBLAST.exe` and the persistent `COBLAST_data` folder beside it.

To install a newer version and keep old databases/results:

1. Close COBLAST completely.
2. Back up the old `COBLAST_data` folder:

```powershell
Copy-Item -Recurse .\COBLAST_data .\COBLAST_data_backup_2026-05-20
```

3. Download and fully extract the new `release` folder.
4. Copy the old `COBLAST_data` folder beside the new `COBLAST.exe`, for example:

```text
C:\COBLAST\COBLAST.exe
C:\COBLAST\COBLAST_data\
```

5. Run this diagnostic from the new folder:

```powershell
.\COBLAST.exe --check-only --skip-smoke --no-browser
```

6. Start the new version:

```powershell
.\COBLAST.exe
```

Keeping `COBLAST_data` should preserve COBLAST-managed databases, the database
registry, and saved result exports. If a database was registered from an
external folder, keep that external FASTA or BLAST database path unchanged, then
use the database-management page to verify it after updating.

To remove an old version, back up `COBLAST_data` if needed, then delete the old
folder containing `COBLAST.exe`. Delete `COBLAST_data` only if the old
databases and results are no longer needed.

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
