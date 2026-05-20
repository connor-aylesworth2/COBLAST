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
