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
