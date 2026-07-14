"""One native folder dialog, shared by the launcher and the Settings page.

The Settings page cannot open a picker on its own — a browser never hands a page
an absolute path, by design. But COBLAST's server *is* the user's machine (Flask
binds 127.0.0.1), so the server opens the dialog and hands the path back to the
page. The launcher reuses the same call for its first-run prompt.

In every caller the dialog is an assist, not the input: the path field stays
typable, so a picker that fails or is unavailable never blocks anyone.
"""

from __future__ import annotations


def ask_directory(title: str, initialdir: str = "", parent=None) -> str | None:
    """Show the native folder dialog; return the chosen path, or None if cancelled.

    Raises RuntimeError when tkinter is missing (headless/stripped build) so
    callers can tell "no dialog available" from "user chose nothing" and point
    the user at the text field instead.
    """
    try:
        import tkinter as tk
        from tkinter import filedialog
    except Exception as exc:  # no tk in this build/session
        raise RuntimeError(f"Folder picker unavailable ({exc}).") from exc

    own_root = parent is None
    root = tk.Tk() if own_root else parent
    try:
        if own_root:
            root.withdraw()
        # Without this the dialog opens behind the browser window (Settings) or
        # the console window (launcher), and looks like nothing happened.
        root.attributes("-topmost", True)
        chosen = filedialog.askdirectory(
            title=title,
            initialdir=initialdir or "",
            mustexist=False,
            parent=root,
        )
    finally:
        if own_root:
            try:
                root.destroy()
            except Exception:
                pass
    return chosen or None
