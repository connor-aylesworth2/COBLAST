"""The first-run prompt must accept a typed path without the native dialog.

Testers lost the folder dialog mid-scroll and had no other way to give a path, so
the window's own field is now the input. This drives that field end to end: build
the real window, press "Use this folder", and check the launcher gets the path.
Never calls filedialog, which is exactly the failure being guarded against.
"""

import pytest

tkinter = pytest.importorskip("tkinter")

import run_COBLAST


def _widgets(parent):
    for child in parent.winfo_children():
        yield child
        yield from _widgets(child)


def test_typed_path_is_accepted_without_the_dialog(monkeypatch, tmp_path):
    try:
        tkinter.Tk().destroy()  # headless CI has no display; _prompt swallows the error
    except tkinter.TclError as exc:
        pytest.skip(f"no display: {exc}")

    data_dir = tmp_path / "COBLAST_data"
    monkeypatch.setenv("COBLAST_DATA_DIR", str(data_dir))  # what the field prefills with

    def press_use_this_folder(self):
        """Stand in for mainloop: click the button a user would click."""
        for widget in _widgets(self):
            if isinstance(widget, tkinter.Button) and widget.cget("text") == "Use this folder":
                widget.invoke()
                return
        raise AssertionError("the prompt has no 'Use this folder' button")

    monkeypatch.setattr(tkinter.Misc, "mainloop", press_use_this_folder)
    chosen = run_COBLAST._prompt_for_data_dir()
    assert chosen == data_dir.resolve()
