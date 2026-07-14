"""The Settings Browse button must never be the only way in.

Cancel and a missing dialog both have to leave the typed path field usable, so
the route reports them instead of writing a bogus path or 500-ing.
"""

import folder_picker
from app import app


def test_browse_returns_picked_path(monkeypatch):
    monkeypatch.setattr(folder_picker, "ask_directory", lambda *a, **kw: r"D:\COBLAST_data")
    payload = app.test_client().post("/settings/browse").get_json()
    assert payload == {"path": r"D:\COBLAST_data"}


def test_browse_cancel_leaves_field_alone(monkeypatch):
    monkeypatch.setattr(folder_picker, "ask_directory", lambda *a, **kw: None)
    assert app.test_client().post("/settings/browse").get_json() == {"path": ""}


def test_browse_without_tkinter_tells_user_to_type(monkeypatch):
    def no_tk(*a, **kw):
        raise RuntimeError("Folder picker unavailable (no tkinter).")

    monkeypatch.setattr(folder_picker, "ask_directory", no_tk)
    payload = app.test_client().post("/settings/browse").get_json()
    assert "Type the folder path" in payload["error"]
