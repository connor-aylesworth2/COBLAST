"""Minimal check that the batch progress endpoint mirrors the in-memory registry."""

from app import app, _batch_progress, _batch_progress_lock


def test_batch_progress_endpoint():
    client = app.test_client()

    # Unknown job ids report zeros so the waiting page stays neutral.
    assert client.get("/batch-progress/nope").get_json() == {
        "done": 0,
        "total": 0,
        "stages": [],
    }

    with _batch_progress_lock:
        _batch_progress["job-1"] = {
            "done": 2,
            "total": 5,
            "stages": {"DB A": {"stage": "Assembling contigs (CAP3)", "since": 1000.0}},
        }
    try:
        payload = client.get("/batch-progress/job-1").get_json()
        assert payload["done"] == 2 and payload["total"] == 5
        assert payload["stages"] == [
            {"label": "DB A", "stage": "Assembling contigs (CAP3)", "since": 1000.0}
        ]
    finally:
        with _batch_progress_lock:
            _batch_progress.pop("job-1", None)


if __name__ == "__main__":
    test_batch_progress_endpoint()
    print("ok")
