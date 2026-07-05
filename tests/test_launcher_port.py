"""Readiness check that gates browser launch on Flask actually listening."""

import socket
import threading
import time

from run_COBLAST import wait_for_port


def test_returns_true_once_something_listens():
    listener = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    listener.bind(("127.0.0.1", 0))
    port = listener.getsockname()[1]
    listener.listen(1)
    try:
        assert wait_for_port(port, timeout=2.0) is True
    finally:
        listener.close()


def test_returns_false_for_dead_port_within_timeout():
    # Grab a free port then close it so nothing is listening there.
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.bind(("127.0.0.1", 0))
    port = sock.getsockname()[1]
    sock.close()

    start = time.monotonic()
    assert wait_for_port(port, timeout=0.5) is False
    assert time.monotonic() - start < 3.0  # honors the timeout, no hang


def test_waits_for_a_late_starting_server():
    port_holder = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    port_holder.bind(("127.0.0.1", 0))
    port = port_holder.getsockname()[1]
    port_holder.close()

    def start_late():
        time.sleep(0.5)
        late = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        late.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        late.bind(("127.0.0.1", port))
        late.listen(1)
        time.sleep(2.0)
        late.close()

    threading.Thread(target=start_late, daemon=True).start()
    assert wait_for_port(port, timeout=5.0) is True
